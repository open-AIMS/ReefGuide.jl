"""Geometry-based assessment methods."""

"""
    assess_reef_site(
        rel_pix::DataFrame,
        rotated::Vector{GI.Wrappers.Polygon},
        suitability_threshold::Float64
    )::Tuple{Float64,Int64,GI.Wrappers.Polygon,Int64}

Assesses the rotations of a search box `geom` for their suitability score (calculated as the
proportion of pixels that meet all specified criteria thresholds).

# Extended help
The scores produced are a proportion of the polygon that are covered by valid pixel points,
relative to the maximum number of points (`max_count`) (0-1). `max_count` is approximate
(determined by user `x_dist` and `y_dist` box dimensions) and doesn't account for buffering
and rotation of the search box. In rare cases scores could be > 1, however returned values
are capped at max 1.

# Arguments
- `rel_pix` : The point data for relevant pixels that are within the search area around a pixel.
- `rotated` : Pre-rotated geometries.
- `suitability_threshold` : Suitability threshold, below which sites are excluded from result sets.

# Returns
- Highest score
- Highest scoring rotation step
- Highest scoring polygon
- Quality control flag for site, indicating if `suitability_threshold` was breached.
"""
function assess_reef_site(
    rel_pix::DataFrame,
    rotated::Vector{GI.Wrappers.Polygon},
    suitability_threshold::Float64
)::Tuple{Float64,Int64,GI.Wrappers.Polygon,Int64}
    # Implementation with pre-rotations
    n_rotations = length(rotated)
    score = @MVector zeros(n_rotations)
    qc_flag = @MVector zeros(Int64, n_rotations)

    for (j, r) in enumerate(rotated)
        score[j] = count(GO.contains.(Ref(r), rel_pix.geometry))

        if score[j] < suitability_threshold
            # Early exit as there's no point in searching further.
            # Changing the rotation is unlikely to improve the score.
            qc_flag[j] = 1
            break
        end

        if score[j] >= (nrow(rel_pix) * 0.98)
            # Found near the best possible score so might as well exit
            break
        end
    end

    best_idx = argmax(score)
    return (
        score[best_idx],
        best_idx,
        rotated[best_idx],
        qc_flag[best_idx]
    )
end

"""
    find_optimal_site_alignment(
        lookup_tbl::DataFrame,
        res::Float64,
        x_dist::Union{Int64,Float64},
        y_dist::Union{Int64,Float64};
        degree_step::Float64=10.0,
        suit_threshold::Float64=0.7
    )::DataFrame

Identify the most suitable site polygons for each pixel in the `search_pixels` DataFrame.
`x_dist` and `y_dist` are x and y lengths of the search polygon in meters. A buffer of the
raster files' resolution is applied to the search box. And angle from a pixel to a reef edge
is identified and used for searching with custom rotation parameters.
Method is currently opperating for CRS in degrees units.

# Arguments
- `lookup_tbl` : DataFrame containing environmental variables for assessment.
- `res` : Resolution of the original raster pixels. Can by found via `abs(step(dims(raster, X)))`.
- `x_dist` : Length of horizontal side of search box (in meters).
- `y_dist` : Length of vertical side of search box (in meters).
- `degree_step` : Degree to perform rotations around identified edge angle.
- `suit_threshold` : Theshold used to skip searching where the proportion of suitable pixels is too low.

# Returns
DataFrame containing highest score, rotation and polygon for each assessment at pixels in indices.
"""
function find_optimal_site_alignment(
    lookup_tbl::DataFrame,
    res::Float64,
    x_dist::Union{Int64,Float64},
    y_dist::Union{Int64,Float64};
    target_crs=EPSG(4326),
    degree_step::Float64=10.0,
    suit_threshold::Float64=0.7
)::DataFrame
    max_count = (
        (x_dist * y_dist) / ceil(degrees_to_meters(res, mean(lookup_tbl.lats)))^2
    )

    # If there are no lookup results - don't try to process polygons
    if nrow(lookup_tbl) == 0
        return DataFrame(;
            score=Float64[],
            orientation=Int64[],
            qc_flag=Int64[],
            geometry=GI.Wrappers.Polygon[]
        )
    end

    search_box = initial_search_box(
        (lookup_tbl.lons[1], lookup_tbl.lats[1]),
        x_dist,
        y_dist,
        target_crs
    )

    # Precalculate rotations
    rotations = 0.0:degree_step:179.0
    rotated_geoms = Vector{GI.Wrappers.Polygon}(undef, length(rotations))
    for (j, r) in enumerate(rotations)
        rotated_geoms[j] = rotate_polygon(search_box, r)
    end

    # Create KD-tree to identify pixels for assessment that are close to each other.
    # We can then apply a heuristic to avoid near-identical assessments.
    @debug "$(now()) : Applying KD-tree filtering on $(nrow(lookup_tbl)) locations"
    inds = Matrix{Float32}([lookup_tbl.lons lookup_tbl.lats]')
    n_locs = size(inds, 2)
    kdtree = KDTree(inds; leafsize=25)

    # Deterministic per-index computation: for each location, find its group of
    # near-identical neighbours and the neighbour closest to the group's center.
    # This step has no cross-iteration dependency, so it can be parallelized freely -
    # unlike the dedup decision below, whose outcome depends on the order locations
    # are visited in.
    group_to_ignore = Vector{Vector{Int64}}(undef, n_locs)
    Threads.@threads for i in 1:n_locs
        coords = inds[:, i]

        # If there are a group of pixels close to each other, only assess the one closest to
        # the center of the group.
        idx, dists = knn(kdtree, coords, 30)  # retrieve 30 closest locations

        # Select current pixel and those ~50% of max shape length away
        x_d = meters_to_degrees(x_dist, coords[2])
        y_d = meters_to_degrees(y_dist, coords[2])
        max_l = max(x_d, y_d) * 0.5
        sel = (dists .<= max_l)

        # Find index of pixel closest to the center of the group
        xs = mean(inds[1, idx[sel]])
        ys = mean(inds[2, idx[sel]])
        closest_idx = argmin(
            map(p2 -> sum((p2 .- (xs, ys)) .^ 2), eachcol(inds[:, idx[sel]]))
        )

        to_keep = idx[sel][closest_idx]
        group_to_ignore[i] = idx[sel][idx[sel] .!= to_keep]
    end

    # Sequential dedup pass: visit locations in a fixed order and accumulate the
    # ignore set. A location already marked ignored by an earlier group is skipped
    # rather than allowed to veto a different group - this is what makes the result
    # independent of `JULIA_NUM_THREADS`/thread scheduling.
    ignored = falses(n_locs)
    for i in 1:n_locs
        ignored[i] && continue
        for j in group_to_ignore[i]
            ignored[j] = true
        end
    end
    ignore_locs = findall(ignored)

    # Create search tree
    tree = STRT.STRtree(lookup_tbl.geometry)

    # Search each location to assess
    assessment_locs = lookup_tbl[Not(ignore_locs), :]
    n_pixels = nrow(assessment_locs)
    @debug "$(now()) : KD-tree filtering - removed $(length(ignore_locs)) near-identical locations, now assessing $(n_pixels) locations"

    best_score = zeros(n_pixels)
    best_poly = Vector{GI.Wrappers.Polygon}(undef, n_pixels)
    best_rotation = zeros(Int64, n_pixels)
    quality_flag = zeros(Int64, n_pixels)

    FLoops.assistant(false)
    @floop for (i, pix) in enumerate(eachrow(assessment_locs))
        lon = pix.lons
        lat = pix.lats

        rotated_copy::Vector{GI.Wrappers.Polygon} =
            move_geom.(
                rotated_geoms,
                Ref((lon, lat))
            )

        # Find pixels within each rotated search area
        in_pixels = unique(reduce(vcat, STRT.query.(Ref(tree), rotated_copy)))
        relevant_pixels = lookup_tbl[in_pixels, :]
        n_matches = nrow(relevant_pixels)

        # Skip if no relevant pixels or if the number of suitable pixels are below required
        # threshold
        if n_matches == 0 || ((n_matches / max_count) < suit_threshold)
            best_score[i] = 0.0
            best_rotation[i] = 0
            best_poly[i] = rotated_copy[1]
            quality_flag[i] = 1
            continue
        end

        best_score[i], best_rotation[i], best_poly[i], quality_flag[i] = assess_reef_site(
            relevant_pixels,
            rotated_copy,
            suit_threshold
        )
    end

    if maximum(best_score) > 0.0
        best_score = min.(best_score ./ max_count, 1.0)
    end

    return DataFrame(;
        score=best_score,
        orientation=best_rotation,
        qc_flag=quality_flag,
        geometry=best_poly
    )
end

"""
    assess_location_quality(
        lookup_tbl::DataFrame,
        assessment_idx::BitVector,
        res::Float64;
        target_crs=EPSG(4326)
    )::Vector{Int8}

Score each pixel based on how many of its immediate neighbours meet criteria.

Note: Score would nominally be out of nine - the total number of cells being the centroid
      +/- 1 pixel in all directions - however it may be slightly more or less due to
      raster resolution and approximations as the area is selected with a square polygon.

# Arguments
- `lookup_tbl` : Lookup table holding all data for the given region
- `assessment_idx` : True/false indicating which entry in the lookup table to assess
- `res` : Raster resolution (assumed to be in decimal degrees)
- `target_crs` : EPSG code of spatial data

# Returns
Indicative percent values (0 - 100) for each assessed pixel.
Results are stored as Int8 to reduce memory use.
"""
function assess_location_quality(
    lookup_tbl::DataFrame,
    assessment_idx::BitVector,
    res::Float64;
    target_crs=EPSG(4326)
)::Vector{Int8}
    # Create square search box
    res = degrees_to_meters(res, lookup_tbl.lats[1])
    x_dist = ceil(Int64, res * 3)  # Search area immediate around target pixel
    y_dist = x_dist
    @debug "Initial search box: $(x_dist)m * $(y_dist)m with res: $(round(res; digits=2))m^2"
    search_box = initial_search_box(
        (lookup_tbl.lons[1], lookup_tbl.lats[1]),
        x_dist * 0.9, # Approach uses the "touches" algorithm to select pixels so reduce
        y_dist * 0.9, # size of polygon to ensure only relevant pixels are "touched"
        target_crs
    )

    max_count = floor(Int64, (x_dist * y_dist) / res^2)

    assessment_locs = lookup_tbl[assessment_idx, :]

    # Create tree specific to search area
    time_taken = @elapsed begin
        tree = STRT.STRtree(assessment_locs.geometry)
    end
    @debug "Took $(time_taken) seconds to create search-specific tree"

    n_pixels_to_assess = nrow(assessment_locs)
    results = zeros(Int8, n_pixels_to_assess)  # Percent value 0 - 100

    @debug "$(now()) : Assessment start"
    Threads.@threads for i in 1:n_pixels_to_assess
        pix = assessment_locs[i, :]

        moved_box::GI.Wrappers.Polygon = move_geom(
            search_box,
            (pix.lons, pix.lats)
        )

        # Find pixels the search area
        pixel_idx = STRT.query(tree, moved_box)
        n_matches = length(pixel_idx)

        # Skip if no relevant pixels or if the number of suitable pixels are below required
        # threshold
        if n_matches == 0
            continue
        end

        # Taking floor to be conservative in estimates
        results[i] = floor(Int8, (n_matches / max_count) * 100)
    end
    @debug "$(now()) : Assessment finished: $(count(results .> 0.0)) / $(n_pixels_to_assess) met criteria"

    if count(results .> 0.0) == 0
        @warn "Empty raster result!"
    end

    # Cap to 0 - 100
    results = clamp!(results, Int8(0), Int8(100))

    return results
end
