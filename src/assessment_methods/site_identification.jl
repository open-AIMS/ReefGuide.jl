"""Methods for identifying potential deployment locations."""

"""
Minimum value a genuinely-assessed pixel can be stored as in the WLC score raster
(see [`assess_region_score`](@ref)). Well below the resolution anyone would read off
a `[0,1]` heatmap, but strictly greater than the raster's `missingval=0`/implicit
sparse-zero background, so a real "worst possible score" pixel stays visually and
numerically distinct from "not part of this region's valid footprint at all".
"""
const SCORE_RASTER_FLOOR = 1.0f-4

"""
    proportion_suitable(
        x::Union{BitMatrix,SparseMatrixCSC{Bool,Int64}}; square_offset::Tuple=(-4, 5)
    )::SparseMatrixCSC{Int8,Int64}

Calculate the the proportion of the subsection that is suitable for deployments.
The `subsection` is the surrounding a rough hectare area centred on each cell of a raster
marked as being suitable according to user-selected criteria.

Cells on the edges of a raster object are assessed using a smaller surrounding area, rather
than shifting the window inward. In usual applications, there will be no target pixel
close to the edge due to the use of buffer areas.

# Arguments
- `x` : Matrix of boolean pixels after filtering with user criteria.
- `square_offset` : The number of pixels +/- around a center "target" pixel to assess as the
                    moving window. Defaults to (-4, 5).
                    Assuming a 10m² pixel, the default `square_offset` resolves to a one
                    hectare area.

# Returns
Matrix of values 0 - 100 indicating the percentage of the area around the target pixel that
meet suitability criteria.
"""
function proportion_suitable(
    x::Union{BitMatrix,SparseMatrixCSC{Bool,Int64}}; square_offset::Tuple=(-4, 5)
)::SparseMatrixCSC{Int8,Int64}
    subsection_dims = size(x)
    target_area = spzeros(Int8, subsection_dims)

    for row_col in findall(x)
        (row, col) = Tuple(row_col)
        x_left = max(col + square_offset[1], 1)
        x_right = min(col + square_offset[2], subsection_dims[2])

        y_top = max(row + square_offset[1], 1)
        y_bottom = min(row + square_offset[2], subsection_dims[1])

        target_area[row, col] = Int8(sum(@views x[y_top:y_bottom, x_left:x_right]))
    end

    return target_area
end

"""
    assess_region(
        params::RegionalAssessmentParameters
    )::Raster

Perform raster suitability assessment based on user-defined criteria.

# Arguments
- params :: RegionalAssessmentParameters

# Returns
Raster, Int8 booleans (0 and 1) indicating which pixel met criteria.
"""
function assess_region(
    params::RegionalAssessmentParameters
)::Raster
    # Make mask of suitable locations
    @debug "$(now()) : Creating mask for region"

    # Estimate size of mask (in MBs)
    lookup_tbl = params.region_data.slope_table
    region_dims = maximum(lookup_tbl.lon_idx), maximum(lookup_tbl.lat_idx)
    mask_size_MB = prod(region_dims) / 1024^2

    # Builds out a set of criteria filters using the regional criteria.
    # NOTE this will only filter over available criteria
    filters = build_criteria_bounds_from_regional_criteria(params.regional_criteria)
    matching_idx = filter_lookup_table_by_criteria(lookup_tbl, filters)

    # Arbitrary threshold - use sparse matrices if likely to exceed memory
    if mask_size_MB < 700
        @debug "$(now()) : Creating mask as a regular matrix (est. size: $(mask_size_MB))"

        indicator = zeros(Int8, region_dims...)
        Threads.@threads for r in eachrow(lookup_tbl[matching_idx, :])
            indicator[r.lon_idx, r.lat_idx] = true
        end
    else
        @debug "$(now()) : Creating mask as a sparse matrix (est. full size: $(mask_size_MB))"

        valid_lons = lookup_tbl[matching_idx, :lon_idx]
        valid_lats = lookup_tbl[matching_idx, :lat_idx]

        indicator = ExtendableSparseMatrix(
            valid_lons, valid_lats, fill(true, count(matching_idx)), region_dims...
        )
    end

    @debug "$(now()) : Returning regional assessment raster"

    return Raster(params.region_data.valid_extent; data=indicator, missingval=0)
end

"""
    assess_region_quality(
        params::RegionalAssessmentParameters
    )::Raster

Perform raster suitability assessment based on user-defined criteria.

# Arguments
- `params` : RegionalAssessmentParameters

# Returns
Raster indicating suitability of target pixel and immediate surrounds.
"""
function assess_region_quality(
    params::RegionalAssessmentParameters
)::Raster
    # Make mask of suitable locations
    @debug "$(now()) : Creating mask for region"

    # Estimate size of mask (in MBs)
    lookup_tbl = params.region_data.slope_table

    # Builds out a set of criteria filters using the regional criteria.
    # NOTE this will only filter over available criteria
    filters = build_criteria_bounds_from_regional_criteria(params.regional_criteria)

    @debug "$(now()) : Applying criteria thresholds to generate mask layer"
    matching_idx = filter_lookup_table_by_criteria(lookup_tbl, filters)

    @debug "$(now()) : Determining quality of each pixel"
    quality_indicator::Vector{Int8} = assess_location_quality(
        lookup_tbl,
        matching_idx,
        params.region_data.res;
        target_crs=params.region_data.crs
    )

    region_dims = size(params.region_data.valid_extent)
    mask_size_MB = prod(region_dims) / 1024^2  # of Int8

    @debug "$(now()) : Marking results in raster"
    # Arbitrary threshold - use sparse matrices if likely to exceed memory
    if mask_size_MB < 700
        @debug "$(now()) : Creating mask as a regular matrix (est. size: $(mask_size_MB))"
        indicator = zeros(Int8, region_dims...)
        valid_entry = findall(matching_idx)

        Threads.@threads for i in 1:length(quality_indicator)
            r = lookup_tbl[valid_entry[i], :]
            indicator[r.lon_idx, r.lat_idx] = quality_indicator[i]
        end
    else
        @debug "$(now()) : Creating mask as a sparse matrix (est. full size: $(mask_size_MB))"

        valid_lons = lookup_tbl[matching_idx, :lon_idx]
        valid_lats = lookup_tbl[matching_idx, :lat_idx]

        indicator = ExtendableSparseMatrix(
            valid_lons, valid_lats, quality_indicator, region_dims...
        )
    end

    @debug "$(now()) : Returning regional quality raster"

    return Raster(params.region_data.valid_extent; data=indicator, missingval=0)
end

"""
    assess_region_score(
        params::RegionalAssessmentParameters
    )::Raster

Continuous MCDA suitability band for a region, in `[0, 1]` where `1` = a perfect
match across all criteria.

Each valid pixel's raw criteria values are mapped to per-criterion `[0,1]`
suitability via the preference function keyed by each criterion's `direction`, then
combined as a per-criterion weighted mean. This is computed alongside
[`assess_region_quality`](@ref) and does not drive pixel selection: the boolean
suitable-pixel set is still the pure AND-gate of
[`filter_lookup_table_by_criteria`](@ref). The score is scattered over every valid
pixel in the region's lookup table, not just the boolean-suitable subset.

# Arguments
- `params` : RegionalAssessmentParameters

# Returns
Raster of `Float32` MCDA suitability scores, `missingval = 0`.
"""
function assess_region_score(
    params::RegionalAssessmentParameters
)::Raster
    @debug "$(now()) : Creating MCDA suitability score band for region"

    lookup_tbl = params.region_data.slope_table

    # Builds out a set of criteria filters using the regional criteria.
    # NOTE this will only score over available criteria
    filters = build_criteria_bounds_from_regional_criteria(params.regional_criteria)

    @debug "$(now()) : Scoring lookup table against criteria preference functions"
    scores::Vector{Float32} = score_lookup_table_by_criteria(lookup_tbl, filters)

    region_dims = size(params.region_data.valid_extent)
    band_size_MB = prod(region_dims) * sizeof(Float32) / 1024^2

    # A row scoring exactly 0.0 (every weighted criterion at its worst) is a real,
    # legitimate result - not the same thing as a pixel outside the region's valid
    # footprint entirely, which is what `missingval=0`/an unset sparse entry means.
    # Both the dense and sparse paths below encode "not assessed" as an implicit/
    # background 0.0f0, so an assessed-but-terrible pixel would otherwise be
    # indistinguishable from background - rendered fully transparent by the
    # frontend's `band == 0 -> transparent` heatmap style (see `heatmapLayerStyle`),
    # looking like a coverage gap rather than a real (if poor) score. Floor stored
    # scores just above zero so only true background ever reads as exactly 0.0; the
    # raw `scores` vector above (e.g. for future non-raster consumers) is untouched.
    raster_scores = max.(scores, SCORE_RASTER_FLOOR)

    @debug "$(now()) : Marking scores in raster"
    # Arbitrary threshold - use sparse matrices if likely to exceed memory
    if band_size_MB < 700
        @debug "$(now()) : Creating score band as a regular matrix (est. size: $(band_size_MB))"
        indicator = zeros(Float32, region_dims...)

        Threads.@threads for i in eachindex(raster_scores)
            indicator[lookup_tbl.lon_idx[i], lookup_tbl.lat_idx[i]] = raster_scores[i]
        end
    else
        @debug "$(now()) : Creating score band as a sparse matrix (est. full size: $(band_size_MB))"

        indicator = ExtendableSparseMatrix(
            lookup_tbl.lon_idx, lookup_tbl.lat_idx, raster_scores, region_dims...
        )
    end

    @debug "$(now()) : Returning regional MCDA suitability score raster"

    return Raster(params.region_data.valid_extent; data=indicator, missingval=0)
end

function assess_sites(params::SuitabilityAssessmentParameters)
    # Convert suitability integer -> float 64
    suitability_threshold = Float64(params.suitability_threshold / 100.0)
    region = params.region

    @debug "$(now()) : Assessing criteria table for $(region)"
    # Get criteria bounds list from criteria
    filters = build_criteria_bounds_from_regional_criteria(params.regional_criteria)

    slope_table = params.region_data.slope_table
    matching_idx::BitVector = filter_lookup_table_by_criteria(
        slope_table,
        filters
    )

    res = params.region_data.res
    @debug "$(now()) : Assessing $(count(matching_idx)) candidate locations in $(region)."
    @debug "Finding optimal site alignment"
    initial_polygons = find_optimal_site_alignment(
        slope_table[matching_idx, :],
        res,
        params.x_dist,
        params.y_dist;
        target_crs=params.region_data.crs,
        degree_step=20.0,
        suit_threshold=suitability_threshold
    )

    return initial_polygons
end
