using ReefGuide
using Test
using Aqua
using DataFrames
using Random
using Arrow
using Statistics
using Tables
using JSON3
import GeoInterface as GI
import GeometryOps as GO

@testset "Aqua" begin
    Aqua.test_undefined_exports(ReefGuide)
    Aqua.test_stale_deps(ReefGuide)
end

@testset "ReefGuide.jl" begin
    # TODO real tests
    @test true
end

@testset "load_target_region error handling" begin
    # A bogus region_id should surface the real KeyError from the dict lookup,
    # not an UndefVarError from referencing the unbound `region_metadata` in the catch block.
    @test_throws KeyError ReefGuide.load_target_region(;
        region_id="bogus", data_source_directory=tempdir()
    )
end

@testset "load_target_region error handling with scope" begin
    # A bogus region_id must fail before `scope` is ever consulted, regardless
    # of which SpatialScope subtype is passed.
    @test_throws KeyError ReefGuide.load_target_region(;
        region_id="bogus",
        data_source_directory=tempdir(),
        scope=ReefGuide.BBoxScope(144.0, -17.0, 146.0, -15.0)
    )
end

@testset "load_target_region backend validation" begin
    @test_throws ArgumentError ReefGuide.load_target_region(;
        region_id="bogus", data_source_directory=tempdir(), backend=:parquets
    )
end

@testset "load_target_region backend=:arrow vs backend=:parquet, same region" begin
    # Real region data from GBR-reef-guidance-assessment (Phase 1 writer output).
    # Cairns-Cooktown is the smallest of the four regions, chosen for speed. This is
    # a smoke test only — full row-for-row Arrow-vs-Parquet equality across all
    # regions is Phase 3's job, not this one.
    data_dir = normpath(
        joinpath(@__DIR__, "..", "..", "GBR-reef-guidance-assessment", "outputs", "MPA")
    )
    @test isdir(data_dir)

    arrow_entry = ReefGuide.load_target_region(;
        region_id="Cairns-Cooktown", data_source_directory=data_dir, backend=:arrow
    )
    parquet_entry = ReefGuide.load_target_region(;
        region_id="Cairns-Cooktown", data_source_directory=data_dir, backend=:parquet
    )

    @test nrow(arrow_entry.slope_table) > 0
    @test nrow(arrow_entry.slope_table) == nrow(parquet_entry.slope_table)
end

@testset "apply_spatial_scope" begin
    lons = [144.0, 144.5, 145.0, 145.5, 146.0]
    lats = [-16.0, -16.2, -16.4, -16.6, -16.8]
    tbl = DataFrame(; lons=lons, lats=lats, value=1:5)

    @testset "BBoxScope: vectorized range predicate, bounds inclusive" begin
        scope = ReefGuide.BBoxScope(144.4, -16.5, 145.1, -16.1)
        result = ReefGuide.apply_spatial_scope(tbl, scope)
        # Rows 2 (144.5,-16.2) and 3 (145.0,-16.4) fall within the box.
        @test sort(result.value) == [2, 3]

        # Un-scoped table is untouched (new DataFrame returned, not a view/mutation).
        @test nrow(tbl) == 5

        # Bounds are inclusive: a bbox tight on a single point still includes it.
        tight_scope = ReefGuide.BBoxScope(lons[1], lats[1], lons[1], lats[1])
        tight_result = ReefGuide.apply_spatial_scope(tbl, tight_scope)
        @test collect(tight_result.value) == [1]
    end

    @testset "PolygonScope: point-in-polygon via GeometryOps.within" begin
        # Triangle covering roughly rows 1-3 (144.0-145.0, -16.0 to -16.4) but not 4-5.
        ring = GI.Wrappers.LinearRing([
            (143.9, -15.9), (145.2, -15.9), (143.9, -16.5), (143.9, -15.9)
        ])
        poly = GI.Wrappers.Polygon([ring])
        scope = ReefGuide.PolygonScope(poly)
        result = ReefGuide.apply_spatial_scope(tbl, scope)
        @test Set(result.value) ⊆ Set([1, 2, 3])
        @test 1 ∈ result.value
    end

    @testset "no scope on load_target_region call path is unaffected" begin
        # Sanity: apply_spatial_scope is never reached when scope is `nothing`
        # (checked via the `!isnothing(scope)` guard in load_target_region), so
        # there's no dispatch method for `Nothing` — confirm that's still true,
        # since a stray `apply_spatial_scope(df, nothing)` method would signal
        # the guard was bypassed somewhere.
        @test !hasmethod(ReefGuide.apply_spatial_scope, Tuple{DataFrame,Nothing})
    end
end

@testset "find_optimal_site_alignment KD-tree filtering is deterministic" begin
    # Regression test for thread-affinity nondeterminism: the same input must produce
    # bit-identical output regardless of how `Threads.@threads` schedules iterations,
    # both across repeated calls and across differing `JULIA_NUM_THREADS` values.
    rng = Random.MersenneTwister(42)
    n = 60
    lons = collect(range(145.0, 145.05; length=n)) .+ 0.0005 .* randn(rng, n)
    lats = collect(range(-16.0, -16.05; length=n)) .+ 0.0005 .* randn(rng, n)
    lookup_tbl = DataFrame(;
        lons=lons, lats=lats, geometry=GI.Wrappers.Point.(lons, lats)
    )

    r1 = ReefGuide.find_optimal_site_alignment(lookup_tbl, 0.0001, 50.0, 50.0)
    r2 = ReefGuide.find_optimal_site_alignment(lookup_tbl, 0.0001, 50.0, 50.0)

    @test r1.score == r2.score
    @test r1.orientation == r2.orientation
    @test r1.qc_flag == r2.qc_flag
    @test GI.coordinates.(r1.geometry) == GI.coordinates.(r2.geometry)
end

"""Axis-aligned square centered at `(cx, cy)` with the given half-width."""
function _square(cx, cy, half)
    ring = GI.Wrappers.LinearRing([
        (cx - half, cy - half),
        (cx + half, cy - half),
        (cx + half, cy + half),
        (cx - half, cy + half),
        (cx - half, cy - half)
    ])
    return GI.Wrappers.Polygon([ring])
end

@testset "assess_reef_site suitability-threshold scale" begin
    # 9 points, all contained in a 10x10 search box centered on them.
    box = _square(2.0, 2.0, 5.0)
    pts = GI.Wrappers.Point[GI.Wrappers.Point(x, y) for x in 1:3, y in 1:3][:]
    rel_pix = DataFrame(;
        geometry=pts, lons=first.(GI.coordinates.(pts)), lats=last.(GI.coordinates.(pts))
    )
    rotated = GI.Wrappers.Polygon[box]

    # raw count is 9; scaled threshold (0.7 * max_count) must be compared on the
    # same raw-count scale, not against the 0-1 fraction directly.
    score, _, _, qc_flag_below = ReefGuide.assess_reef_site(rel_pix, rotated, 20.0, 0.7)
    @test score == 9.0
    @test qc_flag_below == 1  # 9 < 0.7 * 20 = 14 -> flagged

    _, _, _, qc_flag_above = ReefGuide.assess_reef_site(rel_pix, rotated, 5.0, 0.7)
    @test qc_flag_above == 0  # 9 >= 0.7 * 5 = 3.5 -> not flagged

    # Boundary: raw count exactly at the scaled threshold is not flagged (`<`, not `<=`).
    _, _, _, qc_flag_boundary = ReefGuide.assess_reef_site(rel_pix, rotated, 9.0 / 0.7, 0.7)
    @test qc_flag_boundary == 0
end

@testset "filter_sites overlap deduplication" begin
    @testset "two overlapping plus one independent polygon" begin
        a = _square(0.0, 0.0, 1.0)      # spans -1..1
        b = _square(1.5, 0.0, 1.0)      # spans 0.5..2.5 -> overlaps a
        independent = _square(10.0, 10.0, 1.0)

        res_df = DataFrame(;
            score=[5.0, 9.0, 1.0],
            qc_flag=[0, 0, 0],
            geometry=GI.Wrappers.Polygon[a, b, independent]
        )
        result = ReefGuide.filter_sites(res_df)
        @test sort(result.score) == [1.0, 9.0]
    end

    @testset "no-overlap case passes everything through" begin
        a = _square(0.0, 0.0, 1.0)
        b = _square(10.0, 0.0, 1.0)
        res_df = DataFrame(;
            score=[5.0, 9.0], qc_flag=[0, 0], geometry=GI.Wrappers.Polygon[a, b]
        )
        result = ReefGuide.filter_sites(res_df)
        @test sort(result.score) == [5.0, 9.0]
    end

    @testset "three-element overlap chain, A and C not overlapping" begin
        # A(5) <-> B(9) overlap; B(9) <-> C(20) overlap; A and C do not overlap.
        # The global best scorer (C) must survive, and A - which never overlaps
        # C, the polygon that actually survives - must not be collaterally
        # dropped just because it lost a local comparison to B.
        a = _square(0.0, 0.0, 1.0)   # spans -1..1
        b = _square(1.5, 0.0, 1.0)   # spans 0.5..2.5 -> overlaps a
        c = _square(3.0, 0.0, 1.0)   # spans 2..4 -> overlaps b, not a

        res_df = DataFrame(;
            score=[5.0, 9.0, 20.0], qc_flag=[0, 0, 0], geometry=GI.Wrappers.Polygon[a, b, c]
        )
        result = ReefGuide.filter_sites(res_df)
        @test sort(result.score) == [5.0, 20.0]
    end

    @testset "overlapping bboxes but non-overlapping exact geometry are both kept" begin
        # Two triangles whose bounding boxes overlap but whose actual geometry does not -
        # `STRT.query` only tests bbox intersection, so this only passes if the exact
        # `GO.intersects` check is applied on top of it.
        e = GI.Wrappers.Polygon([
            GI.Wrappers.LinearRing([(0.0, 0.0), (4.0, 0.0), (0.0, 2.0), (0.0, 0.0)])
        ])
        f = GI.Wrappers.Polygon([
            GI.Wrappers.LinearRing([(4.0, 4.0), (4.0, 1.0), (1.0, 4.0), (4.0, 4.0)])
        ])
        @test GI.extent(e).X[2] >= GI.extent(f).X[1]  # bboxes do overlap in X
        @test !GO.intersects(e, f)                    # but the triangles themselves don't

        res_df = DataFrame(;
            score=[1.0, 2.0], qc_flag=[0, 0], geometry=GI.Wrappers.Polygon[e, f]
        )
        result = ReefGuide.filter_sites(res_df)
        @test sort(result.score) == [1.0, 2.0]
    end
end

@testset "CriteriaBounds.rule is missing-safe" begin
    # Lookup columns are Union{Missing,T}: a pixel valid under the region's
    # bathymetry-bound mask can lack data for any individual criterion.
    x = Union{Missing,Float32}[0.5, missing, 2.0, 0.0, 1.0]

    cb = ReefGuide.CriteriaBounds("depth", 0.0f0, 1.0f0)
    r = cb.rule(x)
    @test r isa AbstractVector{Bool}          # never Union{Missing,Bool}
    @test r == Bool[1, 0, 0, 1, 1]            # missing fails the check by default

    # A fully-populated (non-nullable) column is unperturbed by the coalesce.
    @test cb.rule(Float32[0.5, 2.0, 0.0]) == Bool[1, 0, 1]

    # An all-missing column (Arrow may hand back eltype Missing) still yields a
    # plain boolean mask, not Union{Missing,Bool}.
    @test cb.rule(fill(missing, 3)) == falses(3)

    # A positive missing_weight lets a no-data pixel through the bound check.
    cb_pass = ReefGuide.CriteriaBounds("depth", 0.0f0, 1.0f0, 1.0f0)
    @test cb_pass.rule(x) == Bool[1, 1, 0, 1, 1]

    # String constructor parses bounds and shares the same missing policy
    @test ReefGuide.CriteriaBounds("depth", "0.0", "1.0").rule(x) == Bool[1, 0, 0, 1, 1]
end

@testset "CriteriaBounds.rule/.score are backend-fidelity-safe (nullable vs plain eltype)" begin
    # Arrow hands back Union{Missing,T} columns even when zero values are actually
    # missing; DuckDB/Parquet round-trips a column with zero missing values as plain
    # T instead (confirmed in the duckdb-switch plan's Phase 0 spike). `rule`/`score`
    # must agree bit-for-bit regardless of which nullability wrapper the loader used.
    cb = ReefGuide.CriteriaBounds("depth", 0.0f0, 1.0f0)
    nullable_col = Union{Missing,Float32}[0.5, 2.0, 0.0, 1.0]
    plain_col = Float32[0.5, 2.0, 0.0, 1.0]

    @test cb.rule(nullable_col) == cb.rule(plain_col)
    @test cb.rule(plain_col) isa AbstractVector{Bool}
    @test cb.score(nullable_col) == cb.score(plain_col)
end

@testset "CriteriaBounds.rule handles the real Turbidity dtype (Union{Missing,UInt16})" begin
    # Turbidity is written as UInt16 in the real production lookup tables (Phase 0
    # spike ground-truth correction), not Float32 like the other continuous criteria.
    cb = ReefGuide.CriteriaBounds("turbidity", 0.0f0, 1000.0f0)

    x = Union{Missing,UInt16}[100, missing, 1247, 0, 500]
    r = cb.rule(x)
    @test r isa AbstractVector{Bool}          # never Union{Missing,Bool}
    @test r == Bool[1, 0, 0, 1, 1]            # missing fails the check by default

    # Plain (non-nullable) UInt16, as Parquet/DuckDB would round-trip a Turbidity
    # column with zero missing values.
    plain = UInt16[100, 1247, 0, 500]
    @test cb.rule(plain) == Bool[1, 0, 1, 1]
    @test cb.score(plain) isa AbstractVector{Float32}
end

@testset "load_scoped_arrow_table vs load_scoped_parquet_table row-for-row equality (Phase 3, duckdb-switch.md)" begin
    # Real regional data, both formats verified fresh (Phase 1/2). Skipped (not
    # failed) if the data directory isn't present on this machine, so the suite
    # still passes on checkouts without the large GBR-reef-guidance-assessment
    # outputs checked out/downloaded.
    mpa_dir = joinpath(
        @__DIR__, "..", "..", "GBR-reef-guidance-assessment", "outputs", "MPA"
    )
    mpa_dir = abspath(mpa_dir)

    region_ids = [
        "Cairns-Cooktown", "Townsville-Whitsunday", "Mackay-Capricorn", "FarNorthern"
    ]

    if !isdir(mpa_dir)
        @warn "Skipping Phase 3 correctness verification: data directory not found" mpa_dir
        @test_skip false
    else
        # Same viewport construction convention as scratch/bbox_spike.jl and
        # scratch/duckdb_spike.jl: median lon/lat center, small (~5-10km) and
        # large (~50km) bbox viewports.
        km_per_deg_lat = 111.32

        function make_bbox(lon_mid, lat_mid, viewport_km)
            km_per_deg_lon = 111.32 * cosd(lat_mid)
            half_span_lat = (viewport_km / 2) / km_per_deg_lat
            half_span_lon = (viewport_km / 2) / km_per_deg_lon
            return ReefGuide.BBoxScope(
                lon_mid - half_span_lon,
                lat_mid - half_span_lat,
                lon_mid + half_span_lon,
                lat_mid + half_span_lat
            )
        end

        function assert_rows_equal(arrow_df, parquet_df; region_id, viewport_km)
            @testset "$(region_id) ($(viewport_km) km viewport)" begin
                # A viewport that happens to select zero rows would make every
                # check below vacuously pass - guard against that explicitly.
                @test nrow(arrow_df) > 0
                @test nrow(arrow_df) == nrow(parquet_df)
                @test Set(names(arrow_df)) == Set(names(parquet_df))

                arrow_keys = sort(collect(zip(arrow_df.lon_idx, arrow_df.lat_idx)))
                parquet_keys = sort(collect(zip(parquet_df.lon_idx, parquet_df.lat_idx)))
                @test arrow_keys == parquet_keys

                # Compare per-column values keyed by (lon_idx, lat_idx) rather than
                # row position, since row order isn't guaranteed identical between
                # the two loaders (Phase 0 findings).
                arrow_sorted = sort(arrow_df, [:lon_idx, :lat_idx])
                parquet_sorted = sort(parquet_df, [:lon_idx, :lat_idx])

                for col in names(arrow_sorted)
                    @test isequal(arrow_sorted[!, col], parquet_sorted[!, col])
                end
            end
        end

        for region_id in region_ids
            arrow_path = joinpath(mpa_dir, "$(region_id)_valid_slopes_lookup.arrow")
            parquet_path = joinpath(mpa_dir, "$(region_id)_valid_slopes_lookup.parquet")

            @test isfile(arrow_path)
            @test isfile(parquet_path)

            # Center the viewport on an actual data point closest to the median
            # lon/lat, rather than the independently-computed per-coordinate
            # median itself - for a non-convex region (e.g. scattered reef
            # clusters), that synthetic point can fall in a gap with no data at
            # all, making the smallest viewport spuriously empty.
            full_arrow_tbl = Arrow.Table(arrow_path)
            lons_full = full_arrow_tbl.lons
            lats_full = full_arrow_tbl.lats
            med_lon = median(lons_full)
            med_lat = median(lats_full)
            nearest_idx = argmin(
                @. (lons_full - med_lon)^2 + (lats_full - med_lat)^2
            )
            lon_mid = lons_full[nearest_idx]
            lat_mid = lats_full[nearest_idx]
            full_arrow_tbl = nothing
            GC.gc()

            for viewport_km in (8.0, 50.0)
                scope = make_bbox(lon_mid, lat_mid, viewport_km)

                arrow_df = ReefGuide.load_scoped_arrow_table(arrow_path, scope)
                parquet_df = ReefGuide.load_scoped_parquet_table(parquet_path, scope)

                assert_rows_equal(arrow_df, parquet_df; region_id, viewport_km)

                arrow_df = nothing
                parquet_df = nothing
                GC.gc()
            end
        end
    end
end

@testset "filter_lookup_table_by_criteria with partial-coverage rows" begin
    lookup = DataFrame(;
        lon_idx=Int32[1, 2, 3, 4],
        lat_idx=Int32[1, 1, 1, 1],
        depth=Union{Missing,Float32}[0.5, missing, 0.7, 0.9],
        slope=Union{Missing,Float32}[10.0, 12.0, missing, 11.0]
    )
    ruleset = ReefGuide.CriteriaBounds[
        ReefGuide.CriteriaBounds("depth", 0.0f0, 1.0f0),
        ReefGuide.CriteriaBounds("slope", 5.0f0, 15.0f0)
    ]

    matches = ReefGuide.filter_lookup_table_by_criteria(lookup, ruleset)
    @test matches isa BitVector                       # declared return type is preserved
    @test matches == BitVector([1, 0, 0, 1])          # rows 2 and 3 excluded via missing

    # Row-indexing patterns the assess_* callers rely on must still work.
    @test nrow(lookup[matches, :]) == 2
    @test count(matches) == 2
end

@testset "derive_criteria_bounds_from_slope_table skips missing values" begin
    region = ReefGuide.RegionMetadata(;
        display_name="Test", id="test-region", available_criteria=["Depth", "Slope"]
    )
    table = DataFrame(;
        Depth=Union{Missing,Float32}[-5.0, missing, -8.0, -3.0],
        Slope=Union{Missing,Float32}[missing, missing, missing, missing]
    )

    bounds = @test_logs (:warn, r"entirely missing") match_mode = :any begin
        ReefGuide.derive_criteria_bounds_from_slope_table(table, region)
    end

    # Depth bounds are the extrema of its non-missing values only.
    @test bounds["Depth"].bounds.min == -8.0f0
    @test bounds["Depth"].bounds.max == -3.0f0
    # A criterion whose column is entirely missing is left out of the result.
    @test !haskey(bounds, "Slope")
end

@testset "an omitted criterion is simply not filtered" begin
    # A BoundedCriteriaDict with no entry for a criterion (as produced for an
    # all-missing column) builds no CriteriaBounds for it, so that column is never
    # filtered - the intended effect of the omission.
    dict = ReefGuide.BoundedCriteriaDict(
        "Depth" => ReefGuide.BoundedCriteria(;
            metadata=ReefGuide.ASSESSMENT_CRITERIA["Depth"],
            bounds=ReefGuide.Bounds(; min=-10.0, max=-2.0)
        )
    )

    filters = ReefGuide.build_criteria_bounds_from_regional_criteria(dict)

    @test [f.name for f in filters] == [:Depth]
end

# MCDA weighted-sum scorer

@testset "MCDA per-direction ramp endpoints + clamping" begin
    lb, ub = 0.0f0, 10.0f0

    cb_low = ReefGuide.CriteriaBounds("x", lb, ub; direction=:lower_is_better)
    s = cb_low.score(Float32[0.0, 10.0, -5.0, 15.0, 5.0])
    @test s[1] === 1.0f0                  # exactly 1 at the lower bound
    @test s[2] === 0.0f0                  # exactly 0 at the upper bound
    @test s[3] === 1.0f0                  # clamped below the ramp
    @test s[4] === 0.0f0                  # clamped above the ramp
    @test isapprox(s[5], 0.5f0)           # linear midpoint

    cb_high = ReefGuide.CriteriaBounds("x", lb, ub; direction=:higher_is_better)
    s = cb_high.score(Float32[0.0, 10.0, -5.0, 15.0, 5.0])
    @test s[1] === 0.0f0                  # exactly 0 at the lower bound
    @test s[2] === 1.0f0                  # exactly 1 at the upper bound
    @test s[3] === 0.0f0                  # clamped below the ramp
    @test s[4] === 1.0f0                  # clamped above the ramp
    @test isapprox(s[5], 0.5f0)
end

@testset "MCDA band triangle" begin
    cb = ReefGuide.CriteriaBounds(
        "Depth", -12.5f0, -2.0f0; direction=:band, band_peak=-5.0f0
    )
    s = cb.score(Float32[-5.0, -12.5, -2.0, -7.0, -4.0, -20.0, 0.0])
    @test s[1] === 1.0f0                       # peak
    @test s[2] === 0.0f0                       # lower bound
    @test s[3] === 0.0f0                       # upper bound
    @test isapprox(s[4], 5.5f0 / 7.5f0)        # linear on the left flank
    @test isapprox(s[5], 2.0f0 / 3.0f0)        # linear on the right flank
    @test s[6] === 0.0f0                       # below the band
    @test s[7] === 0.0f0                       # above the band
end

@testset "MCDA equal-weight mean aggregation across >= 2 criteria" begin
    df = DataFrame(; A=Float32[0.0, 10.0], B=Float32[10.0, 10.0])
    ruleset = ReefGuide.CriteriaBounds[
        ReefGuide.CriteriaBounds("A", 0.0f0, 10.0f0; direction=:higher_is_better),
        ReefGuide.CriteriaBounds("B", 0.0f0, 10.0f0; direction=:higher_is_better)
    ]

    scores = ReefGuide.score_lookup_table_by_criteria(df, ruleset)
    @test scores isa Vector{Float32}
    @test length(scores) == nrow(df)
    @test isapprox(scores[1], 0.5f0)          # (0.0 + 1.0) / 2
    @test isapprox(scores[2], 1.0f0)          # (1.0 + 1.0) / 2
end

@testset "MCDA missing_weight folds a missing cell onto the [0,1] scale" begin
    df = DataFrame(;
        lon_idx=Int32[1, 2],
        lat_idx=Int32[1, 1],
        A=Union{Missing,Float32}[missing, 5.0],
        B=Union{Missing,Float32}[6.0, 6.0]
    )

    # missing_weight = 0.0f0: the missing cell fails the boolean rule and its
    # per-criterion score is 0.
    rs_zero = ReefGuide.CriteriaBounds[
        ReefGuide.CriteriaBounds("A", 0.0f0, 10.0f0; direction=:higher_is_better),
        ReefGuide.CriteriaBounds("B", 0.0f0, 10.0f0; direction=:higher_is_better)
    ]
    matches = ReefGuide.filter_lookup_table_by_criteria(df, rs_zero)
    @test matches isa BitVector
    @test matches == BitVector([0, 1])                     # row 1 excluded via missing
    @test rs_zero[1].score([missing])[1] === 0.0f0         # per-criterion score is 0

    scores_zero = ReefGuide.score_lookup_table_by_criteria(df, rs_zero)
    @test isapprox(scores_zero[1], 0.3f0)                  # (0.0 + 0.6) / 2

    # A mid missing_weight pulls the aggregate for that row up toward 0.5.
    rs_mid = ReefGuide.CriteriaBounds[
        ReefGuide.CriteriaBounds(
            "A", 0.0f0, 10.0f0, 0.5f0; direction=:higher_is_better
        ),
        ReefGuide.CriteriaBounds("B", 0.0f0, 10.0f0; direction=:higher_is_better)
    ]
    scores_mid = ReefGuide.score_lookup_table_by_criteria(df, rs_mid)
    @test isapprox(scores_mid[1], 0.55f0)                  # (0.5 + 0.6) / 2
    @test scores_mid[1] > scores_zero[1]
end

@testset "MCDA boolean path is unchanged (regression) + score band well-formed" begin
    lookup = DataFrame(;
        lon_idx=Int32[1, 2, 3, 4, 5],
        lat_idx=Int32[1, 1, 1, 1, 1],
        Depth=Union{Missing,Float32}[-5.0, -13.0, -3.0, missing, -8.0],
        Slope=Union{Missing,Float32}[10.0, 5.0, 45.0, 20.0, missing]
    )
    ruleset = ReefGuide.CriteriaBounds[
        ReefGuide.CriteriaBounds(
            "Depth", -12.5f0, -2.0f0; direction=:band, band_peak=-5.0f0
        ),
        ReefGuide.CriteriaBounds("Slope", 0.0f0, 40.0f0; direction=:lower_is_better)
    ]

    matches = ReefGuide.filter_lookup_table_by_criteria(lookup, ruleset)
    @test matches isa BitVector
    # Row 1 passes both criteria; rows 2-3 fail on Depth/Slope out of band,
    # rows 4-5 fail on a missing value under the default missing_weight = 0.
    @test matches == BitVector([1, 0, 0, 0, 0])

    scores = ReefGuide.score_lookup_table_by_criteria(lookup, ruleset)
    @test scores isa Vector{Float32}
    @test length(scores) == nrow(lookup)
    @test all(v -> v >= 0.0f0 && v <= 1.0f0, scores)
end

@testset "MCDA config threads Criteria -> BoundedCriteria -> CriteriaBounds" begin
    # Defaults from ASSESSMENT_CRITERIA reach CriteriaBounds via BoundedCriteria.
    dict = ReefGuide.BoundedCriteriaDict(
        "Depth" => ReefGuide.BoundedCriteria(;
            metadata=ReefGuide.ASSESSMENT_CRITERIA["Depth"],
            bounds=ReefGuide.Bounds(; min=-12.5, max=-2.0)
        ),
        "Slope" => ReefGuide.BoundedCriteria(;
            metadata=ReefGuide.ASSESSMENT_CRITERIA["Slope"],
            bounds=ReefGuide.Bounds(; min=0.0, max=40.0)
        )
    )
    filters = ReefGuide.build_criteria_bounds_from_regional_criteria(dict)
    by_name = Dict(f.name => f for f in filters)

    @test by_name[:Depth].direction === :band
    @test by_name[:Depth].band_peak == -5.0f0
    @test by_name[:Depth].missing_weight === 0.0f0
    @test by_name[:Depth].weight === 1.0f0
    @test by_name[:Slope].direction === :lower_is_better
    @test by_name[:Slope].band_peak === nothing

    # A per-request override on BoundedCriteria propagates too.
    dict2 = ReefGuide.BoundedCriteriaDict(
        "Slope" => ReefGuide.BoundedCriteria(;
            metadata=ReefGuide.ASSESSMENT_CRITERIA["Slope"],
            bounds=ReefGuide.Bounds(; min=0.0, max=40.0),
            direction=:band,
            band_peak=15.0,
            missing_weight=0.25f0,
            weight=2.0f0
        )
    )
    f = only(ReefGuide.build_criteria_bounds_from_regional_criteria(dict2))
    @test f.direction === :band
    @test f.band_peak === 15.0f0
    @test f.missing_weight === 0.25f0
    @test f.weight === 2.0f0

    # Overriding direction to :band without an explicit band_peak leaves it
    # `nothing` (metadata default), so the score falls back to the bound midpoint.
    dict3 = ReefGuide.BoundedCriteriaDict(
        "Slope" => ReefGuide.BoundedCriteria(;
            metadata=ReefGuide.ASSESSMENT_CRITERIA["Slope"],
            bounds=ReefGuide.Bounds(; min=0.0, max=40.0),
            direction=:band
        )
    )
    g = only(ReefGuide.build_criteria_bounds_from_regional_criteria(dict3))
    @test g.band_peak === nothing
end

# duckdb-switch plan, Phase 3 — Correctness verification (hard gate, see
# .claude/plans/duckdb-switch.md). Real regional data from
# GBR-reef-guidance-assessment/outputs/MPA (Phase 1 writer output), all 4 regions.
@testset "Phase 3 (duckdb-switch): Arrow vs Parquet correctness verification" begin
    data_dir = normpath(
        joinpath(@__DIR__, "..", "..", "GBR-reef-guidance-assessment", "outputs", "MPA")
    )
    @test isdir(data_dir)

    regions = [
        "Cairns-Cooktown", "Townsville-Whitsunday", "Mackay-Capricorn", "FarNorthern"
    ]

    # Populated by Part 1, reused by Part 2 so each region's real data is only
    # loaded once via each backend (these tables are already-scoped, not full
    # regions, so this stays well under the ~2-3 GB/region full-load floor).
    region_scoped = Dict{String,NamedTuple}()

    @testset "Part 1: PolygonScope row-for-row equality — $(region)" for region in regions
        arrow_path = joinpath(data_dir, "$(region)_valid_slopes_lookup.arrow")
        parquet_path = joinpath(data_dir, "$(region)_valid_slopes_lookup.parquet")
        @test isfile(arrow_path)
        @test isfile(parquet_path)

        # Only the lons/lats columns are decompressed to find the region's
        # median centre, mirroring `_build_scope_mask`'s partial-column access
        # rather than a full-table load.
        tbl = Arrow.Table(arrow_path)
        lons = Tables.getcolumn(tbl, :lons)
        lats = Tables.getcolumn(tbl, :lats)
        med_lon = median(lons)
        med_lat = median(lats)
        # Anchor on an actual data point nearest the median, not the median
        # itself - for a non-convex/scattered region the synthetic median can
        # fall in a real gap, making the polygon spuriously empty (see the
        # identical fix applied to the BBoxScope block above).
        nearest_idx = argmin(@. (lons - med_lon)^2 + (lats - med_lat)^2)
        lon_mid = lons[nearest_idx]
        lat_mid = lats[nearest_idx]
        tbl = nothing
        GC.gc()

        # ~50km-scale right triangle - a "meaningful chunk", same viewport scale
        # as scratch/duckdb_spike.jl's 50km bbox case and the same triangle shape
        # as the existing PolygonScope fixture. Positioned so the anchor point
        # sits at a fixed interior fraction (0.3, 0.3) of the two legs - strictly
        # inside, away from every edge/vertex - rather than at the centroid,
        # which can sit exactly on the hypotenuse for a right isoceles triangle.
        km_per_deg_lat = 111.32
        km_per_deg_lon = 111.32 * cosd(lat_mid)
        leg_lat = 50.0 / km_per_deg_lat
        leg_lon = 50.0 / km_per_deg_lon
        corner_lon = lon_mid - 0.3 * leg_lon
        corner_lat = lat_mid - 0.3 * leg_lat
        ring = GI.Wrappers.LinearRing([
            (corner_lon, corner_lat),
            (corner_lon + leg_lon, corner_lat),
            (corner_lon, corner_lat + leg_lat),
            (corner_lon, corner_lat)
        ])
        poly = GI.Wrappers.Polygon([ring])
        scope = ReefGuide.PolygonScope(poly)

        arrow_scoped = ReefGuide.load_scoped_arrow_table(arrow_path, scope)
        parquet_scoped = ReefGuide.load_scoped_parquet_table(parquet_path, scope)

        @test nrow(arrow_scoped) > 0
        @test nrow(arrow_scoped) == nrow(parquet_scoped)
        @test Set(names(arrow_scoped)) == Set(names(parquet_scoped))

        # Same lon_idx/lat_idx set, order-independent.
        arrow_keys = sort(collect(zip(arrow_scoped.lon_idx, arrow_scoped.lat_idx)))
        parquet_keys = sort(collect(zip(parquet_scoped.lon_idx, parquet_scoped.lat_idx)))
        @test arrow_keys == parquet_keys

        # Align both frames on the same key order, then every column must match
        # exactly, including `missing` in the same positions (`isequal`, not `==`).
        sort!(arrow_scoped, [:lon_idx, :lat_idx])
        sort!(parquet_scoped, [:lon_idx, :lat_idx])
        for col in names(arrow_scoped)
            @test isequal(arrow_scoped[!, col], parquet_scoped[!, col])
        end

        region_scoped[region] = (arrow=arrow_scoped, parquet=parquet_scoped)
    end

    """Narrow a real [min, max] bound to its middle half, so both a pass and a fail case exist."""
    function _narrowed_bounds(lo::Float64, hi::Float64)
        span = hi - lo
        return (Float32(lo + 0.25 * span), Float32(hi - 0.25 * span))
    end

    @testset "Part 2: CriteriaBounds/filter_lookup_table_by_criteria identical across backends — $(region)" for region in
                                                                                                                  regions
        arrow_scoped, parquet_scoped = region_scoped[region]

        bounds_path = joinpath(data_dir, "$(region)_valid_slopes_bounds.json")
        @test isfile(bounds_path)
        bounds_json = JSON3.read(Base.read(bounds_path, String), Dict{String,Any})

        # Depth/Slope mirror the existing synthetic CriteriaBounds/
        # filter_lookup_table_by_criteria tests' logic; Turbidity is included
        # explicitly since it is the one value column with a non-Float32 dtype
        # (Union{Missing,UInt16}) per the duckdb-switch plan's ground truth.
        ruleset = ReefGuide.CriteriaBounds[]
        for crit in ("Depth", "Slope", "Turbidity")
            @test haskey(bounds_json, crit)
            lo, hi = _narrowed_bounds(
                Float64(bounds_json[crit]["min"]), Float64(bounds_json[crit]["max"])
            )
            push!(ruleset, ReefGuide.CriteriaBounds(crit, lo, hi))
        end

        # Per-criterion `rule` must agree bit-for-bit regardless of which
        # backend's nullability wrapper the column came back with.
        for cb in ruleset
            r_arrow = cb.rule(arrow_scoped[!, cb.name])
            r_parquet = cb.rule(parquet_scoped[!, cb.name])
            @test r_arrow isa AbstractVector{Bool}
            @test r_parquet isa AbstractVector{Bool}
            @test r_arrow == r_parquet
        end

        matches_arrow = ReefGuide.filter_lookup_table_by_criteria(arrow_scoped, ruleset)
        matches_parquet = ReefGuide.filter_lookup_table_by_criteria(
            parquet_scoped, ruleset
        )
        @test matches_arrow isa BitVector
        @test matches_parquet isa BitVector
        @test matches_arrow == matches_parquet

        scores_arrow = ReefGuide.score_lookup_table_by_criteria(arrow_scoped, ruleset)
        scores_parquet = ReefGuide.score_lookup_table_by_criteria(parquet_scoped, ruleset)
        @test scores_arrow isa Vector{Float32}
        @test scores_parquet isa Vector{Float32}
        @test scores_arrow == scores_parquet
    end
end
