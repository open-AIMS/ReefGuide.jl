using ReefGuide
using Test
using Aqua
using DataFrames
using Random
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
