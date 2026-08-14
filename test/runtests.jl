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
