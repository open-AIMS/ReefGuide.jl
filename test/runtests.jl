using ReefGuide
using Test
using Aqua
using DataFrames
using Random
import GeoInterface as GI

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

