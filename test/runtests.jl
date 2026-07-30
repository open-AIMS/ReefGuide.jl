using ReefGuide
using Test
using Aqua

@testset "Aqua" begin
    Aqua.test_undefined_exports(ReefGuide)
    Aqua.test_stale_deps(ReefGuide; ignore=[:GeoJSON])
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
