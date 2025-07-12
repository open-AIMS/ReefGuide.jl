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
