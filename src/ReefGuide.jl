module ReefGuide

# System imports
using Base.Threads

# File IO
using Glob, Serialization

# Geospatial
using ArchGDAL, GeoParquet, Rasters

# Collections
using DataFrames, OrderedCollections, SparseArrays

# Multithreading
using FLoops, ThreadsX

# Precompilation
using PrecompileSignatures: @precompile_signatures
using PrecompileTools

# Utilities and helpers for assessments
include("utility/utility.jl")

# Assessment logic
include("assessment_methods/assessment_methods.jl")

export initialise_data,
    assess_sites,
    filter_sites,
    output_geojson,
    assess_region,
    _write_cog,
    _write_tiff

# Auto-generate precompilation signatures for ReefGuide
@precompile_signatures(ReefGuide)

# Force precompilation of methods that slow down initial use.
@setup_workload begin
    @compile_workload begin
        # Enforce precompilation of specific geospatial read/write methods
        GeoParquet.read(joinpath(pkgdir(ReefGuide), "assets", "dummy.parq"))
        GDF.read(joinpath(pkgdir(ReefGuide), "assets", "dummy.gpkg"))

        tmpfile = tempname()
        dims = (X(1.0:1000.0), Y(1:1000.0))
        x = Raster(rand(1000, 1000); dims=dims)
        _write_cog(tmpfile, x; tile_size=(256,), num_threads=1)

        x = Raster(sparse(rand(1000, 1000)); dims=dims)
        _write_cog(tmpfile, x; tile_size=(256,), num_threads=1)

        x2 = ExtendableSparseMatrix{Int8,Int64}(1000, 1000)
        x = Raster(x2; dims=dims)
        _write_cog(tmpfile, x; tile_size=(256,), num_threads=1)
    end
end

end
