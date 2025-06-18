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

end
