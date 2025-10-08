"""Methods for common file I/O."""

"""
    _write_cog(file_path::String, data::Raster, config::Dict)::Nothing

Write out a COG using common options.

# Arguments
- `file_path` : Path to write data out to
- `data` : Raster data to write out
- `tile_size` : Size of tiles to use in the COG (x,y dimension)
- `num_threads` : Number of threads to use for writing
"""
function _write_cog(
    file_path::String, data::Raster; tile_size::Tuple{Integer}, num_threads::Integer
)::Nothing
    Rasters.write(
        file_path,
        data;
        ext=".tiff",
        source="gdal",
        driver="COG",
        options=Dict{String,String}(
            "COMPRESS" => "DEFLATE",
            "LEVEL" => "1",
            "PREDICTOR" => "YES",
            "SPARSE_OK" => "TRUE",
            "OVERVIEW_COUNT" => "5",
            "BIGTIFF" => "IF_SAFER",
            "BLOCKSIZE" => string(first(tile_size)),
            "NUM_THREADS" => string(num_threads)
        ),
        force=true
    )

    return nothing
end

"""
    _write_tiff(file_path::String, data::Raster)::Nothing

Write out a geotiff using common options.

# Arguments
- `file_path` : Path to write data out to
- `data` : Raster data to write out
"""
function _write_tiff(file_path::String, data::Raster)::Nothing
    Rasters.write(
        file_path,
        data;
        ext=".tiff",
        source="gdal",
        driver="gtiff",
        force=true
    )

    return nothing
end
