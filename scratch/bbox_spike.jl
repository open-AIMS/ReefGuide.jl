#!/usr/bin/env julia
# ============================================================================
# THROWAWAY SPIKE SCRIPT — not part of the package source.
#
# Phase 1 de-risking spike for the "fast preview" assessment feature
# (see .claude/plans/2026-07-31_slow_fast_assessment_feature.md).
#
# Question: is "GeoParquet.read the FULL region parquet, then filter in
# memory on :lons/:lats" fast enough to justify a separate "fast" job type,
# or do we need real row-group-level spatial pushdown in Parquet2.jl /
# GeoParquet.jl?
#
# Usage (from the ReefGuide.jl-wt-fast-assessment worktree root):
#   julia --project=. scratch/bbox_spike.jl [path/to/region_valid_slopes_lookup.parq]
#
# Default target file (if no arg given) is the real Cairns-Cooktown region
# lookup parquet from the sibling ReefGuideWorker.jl/data/ directory (real
# production-shaped data, not synthetic) — the smallest of the four regional
# files at ~358 MB.
# ============================================================================

using DataFrames
using GeoParquet
using Statistics
using Printf

const DEFAULT_PARQ = joinpath(
    @__DIR__, "..", "..", "ReefGuideWorker.jl", "data",
    "Cairns-Cooktown_valid_slopes_lookup.parq"
)

parq_path = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_PARQ
parq_path = abspath(parq_path)

println("="^78)
println("BBOX FILTER SPIKE — full load + in-memory filter benchmark")
println("="^78)
println("Target file: $(parq_path)")
println("File size:   $(round(filesize(parq_path) / 1024^2, digits=1)) MB")
println()

# ----------------------------------------------------------------------
# Step 1: Load the FULL region parquet, exactly as load_target_region does
# (ReefGuide.jl/src/utility/regions_criteria_setup.jl:581-583).
# ----------------------------------------------------------------------
println("--- Step 1: GeoParquet.read (full region load) ---")
# Run twice: first to trigger any JIT/precompilation cost, second for the
# "warm" timing that's representative of steady-state worker behaviour.
local slope_table
compile_time = @elapsed begin
    slope_table = GeoParquet.read(parq_path)
end
@printf "First call (incl. compilation): %.3f s\n" compile_time

load_time = @elapsed begin
    global slope_table = GeoParquet.read(parq_path)
end
n_total = nrow(slope_table)
@printf "Warm load time:                 %.3f s\n" load_time
println("Row count (full table):         $(n_total)")
println("Columns: $(names(slope_table))")
println()

# ----------------------------------------------------------------------
# Step 2: Determine :lons/:lats columns (add if not already present, same
# as load_target_region does via add_lat_long_columns_to_dataframe).
# ----------------------------------------------------------------------
has_lonlat = "lons" in names(slope_table) && "lats" in names(slope_table)
println("Has :lons/:lats columns directly: $(has_lonlat)")
if !has_lonlat
    error(
        "This spike assumes :lons/:lats columns are already present (as " *
        "load_target_region expects after add_lat_long_columns_to_dataframe). " *
        "Columns found: $(names(slope_table))"
    )
end

lons = slope_table.lons
lats = slope_table.lats

lon_min, lon_max = extrema(lons)
lat_min, lat_max = extrema(lats)
lon_mid = median(lons)
lat_mid = median(lats)

println()
println("--- Full-table spatial extent ---")
@printf "lon range: [%.5f, %.5f]  (span %.4f deg)\n" lon_min lon_max (lon_max - lon_min)
@printf "lat range: [%.5f, %.5f]  (span %.4f deg)\n" lat_min lat_max (lat_max - lat_min)
@printf "median center: (lon=%.5f, lat=%.5f)\n" lon_mid lat_mid
println()

# ----------------------------------------------------------------------
# Step 3: Construct a realistic "viewport" bbox — a few km across — centred
# on the data's median point (so it's guaranteed to land in a populated
# area rather than an empty corner of the region's bounding envelope).
#
# Degrees-per-km at this latitude:
#   1 deg latitude  ≈ 111.32 km (approx constant everywhere)
#   1 deg longitude ≈ 111.32 * cos(latitude) km
# We use the median latitude of the actual data as our reference latitude.
# ----------------------------------------------------------------------
viewport_km = 5.0  # "a few km across" viewport, e.g. a single reef + margin
km_per_deg_lat = 111.32
km_per_deg_lon = 111.32 * cosd(lat_mid)

half_span_lat = (viewport_km / 2) / km_per_deg_lat
half_span_lon = (viewport_km / 2) / km_per_deg_lon

bbox_lon_min = lon_mid - half_span_lon
bbox_lon_max = lon_mid + half_span_lon
bbox_lat_min = lat_mid - half_span_lat
bbox_lat_max = lat_mid + half_span_lat

println("--- Step 2: constructing a $(viewport_km) km x $(viewport_km) km viewport bbox ---")
@printf "km_per_deg_lon at lat=%.3f: %.3f\n" lat_mid km_per_deg_lon
@printf "bbox lon: [%.6f, %.6f]  (span %.6f deg)\n" bbox_lon_min bbox_lon_max (bbox_lon_max - bbox_lon_min)
@printf "bbox lat: [%.6f, %.6f]  (span %.6f deg)\n" bbox_lat_min bbox_lat_max (bbox_lat_max - bbox_lat_min)
println()

# ----------------------------------------------------------------------
# Step 4: Filter in-memory on :lons/:lats with a simple range predicate,
# exactly the kind of predicate the plan proposes adding to
# load_target_region for the bbox scope case.
# ----------------------------------------------------------------------
println("--- Step 3: in-memory bbox filter (simple range predicate) ---")
filter_time = @elapsed begin
    global viewport_table = filter(
        row -> bbox_lon_min <= row.lons <= bbox_lon_max &&
               bbox_lat_min <= row.lats <= bbox_lat_max,
        slope_table
    )
end
n_filtered = nrow(viewport_table)
@printf "Filter time:                    %.4f s\n" filter_time
println("Row count (viewport, 5km bbox): $(n_filtered)")
@printf "Reduction: %.4f%% of rows retained (%dx smaller)\n" (100 * n_filtered / n_total) (n_total / max(n_filtered, 1))
println()

# Also try a much larger bbox (e.g. 50km, ~ typical zoomed-out viewport) for
# comparison, to see how filter time scales with result-set size (it's the
# predicate scan over ALL rows that dominates cost, not the result size, so
# we'd expect filter time to be roughly constant regardless of bbox size).
viewport_km_large = 50.0
half_span_lat_l = (viewport_km_large / 2) / km_per_deg_lat
half_span_lon_l = (viewport_km_large / 2) / km_per_deg_lon
bbox2 = (
    lon_mid - half_span_lon_l, lon_mid + half_span_lon_l,
    lat_mid - half_span_lat_l, lat_mid + half_span_lat_l
)
filter_time_large = @elapsed begin
    global viewport_table_large = filter(
        row -> bbox2[1] <= row.lons <= bbox2[2] && bbox2[3] <= row.lats <= bbox2[4],
        slope_table
    )
end
n_filtered_large = nrow(viewport_table_large)
println("--- Comparison: $(viewport_km_large) km bbox ---")
@printf "Filter time:                    %.4f s\n" filter_time_large
println("Row count (50km bbox):          $(n_filtered_large)")
println()

# ----------------------------------------------------------------------
# Step 4b: the row-wise `filter(row -> ..., df)` predicate above pays
# per-row DataFrames.jl accessor overhead (DataFrameRow indexing) across
# all 15M rows. A vectorized boolean-mask filter on the raw columns is the
# idiomatic-Julia way to write "a simple range predicate over :lons/:lats"
# and should be dramatically faster — measure it for a fair comparison.
# ----------------------------------------------------------------------
println("--- Step 3b: vectorized boolean-mask filter (idiomatic form) ---")
filter_time_vec = @elapsed begin
    global mask = (bbox_lon_min .<= lons .<= bbox_lon_max) .&
                  (bbox_lat_min .<= lats .<= bbox_lat_max)
    global viewport_table_vec = slope_table[mask, :]
end
n_filtered_vec = nrow(viewport_table_vec)
@printf "Vectorized filter time:         %.4f s\n" filter_time_vec
println("Row count (viewport, 5km bbox): $(n_filtered_vec)")
println()

# ----------------------------------------------------------------------
# Step 5: Summary — total wall time for "load + filter" end to end.
# ----------------------------------------------------------------------
total_time = load_time + filter_time
total_time_vec = load_time + filter_time_vec
println("="^78)
println("SUMMARY")
println("="^78)
@printf "Warm load time (full region):             %.3f s\n" load_time
@printf "Filter time (5km viewport, row-wise):     %.4f s\n" filter_time
@printf "Filter time (5km viewport, vectorized):   %.4f s\n" filter_time_vec
@printf "TOTAL load+filter (row-wise):              %.3f s\n" total_time
@printf "TOTAL load+filter (vectorized):            %.3f s\n" total_time_vec
println("Row count: $(n_total) -> $(n_filtered) ($(round(100*n_filtered/n_total, digits=4))% retained)")
println()

# ----------------------------------------------------------------------
# Step 6 (optional): note on Parquet2.jl / GeoParquet.jl row-group /
# predicate-pushdown API availability. This is investigation only — no
# pushdown is implemented here.
#
# Findings (from reading installed package source under ~/.julia/packages):
#
# - GeoParquet.read(fn; kwargs...) (GeoParquet/src/io.jl) is a thin wrapper:
#     ds = Parquet2.Dataset(fn, kwargs...)
#     df = DataFrame(ds; copycols=false)
#   It always materializes ALL row groups into the DataFrame — kwargs are
#   passed through to Parquet2.Dataset but there is no bbox/predicate kwarg
#   that skips row groups based on spatial extent. GeoParquet.jl writes a
#   whole-file `bbox` in its "geo" metadata (io.jl:32, meta.jl:43) but never
#   reads/uses it to prune anything on the read path — it's metadata-only.
#
# - Parquet2.jl DOES carry per-column-chunk statistics (min/max/n_distinct/
#   n_null) via `use_statistics=true` on `Parquet2.Dataset(...)`, exposed as
#   `VectorWithStatistics` (Parquet2/src/read.jl:57-98, ColumnStatistics).
#   These stats are recorded in a way that follows Parquet's standard
#   "one Statistics struct per column *chunk within a row group*" model
#   (Metadata.thrift), so in principle you COULD open a Dataset, inspect
#   per-row-group min/max on the lons/lats columns, and skip whole row
#   groups that don't intersect the bbox before materializing — but:
#     (a) this requires bypassing GeoParquet.jl's read() entirely and
#         working with Parquet2.Dataset row-group internals directly
#         (no public "select these row groups" API was found — row-group
#         iteration/filtering isn't exposed as a documented DataFrame-
#         producing entry point);
#     (b) it only helps if lons/lats are written with per-row-group
#         min/max stats AND the row groups have spatial locality (i.e. the
#         file isn't just one giant row group, and rows within a row group
#         are roughly spatially clustered rather than randomly ordered) —
#         neither was verified in this spike;
#     (c) no ready-made "spatial predicate pushdown" API exists in either
#         package as of the versions in this Manifest (GeoParquet 0.3.0,
#         Parquet2 per Manifest.toml). Implementing it would be a genuine
#         build task, not a config flag.
# ----------------------------------------------------------------------
println("See script comments (Step 6) for Parquet2.jl/GeoParquet.jl pushdown findings.")
