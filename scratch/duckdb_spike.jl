#!/usr/bin/env julia
# ============================================================================
# THROWAWAY SPIKE SCRIPT — not part of the package source.
#
# Phase 0 de-risking spike for the "DuckDB switch" plan
# (see .claude/plans/duckdb-switch.md, "Phase 0 — Spike").
#
# Question: does a genuine DuckDB query engine (DuckDB.jl / QuackIO.jl) reading
# a Parquet copy of the criteria lookup table round-trip losslessly and give a
# real performance benefit (full-table load and, especially, scoped/bbox load)
# over the current Arrow.jl-based path (`load_scoped_arrow_table` in
# `src/utility/regions_criteria_setup.jl`) — or is any difference within noise?
#
# This is a NEW spike, separate from `scratch/bbox_spike.jl` (which tested
# GeoParquet.jl/Parquet2.jl and found no usable predicate-pushdown API — that
# script is left untouched as the historical record of that finding).
#
# Usage (from the ReefGuide.jl repo root):
#   julia --project=. scratch/duckdb_spike.jl [path/to/region_valid_slopes_lookup.arrow]
#
# Default target file (if no arg given) is the real Cairns-Cooktown region
# lookup Arrow file from the sibling ReefGuideWorker.jl/data/ directory (real
# production-shaped data, not synthetic) — the smallest of the four regional
# files, matching bbox_spike.jl's existing choice for comparability.
#
# The script writes its own scratch Parquet copy (via `QuackIO.write_table`)
# to a temp directory — it does NOT write into ReefGuideWorker.jl/data/, which
# is real production data.
# ============================================================================

using DataFrames
using Arrow
using DuckDB
using QuackIO
using Statistics
using Printf
using ReefGuide

const DEFAULT_ARROW = joinpath(
    @__DIR__, "..", "..", "ReefGuideWorker.jl", "data",
    "Cairns-Cooktown_valid_slopes_lookup.arrow"
)
const N_WARM_TRIALS = 5

arrow_path = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_ARROW
arrow_path = abspath(arrow_path)

scratch_dir = mktempdir()
parquet_path = joinpath(scratch_dir, "Cairns-Cooktown_valid_slopes_lookup.parquet")

println("="^78)
println("DUCKDB/PARQUET SPIKE — round-trip fidelity + full/scoped load benchmarks")
println("="^78)
println("Arrow source file: $(arrow_path)")
println("Arrow file size:   $(round(filesize(arrow_path) / 1024^2, digits=1)) MB")
println("Scratch parquet:   $(parquet_path)")
println()

# ----------------------------------------------------------------------
# Small helpers shared by the benchmark sections below.
# ----------------------------------------------------------------------
"""Run `f()` once (discarded, "incl. compilation") then `n` warm trials, returning the warm elapsed times."""
function timed_trials(f::Function, n::Int; label::String)
    compile_time = @elapsed f()
    @printf "  First call (incl. compilation): %.4f s\n" compile_time
    times = Float64[]
    for _ in 1:n
        push!(times, @elapsed f())
    end
    @printf "  Warm trials (n=%d) %s: %s\n" n label join(round.(times; digits=4), ", ")
    return times
end

function report_stats(label::String, times::Vector{Float64})
    @printf "  %-28s mean=%.4f s  std=%.4f s  min=%.4f s  max=%.4f s\n" label mean(times) std(
        times
    ) minimum(times) maximum(times)
    return (mean=mean(times), std=std(times), min=minimum(times), max=maximum(times))
end

"""Ranges overlap ⇒ difference is noise/negligible, not a legitimate win."""
function verdict(name_a, stats_a, name_b, stats_b)
    overlap = stats_a.max >= stats_b.min && stats_b.max >= stats_a.min
    if overlap
        println(
            "  VERDICT: ranges overlap ([$(round(stats_a.min,digits=4)), $(round(stats_a.max,digits=4))] vs [$(round(stats_b.min,digits=4)), $(round(stats_b.max,digits=4))]) → NEGLIBLE/NOISE, not a legitimate difference."
        )
    else
        faster = stats_a.mean < stats_b.mean ? name_a : name_b
        println(
            "  VERDICT: ranges do NOT overlap ([$(round(stats_a.min,digits=4)), $(round(stats_a.max,digits=4))] vs [$(round(stats_b.min,digits=4)), $(round(stats_b.max,digits=4))]) → LEGITIMATE difference, $(faster) is faster."
        )
    end
    return !overlap
end

# ----------------------------------------------------------------------
# Step 1: Load the real Arrow file, exactly as load_target_region does for
# a full (scope=nothing) load.
# ----------------------------------------------------------------------
println("--- Step 1: load real Arrow file ---")
local arrow_df
compile_time = @elapsed begin
    arrow_df = DataFrame(Arrow.Table(arrow_path))
end
@printf "First call (incl. compilation): %.3f s\n" compile_time
load_time = @elapsed begin
    global arrow_df = DataFrame(Arrow.Table(arrow_path))
end
@printf "Warm load time:                 %.3f s\n" load_time
println("Row count: $(nrow(arrow_df))")
println("Columns:   $(names(arrow_df))")
println()

# ----------------------------------------------------------------------
# Step 2: Write a scratch Parquet copy via QuackIO, then read it back.
# ----------------------------------------------------------------------
println("--- Step 2: QuackIO.write_table (Arrow DataFrame -> scratch Parquet) ---")
write_time = @elapsed QuackIO.write_table(parquet_path, arrow_df; format=:parquet)
@printf "Write time: %.3f s\n" write_time
println("Parquet file size: $(round(filesize(parquet_path) / 1024^2, digits=1)) MB")
println()

println("--- Step 3: QuackIO.read_parquet (scratch Parquet -> DataFrame) ---")
parquet_df = QuackIO.read_parquet(DataFrame, parquet_path)
println("Row count: $(nrow(parquet_df))")
println("Columns:   $(names(parquet_df))")
println()

# ----------------------------------------------------------------------
# Question 1: Round-trip fidelity.
# ----------------------------------------------------------------------
println("="^78)
println("QUESTION 1: ROUND-TRIP FIDELITY")
println("="^78)

same_row_count = nrow(arrow_df) == nrow(parquet_df)
println("Same row count:  $(same_row_count)  ($(nrow(arrow_df)) vs $(nrow(parquet_df)))")

same_columns = Set(names(arrow_df)) == Set(names(parquet_df))
println("Same column set: $(same_columns)")

println("Per-column dtype comparison (Arrow -> Parquet):")
dtypes_match = true
for col in names(arrow_df)
    arrow_type = eltype(arrow_df[!, col])
    parquet_type = eltype(parquet_df[!, col])
    ok = arrow_type == parquet_type
    global dtypes_match &= ok
    @printf "  %-12s %-28s -> %-28s  %s\n" col string(arrow_type) string(parquet_type) (
        ok ? "OK" : "MISMATCH"
    )
end
has_uint16 = any(
    eltype(arrow_df[!, c]) == UInt16 || eltype(parquet_df[!, c]) == UInt16 for
    c in names(arrow_df)
)
println("Any UInt16 column present: $(has_uint16)")
println()

# Row order: is X-major (lon_idx, lat_idx) sort preserved through the
# write/read round trip?
arrow_order_sorted = issorted(collect(zip(arrow_df.lon_idx, arrow_df.lat_idx)))
order_preserved =
    arrow_df.lon_idx == parquet_df.lon_idx && arrow_df.lat_idx == parquet_df.lat_idx
println("Arrow source already X-major (lon_idx,lat_idx) sorted: $(arrow_order_sorted)")
println("Row order preserved Arrow -> Parquet round trip:       $(order_preserved)")

# Spot-check N random rows for exact value equality (including missing
# positions), independent of the row-order question above.
rng_idx = rand(1:nrow(arrow_df), 25)
spot_check_ok = all(
    isequal(arrow_df[i, :], parquet_df[i, :]) for i in rng_idx
)
println("Spot-check (25 random rows, positional) exact equality: $(spot_check_ok)")
println()

# Downstream reliance on row order: consumers (apply_criteria.jl,
# best_fit_polygons.jl, site_identification.jl) all index by the *values* of
# the lon_idx/lat_idx columns (e.g. CartesianIndex.(tmp.lon_idx, tmp.lat_idx),
# indicator[r.lon_idx, r.lat_idx] = ...), never by row position — confirmed by
# grepping regions_criteria_setup.jl's consumers. So even if row order were
# NOT preserved, downstream correctness would be unaffected.
println(
    "Downstream reliance on row order: NONE FOUND. apply_criteria.jl:149-150, " *
    "best_fit_polygons.jl:299-314 and site_identification.jl:68-219 all index " *
    "results by the lon_idx/lat_idx *values* in each row (e.g. " *
    "CartesianIndex.(tmp.lon_idx, tmp.lat_idx), indicator[r.lon_idx, r.lat_idx]), " *
    "never by row position — row order is cosmetic for correctness either way."
)
println()

# ----------------------------------------------------------------------
# Question 2: Full-table load benchmark, Arrow vs DuckDB/QuackIO Parquet.
# ----------------------------------------------------------------------
println("="^78)
println("QUESTION 2: FULL-TABLE LOAD BENCHMARK (N=$(N_WARM_TRIALS) warm trials)")
println("="^78)

println("Arrow (DataFrame(Arrow.Table(path))):")
arrow_full_times = timed_trials(
    () -> DataFrame(Arrow.Table(arrow_path)), N_WARM_TRIALS; label="arrow full load"
)
arrow_full_stats = report_stats("Arrow full load", arrow_full_times)
println()

println("DuckDB/QuackIO (QuackIO.read_parquet(DataFrame, path)):")
parquet_full_times = timed_trials(
    () -> QuackIO.read_parquet(DataFrame, parquet_path),
    N_WARM_TRIALS;
    label="parquet full load"
)
parquet_full_stats = report_stats("Parquet full load", parquet_full_times)
println()
full_load_legitimate = verdict(
    "Arrow", arrow_full_stats, "Parquet/DuckDB", parquet_full_stats
)
println()

# ----------------------------------------------------------------------
# Construct 5km/50km viewport bboxes, same convention as bbox_spike.jl:
# median lon/lat center of the actual data, km-per-degree conversion using
# the median latitude.
# ----------------------------------------------------------------------
lons = arrow_df.lons
lats = arrow_df.lats
lon_mid = median(lons)
lat_mid = median(lats)
km_per_deg_lat = 111.32
km_per_deg_lon = 111.32 * cosd(lat_mid)

function make_bbox(viewport_km::Float64)
    half_span_lat = (viewport_km / 2) / km_per_deg_lat
    half_span_lon = (viewport_km / 2) / km_per_deg_lon
    return (
        lon_mid - half_span_lon,
        lat_mid - half_span_lat,
        lon_mid + half_span_lon,
        lat_mid + half_span_lat
    )
end

println("--- Viewport construction ---")
@printf "median center: (lon=%.5f, lat=%.5f)\n" lon_mid lat_mid
@printf "km_per_deg_lon at lat=%.3f: %.3f\n" lat_mid km_per_deg_lon
println()

# ----------------------------------------------------------------------
# Question 3: Scoped (bbox) load benchmark, 5km and 50km viewports.
#   - Arrow path: ReefGuide.load_scoped_arrow_table (full-decompress-then-mask).
#   - DuckDB path: raw SQL predicate pushdown via read_parquet(...) WHERE ...
# ----------------------------------------------------------------------
println("="^78)
println("QUESTION 3: SCOPED (BBOX) LOAD BENCHMARK (N=$(N_WARM_TRIALS) warm trials)")
println("="^78)

escape_sql(s::AbstractString) = replace(s, "'" => "''")

function duckdb_scoped_load(min_lon, min_lat, max_lon, max_lat)
    con = DuckDB.DB()
    qstr = """
        SELECT * FROM read_parquet('$(escape_sql(parquet_path))')
        WHERE lons BETWEEN $(min_lon) AND $(max_lon)
          AND lats BETWEEN $(min_lat) AND $(max_lat)
    """
    return DataFrame(DuckDB.execute(con, qstr))
end

for viewport_km in (5.0, 50.0)
    min_lon, min_lat, max_lon, max_lat = make_bbox(viewport_km)
    scope = ReefGuide.BBoxScope(min_lon, min_lat, max_lon, max_lat)
    @printf "--- %.0f km viewport --- bbox lon=[%.6f, %.6f] lat=[%.6f, %.6f]\n" viewport_km min_lon max_lon min_lat max_lat

    println("Arrow (ReefGuide.load_scoped_arrow_table):")
    local arrow_n
    arrow_scoped_times = timed_trials(
        () -> (arrow_n = nrow(ReefGuide.load_scoped_arrow_table(arrow_path, scope))),
        N_WARM_TRIALS;
        label="arrow scoped load ($(viewport_km)km)"
    )
    println("  Row count (scoped): $(arrow_n)")
    arrow_scoped_stats = report_stats("Arrow scoped ($(viewport_km)km)", arrow_scoped_times)
    println()

    println("DuckDB (raw SQL predicate pushdown via read_parquet WHERE ...):")
    local duckdb_n
    duckdb_scoped_times = timed_trials(
        () -> (duckdb_n = nrow(duckdb_scoped_load(min_lon, min_lat, max_lon, max_lat))),
        N_WARM_TRIALS;
        label="duckdb scoped load ($(viewport_km)km)"
    )
    println("  Row count (scoped): $(duckdb_n)")
    duckdb_scoped_stats = report_stats(
        "DuckDB scoped ($(viewport_km)km)", duckdb_scoped_times
    )
    println()

    if arrow_n != duckdb_n
        println(
            "  WARNING: row counts differ between Arrow scoped ($(arrow_n)) and DuckDB scoped ($(duckdb_n)) loads for the same bbox!"
        )
    end

    verdict(
        "Arrow scoped ($(viewport_km)km)",
        arrow_scoped_stats,
        "DuckDB scoped ($(viewport_km)km)",
        duckdb_scoped_stats
    )
    println()
end

# ----------------------------------------------------------------------
# Cleanup scratch parquet file/dir (this script must not leave production
# directories touched, and the scratch copy is throwaway).
# ----------------------------------------------------------------------
rm(scratch_dir; recursive=true, force=true)

println("="^78)
println("DONE")
println("="^78)
