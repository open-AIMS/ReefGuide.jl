#!/usr/bin/env julia
# ============================================================================
# THROWAWAY SPIKE SCRIPT — not part of the package source.
#
# Phase 4 benchmark for the "DuckDB switch" plan
# (see .claude/plans/duckdb-switch.md, "Phase 4 — Benchmark & keep/retire decision").
#
# Repeats the Phase 0 methodology (N>=5 warm trials, mean/stddev/min/max,
# non-overlapping-range legitimacy test) across ALL FOUR regions and both
# full-load and scoped-load (5km/50km) cases, for both backends.
#
# Usage (from the ReefGuide.jl repo root):
#   julia --project=. scratch/phase4_benchmark.jl
#
# Data source: GBR-reef-guidance-assessment/outputs/MPA/ (fresh arrow+parquet+bounds
# for all 4 regions, produced in Phase 1).
# ============================================================================

using DataFrames
using Arrow
using DuckDB
using QuackIO
using Statistics
using Printf
using ReefGuide

const N_WARM_TRIALS = 5
const DATA_DIR = abspath(
    joinpath(@__DIR__, "..", "..", "GBR-reef-guidance-assessment", "outputs", "MPA")
)
const REGIONS = [
    "Cairns-Cooktown", "Townsville-Whitsunday", "Mackay-Capricorn", "FarNorthern"
]

function timed_trials(f::Function, n::Int)
    f() # first call incl. compilation, discarded
    times = Float64[]
    for _ in 1:n
        push!(times, @elapsed f())
    end
    return times
end

function stats(times::Vector{Float64})
    return (mean=mean(times), std=std(times), min=minimum(times), max=maximum(times))
end

function verdict(stats_a, stats_b)
    overlap = stats_a.max >= stats_b.min && stats_b.max >= stats_a.min
    return overlap ? "NEGLIGIBLE" : "LEGITIMATE"
end

function escape_sql(s::AbstractString)
    return replace(s, "'" => "''")
end

function duckdb_scoped_load(parquet_path, min_lon, min_lat, max_lon, max_lat)
    con = DuckDB.DB()
    qstr = """
        SELECT * FROM read_parquet('$(escape_sql(parquet_path))')
        WHERE lons BETWEEN $(min_lon) AND $(max_lon)
          AND lats BETWEEN $(min_lat) AND $(max_lat)
    """
    return DataFrame(DuckDB.execute(con, qstr))
end

function make_bbox(lon_mid, lat_mid, viewport_km)
    km_per_deg_lat = 111.32
    km_per_deg_lon = 111.32 * cosd(lat_mid)
    half_lat = (viewport_km / 2) / km_per_deg_lat
    half_lon = (viewport_km / 2) / km_per_deg_lon
    return (
        lon_mid - half_lon, lat_mid - half_lat, lon_mid + half_lon, lat_mid + half_lat
    )
end

results = DataFrame(;
    region=String[],
    load_type=String[],
    backend=String[],
    n_rows=Int[],
    mean_s=Float64[],
    std_s=Float64[],
    min_s=Float64[],
    max_s=Float64[]
)

for region in REGIONS
    arrow_path = joinpath(DATA_DIR, "$(region)_valid_slopes_lookup.arrow")
    parquet_path = joinpath(DATA_DIR, "$(region)_valid_slopes_lookup.parquet")

    if !isfile(arrow_path) || !isfile(parquet_path)
        @warn "Skipping $(region): missing arrow/parquet file" arrow_path parquet_path
        continue
    end

    println("="^78)
    println("REGION: $(region)")
    println("="^78)

    # --- Full load ---
    println("--- Full load ---")
    arrow_full_times = timed_trials(
        () -> DataFrame(Arrow.Table(arrow_path)), N_WARM_TRIALS
    )
    arrow_full_n = nrow(DataFrame(Arrow.Table(arrow_path)))
    arrow_full_stats = stats(arrow_full_times)
    @printf "  Arrow   full load: mean=%.3fs std=%.3fs min=%.3fs max=%.3fs (n=%d)\n" arrow_full_stats.mean arrow_full_stats.std arrow_full_stats.min arrow_full_stats.max arrow_full_n
    push!(
        results,
        (
            region, "full", "arrow", arrow_full_n, arrow_full_stats.mean,
            arrow_full_stats.std, arrow_full_stats.min, arrow_full_stats.max
        )
    )

    parquet_full_times = timed_trials(
        () -> QuackIO.read_parquet(DataFrame, parquet_path), N_WARM_TRIALS
    )
    parquet_full_n = nrow(QuackIO.read_parquet(DataFrame, parquet_path))
    parquet_full_stats = stats(parquet_full_times)
    @printf "  Parquet full load: mean=%.3fs std=%.3fs min=%.3fs max=%.3fs (n=%d)\n" parquet_full_stats.mean parquet_full_stats.std parquet_full_stats.min parquet_full_stats.max parquet_full_n
    push!(
        results,
        (
            region, "full", "parquet", parquet_full_n, parquet_full_stats.mean,
            parquet_full_stats.std, parquet_full_stats.min, parquet_full_stats.max
        )
    )
    println("  Verdict (full load): $(verdict(arrow_full_stats, parquet_full_stats))")

    # --- Scoped loads: 5km, 50km ---
    tbl = Arrow.Table(arrow_path)
    lons = tbl.lons
    lats = tbl.lats
    med_lon = median(lons)
    med_lat = median(lats)
    nearest_idx = argmin(@. (lons - med_lon)^2 + (lats - med_lat)^2)
    lon_mid = lons[nearest_idx]
    lat_mid = lats[nearest_idx]
    tbl = nothing
    GC.gc()

    for viewport_km in (5.0, 50.0)
        min_lon, min_lat, max_lon, max_lat = make_bbox(lon_mid, lat_mid, viewport_km)
        scope = ReefGuide.BBoxScope(min_lon, min_lat, max_lon, max_lat)

        println("--- Scoped load ($(viewport_km) km viewport) ---")
        arrow_scoped_times = timed_trials(
            () -> ReefGuide.load_scoped_arrow_table(arrow_path, scope), N_WARM_TRIALS
        )
        arrow_scoped_n = nrow(ReefGuide.load_scoped_arrow_table(arrow_path, scope))
        arrow_scoped_stats = stats(arrow_scoped_times)
        @printf "  Arrow   scoped(%gkm): mean=%.4fs std=%.4fs min=%.4fs max=%.4fs (n=%d)\n" viewport_km arrow_scoped_stats.mean arrow_scoped_stats.std arrow_scoped_stats.min arrow_scoped_stats.max arrow_scoped_n
        push!(
            results,
            (
                region, "scoped_$(Int(viewport_km))km", "arrow", arrow_scoped_n,
                arrow_scoped_stats.mean, arrow_scoped_stats.std, arrow_scoped_stats.min,
                arrow_scoped_stats.max
            )
        )

        duckdb_scoped_times = timed_trials(
            () -> duckdb_scoped_load(parquet_path, min_lon, min_lat, max_lon, max_lat),
            N_WARM_TRIALS
        )
        duckdb_scoped_n = nrow(
            duckdb_scoped_load(parquet_path, min_lon, min_lat, max_lon, max_lat)
        )
        duckdb_scoped_stats = stats(duckdb_scoped_times)
        @printf "  DuckDB  scoped(%gkm): mean=%.4fs std=%.4fs min=%.4fs max=%.4fs (n=%d)\n" viewport_km duckdb_scoped_stats.mean duckdb_scoped_stats.std duckdb_scoped_stats.min duckdb_scoped_stats.max duckdb_scoped_n
        push!(
            results,
            (
                region, "scoped_$(Int(viewport_km))km", "parquet", duckdb_scoped_n,
                duckdb_scoped_stats.mean, duckdb_scoped_stats.std, duckdb_scoped_stats.min,
                duckdb_scoped_stats.max
            )
        )

        if arrow_scoped_n != duckdb_scoped_n
            println(
                "  WARNING: row count mismatch! arrow=$(arrow_scoped_n) duckdb=$(duckdb_scoped_n)"
            )
        end
        println(
            "  Verdict ($(viewport_km)km scoped): $(verdict(arrow_scoped_stats, duckdb_scoped_stats))"
        )
    end
    println()
end

println("="^78)
println("SUMMARY TABLE")
println("="^78)
show(results; allrows=true, allcols=true)
println()

CSV_PATH = joinpath(@__DIR__, "phase4_benchmark_results.csv")
open(CSV_PATH, "w") do io
    println(io, join(names(results), ","))
    for row in eachrow(results)
        println(io, join(values(row), ","))
    end
end
println()
println("Results written to: $(CSV_PATH)")
