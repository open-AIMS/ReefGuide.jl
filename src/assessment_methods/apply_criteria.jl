"""Methods to filter criteria bounds over rasters and lookup tables"""

"""
    _preference_score(direction, value, lower_bound, upper_bound, band_peak)::Float32

Map a raw criterion `value` to a normalised `[0,1]` suitability score for the given
MCDA preference `direction`. Values outside the ramp domain are clamped.

- `:lower_is_better` — `1` at `lower_bound`, linear down to `0` at `upper_bound`
  (`< lb` → `1`, `> ub` → `0`).
- `:higher_is_better` — `0` at `lower_bound`, linear up to `1` at `upper_bound`
  (`< lb` → `0`, `> ub` → `1`).
- `:band` — triangle: `1` at `band_peak`, linear down to `0` at each of `lower_bound`
  and `upper_bound`; outside `[lb, ub]` → `0`. `band_peak` is expected strictly between
  the bounds; `nothing` falls back to the band midpoint.

A degenerate domain (`lower_bound == upper_bound`, or `band_peak` sitting on a bound) is
guarded against divide-by-zero: it returns `1.0f0` where `value` is in-band, else `0.0f0`.
"""
function _preference_score(
    direction::Symbol,
    value::Real,
    lower_bound::Float32,
    upper_bound::Float32,
    band_peak::OptionalValue{Float32}
)::Float32
    v::Float32 = Float32(value)

    if direction === :band
        peak::Float32 =
            isnothing(band_peak) ? (lower_bound + upper_bound) / 2.0f0 :
            band_peak
        (v < lower_bound || v > upper_bound) && return 0.0f0
        if upper_bound - lower_bound <= 0.0f0
            return 1.0f0
        elseif v <= peak
            return peak - lower_bound <= 0.0f0 ? 1.0f0 :
                   (v - lower_bound) / (peak - lower_bound)
        else
            return upper_bound - peak <= 0.0f0 ? 1.0f0 :
                   (upper_bound - v) / (upper_bound - peak)
        end
    elseif direction === :lower_is_better
        span::Float32 = upper_bound - lower_bound
        span <= 0.0f0 && return v <= lower_bound ? 1.0f0 : 0.0f0
        v <= lower_bound && return 1.0f0
        v >= upper_bound && return 0.0f0
        return (upper_bound - v) / span
    elseif direction === :higher_is_better
        span = upper_bound - lower_bound
        span <= 0.0f0 && return v >= upper_bound ? 1.0f0 : 0.0f0
        v <= lower_bound && return 0.0f0
        v >= upper_bound && return 1.0f0
        return (v - lower_bound) / span
    else
        throw(ArgumentError("Unknown MCDA preference direction: $(direction)"))
    end
end

"""
Bundles a criterion's lookup column name and bounds with two derived functions:
`rule`, the boolean threshold check, and `score`, the normalised `[0,1]` MCDA
suitability map.

Lookup columns are nullable (`Union{Missing,T}`): a pixel valid under the region's
bathymetry-bound mask may still lack data for an individual criterion. `rule` folds
a `missing` to `missing_weight > 0` and always returns a plain boolean mask, so the
default `missing_weight = 0.0f0` excludes no-data pixels. `score` returns `0.0f0`
for a `missing` element; [`score_lookup_table_by_criteria`](@ref) substitutes
`missing_weight` on the `[0,1]` scale when aggregating.
"""
struct CriteriaBounds{F<:Function,G<:Function}
    "The field ID of the criteria"
    name::Symbol
    "Lower bound, inclusive"
    lower_bound::Float32
    "Upper bound, inclusive"
    upper_bound::Float32
    "Score substituted for a `missing` value; `> 0` also passes the boolean rule"
    missing_weight::Float32
    "MCDA preference direction"
    direction::Symbol
    "Value at which suitability peaks when `direction == :band`"
    band_peak::OptionalValue{Float32}
    "Relative weight in the weighted-mean aggregate"
    weight::Float32
    "Takes a value or vector of values, returns whether each matches the criteria"
    rule::F
    "Takes a vector of values, returns the normalised `[0,1]` suitability per element"
    score::G

    function CriteriaBounds(
        name::String, lb::Float32, ub::Float32, missing_weight::Float32=0.0f0;
        direction::Symbol=:higher_is_better,
        band_peak::OptionalValue{Float32}=nothing,
        weight::Float32=1.0f0
    )::CriteriaBounds
        allow_missing::Bool = missing_weight > 0.0f0
        rule = (x) -> coalesce.(lb .<= x .<= ub, allow_missing)
        score =
            (x) -> Float32[
                ismissing(v) ? 0.0f0 : _preference_score(direction, v, lb, ub, band_peak)
                for v in x
            ]
        return new{Function,Function}(
            Symbol(name), lb, ub, missing_weight, direction, band_peak, weight, rule, score
        )
    end

    function CriteriaBounds(
        name::S, lb::S, ub::S, missing_weight::Float32=0.0f0;
        direction::Symbol=:higher_is_better,
        band_peak::OptionalValue{Float32}=nothing,
        weight::Float32=1.0f0
    )::CriteriaBounds where {S<:String}
        return CriteriaBounds(
            name, parse(Float32, lb), parse(Float32, ub), missing_weight;
            direction, band_peak, weight
        )
    end
end

"""
Apply thresholds for each criteria.

# Arguments
- `criteria_stack` : RasterStack of criteria data for a given region
- `lookup` : Lookup dataframe for the region
- `criteria_bounds` : A vector of CriteriaBounds which contains named criteria
  with min/max ranges and a function to apply.

# Returns
BitMatrix indicating locations within desired thresholds
"""
function filter_raster_by_criteria(
    criteria_stack::RasterStack,
    lookup::DataFrame,
    criteria_bounds::Vector{CriteriaBounds}
)::Raster
    # Result store
    data = falses(size(criteria_stack))

    # Apply criteria
    res_lookup = trues(nrow(lookup))
    for filter::CriteriaBounds in criteria_bounds
        res_lookup .= res_lookup .& filter.rule(lookup[!, filter.name])
    end

    tmp = lookup[res_lookup, [:lon_idx, :lat_idx]]
    data[CartesianIndex.(tmp.lon_idx, tmp.lat_idx)] .= true

    res = Raster(criteria_stack.Depth; data=sparse(data), missingval=0)
    return res
end

"""
Filters a lookup table (which contains raster param values too) by building a
bit mask AND'd for all thresholds

# Arguments
- `lookup` : A lookup table to filter
- `ruleset` : Vector of `CriteriaBounds` to filter with

# Returns
BitVector, of pixels that meet criteria
"""
function filter_lookup_table_by_criteria(
    lookup::DataFrame,
    ruleset::Vector{CriteriaBounds}
)::BitVector
    matches::BitVector = fill(true, nrow(lookup))

    for threshold in ruleset
        matches = matches .& threshold.rule(lookup[!, threshold.name])
    end

    return matches
end

"""
Compute the continuous MCDA suitability score for every row of `lookup`.

For each `CriteriaBounds` in `ruleset` the named column is mapped to a normalised
`[0,1]` per-criterion score via `cb.score` (the preference function keyed by
`cb.direction`); a `missing` cell is replaced by `cb.missing_weight` on that same
scale. The per-criterion scores are combined per row as a weighted mean using
`cb.weight` (`sum(weightᵢ * scoreᵢ) / sum(weightᵢ)`), which reduces to the plain
equal-weight mean when all weights are `1.0f0` (the default). `1` = a perfect
match across all criteria.

# Arguments
- `lookup` : A lookup table carrying the criteria value columns
- `ruleset` : Vector of `CriteriaBounds` to score with

# Returns
`Vector{Float32}` of length `nrow(lookup)`. Values land in `[0,1]` unless a
`missing_weight` is deliberately set outside that range.
"""
function score_lookup_table_by_criteria(
    lookup::DataFrame,
    ruleset::Vector{CriteriaBounds}
)::Vector{Float32}
    n = nrow(lookup)
    weighted_sum = zeros(Float32, n)
    weight_total = 0.0f0

    for cb in ruleset
        col = lookup[!, cb.name]
        per_criterion::Vector{Float32} = cb.score(col)
        @inbounds for i in 1:n
            contribution = ismissing(col[i]) ? cb.missing_weight : per_criterion[i]
            weighted_sum[i] += cb.weight * contribution
        end
        weight_total += cb.weight
    end

    weight_total == 0.0f0 && return weighted_sum
    return weighted_sum ./ weight_total
end

"""
    lookup_df_from_raster(raster::Raster, threshold::Union{Int64,Float64})::DataFrame

Build a look up table identifying all pixels in a raster that meet a suitability threshold.

# Arguments
- `raster` : Raster of regional data
- `threshold` : Suitability threshold value (greater or equal than)

# Returns
DataFrame containing indices, lon and lat for each pixel that is intended for further
analysis.
"""
function lookup_df_from_raster(raster::Raster, threshold::Union{Int64,Float64})::DataFrame
    criteria_matches::SparseMatrixCSC{Bool,Int64} = sparse(falses(size(raster)))
    Rasters.read!(raster .>= threshold, criteria_matches)
    indices::Vector{CartesianIndex{2}} = findall(criteria_matches)
    indices_lon::Vector{Float64} = lookup(raster, X)[first.(Tuple.(indices))]
    indices_lat::Vector{Float64} = lookup(raster, Y)[last.(Tuple.(indices))]

    return DataFrame(; indices=indices, lons=indices_lon, lats=indices_lat)
end
