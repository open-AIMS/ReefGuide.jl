
"""
Convert BoundedCriteriaDict to a vector of CriteriaBounds for assessment processing.
Transforms the BoundedCriteriaDict into CriteriaBounds objects that include
evaluation functions. Only includes criteria that are available in the dictionary.

# Arguments
- `bounded_criteria_dict::BoundedCriteriaDict` : Dictionary of bounded criteria to convert

# Returns
Vector of `CriteriaBounds` objects for available criteria.
"""
function build_criteria_bounds_from_regional_criteria(
    bounded_criteria_dict::BoundedCriteriaDict
)::Vector{CriteriaBounds}
    @debug "Converting BoundedCriteriaDict to CriteriaBounds vector"
    criteria_bounds = CriteriaBounds[]

    for (criteria_id, bounded_criteria) in bounded_criteria_dict
        bounds = CriteriaBounds(
            # Field to get in the data
            bounded_criteria.metadata.id,
            # Min/max bounds
            bounded_criteria.bounds.min,
            bounded_criteria.bounds.max
        )
        push!(criteria_bounds, bounds)
        @debug "Added criteria bounds" criteria_id = criteria_id min_val =
            bounded_criteria.bounds.min max_val = bounded_criteria.bounds.max
    end

    @debug "Built CriteriaBounds vector" total_criteria = length(criteria_bounds) criteria_ids = [
        cb.name for cb in criteria_bounds
    ]
    return criteria_bounds
end

"""
Convert a min/max tuple to a Bounds struct.

# Arguments
- `min_max::Tuple{Number,Number}` : Tuple containing (minimum, maximum) values

# Returns
`Bounds` struct with converted float values.
"""
function bounds_from_tuple(min_max::Tuple{Number,Number})::Bounds
    return Bounds(; min=min_max[1], max=min_max[2])
end

"""
Generate the filename for slope lookup data for a given region.

# Arguments
- `region::RegionMetadata` : Region metadata containing ID

# Returns
String filename in format "{region_id}_slope_lookup.arrow"
"""
function get_slope_filename(region::RegionMetadata)::String
    filename = "$(region.id)$(SLOPES_LOOKUP_SUFFIX)"
    @debug "Generated slope filename" region_id = region.id filename
    return filename
end

function dms_to_decimal(d, m, s, o)
    # 40°26'46"N

    decimal_degrees = d + (m / 60) + (s / 3600)

    if lowercase(o) == "s" || lowercase(o) == "w"
        # South, so invert
        decimal_degrees = -decimal_degrees
    end

    return decimal_degrees
end
function dms_to_decimal(dms::String)
    expanded = split(dms)
    if length(expanded) == 4
        if lowercase(expanded[1]) in ["n", "e"]
            o, d, m, s = expanded
        elseif lowercase(expanded[1]) in ["s", "w"]
            o, d, m, s = expanded
        else
            throw(ArgumentError("Unknown format!"))
        end

        return dms_to_decimal(d, m, s, o)
    end

    mat = match(r"(\d+)°(\d+)'(\d+)\"([NSEW])", dms)

    if !isnothing(mat)
        d = parse(Int, mat.captures[1])
        m = parse(Int, mat.captures[2])
        s = parse(Int, mat.captures[3])
        orient = mat.captures[4][1]
    else
        mat = match(r"([NSEW])(\d+)°(\d+)'(\d+)\"", dms)
        orient = mat.captures[1][1]
        d = parse(Int, mat.captures[2])
        m = parse(Int, mat.captures[3])
        s = parse(Int, mat.captures[4])
    end

    return dms_to_decimal(d, m, s, orient)
end