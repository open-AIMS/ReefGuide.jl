# =============================================================================
# Assessment Parameters Constants
# =============================================================================

const DEFAULT_SUITABILITY_THRESHOLD::Int64 = 80

# =============================================================================
# Assessment Parameters Data Structures
# =============================================================================

"""
Payload for regional assessment actions - this includes all merged bounds and
regional data.

# Fields
- `region::String` : The region that is being assessed
- `regional_criteria::BoundedCriteriaDict` : The criteria to assess, including user provided bounds
- `region_data::RegionalDataEntry` : The data to consider for this region
"""
struct RegionalAssessmentParameters
    region::String
    regional_criteria::BoundedCriteriaDict
    region_data::RegionalDataEntry

    function RegionalAssessmentParameters(;
        region::String,
        regional_criteria::BoundedCriteriaDict,
        region_data::RegionalDataEntry
    )
        return new(region, regional_criteria, region_data)
    end
end

"""
Payload for suitability assessment actions - this includes all merged bounds and
regional data plus spatial dimensions.

# Fields
- `region::String` : The region that is being assessed
- `regional_criteria::BoundedCriteriaDict` : The criteria to assess, including user provided bounds
- `region_data::RegionalDataEntry` : The data to consider for this region
- `suitability_threshold::Int64` : The cutoff to consider a site suitable
- `x_dist::Int64` : X dimension of polygon (metres)
- `y_dist::Int64` : Y dimension of polygon (metres)
"""
struct SuitabilityAssessmentParameters
    # Regional criteria
    region::String
    regional_criteria::BoundedCriteriaDict
    region_data::RegionalDataEntry
    suitability_threshold::Int64

    # Additional criteria
    x_dist::Int64
    y_dist::Int64

    function SuitabilityAssessmentParameters(;
        region::String,
        regional_criteria::BoundedCriteriaDict,
        region_data::RegionalDataEntry,
        suitability_threshold::Int64,
        x_dist::Int64,
        y_dist::Int64
    )
        @debug "Created SuitabilityAssessmentParameters" region suitability_threshold x_dist y_dist
        return new(
            region, regional_criteria, region_data, suitability_threshold, x_dist, y_dist
        )
    end
end
