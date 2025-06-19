# ReefGuide.jl

A Julia package for coral reef site assessment and suitability analysis. ReefGuide provides tools for evaluating potential deployment locations on coral reefs based on multiple environmental criteria including depth, slope, turbidity, wave conditions, and rugosity.

## Prerequisites

- Git
- Julia v1.11 or later (managed via juliaup)

## Installation & Setup

1. **Install Julia using juliaup**

   ```bash
   # Install juliaup (Julia version manager)
   curl -fsSL https://install.julialang.org | sh

   # Install and set Julia v1.11 as default
   juliaup add 1.11
   juliaup default 1.11
   ```

2. **Clone and setup the project**

   ```bash
   git clone <repository-url>
   cd ReefGuide.jl
   julia --project=. -e "using Pkg; Pkg.instantiate()"
   ```

3. **Development environment (recommended)**
   ```bash
   cd sandbox
   ./init.sh    # Sets up local project with Revise.jl for development
   ./start.sh   # Starts Julia with hot-reloading enabled
   ```

## Quick Start

### 1. Initialize Data

First, specify your data directory and load the regional datasets. For the data spec, see [spec](#data-specification).

```julia
using ReefGuide

# Point to your data directory containing regional GeoTIFF and Parquet files
data_dir = "/path/to/your/data"
regional_data = initialise_data(data_dir)
```

The `initialise_data()` function loads:

- Regional raster stacks (depth, slope, turbidity, wave data, etc.)
- Slope lookup tables with valid reef coordinates
- Canonical reef outline geometries
- Computed criteria bounds for each region

### 2. Regional Assessment

Perform raster-based suitability assessment over entire regions:

```julia
# Configure assessment parameters
params = RegionalAssessmentParameters(
    region="Townsville-Whitsunday",
    regional_criteria=criteria_dict,  # User-defined bounds for each criteria
    region_data=regional_data.regions["Townsville-Whitsunday"]
)

# Generate suitability raster
suitability_raster = assess_region(params)
```

### 3. Site Assessment

Identify optimal deployment sites with custom polygon fitting:

```julia
# Configure site assessment parameters
params = SuitabilityAssessmentParameters(
    region="Townsville-Whitsunday",
    regional_criteria=criteria_dict,
    region_data=regional_data.regions["Townsville-Whitsunday"],
    suitability_threshold=80,  # Minimum suitability percentage
    x_dist=100,  # Polygon width in meters
    y_dist=100   # Polygon height in meters
)

# Find optimal sites
sites = assess_sites(params)
filtered_sites = filter_sites(sites)
```

### 4. Export Results

```julia
# Export sites as GeoJSON
output_geojson("results/optimal_sites.geojson", filtered_sites)

# Export raster as GeoTIFF
_write_tiff("results/suitability.tif", suitability_raster)
```

## Available Regions

- **Townsville/Whitsunday Management Area** (`Townsville-Whitsunday`)
- **Cairns/Cooktown Management Area** (`Cairns-Cooktown`)
- **Mackay/Capricorn Management Area** (`Mackay-Capricorn`)
- **Far Northern Management Area** (`FarNorthern`)

## Assessment Criteria

The package supports assessment based on:

- **Depth** - Water depth from mean astronomical tide
- **Slope** - Reef slope angle in degrees
- **Turbidity** - Water clarity (Secchi depth)
- **Wave Height** - Significant wave height (90th percentile)
- **Wave Period** - Time between waves (90th percentile)
- **Rugosity** - Sea floor roughness (Townsville region only)

## Key Functions

- `initialise_data(data_dir)` - Load all regional data
- `assess_region(params)` - Regional suitability assessment
- `assess_sites(params)` - Site-specific polygon assessment
- `filter_sites(results)` - Remove overlapping/low-quality sites
- `output_geojson(path, data)` - Export results to GeoJSON

## Development and methodology

The package uses:

- **Rasters.jl** for geospatial raster operations
- **GeoDataFrames.jl** for vector geometry handling
- **ArchGDAL.jl** for GDAL integration
- **GeometryOps.jl** for geometric computations
- **FLoops.jl** for parallel processing

For detailed scientific methodology and parameter guidance, consult with the project team.

## Data specification

### Summary Table

| File Type               | Format     | Pattern                             | Regions                  | Purpose                                                                       |
| ----------------------- | ---------- | ----------------------------------- | ------------------------ | ----------------------------------------------------------------------------- |
| **Slope Lookup Table**  | Parquet    | `{region}_valid_slopes_lookup.parq` | All (4 files)            | Point data with coordinates and environmental values for valid reef locations |
| **Valid Extent Raster** | GeoTIFF    | `{region}_valid_slopes.tif`         | All (4 files)            | Binary mask indicating valid reef slope areas                                 |
| **Depth Raster**        | GeoTIFF    | `{region}*_bathy.tif`               | All (4 files)            | Bathymetry data from Mean Astronomical Tide                                   |
| **Slope Raster**        | GeoTIFF    | `{region}*_slope.tif`               | All (4 files)            | Reef slope angles in degrees                                                  |
| **Turbidity Raster**    | GeoTIFF    | `{region}*_turbid.tif`              | All (4 files)            | Water clarity (Secchi depth in meters)                                        |
| **Wave Height Raster**  | GeoTIFF    | `{region}*_waves_Hs.tif`            | All (4 files)            | Significant wave height, 90th percentile (meters)                             |
| **Wave Period Raster**  | GeoTIFF    | `{region}*_waves_Tp.tif`            | All (4 files)            | Wave period, 90th percentile (seconds)                                        |
| **Rugosity Raster**     | GeoTIFF    | `{region}*_rugosity.tif`            | Townsville only (1 file) | Sea floor roughness (standard deviation)                                      |
| **Reef Outlines**       | GeoPackage | `rrap_canonical_outlines.gpkg`      | Global (1 file)          | Canonical reef boundary geometries                                            |

### Detailed File Specifications

#### Slope Lookup Tables (`{region}_valid_slopes_lookup.parq`)

- **Format**: Apache Parquet with spatial geometry
- **Required Columns**:
  - `geometry` - Point geometries for valid reef locations
  - `lon_idx`, `lat_idx` - Raster grid indices
  - `Depth` - Bathymetry values (meters)
  - `Slope` - Slope angles (degrees)
  - `Turbidity` - Secchi depth (meters)
  - `WavesHs` - Wave height (meters)
  - `WavesTp` - Wave period (seconds)
  - `Rugosity` - Sea floor roughness (Townsville only)
- **Auto-generated**: `lons`, `lats` - Decimal degree coordinates (added by code)
- **Purpose**: Pre-filtered point dataset containing environmental values for computationally efficient site assessment

#### Valid Extent Rasters (`{region}_valid_slopes.tif`)

- **Format**: GeoTIFF raster
- **Data Type**: Binary (0/1 or boolean)
- **CRS**: Same as regional raster stack
- **Purpose**: Spatial mask defining areas with valid reef slopes suitable for assessment

#### Environmental Criteria Rasters (`{region}*_{criteria}.tif`)

- **Format**: GeoTIFF raster
- **Data Type**: Float32 or appropriate numeric type
- **CRS**: Consistent across all rasters in a region
- **Units**: As specified in criteria metadata
- **Purpose**: Gridded environmental data for spatial analysis and visualization

#### Reef Outlines (`rrap_canonical_outlines.gpkg`)

- **Format**: GeoPackage vector file
- **Geometry Type**: Polygon/MultiPolygon
- **Required Fields**: Geometry column with reef boundary polygons
- **CRS**: Geographic coordinate system (typically EPSG:4326)
- **Purpose**: Reference reef boundaries for spatial analysis and edge detection

### File Tree Structure

```
data/
├── rrap_canonical_outlines.gpkg                    # Global reef outlines
├── Townsville-Whitsunday_valid_slopes_lookup.parq  # Slope lookup table
├── Townsville-Whitsunday_valid_slopes.tif          # Valid extent mask
├── Townsville-Whitsunday_*_bathy.tif               # Depth raster
├── Townsville-Whitsunday_*_slope.tif               # Slope raster
├── Townsville-Whitsunday_*_turbid.tif              # Turbidity raster
├── Townsville-Whitsunday_*_waves_Hs.tif            # Wave height raster
├── Townsville-Whitsunday_*_waves_Tp.tif            # Wave period raster
├── Townsville-Whitsunday_*_rugosity.tif            # Rugosity raster (unique to Townsville)
├── Cairns-Cooktown_valid_slopes_lookup.parq        # Slope lookup table
├── Cairns-Cooktown_valid_slopes.tif                # Valid extent mask
├── Cairns-Cooktown_*_bathy.tif                     # Depth raster
├── Cairns-Cooktown_*_slope.tif                     # Slope raster
├── Cairns-Cooktown_*_turbid.tif                    # Turbidity raster
├── Cairns-Cooktown_*_waves_Hs.tif                  # Wave height raster
├── Cairns-Cooktown_*_waves_Tp.tif                  # Wave period raster
├── Mackay-Capricorn_valid_slopes_lookup.parq       # Slope lookup table
├── Mackay-Capricorn_valid_slopes.tif               # Valid extent mask
├── Mackay-Capricorn_*_bathy.tif                    # Depth raster
├── Mackay-Capricorn_*_slope.tif                    # Slope raster
├── Mackay-Capricorn_*_turbid.tif                   # Turbidity raster
├── Mackay-Capricorn_*_waves_Hs.tif                 # Wave height raster
├── Mackay-Capricorn_*_waves_Tp.tif                 # Wave period raster
├── FarNorthern_valid_slopes_lookup.parq            # Slope lookup table
├── FarNorthern_valid_slopes.tif                    # Valid extent mask
├── FarNorthern_*_bathy.tif                         # Depth raster
├── FarNorthern_*_slope.tif                         # Slope raster
├── FarNorthern_*_turbid.tif                        # Turbidity raster
├── FarNorthern_*_waves_Hs.tif                      # Wave height raster
└── FarNorthern_*_waves_Tp.tif                      # Wave period raster
```

**Total Files Required**: 30 files

- 1 global reef outlines file
- 4 regions × 7 files each (5 base criteria + lookup + extent)
- 1 additional rugosity file for Townsville

**Note**: The `*` in filenames represents potential additional naming components that may exist between the region ID and the criteria suffix, as the code uses glob pattern matching to find files.
