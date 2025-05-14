# firedpy-layer

firedpy-layer is a Python-based tool for analyzing global fire data. It processes fire event data from GeoPackage files, creates a global grid of fire frequency and burned area, and applies a land cover mask to exclude urban areas and croplands from the analysis.

## Overview

firedpy-layer performs the following operations:

1. Creates a land cover mask from PROBAV land cover data, marking urban areas and croplands as 0 and other land cover types as 1
2. Extends the land cover mask to include polar regions
3. Processes fire data from GeoPackage files to calculate fire counts and burned area per grid cell
4. Applies the land cover mask to the fire data to exclude urban areas and croplands
5. Creates visualizations of the masked fire data

## Requirements

- Python 3.x
- Required Python packages:
  - numpy
  - xarray
  - rasterio
  - geopandas
  - matplotlib
  - cartopy
  - netCDF4

## Data Requirements

1. **Fire Data**: GeoPackage (.gpkg) files containing fire event data with geometries and burned area information
2. **Land Cover Data**: PROBAV land cover data in GeoTIFF format

## Directory Structure

- `data/`: Directory for input GeoPackage files containing fire event data
- `modis/`: Directory for PROBAV land cover data
- `output/`: Directory where intermediate NetCDF files are saved
- `firedpy_layer_extended/`: Directory where masked fire data is saved
- `firedpy_layer_extended/plots/`: Directory where visualizations are saved

## Usage

### Data Preparation

1. Place your GeoPackage (.gpkg) files containing fire event data in the `data/` directory
2. Place the PROBAV land cover GeoTIFF file in the `modis/` directory with the name `PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif`

### Running the Analysis

Run the following scripts in sequence:

```bash
python create_landcover_mask.py
python extend_landcover_mask.py
python global_grid_analysis.py
python multiply_nc_files_extended.py
python plot_nc_files.py
```

Alternatively, you can use the provided shell script:

```bash
./run.sh
```

### Output

The analysis produces the following outputs:

1. **NetCDF Files**:
   - `output/fires_per_year.nc`: Annual average number of fires per grid cell
   - `output/burned_area_per_year.nc`: Annual average burned area per grid cell
   - `firedpy_layer_extended/masked_fires_per_year.nc`: Masked annual average number of fires (excluding urban areas and croplands)
   - `firedpy_layer_extended/masked_burned_area_per_year.nc`: Masked annual average burned area (excluding urban areas and croplands)

2. **GeoTIFF Files**:
   - `output/fire_counts.tif`: Total number of fires per grid cell
   - `output/annual_burned_area.tif`: Annual average burned area per grid cell

3. **Visualizations**:
   - `output/fire_counts_map.png`: Global map of total fire counts
   - `output/annual_burned_area_map.png`: Global map of annual burned area
   - `firedpy_layer_extended/plots/masked_fires_per_year.png`: Global map of masked annual fire frequency
   - `firedpy_layer_extended/plots/masked_burned_area_per_year.png`: Global map of masked annual burned area

## Notes

- The analysis uses a 0.5-degree global grid
- The time period covered is 2000-2021 (22 years)
- Urban areas and croplands (values 40 and 50 in the PROBAV land cover data) are excluded from the final analysis
