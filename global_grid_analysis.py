#!/usr/bin/env python3
"""
Global Grid Analysis for Fire Data

This script:
1. Reads all .gpkg files in the data directory
2. Creates a 0.5 degree global grid
3. Calculates the total number of fires per grid cell
4. Calculates the burned area per grid cell divided by 22 years
5. Saves the results to GeoTIFF files
"""

import os
import sys
import glob
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import rasterio
from rasterio.transform import from_origin
import xarray as xr
import warnings
warnings.filterwarnings('ignore')

# Constants
CELL_SIZE = 0.5  # Grid cell size in degrees
YEARS = 22  # Number of years (2000-2021)
GLOBAL_BOUNDS = (-180, -90, 180, 90)  # (xmin, ymin, xmax, ymax)

# Directories
DATA_DIR = "data"
OUTPUT_DIR = "output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def create_global_grid():
    """Create a global grid with 0.5 degree cell size."""
    xmin, ymin, xmax, ymax = GLOBAL_BOUNDS
    
    # Calculate the number of cells in each dimension
    nx = int((xmax - xmin) / CELL_SIZE)
    ny = int((ymax - ymin) / CELL_SIZE)
    
    print(f"Creating global grid with dimensions: {nx} x {ny} cells")
    
    # Initialize arrays for fire counts and burned area
    fire_counts = np.zeros((ny, nx), dtype=np.int32)
    burned_area = np.zeros((ny, nx), dtype=np.float32)
    
    return fire_counts, burned_area

def process_gpkg_files(data_dir, fire_counts, burned_area=None):
    """Process all .gpkg files in the data directory and update the global grid arrays."""
    # Get all .gpkg files
    gpkg_files = glob.glob(os.path.join(data_dir, "*.gpkg"))
    
    if not gpkg_files:
        print(f"No .gpkg files found in {data_dir}")
        return
    
    print(f"Found {len(gpkg_files)} .gpkg files to process")
    
    # Process each file
    for file_idx, file_path in enumerate(gpkg_files):
        print(f"Processing file {file_idx+1}/{len(gpkg_files)}: {os.path.basename(file_path)}")
        
        try:
            # Read the GeoPackage file
            gdf = gpd.read_file(file_path)
            print(f"  Loaded file with {len(gdf)} events")
            
            # Reproject to WGS84 if needed
            if gdf.crs and gdf.crs != "EPSG:4326":
                print(f"  Reprojecting from {gdf.crs} to EPSG:4326...")
                gdf = gdf.to_crs("EPSG:4326")
            
            # Process each fire event
            print(f"  Assigning events to grid cells...")
            processed = 0
            
            for idx, row in gdf.iterrows():
                processed += 1
                if processed % 100000 == 0:
                    print(f"    Processed {processed}/{len(gdf)} events...")
                
                try:
                    # Get the centroid of the geometry
                    geom = row.geometry
                    if geom is None or geom.is_empty:
                        continue
                    
                    centroid = geom.centroid
                    x, y = centroid.x, centroid.y
                    
                    # Skip if outside global bounds
                    if (x < GLOBAL_BOUNDS[0] or x >= GLOBAL_BOUNDS[2] or 
                        y < GLOBAL_BOUNDS[1] or y >= GLOBAL_BOUNDS[3]):
                        continue
                    
                    # Calculate the grid indices
                    i = int((x - GLOBAL_BOUNDS[0]) / CELL_SIZE)
                    j = int((y - GLOBAL_BOUNDS[1]) / CELL_SIZE)
                    
                    # Ensure indices are within bounds
                    if (0 <= i < fire_counts.shape[1] and 0 <= j < fire_counts.shape[0]):
                        # Update fire count
                        fire_counts[j, i] += 1
                        
                        # Update burned area if available
                        if 'tot_ar_km2' in row and not np.isnan(row['tot_ar_km2']):
                            burned_area[j, i] += row['tot_ar_km2']
                
                except Exception as e:
                    # Skip problematic geometries but continue processing
                    continue
            
            print(f"  Completed processing {processed} events from {os.path.basename(file_path)}")
            
        except Exception as e:
            print(f"  Error processing file {file_path}: {e}")
            continue
    
    return fire_counts, burned_area

def save_results(fire_counts, burned_area):
    """Save the results to NetCDF files and create visualizations."""
    xmin, ymin, xmax, ymax = GLOBAL_BOUNDS
    
    # Calculate annual average burned area
    annual_burned_area = burned_area / YEARS
    
    # Create coordinate arrays
    ny, nx = fire_counts.shape
    lons = np.linspace(xmin + CELL_SIZE/2, xmax - CELL_SIZE/2, nx)
    lats = np.linspace(ymin + CELL_SIZE/2, ymax - CELL_SIZE/2, ny)
    
    # Create xarray DataArrays
    # 1. Fires per year
    fires_per_year = fire_counts / YEARS  # Convert total fires to annual average
    fires_da = xr.DataArray(
        fires_per_year,
        coords=[('latitude', lats), ('longitude', lons)],
        attrs={
            'long_name': 'Annual average number of fires',
            'units': 'fires/year',
            'description': f'Average number of fires per year over {YEARS} years (2000-2021)'
        }
    )
    
    # 2. Burned area per year
    burned_area_da = xr.DataArray(
        annual_burned_area,
        coords=[('latitude', lats), ('longitude', lons)],
        attrs={
            'long_name': 'Annual burned area',
            'units': 'km²/year',
            'description': f'Average burned area per year over {YEARS} years (2000-2021)'
        }
    )
    
    # Create datasets
    fires_ds = xr.Dataset({'fires_per_year': fires_da})
    burned_area_ds = xr.Dataset({'burned_area_per_year': burned_area_da})
    
    # Add global attributes
    for ds in [fires_ds, burned_area_ds]:
        ds.attrs['title'] = 'Global Fire Analysis'
        ds.attrs['source'] = 'FireDPy analysis'
        ds.attrs['references'] = 'Generated from MODIS fire data'
        ds.attrs['cell_size_degrees'] = CELL_SIZE
        ds.attrs['time_period'] = '2000-2021'
        ds.attrs['created'] = np.datetime64('now').astype(str)
    
    # Save to NetCDF files
    fires_nc_file = os.path.join(OUTPUT_DIR, "fires_per_year.nc")
    burned_area_nc_file = os.path.join(OUTPUT_DIR, "burned_area_per_year.nc")
    
    fires_ds.to_netcdf(fires_nc_file)
    burned_area_ds.to_netcdf(burned_area_nc_file)
    
    print(f"Saved fires per year to {fires_nc_file}")
    print(f"Saved burned area per year to {burned_area_nc_file}")
    
    # Also save as GeoTIFF for backward compatibility
    transform = from_origin(xmin, ymax, CELL_SIZE, CELL_SIZE)
    
    # Save fire counts as GeoTIFF
    fire_counts_file = os.path.join(OUTPUT_DIR, "fire_counts.tif")
    with rasterio.open(
        fire_counts_file,
        'w',
        driver='GTiff',
        height=fire_counts.shape[0],
        width=fire_counts.shape[1],
        count=1,
        dtype=fire_counts.dtype,
        crs='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs',
        transform=transform,
    ) as dst:
        dst.write(fire_counts, 1)
    
    # Save annual burned area as GeoTIFF
    annual_ba_file = os.path.join(OUTPUT_DIR, "annual_burned_area.tif")
    with rasterio.open(
        annual_ba_file,
        'w',
        driver='GTiff',
        height=annual_burned_area.shape[0],
        width=annual_burned_area.shape[1],
        count=1,
        dtype=annual_burned_area.dtype,
        crs='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs',
        transform=transform,
    ) as dst:
        dst.write(annual_burned_area, 1)
    
    # Create global visualizations
    plt.figure(figsize=(15, 8))
    plt.imshow(fire_counts, origin='upper', cmap='viridis', 
               extent=[xmin, xmax, ymin, ymax])
    plt.colorbar(label='Number of Fires')
    plt.title('Total Number of Fires (2000-2021)')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "fire_counts_map.png"), dpi=300)
    
    plt.figure(figsize=(15, 8))
    plt.imshow(annual_burned_area, origin='upper', cmap='hot_r', 
               extent=[xmin, xmax, ymin, ymax])
    plt.colorbar(label='Annual Burned Area (km²/year)')
    plt.title('Annual Burned Area (km²/year)')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "annual_burned_area_map.png"), dpi=300)
    
    print("Saved global visualizations")

def main():
    # Create global grid
    fire_counts, burned_area = create_global_grid()
    
    # Process all .gpkg files
    process_gpkg_files(DATA_DIR, fire_counts, burned_area)
    
    # Save results
    save_results(fire_counts, burned_area)
    
    # Print summary statistics
    total_fires = np.sum(fire_counts)
    total_burned_area = np.sum(burned_area)
    annual_burned_area = total_burned_area / YEARS
    nonzero_cells = np.count_nonzero(fire_counts)
    
    print("\n--- Summary Statistics ---")
    print(f"Total number of fires: {total_fires}")
    print(f"Total burned area: {total_burned_area:.2f} km²")
    print(f"Annual burned area: {annual_burned_area:.2f} km²/year")
    print(f"Number of 0.5-degree cells with fire events: {nonzero_cells}")
    print(f"Percentage of global cells with fire events: {nonzero_cells / (fire_counts.shape[0] * fire_counts.shape[1]) * 100:.2f}%")
    
    print("\nAnalysis complete.")

if __name__ == "__main__":
    main()
