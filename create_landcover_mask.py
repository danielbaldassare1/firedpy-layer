#!/usr/bin/env python3
"""
Create a global land cover mask NetCDF file from PROBAV land cover data.
The output will have a resolution of 0.5 degrees, where:
- 0 = urban land cover or croplands (values 40 and 50 in the original data)
- 1 = all other valid land cover types
"""

import os
import numpy as np
import rasterio
import xarray as xr
from rasterio.transform import rowcol

def create_landcover_mask(input_tif, output_nc, resolution=0.5):
    """
    Create a global land cover mask NetCDF file from PROBAV land cover data.
    Uses undersampling approach for efficiency.
    
    Parameters:
    -----------
    input_tif : str
        Path to the input TIFF file
    output_nc : str
        Path to the output NetCDF file
    resolution : float
        Resolution of the output grid in degrees (default: 0.5)
    """
    print(f"Reading input file: {input_tif}")
    
    # Read the input TIFF file
    with rasterio.open(input_tif) as src:
        # Get the NoData value
        nodata = src.nodata
        print(f"NoData value: {nodata}")
        
        # Define global grid bounds
        xmin, ymin, xmax, ymax = -180.0, -90.0, 180.0, 90.0
        
        # Calculate the number of pixels in each dimension
        width = int((xmax - xmin) / resolution)
        height = int((ymax - ymin) / resolution)
        
        print(f"Output dimensions: {width} x {height} pixels")
        
        # Create latitude and longitude arrays for the output grid
        lons = np.linspace(xmin + resolution/2, xmax - resolution/2, width)
        lats = np.linspace(ymax - resolution/2, ymin + resolution/2, height)
        
        # Create an empty array for the output data
        mask = np.ones((height, width), dtype=np.float32)
        
        # Create a meshgrid of the output coordinates
        lon_grid, lat_grid = np.meshgrid(lons, lats)
        
        print("Sampling points from the original raster...")
        
        # Convert the output grid coordinates to pixel coordinates in the input raster
        rows, cols = rowcol(src.transform, lon_grid.flatten(), lat_grid.flatten())
        
        # Create a mask for valid pixel coordinates
        valid_pixels = (
            (rows >= 0) & (rows < src.height) & 
            (cols >= 0) & (cols < src.width)
        )
        
        # Read the values at the valid pixel coordinates
        values = np.full(len(rows), nodata, dtype=np.uint8)
        if np.any(valid_pixels):
            # Only read values for valid pixel coordinates
            valid_rows = rows[valid_pixels]
            valid_cols = cols[valid_pixels]
            
            # Read the values in chunks to avoid memory issues
            chunk_size = 1000000  # Adjust based on available memory
            for i in range(0, len(valid_rows), chunk_size):
                chunk_rows = valid_rows[i:i+chunk_size]
                chunk_cols = valid_cols[i:i+chunk_size]
                
                # Read the values for this chunk
                chunk_values = src.read(1, window=((min(chunk_rows), max(chunk_rows) + 1),
                                                  (min(chunk_cols), max(chunk_cols) + 1)))
                
                # Extract the values at the specific coordinates
                for j, (r, c) in enumerate(zip(chunk_rows - min(chunk_rows), 
                                              chunk_cols - min(chunk_cols))):
                    if 0 <= r < chunk_values.shape[0] and 0 <= c < chunk_values.shape[1]:
                        values[valid_pixels.nonzero()[0][i+j]] = chunk_values[r, c]
        
        # Reshape the values to match the output grid
        values = values.reshape(height, width)
        
        print("Creating mask...")
        
        # Set urban and cropland areas to 0
        mask[(values == 40) | (values == 50)] = 0
        
        # Set NoData areas to NaN
        mask[values == nodata] = np.nan
        
        # Create an xarray Dataset
        print("Creating xarray Dataset...")
        ds = xr.Dataset(
            data_vars={
                "landcover_mask": (["lat", "lon"], mask),
            },
            coords={
                "lat": lats,
                "lon": lons,
            },
            attrs={
                "description": "Global land cover mask at 0.5 degree resolution",
                "source": os.path.basename(input_tif),
                "legend": "0 = urban land cover or croplands, 1 = all other valid land cover types",
                "method": "Undersampling from original PROBAV land cover data",
            }
        )
        
        # Add variable attributes
        ds["landcover_mask"].attrs = {
            "long_name": "Land cover mask",
            "units": "1",
            "description": "0 = urban land cover or croplands, 1 = all other valid land cover types",
        }
        
        # Save to NetCDF
        print(f"Saving to NetCDF: {output_nc}")
        ds.to_netcdf(output_nc)
        
        print("Done!")

if __name__ == "__main__":
    # Define input and output paths
    input_tif = "modis/PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif"
    output_nc = "modis/global_landcover_mask_0.5deg.nc"
    
    # Create the landcover mask
    create_landcover_mask(input_tif, output_nc, resolution=0.5)
