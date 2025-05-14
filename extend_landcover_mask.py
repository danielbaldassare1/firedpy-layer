#!/usr/bin/env python3
"""
Extend the existing landcover mask to include polar regions.
This script reads the existing landcover mask NetCDF file, fills NaN values with 1,
and saves the result to a new NetCDF file.
"""

import xarray as xr
import numpy as np

def extend_landcover_mask(input_nc, output_nc):
    """
    Extend the existing landcover mask to include polar regions.
    
    Parameters:
    -----------
    input_nc : str
        Path to the input NetCDF file
    output_nc : str
        Path to the output NetCDF file
    """
    print(f"Reading input file: {input_nc}")
    
    # Read the input NetCDF file
    ds = xr.open_dataset(input_nc)
    
    # Get the landcover mask
    mask = ds["landcover_mask"].values
    
    # Find the valid latitude range
    valid_mask = ~np.isnan(mask)
    valid_rows = np.where(np.any(valid_mask, axis=1))[0]
    min_lat_idx = valid_rows.min()
    max_lat_idx = valid_rows.max()
    
    print(f"Valid latitude range: {ds.lat.values[min_lat_idx]} to {ds.lat.values[max_lat_idx]}")
    print(f"Valid latitude indices: {min_lat_idx} to {max_lat_idx} out of {len(ds.lat)}")
    
    # Fill NaN values with 1
    extended_mask = np.nan_to_num(mask, nan=1.0)
    
    print(f"Filled {np.isnan(mask).sum()} NaN values with 1")
    
    # Create a new dataset with the extended mask
    extended_ds = xr.Dataset(
        data_vars={
            "landcover_mask": (["lat", "lon"], extended_mask),
        },
        coords={
            "lat": ds.lat,
            "lon": ds.lon,
        },
        attrs=ds.attrs
    )
    
    # Update attributes
    extended_ds.attrs["polar_regions_filled"] = "True"
    extended_ds["landcover_mask"].attrs = ds["landcover_mask"].attrs
    extended_ds["landcover_mask"].attrs["polar_regions_filled"] = "True"
    
    # Save to NetCDF
    print(f"Saving to NetCDF: {output_nc}")
    extended_ds.to_netcdf(output_nc)
    
    print("Done!")

if __name__ == "__main__":
    # Define input and output paths
    input_nc = "modis/global_landcover_mask_0.5deg.nc"
    output_nc = "modis/global_landcover_mask_0.5deg_extended.nc"
    
    # Extend the landcover mask
    extend_landcover_mask(input_nc, output_nc)
