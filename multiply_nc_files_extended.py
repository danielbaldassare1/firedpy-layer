#!/usr/bin/env python3
"""
Script to multiply a landcover mask with NetCDF files.
This script multiplies the extended landcover mask file with each .nc file in the output directory
and saves the results in the firedpy_layer directory.
"""

import os
import xarray as xr
import glob
import numpy as np

def multiply_nc_files(mask_file, input_dir, output_dir):
    """
    Multiply each .nc file in input_dir with the mask_file and save to output_dir.
    Memory-efficient version that processes one variable at a time.
    
    Parameters:
    -----------
    mask_file : str
        Path to the mask NetCDF file
    input_dir : str
        Directory containing input NetCDF files
    output_dir : str
        Directory where output files will be saved
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Get list of .nc files in input directory
    input_files = glob.glob(os.path.join(input_dir, "*.nc"))
    print(f"Found {len(input_files)} .nc files in {input_dir}")
    
    # Process each input file
    for input_file in input_files:
        filename = os.path.basename(input_file)
        output_file = os.path.join(output_dir, f"masked_{filename}")
        
        print(f"Processing: {filename}")
        
        # Get information about the input file without loading all data
        with xr.open_dataset(input_file) as input_ds:
            # Create a new dataset with the same coordinates and attributes
            coords_dict = {name: coord for name, coord in input_ds.coords.items()}
            attrs_dict = input_ds.attrs.copy()
            var_names = list(input_ds.data_vars)
            var_attrs = {name: input_ds[name].attrs.copy() for name in var_names}
            var_dims = {name: input_ds[name].dims for name in var_names}
            
            # Close the dataset to free memory
            input_ds.close()
        
        # Process one variable at a time
        for var_name in var_names:
            print(f"  Processing variable: {var_name}")
            
            # Open only the specific variable from the input file
            with xr.open_dataset(input_file) as input_ds:
                var_data = input_ds[var_name]
                
                # Open the mask file
                with xr.open_dataset(mask_file) as mask_ds:
                    # Get the first variable from the mask dataset
                    mask_var_name = list(mask_ds.data_vars)[0]
                    mask_data = mask_ds[mask_var_name]
                    
                    # Multiply the variable with the mask
                    # Convert to numpy arrays for more efficient multiplication
                    result = var_data.values * mask_data.values
                    
                    # If this is the first variable, create the output file
                    if not os.path.exists(output_file):
                        # Create a new dataset with the same coordinates
                        result_ds = xr.Dataset(coords=coords_dict, attrs=attrs_dict)
                        
                        # Add the variable to the dataset
                        result_ds[var_name] = (var_dims[var_name], result)
                        result_ds[var_name].attrs = var_attrs[var_name]
                        
                        # Save to netCDF
                        result_ds.to_netcdf(output_file, mode='w')
                        result_ds.close()
                    else:
                        # Append to existing file
                        with xr.open_dataset(output_file) as result_ds:
                            # Create a new dataset with just this variable
                            temp_ds = xr.Dataset(coords=coords_dict)
                            temp_ds[var_name] = (var_dims[var_name], result)
                            temp_ds[var_name].attrs = var_attrs[var_name]
                            
                            # Save to a temporary file
                            temp_file = output_file + '.temp'
                            temp_ds.to_netcdf(temp_file)
                            
                            # Close datasets
                            temp_ds.close()
                            result_ds.close()
                        
                        # Merge the temporary file with the output file
                        with xr.open_dataset(output_file) as ds1, xr.open_dataset(temp_file) as ds2:
                            merged = xr.merge([ds1, ds2])
                            merged.to_netcdf(output_file + '.merged')
                            ds1.close()
                            ds2.close()
                            merged.close()
                        
                        # Replace the output file with the merged file
                        os.replace(output_file + '.merged', output_file)
                        os.remove(temp_file)
                    
                    # Free memory
                    del result
                    
                # Close datasets
                mask_ds.close()
                input_ds.close()
        
        print(f"Completed processing: {filename}")
    
    print("Processing complete!")

if __name__ == "__main__":
    # Define paths
    mask_file = "modis/global_landcover_mask_0.5deg_extended.nc"
    input_dir = "output"
    output_dir = "firedpy_layer_extended"
    
    # Run the multiplication
    multiply_nc_files(mask_file, input_dir, output_dir)
