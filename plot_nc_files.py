#!/usr/bin/env python3
"""
Script to create global plots from NetCDF files in the firedpy_layer directory.
This script loads each .nc file, creates a global map plot, and saves it as a PNG image.
"""

import os
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as ticker
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from pathlib import Path

def plot_nc_file(nc_file, output_dir):
    """
    Create a global plot from a NetCDF file.
    
    Parameters:
    -----------
    nc_file : str
        Path to the NetCDF file
    output_dir : str
        Directory where output plots will be saved
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Get filename without extension for plot title and output filename
    filename = os.path.basename(nc_file)
    name_without_ext = os.path.splitext(filename)[0]
    
    print(f"Creating plot for: {filename}")
    
    # Load the NetCDF file
    ds = xr.open_dataset(nc_file)
    
    # Get the first data variable
    var_name = list(ds.data_vars)[0]
    data = ds[var_name]
    
    # Get the dimensions
    dims = data.dims
    
    # Check if we have latitude and longitude dimensions
    lat_dim = next((dim for dim in dims if 'lat' in dim.lower()), None)
    lon_dim = next((dim for dim in dims if 'lon' in dim.lower()), None)
    
    if lat_dim is None or lon_dim is None:
        print(f"Warning: Could not identify lat/lon dimensions in {filename}")
        print(f"Available dimensions: {dims}")
        ds.close()
        return
    
    # Get latitude and longitude values
    lats = ds[lat_dim].values
    lons = ds[lon_dim].values
    
    # Create a figure with a map projection that better shows polar regions
    plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=ccrs.Mollweide())
    
    # Add coastlines, borders, and grid lines
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    
    # Add gridlines without any labels
    gl = ax.gridlines(draw_labels=False, linewidth=0.5, color='gray', alpha=0.5, linestyle='--',
                     xlocs=range(-180, 181, 60), ylocs=range(-90, 91, 30))
    
    # Set global extent
    ax.set_global()
    
    # Add a note about the projection
    plt.figtext(0.01, 0.01, "Mollweide projection", fontsize=8, color='gray')
    
    # Create a colormap based on the data
    data_values = data.values
    data_min = np.nanmin(data_values)
    data_max = np.nanmax(data_values)
    
    # Use white to dark red color scheme for all plots
    if data_max > 0 and data_min >= 0 and data_max / max(data_min, 1e-10) > 100:
        # Use log scale for data with large range
        # Adjust the minimum value to ensure good contrast
        vmin = max(data_min, data_max / 1000)  # Ensure at least 3 orders of magnitude range
        norm = colors.LogNorm(vmin=vmin, vmax=data_max)
    else:
        # Use linear scale for data with smaller range
        norm = colors.Normalize(vmin=data_min, vmax=data_max)
    
    # Use white to dark red color scheme
    cmap = plt.colormaps['Reds']
    
    # Plot the data
    # For 2D data (lat, lon)
    if len(data.shape) == 2:
        im = ax.pcolormesh(lons, lats, data.values, transform=ccrs.PlateCarree(), 
                          cmap=cmap, norm=norm, shading='auto')
    # For 3D data (time, lat, lon) - plot the first time step
    elif len(data.shape) == 3:
        im = ax.pcolormesh(lons, lats, data.values[0], transform=ccrs.PlateCarree(), 
                          cmap=cmap, norm=norm, shading='auto')
        print(f"Note: Plotting first time step only for {filename}")
    else:
        print(f"Warning: Unexpected data shape {data.shape} in {filename}")
        ds.close()
        return
    
    # Add colorbar with more ticks for better reference
    if isinstance(norm, colors.LogNorm):
        # For logarithmic scale, use logarithmically spaced ticks
        tick_locator = ticker.LogLocator(base=10)
        cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, shrink=0.8, 
                           ticks=tick_locator)
    else:
        # For linear scale, use more ticks than default
        cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, shrink=0.8)
    
    # Set publication-quality colorbar label
    if 'fires' in var_name.lower():
        cbar.set_label("Fire frequency (fires/year)")
    elif 'burned' in var_name.lower():
        cbar.set_label("Burned area (km²/year)")
    else:
        cbar.set_label(f"{var_name} ({data.units if hasattr(data, 'units') else ''})")
    
    # Add precise, scientific title
    if 'fires' in var_name.lower():
        plt.title(f"Annual Fire Frequency (2000-2021)")
    elif 'burned' in var_name.lower():
        plt.title(f"Annual Burned Area (2000-2021)")
    else:
        plt.title(f"{var_name}")
    
    # Save the plot
    output_file = os.path.join(output_dir, f"{name_without_ext}.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved plot to: {output_file}")
    
    # Close the figure and dataset
    plt.close()
    ds.close()

def plot_all_nc_files(input_dir, output_dir):
    """
    Create plots for all NetCDF files in the input directory.
    
    Parameters:
    -----------
    input_dir : str
        Directory containing NetCDF files
    output_dir : str
        Directory where output plots will be saved
    """
    # Get all .nc files in the input directory
    nc_files = list(Path(input_dir).glob("*.nc"))
    
    if not nc_files:
        print(f"No .nc files found in {input_dir}")
        return
    
    print(f"Found {len(nc_files)} .nc files in {input_dir}")
    
    # Process each file
    for nc_file in nc_files:
        plot_nc_file(str(nc_file), output_dir)
    
    print("All plots created successfully!")

if __name__ == "__main__":
    # Define paths
    input_dir = "firedpy_layer_extended"
    output_dir = "firedpy_layer_extended/plots"
    
    # Create plots for all .nc files
    plot_all_nc_files(input_dir, output_dir)
