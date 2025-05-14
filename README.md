# FireDPy

This package produces global geospatial layers of burned area and number of fires per year at a 0.5 degree resolution.
To run this analysis, download firedpy data and MODIS data, run scripts sequentially.

1. Put GPKG files in `data/`
2. Put PROBAV tif in `modis/`
3. Run:
```
python create_landcover_mask.py
python extend_landcover_mask.py
python global_grid_analysis.py
python multiply_nc_files_extended.py
python plot_nc_files.py
```
