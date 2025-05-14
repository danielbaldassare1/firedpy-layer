# FireDPy

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
