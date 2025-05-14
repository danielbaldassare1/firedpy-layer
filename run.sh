#!/bin/bash
python create_landcover_mask.py
python extend_landcover_mask.py
python global_grid_analysis.py
python multiply_nc_files_extended.py
python plot_nc_files.py
