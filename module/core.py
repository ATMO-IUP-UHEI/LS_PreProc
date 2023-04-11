import numpy as np
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
import sys

import functions.atm as atm
import functions.l1b as l1b
import functions.aster as aster
import functions.era5 as era5
import functions.cams as cams

import cfgrib

import time

# prepare atm file
atm_data = atm.prepare_atm()

# gather (lat, lon) grid from l1b file
# gather time information from l1b file
l1b_file = sys.argv[1]
l1b_name = l1b_file.split("/")[-1]
if not l1b_name.startswith("L1B_"):
    sys.exit("L1B file name has to start with L1B_")
l1b.import_l1b_grid_example(atm_data, l1b_file)

# gather surface elevation from aster file
# the spatial grid is interpolated onto the l1b spatial grid
atm_data = aster.main(atm_data, "data/external/aster/download")

# era5 has highest resolution for pressure grid (137 model layers)
# the spatial grid is interpolated onto the l1b spatial grid
# the vertical pressure grid is gathered from era5
# the temporal grid is interpolated onto the l1b temporal grid
atm_data = era5.main(atm_data, "data/external/era5/download")

# cams has vertical trace gas profiles
# the spatial grid is interpolated onto the l1b spatial grid
# the vertical pressure grid is interpolated onto the era5 vertical pressure grid
# the temporal grid is interpolated onto the l1b temporal grid
atm_data = cams.main(atm_data, "data/external/cams/download")

atm_name = "ATM_" + l1b_name.strip("L1B_")
print(f"finished atm_data = {atm_data}")
print(f"saving to {atm_name}")
atm_data.to_netcdf(f"data/output/{atm_name}", mode = "w", format = "NETCDF4")
