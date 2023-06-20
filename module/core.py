import numpy as np
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
import sys
import os
import cfgrib
import time

import functions.atm as atm
import functions.l1b as l1b
import functions.aster as aster
import functions.era5 as era5
import functions.cams as cams

def main():
    # prepare atm file
    atm_data = atm.prepare_atm()

    # l1b_file should be something like "SYNTH_SPECTRA/L1B_scenario.nc"
    l1b_file = sys.argv[1]
    l1b_name = l1b_file.split("/")[-1]
    if not l1b_name.startswith("L1B_"):
        sys.exit("L1B file name needs to start with L1B_")
    prefix, scenario_name = l1b_name[:3], l1b_name[4:]

    # gather (lat, lon) grid from l1b file
    # gather time information from l1b file
    l1b.import_l1b_grid_example(atm_data, l1b_file)

    # path at which external raw data is stored, so that it can be imported here
    data_path = os.path.join(
        "/", "Users", "leonie", "phd", "RemoTeC", "preprocess", "ATM_preproc", "module", "data", "external"
    )

    # gather surface elevation from aster file
    # the spatial grid is interpolated onto the l1b spatial grid
    aster_path = os.path.join(data_path, "aster", "download")
    if not os.path.exists(aster_path):
        os.makedirs(aster_path)
    atm_data = aster.main(atm_data, aster_path)

    # era5 has highest resolution for pressure grid (137 model layers)
    # the spatial grid is interpolated onto the l1b spatial grid
    # the vertical pressure grid is gathered from era5
    # the temporal grid is interpolated onto the l1b temporal grid
    era5_path = os.path.join(data_path, "era5", "download")
    if not os.path.exists(era5_path):
        os.makedirs(era5_path)
    atm_data = era5.main(atm_data, era5_path)

    # cams has vertical trace gas profiles
    # the spatial grid is interpolated onto the l1b spatial grid
    # the vertical pressure grid is interpolated onto the era5 vertical pressure grid
    # the temporal grid is interpolated onto the l1b temporal grid
    cams_path = os.path.join(data_path, "cams", "download")
    if not os.path.exists(cams_path):
        os.makedirs(cams_path)
    atm_data = cams.main(atm_data, cams_path)

    print(f"finished atm_data = {atm_data}")
    atm.write_output(atm_data, scenario_name)

if __name__ == "__main__":
    main()
