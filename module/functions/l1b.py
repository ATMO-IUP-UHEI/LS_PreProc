import xarray as xr

def import_l1b_grid_example(atm_data, l1b_path):
    l1b_data = xr.open_dataset(l1b_path)

    atm_data["datetime"] = l1b_data.datetime
    atm_data["latitude"] = l1b_data.latitude
    atm_data["longitude"] = l1b_data.longitude

    l1b_data.close()
################################################################################
