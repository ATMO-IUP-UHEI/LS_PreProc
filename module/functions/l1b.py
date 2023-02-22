import xarray as xr

def import_l1b_grid_example(atm_data, l1b_path):
    l1b_data = xr.open_dataset(l1b_path)

    atm_data["datetime"] = l1b_data.datetime.values
    atm_data.datetime.attrs = {}
    atm_data.datetime.attrs["standard_name"] = "date and time of measurement"

    nline = len(l1b_data.line)
    nsample = len(l1b_data.sample)

    atm_data["latitude"] = (("line", "sample"), l1b_data.latitude.values)
    atm_data.latitude.attrs["standard_name"] = "latitude"
    atm_data.latitude.attrs["units"] = "degrees north"

    atm_data["longitude"] = (("line", "sample"), l1b_data.longitude.values)
    atm_data.longitude.attrs["standard_name"] = "longitude"
    atm_data.longitude.attrs["units"] = "degrees east"

    l1b_data.close()
