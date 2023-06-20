import os

def prepare_atm():
    import xarray as xr

    atm_data = xr.Dataset()

    return atm_data



def write_output(atm_data, scenario_name):
    atm_path = os.path.join("SYNTH_SPECTRA", f"ATM_{scenario_name}")

    # Don't write _FillValue = NaN into netcdf file
    for var in atm_data.data_vars:
        atm_data[var].encoding.update({"_FillValue": None})

    atm_data.to_netcdf(atm_path, mode="w", format="NETCDF4")
