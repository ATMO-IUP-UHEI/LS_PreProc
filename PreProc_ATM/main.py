import configparser
import sys
import os
import xarray as xr
from datetime import datetime
import time
from cftime import date2num

import functions.aster as aster
import functions.era5 as era5
import functions.egg4 as egg4


def main():
    config_file = f"{os.path.dirname(__file__)}/config/config.ini"
    if not os.path.isfile(config_file):
        sys.exit("no config file provided")
    config = configparser.ConfigParser()
    config.read(config_file)

    general = config["general"]
    sources = config["sources"]
    general["auxiliary_path"] = \
        f"{os.path.dirname(__file__)}/data"

    sources = generate_source_list(sources)
    dims = generate_dims()

    l1b = xr.open_dataset("DATA_IN/L1B_DATA.nc")
    atm = create_atm(l1b, dims)

    atm = get_data(atm, general, sources, dims)
    atm = set_attributes(atm, dims)
    write_data(atm)

    return


def generate_source_list(variables):
    sources = {}

    for var in variables:
        source = variables[var]

        if source not in sources:
            sources[source] = []

        sources[source].append(var)

    return sources


def generate_dims():
    dims = {
        "x": "line",
        "y": "frame",
        "z": "level"
    }
    return dims


def create_atm(l1b, dims):
    atm = xr.Dataset()

    # import time as np.datetime64 object so that interpolation can work.
    # later convert to float32
    atm["time"] = l1b.time
    atm["latitude"] = l1b.latitude
    atm["longitude"] = l1b.longitude

    return atm


def get_data(atm, general, sources, dims):
    if "aster" in sources:
        # TODO:
        # make sure this only gets the variables that are actually specified
        atm = aster.main(
            atm,
            general["aster_path"],
            dims
        )
    if "era5" in sources:
        # TODO:
        # make sure this only gets the variables that are actually specified
        atm = era5.main(
            atm,
            general["auxiliary_path"],
            general["era5_path"],
            general["nlevel"],
            dims
        )
    if "egg4" in sources:
        # TODO:
        # make sure this only gets the variables that are actually specified
        atm = egg4.main(
            atm,
            general["egg4_path"],
            dims
        )

    return atm


def set_attributes(atm, dims):
    # convert np.datetime64 into float32
    time, time_units_string = get_time(atm)
    atm["time"] = xr.DataArray(
        data=time,
        dims=(dims["y"]),
    ).astype("float32")
    atm.time.attrs["long_name"] = "date and time in UTC"
    atm.time.attrs["units"] = time_units_string

    assert "time" in atm.data_vars
    assert atm.time.dims == (dims["y"],)
    assert atm.time.dtype == "float32"
    assert atm.time.attrs["units"].startswith("seconds since ")
    atm.time.attrs["long_name"] = "date and time in UTC"

    assert "latitude" in atm.data_vars
    assert atm.latitude.dims == (dims["y"], dims["x"])
    assert atm.latitude.dtype == "float32"
    atm.latitude.attrs["long_name"] = "latitude at pixel center"
    atm.latitude.attrs["units"] = "degrees north"

    assert "longitude" in atm.data_vars
    assert atm.longitude.dims == (dims["y"], dims["x"])
    assert atm.longitude.dtype == "float32"
    atm.longitude.attrs["long_name"] = "longitude at pixel center"
    atm.longitude.attrs["units"] = "degrees east"

    assert "pressure" in atm.data_vars
    assert atm.pressure.dims == (dims["y"], dims["x"], dims["z"])
    assert atm.pressure.dtype == "float32"
    atm.pressure.attrs["long_name"] = "pressure"
    atm.pressure.attrs["units"] = "Pa"

    assert "temperature" in atm.data_vars
    assert atm.temperature.dims == (dims["y"], dims["x"], dims["z"])
    assert atm.temperature.dtype == "float32"
    atm.temperature.attrs["long_name"] = "temperature"
    atm.temperature.attrs["units"] = "K"

    assert "geometric_altitude" in atm.data_vars
    assert atm.geometric_altitude.dims == (dims["y"], dims["x"], dims["z"])
    assert atm.geometric_altitude.dtype == "float32"
    atm.geometric_altitude.attrs["long_name"] = \
        "geometric altitude above sea level"
    atm.geometric_altitude.attrs["units"] = "m"

    assert "h2o" in atm.data_vars
    assert atm.h2o.dims == (dims["y"], dims["x"], dims["z"])
    assert atm.h2o.dtype == "float32"
    atm.h2o.attrs["long_name"] = "H2O mole fraction"
    atm.h2o.attrs["units"] = "mol mol-1"

    assert "co2" in atm.data_vars
    assert atm.co2.dims == (dims["y"], dims["x"], dims["z"])
    assert atm.co2.dtype == "float32"
    atm.co2.attrs["long_name"] = "CO2 mole fraction"
    atm.co2.attrs["units"] = "mol mol-1"

    assert "ch4" in atm.data_vars
    assert atm.ch4.dims == (dims["y"], dims["x"], dims["z"])
    assert atm.ch4.dtype == "float32"
    atm.ch4.attrs["long_name"] = "CH4 mole fraction"
    atm.ch4.attrs["units"] = "mol mol-1"

    assert "wind_u_component" in atm.data_vars
    assert atm.wind_u_component.dims == (dims["y"], dims["x"], dims["z"])
    assert atm.wind_u_component.dtype == "float32"
    atm.wind_u_component.attrs["long_name"] = "u-component of wind speed"
    atm.wind_u_component.attrs["units"] = "m s-1"

    assert "wind_v_component" in atm.data_vars
    assert atm.wind_v_component.dims == (dims["y"], dims["x"], dims["z"])
    assert atm.wind_v_component.dtype == "float32"
    atm.wind_v_component.attrs["long_name"] = "v-component of wind speed"
    atm.wind_v_component.attrs["units"] = "m s-1"

    return atm


def get_time(atm):
    # convert numpy.datetime64 into datetime object and then into string
    reference_time = str(min(atm.time.values))[:-3]
    time_units_string = f"seconds since {reference_time}"

    time = [
        datetime.strptime(str(t)[:-3], "%Y-%m-%dT%H:%M:%S.%f")
        for t in atm.time.values
    ]

    time = date2num(time, time_units_string)

    return time, time_units_string


def write_data(atm):
    output = "DATA_IN/ATM_DATA.nc"

    history_string = \
        "Created using the ATM preprocessor for RemoTeC."

    atm.attrs["history"] = history_string

    for var in atm.data_vars:
        atm[var].encoding.update({"_FillValue": None})

    atm.to_netcdf(
        output,
        mode="w",
        format="NETCDF4"
    )


if __name__ == "__main__":
    start_time = time.time()
    main()
    stop_time = time.time()
    print(f"Total runtime: {stop_time - start_time:.2f} s.")
