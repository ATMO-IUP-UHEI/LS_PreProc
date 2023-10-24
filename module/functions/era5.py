import numpy as np
import xarray as xr
import sys
import os
from datetime import datetime

import functions.constants as constants


def main(atm, era5_folder_path, nlevel, dims):
    print(datetime.now(), "getting era5 data")
    # get era5_data from the folder
    # name(s) of the required file will be constructed from the information in
    # atm
    # era5_data comes in multilayer (ml) and surface (sfc) files
    era5_ml = get_era5(atm, era5_folder_path, "ml")
    era5_sfc = get_era5(atm, era5_folder_path, "sfc")

    # drop unnecessary variables from era5_data and rename needed ones
    era5_ml = prepare_era5(era5_ml, "ml")
    era5_sfc = prepare_era5(era5_sfc, "sfc")

    # merge multilevel and surface datasets
    era5 = merge_era5(era5_ml, era5_sfc, dims)

    # calculate pressure grid on era5_levels
    era5 = calculate_era5_pressure(era5, era5_folder_path)

    # calculate geometric altitude
    era5 = calculate_era5_geometric_altitude(atm, era5)

    # calculate h2o mole fraction
    era5 = calculate_era5_h2o_mole_fraction(era5)

    # interpolate era5 to atm time and lat/lon grid
    era5 = interpolate_era5_onto_atm(atm, era5, dims)

    # correct each pixel for the error in the low-resolution surface elevation
    # using high resolution surface elevation from atm
    era5 = correct_surface_elevation(atm, era5)

    # interpolate era5 pressure grid onto new pressure grid
    import time
    start = time.time()
    era5 = interpolate_era5_onto_pressure_grid(atm, era5, nlevel, dims)
    end = time.time()
    print(f"inefficient part of era5: {end - start} s.")

    atm["pressure"] = xr.DataArray(
        data=era5.pressure.values,
        dims=(dims["y"], dims["x"], dims["z"])
    ).astype("float32")
    atm.pressure.attrs["source"] = "era5"

    atm["temperature"] = xr.DataArray(
        data=era5.temperature.values,
        dims=(dims["y"], dims["x"], dims["z"])
    ).astype("float32")
    atm.temperature.attrs["source"] = "era5"

    atm["geometric_altitude"] = xr.DataArray(
        data=era5.geometric_altitude.values,
        dims=(dims["y"], dims["x"], dims["z"])
    ).astype("float32")
    atm.geometric_altitude.attrs["source"] = \
        f"era5. surface elevation: {atm.surface_elevation.attrs['source']}"

    atm = atm.drop("surface_elevation")

    atm["h2o"] = xr.DataArray(
        data=era5.h2o.values,
        dims=(dims["y"], dims["x"], dims["z"])
    ).astype("float32")
    atm.h2o.attrs["source"] = "era5"

    return atm


def get_era5(atm, era5_folder_path, era5_type):
    file_name_list = generate_file_list_from_atm(atm, era5_type)

    era5_data_set = False
    for file_name in file_name_list:
        era5 = xr.open_dataset(
            f"{era5_folder_path}/{file_name}", engine="cfgrib"
        )

        if not era5_data_set:
            merged = era5
            era5_data_set = True

        era5 = xr.merge([era5, merged])
        merged.close()
    return era5


def generate_file_list_from_atm(atm, era5_type):
    yyyymmdd = "".join(str(min(atm.time).values)[:10].split("-"))

    file_name_list = [f"era5_{era5_type}_{yyyymmdd}.grb"]

    return file_name_list


def prepare_era5(era5, era5_type):
    if era5_type == "ml":
        era5 = era5.drop(["step", "valid_time"])
        era5 = era5.rename({"hybrid": "era5_level"})
        era5 = era5.drop_vars("era5_level")
        era5 = era5.rename({"t": "temperature", "q": "specific_humidity"})

        era5 = era5.assign(temperature=era5.temperature.astype(np.float32))
        era5 = era5.assign(
            specific_humidity=era5.specific_humidity.astype(np.float32)
        )

    if era5_type == "sfc":
        era5 = era5.drop(["number", "step", "surface", "valid_time"])
        era5 = era5.rename(
            {"z": "surface_geopotential", "sp": "surface_pressure"}
        )

        era5 = era5.assign(
            surface_geopotential=era5.surface_geopotential.astype(np.float32)
        )
        era5 = era5.assign(
            surface_pressure=era5.surface_pressure.astype(np.float32)
        )
        era5 = era5.assign(u10=era5.u10.astype(np.float32))
        era5 = era5.assign(v10=era5.v10.astype(np.float32))

    return era5


def merge_era5(era5_ml, era5_sfc, dims):
    # latitude and longitude grids are almost the same but not quite
    # this leads to problems when merging.
    # Don't want to interpolate onto atm grid yet, as it is very fine, meaning
    # we have to do a lot of calculations. Much better to do them first, and
    # then interpolate onto the fine grid.
    # Therefore, we will here interpolate on a dummy grid, which is almost like
    # the grids for the ml and sfc data

    min_latitude = max(min(era5_ml.latitude), min(era5_sfc.latitude))
    max_latitude = min(max(era5_ml.latitude), max(era5_sfc.latitude))
    Nlat = len(era5_ml.latitude.values)

    min_longitude = max(min(era5_ml.longitude), min(era5_sfc.longitude))
    max_longitude = min(max(era5_ml.longitude), max(era5_sfc.longitude))
    Nlon = len(era5_ml.longitude.values)

    new_latitude = np.linspace(min_latitude, max_latitude, Nlat)
    new_longitude = np.linspace(min_longitude, max_longitude, Nlon)

    era5_ml = era5_ml.interp(
        latitude=new_latitude,
        longitude=new_longitude
    )

    era5_sfc = era5_sfc.interp(
        latitude=new_latitude,
        longitude=new_longitude
    )

    era5 = xr.merge([era5_ml, era5_sfc])

    return era5


def calculate_era5_pressure(era5, era5_folder_path):
    # get standard ecmwd atmosphere parameters from file downloaded from
    # https://confluence.ecmwf.int/display/UDOC/L137+model+level+definitions
    try:
        auxiliary_path = os.path.join(era5_folder_path, "..", "auxiliary")
        n, a, b, ph, pf, gpa, gma, t, rho = np.genfromtxt(
            os.path.join(auxiliary_path, "standard_atmosphere.csv"),
            delimiter=",", skip_header=2, unpack=True)
    except FileNotFoundError:
        sys.exit("standard_atmosphere.csv missing. Download from "
                 + "https://confluence.ecmwf.int/display/UDOC/L137+model+"
                 + "level+definitions and put it into "
                 + "data/external/era5/auxiliary/")

    # pressure
    # pressure for layers is given at lower boundary of the model level
    # pressure calculated using a and b parameters and the following manual:
    # https://www.ecmwf.int/sites/default/files/elibrary/2015/9210-part-iii-dynamics-and-numerical-procedures.pdf
    pressure = a[None, :, None, None] + b[None, :, None, None] \
        * era5.surface_pressure.values[:, None, :, :]

    era5["pressure"] = xr.DataArray(
        data=pressure,
        dims=("time", "era5_level", "latitude", "longitude")
    )

    era5 = era5.drop("surface_pressure")

    return era5


def calculate_era5_geometric_altitude(atm, era5):

    # get surface elevation from surface geopotential
    surface_elevation = \
        era5.surface_geopotential / constants.gravitational_acceleration()

    geometric_altitude = \
        surface_elevation \
        + constants.scale_height \
        * np.log(era5.pressure[:, -1, :, :]/era5.pressure[:, :, :, :])

    geometric_altitude = geometric_altitude.transpose(
        "time", "era5_level", "latitude", "longitude"
    )

    era5["geometric_altitude"] = xr.DataArray(
        data=geometric_altitude,
        dims=("time", "era5_level", "latitude", "longitude")
    )

    era5 = era5.drop("surface_geopotential")

    return era5


def calculate_era5_h2o_mole_fraction(era5):
    # water vapor mole fraction
    # specific humidity is converted into h2o mole fraction using
    # s = m_H2O / m_HUM (s = specific humidity)
    # m_X = N_X * m_X
    # m_HUM = m_DRY + m_H2O
    # to find: h2o_mole_fraction = N_H2O / N_DRY = s/(1-s) * M_DRY/M_H2O
    h2o = era5.specific_humidity \
        / (1 - era5.specific_humidity) \
        * constants.molar_mass_dry_air \
        / constants.molar_mass_h2o

    era5["h2o"] = xr.DataArray(
        data=h2o,
        dims=("time", "era5_level", "latitude", "longitude")
    )

    era5 = era5.drop("specific_humidity")

    return era5


def interpolate_era5_onto_atm(atm, era5, dims):
    era5 = era5.interp(
        time=atm.time,
        latitude=atm.latitude,
        longitude=atm.longitude
    )

    era5 = era5.transpose(dims["y"], dims["x"], "era5_level")

    return era5


def correct_surface_elevation(atm, era5):
    above_ground = era5.geometric_altitude >= atm.surface_elevation
    era5["pressure"] = era5.pressure.where(above_ground)
    era5["temperature"] = era5.temperature.where(above_ground)
    era5["geometric_altitude"] = era5.geometric_altitude.where(above_ground)
    era5["h2o"] = era5.h2o.where(above_ground)

    ground_level = \
        era5.geometric_altitude.argmin(dim="era5_level", skipna=True)

    elevation_difference = \
        era5.geometric_altitude[..., ground_level] - atm.surface_elevation

    scale_height = \
        era5.temperature[..., ground_level] \
        * constants.gas_constant \
        / constants.molar_mass_dry_air \
        / constants.gravitational_acceleration()

    era5.geometric_altitude[..., ground_level] = \
        atm.surface_elevation

    era5.pressure[..., ground_level] = \
        era5.pressure[..., ground_level] \
        * np.exp(elevation_difference / scale_height)

    era5.temperature[..., ground_level] = \
        era5.temperature[..., ground_level] \
        - constants.lapse_rate_lower_troposphere \
        * elevation_difference

    # era5.specific_humidity[..., ground_level] is assumed to have a constant
    # profile, no change necessary

    return era5


def interpolate_era5_onto_pressure_grid(atm, era5, nlevel, dims):
    nlevel = int(nlevel)

    old_level_grid = np.tile(
        era5.era5_level.values, (len(atm[dims["y"]]), len(atm[dims["x"]]), 1)
    )
    old_pressure_grid = era5.pressure.values

    new_min_pressure = \
        era5.pressure.min(dim="era5_level", skipna=True)
    new_max_pressure = \
        era5.pressure.max(dim="era5_level", skipna=True)

    new_pressure_grid = np.linspace(new_min_pressure, new_max_pressure, nlevel)
    new_pressure_grid = np.moveaxis(new_pressure_grid, [0, 1, 2], [2, 0, 1])

    # something like this would be cool but doesn't exist for 3D grids
    # new_level_grid = \
    #     np.interp(new_pressure_grid, old_pressure_grid, old_level_grid)
    Nx = len(atm[dims["x"]])
    Ny = len(atm[dims["y"]])
    Nz = nlevel

    new_level_grid = np.empty(shape=new_pressure_grid.shape)
    for y in range(Ny):
        for x in range(Nx):
            new_level_grid[y, x, :] = np.interp(
                new_pressure_grid[y, x, :],
                old_pressure_grid[y, x, :],
                old_level_grid[y, x, :]
            )

    new_level_grid = xr.DataArray(
        data=new_level_grid,
        dims=(dims["y"], dims["x"], dims["z"])
    ).astype("float32")

    # something like this would be cool but doesn't exist for 3D grids
    # era5 = era5.interp(era5_level=new_level_grid)
    # instead of doing a 3D interpolation, it tries to split the era5_level
    # dimension into 3 coordinates, not what we want.
    # Finding an efficient way to do this will save the bulk of the calculation
    # time, sadly. This is super slow.
    for variable in ["pressure", "temperature", "geometric_altitude", "h2o"]:
        variable_array = np.empty(shape=(Ny, Nx, Nz), dtype="float32")
        for y in range(Ny):
            for x in range(Nx):
                data = era5[variable].isel(frame=y, line=x)
                level = new_level_grid[y, x, :]
                variable_array[y, x, :] = data.interp(era5_level=level)

        era5 = era5.drop(variable)

        era5[variable] = xr.DataArray(
            data=variable_array,
            dims=(dims["y"], dims["x"], dims["z"])
        ).astype("float32")

        variable_array = None

    return era5
