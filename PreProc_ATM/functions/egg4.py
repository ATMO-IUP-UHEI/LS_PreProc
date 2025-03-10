import xarray as xr
import numpy as np
from datetime import datetime

import functions.constants as constants


def main(atm, egg4_folder_path, dims):
    print(datetime.now(), "getting egg4 data")
    # get egg4 from the folder.
    # name(s) of the required file will be constructed from the information in
    # atm
    egg4 = get_egg4(atm, egg4_folder_path)

    # drop unnecessary variables from egg4 and rename needed ones
    egg4 = prepare_egg4(egg4)

    # if no data was available for desired year, previous years can be
    # downloaded. This script handles this case as an input.
    egg4, year_differs, used_year = correct_year(atm, egg4)

    # calculate egg4 pressure
    egg4 = calculate_egg4_pressure(egg4)

    # calculate co2 mole fraction
    egg4 = calculate_egg4_co2_mole_fraction(egg4)

    # calculate ch4 mole fraction
    egg4 = calculate_egg4_ch4_mole_fraction(egg4)

    # interpolate egg4 to atm time and lat/lon grid
    egg4 = interpolate_egg4_onto_atm(atm, egg4, dims)

    # interpolate egg4 to atm pressure grid
    egg4 = interpolate_egg4_onto_pressure_grid(atm, egg4, dims)

    atm["co2"] = xr.DataArray(
        data=egg4.co2.values,
        dims=(dims["y"], dims["x"], dims["z"])
    ).astype("float32")
    atm.co2.attrs["source"] = "egg4"
    if year_differs:
        atm.co2.attrs["source"] += f" for year {used_year}"

    atm["ch4"] = xr.DataArray(
        data=egg4.ch4.values,
        dims=(dims["y"], dims["x"], dims["z"])
    ).astype("float32")
    atm.ch4.attrs["source"] = "egg4"
    if year_differs:
        atm.ch4.attrs["source"] += f" for year {used_year}"

    return atm


def get_egg4(atm, egg4_folder_path):
    file_name_list = generate_file_list_from_atm_data(atm)

    egg4_is_set = False
    for file_name in file_name_list:
        data = xr.open_dataset(
            f"{egg4_folder_path}/{file_name}", engine="cfgrib"
        )

        if not egg4_is_set:
            egg4 = data
            egg4_is_set = True
            continue

        egg4 = xr.concat([egg4, data], dim="time")
        data.close()

    return egg4


def generate_file_list_from_atm_data(atm):
    yyyymmdd = "".join(str(min(atm.time).values)[:10].split("-"))

    file_name_list = [f"egg4_{yyyymmdd}.grb"]

    return file_name_list


def prepare_egg4(egg4):
    # prepare egg4
    egg4 = egg4.drop(["step", "time", "number"])
    egg4 = egg4.rename({"valid_time": "time", "step": "time"})
    egg4 = egg4.rename({"isobaricInhPa": "pressure"})

    return egg4


def correct_year(atm, egg4):
    # if year in atm differs from year in egg4 it means that no data was
    # available for the desired year. Instead, data for a preceding year was
    # downloaded and this needs to be accounted for. We will just act as if
    # the data is from the desired year, but write a little annotation into
    # the data's attributes.
    year_differs = False
    atm_time = str(atm.time.values[0])
    egg4_time = str(egg4.time.values[0])

    desired_year = int(atm_time[:4])
    used_year = int(egg4_time[:4])

    if desired_year == used_year:
        return egg4, year_differs, None

    year_differs = True
    for index, time in enumerate(egg4.time.values):
        # numpy datetime64 objects are hard to work with. Convert them into
        # regular datetime objects, do changes, then convert back.
        epoch = np.datetime64("1970-01-01T00:00:00Z")
        # convert to timestamp:
        ts = (time - epoch) / np.timedelta64(1, "s")
        # timestamp to standard utctime
        dt = datetime.utcfromtimestamp(ts)
        # replace year
        dt = dt.replace(year=desired_year)
        # utctime to numpy datetime64
        egg4.time.values[index] = np.datetime64(dt)

    return egg4, year_differs, used_year


def calculate_egg4_pressure(egg4):
    # convert from hPa to Pa
    egg4.pressure.values[:] *= 100

    return egg4


def calculate_egg4_co2_mole_fraction(egg4):
    # write mole fraction and not mass fraction.
    # use m_x = N_X * M_X
    egg4.co2.values = \
        egg4.co2 * constants.molar_mass_dry_air / constants.molar_mass_co2

    return egg4


def calculate_egg4_ch4_mole_fraction(egg4):
    # write mole fraction and not mass fraction.
    # use m_x = N_X * M_X
    egg4.ch4.values = \
        egg4.ch4 * constants.molar_mass_dry_air / constants.molar_mass_ch4

    return egg4


def interpolate_egg4_onto_atm(atm, egg4, dims):
    egg4 = egg4.interp(
        time=atm.time,
        latitude=atm.latitude,
        longitude=atm.longitude
    )

    egg4 = egg4.transpose(dims["y"], dims["x"], "pressure")

    return egg4


def interpolate_egg4_onto_pressure_grid(atm, egg4, dims):
    new_pressure_grid = atm.pressure
    old_pressure_grid = egg4.pressure

    # perform interpolation
    egg4_interpolated = xr.apply_ufunc(
        interpolate_pressures,
        new_pressure_grid, old_pressure_grid, egg4,
        input_core_dims=[["level"], ["pressure"], ["pressure"]],
        output_core_dims=[["level"]],
        vectorize=True
    )

    return egg4_interpolated


def interpolate_pressures(new_pressure_array, old_pressure_array, egg4):
    # If the atm pressure grid doesn't surround the egg4
    # pressure grid, there will be nan values at the boundaries.
    # Extrapolate vertical profiles constantly upwards and
    # downwards by overwriting all nan values with nearest
    # neighbor.

    # Perform interpolation.
    # interpolation has problems if the order of the pressure
    # grids does not align, therefore flip.
    egg4_interpolated = np.interp(
        new_pressure_array, old_pressure_array[::-1], egg4[::-1],
        left=np.nan, right=np.nan
    )

    # If old grid doesn't surround new grid, we have nans.
    # extrapolate them constantly
    valid_mask = ~np.isnan(egg4_interpolated)
    if False in valid_mask:
        left = np.where(valid_mask)[0][0]
        right = np.where(valid_mask)[0][-1]

        egg4_interpolated[:left] = egg4_interpolated[left]
        egg4_interpolated[right:] = egg4_interpolated[right]

    return egg4_interpolated
