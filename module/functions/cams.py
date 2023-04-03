import xarray as xr
import cfgrib
import numpy as np

import functions.constants as constants

def main(atm_data, cams_folder_path):
    print("getting cams data")

    # get cams_data from the folder.
    # name(s) of the required file will be constructed from the information in atm_data
    cams_data = get_cams_data(atm_data, cams_folder_path)

    # drop unnecessary variables from cams_data and rename needed ones
    cams_data = prepare_cams_data(cams_data)

    # prepare variables for atm_data
    atm_data = prepare_atm_data(atm_data)

    # interpolation begins
    interpolated_full = cams_data.interp(latitude = atm_data.latitude, longitude = atm_data.longitude, datetime = atm_data.datetime, pressure = atm_data.pressure)

    nline = len(atm_data.line.values)
    nsample = len(atm_data.sample.values)

    for line in range(nline):
        for sample in range(nsample):
            interpolated = interpolated_full.isel(line = line, sample = sample)

            # if the vertical pressure profile of cams does not surround the vertical pressure profile of
            # the atm_data, then the trace gas profiles will be cut off. They are extended constantly upwards
            # and downwards by replacing all nan values at the array boundaries by their nearest neighbor non-nan
            # value
            # write molar mixing ratio to atm_data, not mass fraction
            # use:
            #     m_x = N_X * M_X

            # co2
            non_nan_indices = np.where(~np.isnan(interpolated.co2_mass_fraction.values))[0]
            first_non_nan = non_nan_indices[0]
            last_non_nan = non_nan_indices[-1]
            interpolated.co2_mass_fraction.values[:first_non_nan] = interpolated.co2_mass_fraction.values[first_non_nan]
            interpolated.co2_mass_fraction.values[last_non_nan:] = interpolated.co2_mass_fraction.values[last_non_nan]
            atm_data.co2[:, line, sample] = interpolated.co2_mass_fraction.values * constants.molar_mass_dry_air / constants.molar_mass_co2

            # ch4
            non_nan_indices = np.where(~np.isnan(interpolated.ch4_mass_fraction.values))[0]
            first_non_nan = non_nan_indices[0]
            last_non_nan = non_nan_indices[-1]
            interpolated.ch4_mass_fraction.values[:first_non_nan] = interpolated.ch4_mass_fraction.values[first_non_nan]
            interpolated.ch4_mass_fraction.values[last_non_nan:] = interpolated.ch4_mass_fraction.values[last_non_nan]
            atm_data.ch4[:, line, sample] = interpolated.ch4_mass_fraction.values * constants.molar_mass_dry_air / constants.molar_mass_ch4

    cams_data.close()

    return atm_data



def get_cams_data(atm_data, cams_folder_path):
    file_name_list = generate_file_list_from_atm_data(atm_data)

    cams_data_set = False
    for file_name in file_name_list:
        # loop can theoretically handle multiple datetime entries, but this is not implemented
        # in the rest of the scripts
        data = xr.open_dataset(f"{cams_folder_path}/{file_name}", engine = "cfgrib")

        if not cams_data_set:
            cams_data = data
            cams_data_set = True
            continue

        cams_data = xr.concat([cams_data, data], dim = "time")
        data.close()

    return cams_data



def generate_file_list_from_atm_data(atm_data):
    datetime = str(atm_data.datetime.values)
    yyyymmdd = "".join(datetime[:10].split("-"))

    file_name_list = [f"egg4_{yyyymmdd}.grb"]

    return file_name_list



def prepare_cams_data(cams_data):
    # prepare cams_data
    cams_data = cams_data.drop(["step", "time", "number"])
    cams_data = cams_data.rename({"valid_time": "datetime", "step": "datetime"})
    cams_data = cams_data.rename({"isobaricInhPa": "pressure"})
    cams_data = cams_data.rename({"co2": "co2_mass_fraction", "ch4": "ch4_mass_fraction"})

    cams_data.latitude.attrs["standard_name"] = "latitude"
    cams_data.latitude.attrs["units"] = "degrees north"

    cams_data.longitude.attrs["standard_name"] = "longitude"
    cams_data.longitude.attrs["units"] = "degrees east"

    cams_data.pressure.attrs["standard_name"] = "pressure"
    cams_data.pressure.values[:] *= 100 # convert from hPa to Pa = kg m-1 s-2
    cams_data.pressure.attrs["units"] = "kg m-1 s-2"

    cams_data.co2_mass_fraction.attrs["standard_name"] = "CO2 mass fraction"
    cams_data.co2_mass_fraction.attrs["units"] = "kg kg-1"

    cams_data.ch4_mass_fraction.attrs["standard_name"] = "CH4 mass fraction"
    cams_data.ch4_mass_fraction.attrs["units"] = "kg kg-1"

    return cams_data



def prepare_atm_data(atm_data):
    # prepare atm_data
    nlevel = len(atm_data.level)
    nline = len(atm_data.line)
    nsample = len(atm_data.sample)

    atm_data["co2"] = (("level", "line", "sample"), np.empty(shape = (nlevel, nline, nsample)))
    atm_data.co2.attrs["standard_name"] = "CO2 mole fraction"
    atm_data.co2.attrs["units"] = "mol mol-1"

    atm_data["ch4"] = (("level", "line", "sample"), np.empty(shape = (nlevel, nline, nsample)))
    atm_data.ch4.attrs["standard_name"] = "CH4 mole fraction"
    atm_data.ch4.attrs["units"] = "mol mol-1"

    return atm_data


if __name__ == "__main__":
    main(atm_data, cams_folder_path)
