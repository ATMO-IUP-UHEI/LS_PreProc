import xarray as xr
import cfgrib
import numpy as np

import functions.constants as constants

def import_cams(atm_data, cams_path):
    # get cams_data
    cams_date_list = ["20190806"]
    cams_data_set = False
    for cams_date in cams_date_list:
        data = xr.open_dataset(f"{cams_path}/cams_rea{cams_date}.grb", engine = "cfgrib")

        if not cams_data_set:
            cams_data = data
            cams_data_set = True
            continue

        cams_data = xr.concat([cams_data, data], dim = "time")
        data.close()

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

    # prepare atm_data
    nlevel = len(atm_data.level)
    nline = len(atm_data.line)
    nsample = len(atm_data.sample)

    atm_data["co2"] = (("level", "line", "sample"), np.empty(shape = (nlevel, nline, nsample)))
    atm_data.co2.attrs["standard_name"] = "CO2 molar mixing ratio"
    atm_data.co2.attrs["units"] = "mol mol-1"

    atm_data["ch4"] = (("level", "line", "sample"), np.empty(shape = (nlevel, nline, nsample)))
    atm_data.ch4.attrs["standard_name"] = "CH4 molar mixing ratio"
    atm_data.ch4.attrs["units"] = "mol mol-1"

    # interpolation begins
    datetime = atm_data.datetime.values
    for line in range(nline):
        for sample in range(nsample):
            latitude = atm_data.latitude.values[line, sample]
            longitude = atm_data.longitude.values[line, sample]
            pressure = atm_data.pressure.values[:, line, sample]

            # interpolate
            interpolated = cams_data.interp(latitude = latitude, longitude = longitude, datetime = datetime, pressure = pressure)

            # if the vertical pressure profile of cams does not surround the vertical pressure profile of
            # the atm_data, then the trace gas profiles will be cut off. They are extended constantly upwards
            # and downwards by replacing all nan values at the array boundaries by their nearest neighbor non-nan
            # value

            # co2
            non_nan_indices = np.where(~np.isnan(interpolated.co2_mass_fraction.values))[0]
            first_non_nan = non_nan_indices[0]
            last_non_nan = non_nan_indices[-1]
            interpolated.co2_mass_fraction.values[:first_non_nan] = interpolated.co2_mass_fraction.values[first_non_nan]
            interpolated.co2_mass_fraction.values[last_non_nan:] = interpolated.co2_mass_fraction.values[last_non_nan]

            # ch4
            non_nan_indices = np.where(~np.isnan(interpolated.ch4_mass_fraction.values))[0]
            first_non_nan = non_nan_indices[0]
            last_non_nan = non_nan_indices[-1]
            interpolated.ch4_mass_fraction.values[:first_non_nan] = interpolated.ch4_mass_fraction.values[first_non_nan]
            interpolated.ch4_mass_fraction.values[last_non_nan:] = interpolated.ch4_mass_fraction.values[last_non_nan]

            # write to atm_data
            # instead of mass fraction, molar mixing ratio is written to atm_data.
            # use:
            # m_x = N_X * M_X
            atm_data.co2[:, line, sample] = interpolated.co2_mass_fraction.values * constants.molar_mass_dry_air / constants.molar_mass_co2
            atm_data.ch4[:, line, sample] = interpolated.ch4_mass_fraction.values * constants.molar_mass_dry_air / constants.molar_mass_ch4

    cams_data.close()
