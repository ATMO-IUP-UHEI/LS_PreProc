import numpy as np
import xarray as xr
import cfgrib

import functions.constants as const

def import_era5(atm_data, era5_folder_path):
    # era5 data comes in hourly files. Every time step includes a multi-layer dataset and a surface dataset.
    # Here, the time steps surrounding the measurement are loaded in and combined into one dataset.
    era5_datetime_list = ["2019080614", "2019080615"]
    era5_data_set = False
    for era5_datetime in era5_datetime_list:
        ml = xr.open_dataset(f"{era5_folder_path}/ml{era5_datetime}.grb", engine = "cfgrib")
        sfc = xr.open_dataset(f"{era5_folder_path}/sfc{era5_datetime}.grb", engine = "cfgrib")
        merged = xr.merge([ml, sfc])
        ml.close()
        sfc.close()

        if not era5_data_set:
            era5_data = merged
            era5_data_set = True
            continue

        era5_data = xr.concat([era5_data, merged], dim = "time")
        merged.close()

    ndatetime = len(era5_data.datetime)
    nhybrid = len(era5_data.hybrid)
    nlatitude = len(era5_data.latitude)
    nlongitude = len(era5_data.longitude)

    # pressure for layers is given at lower boundary of the model level
    # pressure calculated using a and b parameters and the following manual:
    # https://www.ecmwf.int/sites/default/files/elibrary/2015/9210-part-iii-dynamics-and-numerical-procedures.pdf
    # a and b parameters downloaded from https://confluence.ecmwf.int/display/UDOC/L137+model+level+definitions
    era5_data["pressure"] = (("datetime", "hybrid", "latitude", "longitude"), np.empty(shape = (ndatetime, nhybrid, nlatitude, nlongitude)))
    n, a, b, ph, pf, gpa, gma, t, rho = np.genfromtxt("data/era5/standard_atmosphere.csv", delimiter = ",", skip_header = 2, unpack = True)
    for hybrid in range(nhybrid):
        era5_data.pressure.values[:, hybrid, :, :] = a[hybrid] + b[hybrid] * era5_data.sp.values[:, :, :]

    era5_data["geometric_altitude"] = (("hybrid"), gma)
    era5_data.geometric_altitude.attrs["standard_name"] = "geometric altitude"
    era5_data.geometric_altitude.attrs["units"] = "m"

    # cleaning up the dataset. Renaming variables and giving them attributes so they are more readable
    era5_data = era5_data.drop(["number", "step", "time", "surface", "sp"])
    era5_data = era5_data.rename({"hybrid": "model_level", "valid_time": "datetime", "time": "datetime"})
    era5_data = era5_data.rename({"z": "surface_geopotential", "t": "temperature", "q": "specific_humidity"})

    era5_data.latitude.attrs["standard_name"] = "latitude"
    era5_data.latitude.attrs["units"] = "degrees north"

    era5_data.longitude.attrs["standard_name"] = "longitude"
    era5_data.longitude.attrs["units"] = "degrees east"

    era5_data.pressure.attrs["standard_name"] = "pressure"
    era5_data.pressure.attrs["units"] = "kg m-1 s-2"

    era5_data.surface_geopotential.attrs["standard_name"] = "geopotential at surface"
    era5_data.surface_geopotential.attrs["units"] = "m2 s-2"

    era5_data.temperature.attrs["standard_name"] = "temperature"
    era5_data.temperature.attrs["units"] = "K"

    era5_data.specific_humidity.attrs["standard_name"] = "specific humidity"
    era5_data.specific_humidity.attrs["units"] = "kg kg-1"

    # interpolate era5 data onto l1b grid and time and write to atm_data
    # atm_data will get its own pressure grid
    nlevel = 60
    nline = len(atm_data.line)
    nsample = len(atm_data.sample)

    atm_data["geometric_altitude"] = (("level", "line", "sample"), np.empty(shape = (nlevel, nline, nsample)))
    atm_data.geometric_altitude.attrs["standard_name"] = "geometric altitude"
    atm_data.geometric_altitude.attrs["units"] = "m"

    atm_data["pressure"] = (("level", "line", "sample"), np.empty(shape = (nlevel, nline, nsample)))
    atm_data.pressure.attrs["standard_name"] = "pressure"
    atm_data.pressure.attrs["units"] = "kg m-1 s-2"

    atm_data["temperature"] = (("level", "line", "sample"), np.empty(shape = (nlevel, nline, nsample)))
    atm_data.temperature.attrs["standard_name"] = "temperature"
    atm_data.temperature.attrs["units"] = "K"

    atm_data["h2o"] = (("level", "line", "sample"), np.empty(shape = (nlevel, nline, nsample)))
    atm_data.h2o.attrs["standard_name"] = "molar mixing ratio of H2O"
    atm_data.h2o.attrs["units"] = "mol mol-1"

    # interpolation begins
    datetime = atm_data.datetime.values
    for line in range(nline):
        for sample in range(nsample):
            latitude = atm_data.latitude.values[line, sample]
            longitude = atm_data.longitude.values[line, sample]

            # interpolate
            interpolated = era5_data.interp(latitude = latitude, longitude = longitude, datetime = datetime)

            # correct surface values
            # era5 has a worse spatial resolution than other digital elevation models
            # surface elevation will therefore not be correct. When interpolating to a different surface elevation the bottom value of the data variables
            # must be corrected to account for this fact.
            dem_surface_elevation = atm_data.surface_elevation[line, sample].values # digital elevation model height
            era5_surface_elevation = interpolated.geometric_altitude.values[-1]

            # compare surface elevation
            surface_altitude_difference = era5_surface_elevation - dem_surface_elevation

            # if era5 < dem, the profile must be cut off at the bottom, such that the lowest level of the altitude profile
            # is above the ground
            if surface_altitude_difference < 0:
                interpolated = interpolated.where(interpolated.geometric_altitude >= dem_surface_elevation, drop = True)
                # recalculate surface elevation difference between dem and new cut off era5 dataset
                era5_surface_elevation = interpolated.geometric_altitude.values[-1]
                surface_altitude_difference = era5_surface_elevation - dem_surface_elevation

            # the lowest level of the altitude profile is just above the ground
            # data profiles must be extended downwards to touch the surface
            scale_height = interpolated.temperature.values[-1] * const.constants["gas_constant"] / const.constants["molar_mass_dry_air"] / const.gravitational_acceleration()

            # surface altitude is set to digital elevation model altitude
            interpolated.geometric_altitude.values[-1] = dem_surface_elevation

            # pressure is extended according to barometric height formula
            interpolated.pressure.values[-1] = interpolated.pressure.values[-1] * np.exp(surface_altitude_difference / scale_height)
            
            # temperature is extended according to adiabatic lapse rate
            interpolated.temperature.values[-1] = interpolated.temperature.values[-1] - const.constants["lapse_rate_lower_troposphere"] * surface_altitude_difference

            # specific humidity is extended downwards constantly (no change necessary, because value at index -1 is already correct)

            # interpolated to nlevel levels between surface and highest era5 level
            level = np.linspace(min(interpolated.model_level), max(interpolated.model_level), nlevel)
            interpolated = interpolated.interp(model_level = level)

            # write to atm_data
            atm_data.geometric_altitude.values[:, line, sample] = interpolated.geometric_altitude.values
            atm_data.pressure.values[:, line, sample] = interpolated.pressure.values
            atm_data.temperature.values[:, line, sample] = interpolated.temperature.values
            # instead of specific humidity, h2o molar mixing ratio is written to atm_data
            # use the following equations:
            # s = m_H2O / m_HUM (s = specific humidity)
            # m_X = N_X * m_X
            # m_HUM = m_DRY + m_H2O
            # to find: molar_mixing_ratio = N_H2O / N_DRY = s/(1-s) * M_DRY/M_H2O
            s = interpolated.specific_humidity.values
            h2o = s / (1 - s) * const.constants["molar_mass_dry_air"] / const.constants["molar_mass_h2o"]
            atm_data.h2o[:, line, sample] = h2o

    era5_data.close()
################################################################################
