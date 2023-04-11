import numpy as np
import xarray as xr
import cfgrib

import functions.constants as constants

def main(atm_data, era5_folder_path):
    print("getting era5 data")
    # get era5_data from the folder
    # name(s) of the required file will be constructed from the information in atm_data
    # era5_data comes in multilayer (ml) and surface (sfc) files
    era5_ml_data = get_era5_data(atm_data, era5_folder_path, "ml")
    era5_sfc_data = get_era5_data(atm_data, era5_folder_path, "sfc")

    # drop unnecessary variables from era5_data and rename needed ones
    era5_ml_data = prepare_era5_data(era5_ml_data, "ml")
    era5_sfc_data = prepare_era5_data(era5_sfc_data, "sfc")

    # prepare variables for atm_data
    atm_data = prepare_atm_data(atm_data)

    # interpolate to datetime and lat/lon grid of measurement
    interpolated = interpolate_era5_to_atm(era5_ml_data, era5_sfc_data, atm_data)

    # calculate quantities of interest
    interpolated = calculate_quantities_of_interest(interpolated)

    # correct for high-resolution surface elevation
    interpolated = correct_elevation(interpolated, atm_data)

    # interpolate linearly with pressure onto vertical atm_data grid
    atm_data = write_to_atm_data(interpolated, atm_data)

    return atm_data

    # geometric altitude (needs geopotential, which in turn needs virtual temperature, which in turn needs mixing ratio, which in turn needs dew point temperature, which in turn needs relative humidity, which in turn needs mixing ratio, which in turn needs dew point (something's wrong, lmao))
    # virtual temperature
    interpolated["virtual_temperature"] = (("level", "line", "sample"), np.empty(shape = (nlevel, nline, nsample)))
    interpolated.virtual_temperature.attrs["standard_name"] = "virtual temperature"
    interpolated.virtual_temperature.attrs["units"] = "K"
    # virtual temperature is approximately temperature. for now just set it equal and worry about it later
    interpolated.virtual_temperature.values = interpolated.temperature.values

    interpolated["geometric_altitude"] = (("level", "line", "sample"), np.empty(shape = (nlevel, nline, nsample)))
    interpolated.geometric_altitude.attrs["standard_name"] = "geometric altitude"
    interpolated.geometric_altitude.attrs["units"] = "m"

    # print("################################################################################")
    # # geopotential calculated from the same manual, formula 2.22
    # era5_data["geopotential"] = (("datetime", "model_level", "latitude", "longitude"), np.empty(shape = (ndatetime, nmodel_level, nlatitude, nlongitude)))
    # era5_data.geopotential.attrs["standard_name"] = "geopotential"
    # era5_data.geopotential.attrs["units"] = "m2 s-2"

    # era5_data["virtual_temperature"] = (("datetime", "model_level", "latitude", "longitude"), np.empty(shape = (ndatetime, nmodel_level, nlatitude, nlongitude)))
    # era5_data.virtual_temperature.attrs["standard_name"] = "virtual temperature"
    # era5_data.virtual_temperature.attrs["units"] = "K"
    # for model_level in range(nmodel_level):
    # #     https://en.wikipedia.org/wiki/Virtual_temperature
    #     era5_data.virtual_temperature[:, model_level, :, :] = era5_data.temperature[:, model_level, :, :] * (w + e)/(e * (1 + w))


    # era5_data["tmp"] = (("datetime", "model_level", "latitude", "longitude"), np.empty(shape = (ndatetime, nmodel_level, nlatitude, nlongitude)))
    # for i in range(nmodel_level-1): # bottom level is not changed, geopotential will be set to surface geopotential, therefore iterate only to nmodel_level-1
    #     era5_data.tmp[:, i, :, :] = np.log(era5_data.pressure[:, i+1, :, :]/era5_data.pressure[:, i, :, :])
    # era5_data.tmp[:, -1, :, :] = 0
    # print(era5_data.tmp[0, :, 0, 0])
    # #for i in range(nmodel_level):
    #     #era5_data.tmp[:, i, :, :] = era5_data.tmp[:, i:, :, :].sum(dim = "model_level")
    # print(era5_data.tmp[0, :, 0, 0])

    # for model_level in range(nmodel_level):
    #     era5_data.geopotential.values[:, model_level, :, :] = era5_data.surface_geopotential.values[:, :, :] + R * era5_data.virtual_temperature[:, :model_level+1, :, :] * era5_data.tmp[:, :model_level+1, :, :].sum(dim = "model_level")

    # era5_data["geometric_altitude"] = (("datetime", "model_level", "latitude", "longitude"), era5_data.geopotential.values/constants.gravitational_acceleration())
    # era5_data.geometric_altitude.attrs["standard_name"] = "geometric altitude"
    # era5_data.geometric_altitude.attrs["units"] = "m"

    # print("TESTING TESTING")
    # #print(era5_data.geopotential[0, :, 1, 1].values)
    # #print(era5_data.surface_geopotential[0, 1, 1].values)
    # #print(era5_data.geometric_altitude[0, :, 1, 1].values)
    # print("################################################################################")
    # import sys
    # sys.exit()



def get_era5_data(atm_data, era5_folder_path, era5_type):
    file_name_list = generate_file_list_from_atm_data(atm_data, era5_type)

    era5_data_set = False
    for file_name in file_name_list:
        era5_data = xr.open_dataset(f"{era5_folder_path}/{file_name}", engine = "cfgrib")

        if not era5_data_set:
            merged = era5_data
            era5_data_set = True

        era5_data = xr.merge([era5_data, merged])
        merged.close()
    return era5_data



def generate_file_list_from_atm_data(atm_data, era5_type):
    yyyymmdd = "".join((atm_data.attrs["ISO 8601 datetime"].split("T")[0]).split("-"))

    file_name_list = [f"era5_{era5_type}_{yyyymmdd}.grb"]

    return file_name_list



def prepare_era5_data(era5_data, era5_type):
    era5_data.latitude.attrs["standard_name"] = "latitude"
    era5_data.latitude.attrs["units"] = "degrees north"

    era5_data.longitude.attrs["standard_name"] = "longitude"
    era5_data.longitude.attrs["units"] = "degrees east"

    if era5_type == "ml":
        era5_data = era5_data.drop(["step", "valid_time"])
        era5_data = era5_data.rename({"hybrid": "level", "time": "datetime"})
        era5_data = era5_data.drop_vars("level")
        era5_data = era5_data.rename({"t": "temperature", "q": "specific_humidity"})

        era5_data.temperature.attrs["standard_name"] = "temperature"
        era5_data.temperature.attrs["units"] = "K"
        era5_data = era5_data.assign(temperature = era5_data.temperature.astype(np.float64))

        era5_data.specific_humidity.attrs["standard_name"] = "specific humidity"
        era5_data.specific_humidity.attrs["units"] = "kg kg-1"
        era5_data = era5_data.assign(specific_humidity = era5_data.specific_humidity.astype(np.float64))

    if era5_type == "sfc":
        era5_data = era5_data.drop(["number", "step", "surface", "valid_time"])
        era5_data = era5_data.rename({"time": "datetime"})
        era5_data = era5_data.rename({"z": "surface_geopotential", "sp": "surface_pressure"})

        era5_data.surface_geopotential.attrs["standard_name"] = "geopotential at surface"
        era5_data.surface_geopotential.attrs["units"] = "m2 s-2"
        era5_data = era5_data.assign(surface_geopotential = era5_data.surface_geopotential.astype(np.float64))

        era5_data.surface_pressure.attrs["standard_name"] = "surface pressure"
        era5_data.surface_pressure.attrs["units"] = "kg m-1 s-2"
        era5_data = era5_data.assign(surface_pressure = era5_data.surface_pressure.astype(np.float64))

        era5_data = era5_data.assign(u10 = era5_data.u10.astype(np.float64))
        era5_data = era5_data.assign(v10 = era5_data.v10.astype(np.float64))

    return era5_data



def prepare_atm_data(atm_data):
    nlevel = 60
    nline = len(atm_data.line)
    nsample = len(atm_data.sample)

    atm_data["pressure"] = (("level", "line", "sample"), np.empty(shape = (nlevel, nline, nsample)))
    atm_data.pressure.attrs["standard_name"] = "pressure"
    atm_data.pressure.attrs["units"] = "kg m-1 s-2"

    atm_data["temperature"] = (("level", "line", "sample"), np.empty(shape = (nlevel, nline, nsample)))
    atm_data.temperature.attrs["standard_name"] = "temperature"
    atm_data.temperature.attrs["units"] = "K"

    atm_data["geometric_altitude"] = (("level", "line", "sample"), np.empty(shape = (nlevel, nline, nsample)))
    atm_data.geometric_altitude.attrs["standard_name"] = "geometric altitude"
    atm_data.geometric_altitude.attrs["units"] = "m"

    atm_data["h2o"] = (("level", "line", "sample"), np.empty(shape = (nlevel, nline, nsample)))
    atm_data.h2o.attrs["standard_name"] = "H2O mole fraction"
    atm_data.h2o.attrs["units"] = "mol mol-1"

    return atm_data



def interpolate_era5_to_atm(era5_ml_data, era5_sfc_data, atm_data):
    interpolated_ml = era5_ml_data.interp(datetime = np.datetime64(atm_data.attrs["ISO 8601 datetime"]), latitude = atm_data.latitude, longitude = atm_data.longitude)
    interpolated_sfc = era5_sfc_data.interp(datetime = np.datetime64(atm_data.attrs["ISO 8601 datetime"]), latitude = atm_data.latitude, longitude = atm_data.longitude)
    interpolated = xr.merge([interpolated_ml, interpolated_sfc])

    return interpolated



def calculate_quantities_of_interest(interpolated):
    nlevel = len(interpolated.level)
    nline = len(interpolated.line)
    nsample = len(interpolated.sample)

    # get standard ecmwd atmosphere parameters from file downloaded from https://confluence.ecmwf.int/display/UDOC/L137+model+level+definitions
    try:
        n, a, b, ph, pf, gpa, gma, t, rho = np.genfromtxt("data/external/era5/auxiliary/standard_atmosphere.csv", delimiter = ",", skip_header = 2, unpack = True)
    except:
        sys.exit("standard_atmosphere.csv missing. Download from https://confluence.ecmwf.int/display/UDOC/L137+model+level+definitions and put it into data/external/era5/auxiliary/")

    # pressure
    # pressure for layers is given at lower boundary of the model level
    # pressure calculated using a and b parameters and the following manual:
    # https://www.ecmwf.int/sites/default/files/elibrary/2015/9210-part-iii-dynamics-and-numerical-procedures.pdf
    interpolated["pressure"] = (("level", "line", "sample"), np.empty(shape = (nlevel, nline, nsample)))
    interpolated.pressure.attrs["standard_name"] = "pressure"
    interpolated.pressure.attrs["units"] = "kg m-1 s-2"

    pressures = np.tile(interpolated.surface_pressure, (nlevel, 1, 1))
    pressures = np.swapaxes(pressures, 0, 2)
    pressures = a + b * pressures
    pressures = np.swapaxes(pressures, 0, 2)
    interpolated.pressure.values = pressures

    # water vapor mole fraction
    # specific humidity is converted into h2o mole fraction using
    # s = m_H2O / m_HUM (s = specific humidity)
    # m_X = N_X * m_X
    # m_HUM = m_DRY + m_H2O
    # to find: h2o_mole_fraction = N_H2O / N_DRY = s/(1-s) * M_DRY/M_H2O
    interpolated["h2o"] = (("level", "line", "sample"), np.empty(shape = (nlevel, nline, nsample)))
    interpolated.h2o.attrs["standard_name"] = "H2O mole fraction"
    interpolated.h2o.attrs["units"] = "mol mol-1"

    h2o = interpolated.specific_humidity / (1 - interpolated.specific_humidity) * constants.molar_mass_dry_air / constants.molar_mass_h2o
    interpolated.h2o.values = h2o

    # geometric altitude
    interpolated["geometric_altitude"] = (("level", "line", "sample"), np.empty(shape = (nlevel, nline, nsample)))
    interpolated.geometric_altitude.attrs["standard_name"] = "geometric altitude"
    interpolated.geometric_altitude.attrs["units"] = "m"

    print("TODO this is only for the standard atmosphere, you need to calculate geometric altitude using the geopotential and gravitational constant")
    gma = np.tile(gma, (nline, nsample, 1))
    gma = np.swapaxes(gma, 0, 2)
    gma = np.swapaxes(gma, 1, 2)

    interpolated.geometric_altitude.values = gma

    return interpolated



def correct_elevation(interpolated, atm_data):
    # operations need to be performed on each line, sample pixel of the interpolated dataset
    # this is done here by grouping the dataset by location. Since grouping by multiple
    # dimensions is not supported by xarray, yet, line and sample are stacked into a new
    # dimension called location, first. After modifying the dataset it is unstacked again
    interpolated = interpolated.stack(location=["line", "sample"])
    interpolated = interpolated.groupby("location")
    interpolated = interpolated.apply(correct_elevation_for_pixel, args=[atm_data])
    interpolated = interpolated.unstack("location")
    interpolated = interpolated.drop_vars(["line", "sample"])

    return interpolated



def correct_elevation_for_pixel(interpolated_pixel, atm_data):
    nlevel = len(interpolated_pixel.level)

    # drop all altitudes from era5 that are below the dem surface elevation
    dem_surface_elevation = atm_data.isel(line=interpolated_pixel.line, sample=interpolated_pixel.sample).surface_elevation
    interpolated_pixel = interpolated_pixel.where(interpolated_pixel.geometric_altitude >= dem_surface_elevation)

    # lowest value has to be set to surface level
    # if no nan value exists, replace bottom level with ground level
    # if at least one nan value exists at and above the bottom level,
    # replace first nan value with ground level
    if np.isnan(interpolated_pixel.pressure[-1]):
        ground_level = np.nanargmax(interpolated_pixel.pressure) + 1
    else:
        ground_level = nlevel - 1

    # difference between digital elevation model surface elevation
    # and the lowest elevation given in the era5 model above the ground
    # the lowest level of the altitude profile is just above the ground
    # data profiles must be extended downwards to touch the surface
    era5_lowest_elevation = interpolated_pixel.geometric_altitude[ground_level-1]
    elevation_difference = era5_lowest_elevation - dem_surface_elevation

    # surface altitude is set to digital elevation model altitude
    interpolated_pixel.geometric_altitude[ground_level] = dem_surface_elevation

    # pressure is extended according to barometric height formula
    scale_height = interpolated_pixel.temperature[ground_level] * constants.gas_constant / constants.molar_mass_dry_air / constants.gravitational_acceleration()
    interpolated_pixel.pressure[ground_level] = interpolated_pixel.pressure[ground_level-1] * np.exp(elevation_difference / scale_height)

    # temperature is extended according to adiabatic lapse rate
    interpolated_pixel.temperature[ground_level] = interpolated_pixel.temperature[ground_level-1] - constants.lapse_rate_lower_troposphere * elevation_difference

    # specific humidity is extended downwards constantly (no change necessary, because value at index -1 is already correct)

    return interpolated_pixel



def write_to_atm_data(interpolated, atm_data):
    # operations need to be performed on each line, sample pixel of the interpolated dataset
    # this is done here by grouping the dataset by location. Since grouping by multiple
    # dimensions is not supported by xarray, yet, line and sample are stacked into a new
    # dimension called location, first. After modifying the dataset it is unstacked again
    interpolated = interpolated.stack(location=["line", "sample"])
    interpolated = interpolated.groupby("location")
    interpolated = interpolated.apply(interpolate_pixel_to_atm_pressure, args=[atm_data])
    interpolated = interpolated.unstack("location")
    interpolated = interpolated.drop_vars(["line", "sample"])

    atm_data.temperature.values = interpolated.temperature.values
    atm_data.pressure.values = interpolated.pressure.values
    atm_data.geometric_altitude.values = interpolated.geometric_altitude.values
    atm_data.h2o.values = interpolated.h2o.values

    return atm_data



def interpolate_pixel_to_atm_pressure(interpolated_pixel, atm_data):
    atm_nlevel = 60
    min_pressure = np.nanmin(interpolated_pixel.pressure)
    min_pressure = np.nanmax([2, min_pressure])
    max_pressure = np.nanmax(interpolated_pixel.pressure)
    level_grid = interpolated_pixel.level.values
    pressure_grid = interpolated_pixel.pressure.values

    # interpolate onto grid linear in pressure (there has to be a built in function for this)
    new_pressure_grid = np.linspace(min_pressure, max_pressure, atm_nlevel)
    new_level_grid = np.interp(new_pressure_grid, pressure_grid, level_grid)
    interpolated_pixel = interpolated_pixel.interp(level=new_level_grid)

    interpolated_pixel = interpolated_pixel.drop_vars(["level"])

    return interpolated_pixel



if __name__ == "__main__":
    main(atm_data, era5_folder_path)
