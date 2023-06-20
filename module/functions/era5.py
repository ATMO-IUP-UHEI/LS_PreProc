import numpy as np
import xarray as xr
import cfgrib
import sys
import os

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
    interpolated = calculate_quantities_of_interest(interpolated, era5_folder_path, atm_data)

    # correct for high-resolution surface elevation
    interpolated = correct_elevation(interpolated, atm_data)

    # interpolate linearly with pressure onto vertical atm_data grid
    interpolated = interpolate_interpolated_to_target_pressure_grid(interpolated, atm_data)

    # write to atm_data
    atm_data = write_interpolated_to_atm_data(interpolated, atm_data)

    return atm_data



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

    atm_data["pressure"] = (("line", "sample", "level"), np.empty(shape = (nline, nsample, nlevel)))
    atm_data.pressure.attrs["standard_name"] = "pressure"
    atm_data.pressure.attrs["units"] = "kg m-1 s-2"

    atm_data["temperature"] = (("line", "sample", "level"), np.empty(shape = (nline, nsample, nlevel)))
    atm_data.temperature.attrs["standard_name"] = "temperature"
    atm_data.temperature.attrs["units"] = "K"

    atm_data["geometric_altitude"] = (("line", "sample", "level"), np.empty(shape = (nline, nsample, nlevel)))
    atm_data.geometric_altitude.attrs["standard_name"] = "geometric altitude"
    atm_data.geometric_altitude.attrs["units"] = "m"

    atm_data["h2o"] = (("line", "sample", "level"), np.empty(shape = (nline, nsample, nlevel)))
    atm_data.h2o.attrs["standard_name"] = "H2O mole fraction"
    atm_data.h2o.attrs["units"] = "mol mol-1"

    return atm_data



def interpolate_era5_to_atm(era5_ml_data, era5_sfc_data, atm_data):
    interpolated_ml = era5_ml_data.interp(datetime = np.datetime64(atm_data.attrs["ISO 8601 datetime"]), latitude = atm_data.latitude, longitude = atm_data.longitude)
    interpolated_sfc = era5_sfc_data.interp(datetime = np.datetime64(atm_data.attrs["ISO 8601 datetime"]), latitude = atm_data.latitude, longitude = atm_data.longitude)
    interpolated = xr.merge([interpolated_ml, interpolated_sfc])

    interpolated = interpolated.transpose("line", "sample", "level")

    return interpolated



def calculate_quantities_of_interest(interpolated, era5_folder_path, atm_data):
    nline = len(interpolated.line)
    nsample = len(interpolated.sample)
    nlevel = len(interpolated.level)

    # get standard ecmwd atmosphere parameters from file downloaded from https://confluence.ecmwf.int/display/UDOC/L137+model+level+definitions
    try:
        auxiliary_path = os.path.join(era5_folder_path, "..", "auxiliary")
        n, a, b, ph, pf, gpa, gma, t, rho = np.genfromtxt(
            os.path.join(auxiliary_path, "standard_atmosphere.csv"),
            delimiter = ",", skip_header = 2, unpack = True)
    except:
        sys.exit("standard_atmosphere.csv missing. Download from https://confluence.ecmwf.int/display/UDOC/L137+model+level+definitions and put it into data/external/era5/auxiliary/")

    # pressure
    # pressure for layers is given at lower boundary of the model level
    # pressure calculated using a and b parameters and the following manual:
    # https://www.ecmwf.int/sites/default/files/elibrary/2015/9210-part-iii-dynamics-and-numerical-procedures.pdf
    interpolated["pressure"] = (("line", "sample", "level"), np.empty(shape = (nline, nsample, nlevel)))
    interpolated.pressure.attrs["standard_name"] = "pressure"
    interpolated.pressure.attrs["units"] = "kg m-1 s-2"

    pressures = a[None, None, :] + b[None, None, :] * interpolated.surface_pressure.values[:, :, None]
    interpolated.pressure.values = pressures

    # water vapor mole fraction
    # specific humidity is converted into h2o mole fraction using
    # s = m_H2O / m_HUM (s = specific humidity)
    # m_X = N_X * m_X
    # m_HUM = m_DRY + m_H2O
    # to find: h2o_mole_fraction = N_H2O / N_DRY = s/(1-s) * M_DRY/M_H2O
    interpolated["h2o"] = (("line", "sample", "level"), np.empty(shape = (nline, nsample, nlevel)))
    interpolated.h2o.attrs["standard_name"] = "H2O mole fraction"
    interpolated.h2o.attrs["units"] = "mol mol-1"

    h2o = interpolated.specific_humidity / (1 - interpolated.specific_humidity) * constants.molar_mass_dry_air / constants.molar_mass_h2o
    interpolated.h2o.values = h2o

    # geometric altitude
    interpolated["geometric_altitude"] = (("line", "sample", "level"), np.empty(shape = (nline, nsample, nlevel)))
    interpolated.geometric_altitude.attrs["standard_name"] = "geometric altitude"
    interpolated.geometric_altitude.attrs["units"] = "m"

    # TODO: what happens here?
    print("todo: is low resolution surface elevation actually compared to high resolution surface elevation later?")
    gma = atm_data.surface_elevation.values[:, :, None] + constants.scale_height * np.log(interpolated.pressure[..., -1]/interpolated.pressure[...])
    interpolated.geometric_altitude.values = gma

    return interpolated



def correct_elevation(interpolated, atm_data):

    interpolated.geometric_altitude[0, 0].values += -20
    interpolated.geometric_altitude[0, 1].values += +20

    interpolated = interpolated.where(interpolated.geometric_altitude >= atm_data.surface_elevation)

    ground_level = interpolated.geometric_altitude.argmin(dim="level", skipna=True)

    elevation_difference = interpolated.geometric_altitude[..., ground_level] - atm_data.surface_elevation

    interpolated.geometric_altitude[..., ground_level] = atm_data.surface_elevation

    scale_height = interpolated.temperature[..., ground_level] * constants.gas_constant / constants.molar_mass_dry_air / constants.gravitational_acceleration()
    interpolated.pressure[..., ground_level] = interpolated.pressure[..., ground_level] * np.exp(elevation_difference / scale_height)

    interpolated.temperature[..., ground_level] = interpolated.temperature[..., ground_level] - constants.lapse_rate_lower_troposphere * elevation_difference

    # interpolated.specific_humidity[..., ground_level] is assumed to have a constant profile, no change necessary

    return interpolated



def correct_elevation_old(interpolated, atm_data):
    # operations need to be performed on each line, sample pixel of the interpolated dataset
    # this is done here by grouping the dataset by location. Since grouping by multiple
    # dimensions is not supported by xarray, yet, line and sample are stacked into a new
    # dimension called location, first. After modifying the dataset it is unstacked again
    interpolated = interpolated.stack(location=["line", "sample"])
    interpolated = interpolated.groupby("location")
    interpolated = interpolated.apply(correct_elevation_for_pixel, args=[atm_data])
    interpolated = interpolated.unstack("location")
    interpolated = interpolated.drop_vars(["line", "sample"])
    interpolated = interpolated.transpose("line", "sample", "level")

    return interpolated



def correct_elevation_for_pixel(interpolated_pixel, atm_data):
    nlevel = len(interpolated_pixel.level)

    # drop all altitudes from era5 that are below the dem surface elevation
    dem_surface_elevation = atm_data.isel(line=interpolated_pixel.line, sample=interpolated_pixel.sample).surface_elevation
    interpolated_pixel = interpolated_pixel.where(interpolated_pixel.geometric_altitude >= dem_surface_elevation)

    # lowest value has to be set to surface level
    # if no nan value exists, replace bottom level with ground level
    # if at least one nan value exists at and above the bottom level,
    # replace first non nan value with ground level
    ground_level = np.nanargmax(interpolated_pixel.pressure)

    # difference between digital elevation model surface elevation
    # and the lowest elevation given in the era5 model above the ground
    # the lowest level of the altitude profile is just above the ground
    # data profiles must be extended downwards to touch the surface
    era5_lowest_elevation = interpolated_pixel.geometric_altitude[ground_level]
    elevation_difference = era5_lowest_elevation - dem_surface_elevation

    # surface altitude is set to digital elevation model altitude
    interpolated_pixel.geometric_altitude[ground_level] = dem_surface_elevation

    # pressure is extended according to barometric height formula
    scale_height = interpolated_pixel.temperature[ground_level] * constants.gas_constant / constants.molar_mass_dry_air / constants.gravitational_acceleration()
    interpolated_pixel.pressure[ground_level] = interpolated_pixel.pressure[ground_level] * np.exp(elevation_difference / scale_height)

    # temperature is extended according to adiabatic lapse rate
    interpolated_pixel.temperature[ground_level] = interpolated_pixel.temperature[ground_level] - constants.lapse_rate_lower_troposphere * elevation_difference

    # specific humidity is extended downwards constantly (no change necessary)

    return interpolated_pixel



def interpolate_interpolated_to_target_pressure_grid(interpolated, atm_data):
    min_pressure_grid = interpolated.pressure.min(dim="level", skipna=True)
    min_pressure_grid = min_pressure_grid.where(min_pressure_grid > 2, other=2)
    max_pressure_grid = interpolated.pressure.max(dim="level", skipna=True)

    pressure_grid = np.linspace(min_pressure_grid, max_pressure_grid, len(atm_data.level))
    pressure_grid = pressure_grid.transpose([1, 2, 0])
    pressure_grid = xr.DataArray(
        pressure_grid,
        dims=["line", "sample", "level_target"],
        coords={"line": min_pressure_grid.line.values, "sample": min_pressure_grid.sample.values}
    )

    import xgcm

    interpolated = interpolated.assign_coords({
        "line": ("line", range(len(interpolated.line.values))),
        "sample": ("sample", range(len(interpolated.sample.values))),
        "level": ("level", range(len(interpolated.level.values)))})

    interpolated.line.attrs["axis"] = "X"
    interpolated.sample.attrs["axis"] = "Y"
    interpolated.level.attrs["axis"] = "Z"

    interpolated = interpolated.chunk({"line": 250, "sample": 250, "level": -1})

    grid = xgcm.Grid(interpolated, periodic=False)
    from functools import partial
    
    def transformer(data, target_levels, target_data):
        data = xr.DataArray(
            data,
            dims=["level"],
            coords={"level": ("level", range(len(data)))})
        data.level.attrs["axis"] = "Z"

        target_data = xr.DataArray(
            target_data,
            dims=["level"],
            coords={"level": ("level", range(len(target_data)))})
        target_data.level.attrs["axis"] = "Z"

        return grid.transform(data, 'Z', target_levels, target_data=target_data, method='log')

    for var in set(interpolated.data_vars) - {"pressure"}:
        interpolated[var].persist()
        if "level" in interpolated[var].dims:
            trans_var = xr.apply_ufunc(
                transformer,
                *[interpolated[var], pressure_grid, interpolated.pressure],
                input_core_dims=[['level'], ['level_target'], ['level']],
                output_core_dims=[['level_target']],
                dask="parallelized",
                vectorize=True
            )

            interpolated[var] = trans_var

    interpolated["pressure"] = pressure_grid

    sys.exit(interpolated)



    return interpolated



def interpolate_interpolated_to_target_pressure_grid_toffer(interpolated, atm_data):
    min_pressure_grid = interpolated.pressure.min(dim="level", skipna=True)
    min_pressure_grid = min_pressure_grid.where(min_pressure_grid > 2, other=2)
    max_pressure_grid = interpolated.pressure.max(dim="level", skipna=True)

    pressure_grid = np.linspace(min_pressure_grid, max_pressure_grid, len(atm_data.level))
    pressure_grid = pressure_grid.transpose([1, 2, 0])
    pressure_grid = xr.DataArray(
        pressure_grid,
        dims=["line", "sample", "level_target"],
        coords={"line": min_pressure_grid.line.values, "sample": min_pressure_grid.sample.values}
    )
    
    # interpolated.pressure: line, sample, level   3, 4, 137
    # pressure_grid: line, sample, level_target    3, 4, 60

    # # with dask
    # pressure_difference = pressure_grid.chunk({"line": 10, "sample": 10}) - interpolated.pressure.chunk({"line": 10, "sample": 10})
    # print("compute")
    # import dask.distributed
    # client = dask.distributed.Client()
    # print(client)
    # client.compute(pressure_difference)
    # pressure_difference = pressure_difference.compute(parallel=False)
    # print("done computing")

    # without dask
    pressure_difference = pressure_grid - interpolated.pressure

    transition_matrix = xr.zeros_like(pressure_difference)

    # TODO: document what this does :)
    level_index_above_level_target = pressure_difference.where(pressure_difference>=0).argmin(dim="level", skipna=True)
    level_index_below_level_target = pressure_difference.where(pressure_difference<=0).argmax(dim="level", skipna=True)
    pressure_diff_at_upper_level = pressure_difference[..., level_index_above_level_target]
    pressure_diff_at_lower_level = pressure_difference[..., level_index_below_level_target]
    pressure_diff_between_levels = pressure_diff_at_upper_level - pressure_diff_at_lower_level
    upper_weight = np.abs(pressure_diff_at_lower_level) / pressure_diff_between_levels
    lower_weight = pressure_diff_at_upper_level / pressure_diff_between_levels
    # handle case where the pressure on a level equals the pressure on a target_level, therefore pressure_diff is zero
    upper_weight = upper_weight.fillna(0)
    lower_weight = lower_weight.fillna(1)
    transition_matrix[..., level_index_above_level_target] = upper_weight
    transition_matrix[..., level_index_below_level_target] = lower_weight

    magic_interpolation_constant = 42
    for variable in interpolated.data_vars:
        if "level" in interpolated[variable].dims:
            interpolated[f"{variable}_target"] = interpolated.fillna(magic_interpolation_constant)[variable].dot(transition_matrix, dims="level")
            interpolated = interpolated.drop_vars(variable).rename({f"{variable}_target": variable})

    interpolated = interpolated.rename_dims(level_target="level")

    print(interpolated)
    sys.exit()

    return interpolated



def write_interpolated_to_atm_data(interpolated, atm_data):
    atm_data.temperature.values = interpolated.temperature.values
    atm_data.pressure.values = interpolated.pressure.values
    atm_data.geometric_altitude.values = interpolated.geometric_altitude.values
    atm_data.h2o.values = interpolated.h2o.values

    return atm_data



def write_to_atm_data_old(interpolated, atm_data):
    # operations need to be performed on each line, sample pixel of the interpolated dataset
    # this is done here by grouping the dataset by location. Since grouping by multiple
    # dimensions is not supported by xarray, yet, line and sample are stacked into a new
    # dimension called location, first. After modifying the dataset it is unstacked again
    interpolated = interpolated.stack(location=["line", "sample"])
    interpolated = interpolated.groupby("location")
    interpolated = interpolated.apply(interpolate_pixel_to_atm_pressure, args=[atm_data])
    interpolated = interpolated.unstack("location")
    interpolated = interpolated.drop_vars(["line", "sample"])
    interpolated = interpolated.transpose("line", "sample", "level")

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
