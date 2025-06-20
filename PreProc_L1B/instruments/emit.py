from datetime import datetime
from datetime import timedelta
from cftime import date2num
import netCDF4 as nc
import xarray as xr
import numpy as np
import sys
import os

def import_data(config, dims):
    
    # Read config file
    inp_path = config.get("config","inp_path")
    tile_list = config.get("config","L1B").split()
    snr = config.get("config","SNR")
    
    # Loop over input files to be concatenated
    input_file_list = []
    for i, tile in enumerate(tile_list):

        tile_path = os.path.join(inp_path, tile)
        tile_root, tile_band = import_tile(tile_path, snr, dims)
        input_file_list.append(tile)
        
        if i == 0:
            root = tile_root
            band = tile_band
        else:
            root = xr.concat(
                (root, tile_root), dim="frame", data_vars="minimal")
            band = xr.concat(
                (band, tile_band), dim="frame", data_vars="minimal")
    
    return root, [band], input_file_list

def import_tile(l1bpath, snr, dims):
    
    if 'OBS' in l1bpath:
        l1bobs = xr.open_datatree(l1bpath)
        l1brad = xr.open_datatree(l1bpath.replace("_OBS_", "_RAD_"))
    elif 'RAD' in l1bpath:
        l1brad = xr.open_datatree(l1bpath)
        l1bobs = xr.open_datatree(l1bpath.replace("_RAD_", "_OBS_"))
    else:
        sys.exit("L1B file must contain either 'OBS' or 'RAD' in the filename.")

    names=list(l1bobs["sensor_band_parameters"]["observation_bands"].values)
    
    data = xr.Dataset()
    band1_data = xr.Dataset()
    
    time, time_units_string = get_time(l1brad)
    data["time"] = xr.DataArray(
        data=time,
        dims=(dims["y"]),
    ).astype("float32")
    data.time.attrs["units"] = time_units_string
    
    data["latitude"] = xr.DataArray(
        data=l1brad["location"]["lat"],
        dims=(dims["y"], dims["x"]),
    ).astype("float32")
    
    data["longitude"] = xr.DataArray(
        data=l1brad["location"]["lon"],
        dims=(dims["y"], dims["x"]),
    ).astype("float32")

    data["elevation"] = xr.DataArray(
        data=l1brad["location"]["elev"],
        dims=(dims["y"], dims["x"]),
    ).astype("float32")
   
    idx = names.index("To-sensor azimuth (0 to 360 degrees CW from N)")
    data["viewing_azimuth_angle"] = xr.DataArray(
        data=l1bobs["obs"][:,:,idx],
        dims=(dims["y"], dims["x"]),
    ).astype("float32")
    
    idx = names.index("To-sensor zenith (0 to 90 degrees from zenith)")
    data["viewing_zenith_angle"] = xr.DataArray(
        data=l1bobs["obs"][:,:,idx],
        dims=(dims["y"], dims["x"]),
    ).astype("float32")
    
    idx = names.index("To-sun azimuth (0 to 360 degrees CW from N)")
    data["solar_azimuth_angle"] = xr.DataArray(
        data=l1bobs["obs"][:,:,idx],
        dims=(dims["y"], dims["x"]),
    ).astype("float32")

    idx = names.index("To-sun zenith (0 to 90 degrees from zenith)")
    data["solar_zenith_angle"] = xr.DataArray(
        data=l1bobs["obs"][:,:,idx],
        dims=(dims["y"], dims["x"]),
    ).astype("float32")
    
    band1_data["wavelength"] = xr.DataArray(
        data=l1brad["sensor_band_parameters"]["wavelengths"],
        dims=(dims["z"]),
    ).astype("float32")

    band1_data["fwhm"] = xr.DataArray(
        data=l1brad["sensor_band_parameters"]["fwhm"],
        dims=(dims["z"]),
    ).astype("float32")

    radiance, radiance_noise = get_radiance(l1brad, snr)
    band1_data["radiance"] = xr.DataArray(
        data=radiance,
        dims=(dims["y"], dims["x"],dims["z"]),
    ).astype("float32")

    band1_data["radiance_noise"] = xr.DataArray(
        data=radiance_noise,
        dims=(dims["y"], dims["x"],dims["z"]),
    ).astype("float32")

    band1_data["flat_field_update"] = xr.DataArray(
        data=l1brad["flat_field_update"],
        dims=(dims["x"],dims["z"]),
    ).astype("float32")
    
    return data, band1_data

def get_time(l1b):
    # Calculate the  observation times for the individual frames using start and
    # stop times of the tile and linearly interpolate between all
    # downtrack pixels. Write the start time into units string.
    
    start_time = l1b.ds.attrs.get("time_coverage_start")
    stop_time  = l1b.ds.attrs.get("time_coverage_end")
    start_time = datetime.strptime(start_time, "%Y-%m-%dT%H:%M:%S%z")
    stop_time  = datetime.strptime(stop_time, "%Y-%m-%dT%H:%M:%S%z")
    Ntemporal  = l1b.sizes["downtrack"]
    
    reference_time = datetime.strftime(start_time, "%Y-%m-%d %H:%M:%S%z")
    time_units_string = f"seconds since {reference_time}"
    
    time = []
    time_diff_total = \
        timedelta(seconds=(stop_time - start_time).total_seconds())
    
    for i in range(Ntemporal):
        time_measurement = start_time + i/(Ntemporal - 1) * time_diff_total
        time.append(time_measurement)

    time = date2num(time, time_units_string)
    
    return time, time_units_string


def get_radiance(l1b, snr):
    # Get radiance and wavelength from the L1B file.
    # Convert radiance from microW cm-2 sr-1 nm-1 to photons s-1 cm-2 sr-1 nm-1.
    # Calculate radiance noise under shot-noise assumption based on
    # SNR for the brightest measurement.

    wavelength = l1b["sensor_band_parameters"]["wavelengths"]
    radiance_microWatts   = l1b["radiance"] # in microW cm-2 sr-1 nm-1
    
    # convert units
    # microW cm-2 sr-1 nm-1 -> photons s-1 cm-2 sr-1 nm-1
    radiance_photons = radiance_microWatts * 1e-6 * wavelength * 1e-9 / (6.62607015e-34 * 299792458)
    
    # Calculate radiance noise under shot-noise assumption based on
    # SNR for the brightest measurement:
    #
    # noise_pixel = sqrt(rad_pixel/rad_maxI) * noise_maxI
    # noise_maxI  = rad_maxI / SNR
    # ============================
    # noise_pixel = sqrt(rad_pixel) * sqrt(rad_maxI) / SNR
    #
    # For SNR, see Thompson et al., RSE, 2024,
    # https://doi.org/10.1016/j.rse.2023.113986).
    # SNR is user-provided in the config file.
    radiance_masked = radiance_photons.where(radiance_photons >= 0, other=np.nan)
    radiance_noise = np.sqrt(radiance_masked) * np.sqrt(radiance_masked.max()) / float(snr)
   
    # Replace small noise values by noise floor
    threshold = 0.1 * radiance_masked.max() / float(snr)
    radiance_noise = radiance_noise.fillna(threshold)
    radiance_noise = radiance_noise.where(radiance_noise >= threshold, other=threshold)
    
    return radiance_photons, radiance_noise

