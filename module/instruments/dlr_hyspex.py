import numpy as np
import xarray as xr
from datetime import datetime
from cftime import date2num
import os


def import_data(config, dims):
    input_data = xr.open_dataset(os.path.join(config["path"], config["l1b"]))

    dlr_hyspex_data = xr.Dataset()

    time, time_units_string = get_time(input_data)
    dlr_hyspex_data["time"] = xr.DataArray(
        data=time,
        dims=(dims["y"]),
    ).astype("float32")
    dlr_hyspex_data.time.attrs["units"] = time_units_string

    dlr_hyspex_data["latitude"] = xr.DataArray(
        data=get_latitude(input_data),
        dims=(dims["y"], dims["x"]),
    ).astype("float32")

    dlr_hyspex_data["longitude"] = xr.DataArray(
        data=get_longitude(input_data),
        dims=(dims["y"], dims["x"]),
    ).astype("float32")

    dlr_hyspex_data["solar_zenith_angle"] = xr.DataArray(
        data=get_sza(input_data),
        dims=(dims["y"], dims["x"]),
    ).astype("float32")

    dlr_hyspex_data["solar_azimuth_angle"] = xr.DataArray(
        data=get_saa(input_data),
        dims=(dims["y"], dims["x"])
    ).astype("float32")

    dlr_hyspex_data["viewing_zenith_angle"] = xr.DataArray(
        data=get_vza(input_data),
        dims=(dims["y"], dims["x"]),
    ).astype("float32")

    dlr_hyspex_data["viewing_azimuth_angle"] = xr.DataArray(
        data=get_vaa(input_data),
        dims=(dims["y"], dims["x"]),
    ).astype("float32")

    dlr_hyspex_data["observer_altitude"] = xr.DataArray(
        data=get_z(input_data),
        dims=(dims["y"], dims["x"]),
    ).astype("float32")

    band1_data = xr.Dataset()

    wavelength, radiance, radiance_noise = get_spectrum(input_data)

    band1_data["wavelength"] = xr.DataArray(
        data=wavelength,
        dims=(dims["z"]),
    ).astype("float32")

    band1_data["radiance"] = xr.DataArray(
        data=radiance,
        dims=(dims["y"], dims["x"], dims["z"]),
    ).astype("float32")

    band1_data["radiance_noise"] = xr.DataArray(
        data=radiance_noise,
        dims=(dims["y"], dims["x"], dims["z"]),
    ).astype("float32")

    print("TODO LS: Maybe read surface elevation data?")

    band_list = [band1_data]
    input_file_list = [config["l1b"]]

    return dlr_hyspex_data, band_list, input_file_list


def get_time(input_data):
    time = input_data.time.values
    start_time = np.datetime_as_string(time[0])[:-3]
    start_time = datetime.strptime(start_time, "%Y-%m-%dT%H:%M:%S.%f")

    reference_time = datetime.strftime(start_time, "%Y-%m-%d %H:%M:%S.%f")
    time_units_string = f"seconds since {reference_time}"

    out_time = []
    for measurement_time in time:
        measurement_time = np.datetime_as_string(measurement_time)[:-3]
        measurement_time = \
            datetime.strptime(measurement_time, "%Y-%m-%dT%H:%M:%S.%f")
        out_time.append(measurement_time)

    out_time = date2num(out_time, time_units_string)

    return out_time, time_units_string


def get_latitude(input_data):
    return input_data.pixel_center_lat.values


def get_longitude(input_data):
    return input_data.pixel_center_lon.values


def get_sza(input_data):
    # sza is only provided along time axis.
    # RemoTeC needs it for each pixel. Copy the values along the spatial axis.
    Ntemporal = len(input_data.time)
    Nspatial = len(input_data.across_track_pixel_index)
    sza = np.array(input_data.solar_zenith_ang.values)
    sza = np.broadcast_to(sza, shape=(Nspatial, Ntemporal)).transpose()
    return sza


def get_saa(input_data):
    # saa is only provided along time axis.
    # RemoTeC needs it for each pixel. Copy the values along the spatial axis.
    Ntemporal = len(input_data.time)
    Nspatial = len(input_data.across_track_pixel_index)
    saa = np.array(input_data.solar_azimuth_ang.values)
    saa = np.broadcast_to(saa, shape=(Nspatial, Ntemporal)).transpose()
    return saa


def get_vza(input_data):
    return input_data.pixel_zenith_ang.values


def get_vaa(input_data):
    return input_data.pixel_azimuth_ang.values


def get_z(input_data):
    # height above reference ellipsoid
    hae = input_data.sensor_height
    # surface height above reference ellipsoid
    sfc = input_data.pixel_center_height
    # sensor height above surface
    z = hae - sfc
    return z.values


def get_spectrum(input_data):
    # Wavelength in nm
    wavelength = input_data.wavelength.mean(dim="across_track_pixel_index")

    # Radiance in mW m-2 sr-1 nm-1
    radiance = input_data.radiance

    # Convert units
    # mW m-2 sr-1 nm-1 -> photons s-1 cm-2 sr-1 nm-1
    planck_constant = 6.62607015e-34  # J s
    light_speed = 299792458  # m s-1
    radiance = \
        radiance * 1e-4 * wavelength * 1e-9 \
        / planck_constant / light_speed

    # Calculate radiance noise from radiance using signal-to-noise ratio
    print("TODO LS: I AM MAKING UP AN SNR FOR NOW. CHANGE THIS !!!!")
    snr = 100/1
    radiance_noise = radiance / snr

    return wavelength.values, radiance.values, radiance_noise.values
