import h5py
import zipfile
from io import BytesIO
import xarray as xr
import numpy as np
import sys
import os
from datetime import datetime
from datetime import timedelta

def import_data(config, dims):
    l1b = get_data(os.path.join(config["path"], config["l1b"]))
    l2b = get_data(os.path.join(config["path"], config["l2b"]))

    prisma_data = xr.Dataset()

    prisma_data["time"] = (
        (dims["y"]), get_time(l1b)
    )

    prisma_data["latitude"] = (
        (dims["y"], dims["x"]), get_latitude(l2b)
    )

    prisma_data["longitude"] = (
        (dims["y"], dims["x"]), get_longitude(l2b)
    )

    prisma_data["solar_zenith_angle"] = (
        (dims["y"], dims["x"]), get_sza(l2b)
    )

    prisma_data["viewing_zenith_angle"] = (
        (dims["y"], dims["x"]), get_vza(l2b)
    )

    print("TODO LS: Maybe get surface elevation?")
    # for surface elevation maybe check l2 product 8.2.1 HEIGHT_OFF ??
    # but maybe also page before or after, not exactly sure. This should be
    # provided starting from the L2B product.

    band1_data = xr.Dataset()

    wavelength, radiance, radiance_noise = get_spectrum(l1b, "swir")

    band1_data["wavelength"] = (
        (dims["z"]), wavelength
    )

    band1_data["radiance"] = (
        (dims["y"], dims["x"], dims["z"]), radiance
    )

    band1_data["radiance_noise"] = (
        (dims["y"], dims["x"], dims["z"]), radiance_noise
    )

    band2_data = xr.Dataset()

    wavelength, radiance, radiance_noise = get_spectrum(l1b, "vnir")

    band2_data["wavelength"] = (
        (dims["z"]), wavelength
    )

    band2_data["radiance"] = (
        (dims["y"], dims["x"], dims["z"]), radiance
    )

    band2_data["radiance_noise"] = (
        (dims["y"], dims["x"], dims["z"]), radiance_noise
    )

    band_list = [band1_data, band2_data]
    input_file_list = [config["l1b"], config["l2b"]]

    return prisma_data, band_list, input_file_list


def get_data(zip_file_path):
    zip_contents = zipfile.ZipFile(zip_file_path)
    # assume .zip folder contains only the he5 file and nothing else
    he5_contents = zip_contents.open(zip_contents.namelist()[0]).read()
    data = h5py.File(BytesIO(he5_contents), mode="r")
    return data


def get_time(l1b):
    # check documentation p.130 at the very bottom. This field should give
    # the UTC time in MJD2000 Decimal days format, but online converters
    # get a date in 1987 or something like that. Something seems to be wrong.
    # Instead, calculate the times for the individual frames using start and
    # stop times of the measurements and linearly interpolate between all
    # along-track pixels.
    # You could get size of (temporal) using the time field and get its size
    # even if you can't use the values.
    # time = l1b["/HDFEOS/SWATHS/PRS_L1_HCO/Geolocation Fields/Time"]

    start_time = l1b.attrs["Product_StartTime"].decode("utf-8")
    stop_time = l1b.attrs["Product_StopTime"].decode("utf-8")
    start_time = datetime.strptime(start_time, "%Y-%m-%dT%H:%M:%S.%f")
    stop_time = datetime.strptime(stop_time, "%Y-%m-%dT%H:%M:%S.%f")
    Ntemporal = len(l1b["/HDFEOS/SWATHS/PRS_L1_HCO/Geolocation Fields/Time"])

    time_diff = timedelta(seconds=(stop_time - start_time).total_seconds())

    time = []
    for i in range(Ntemporal):
        timestamp = start_time + i/(Ntemporal - 1) * time_diff
        time_string = datetime.strftime(timestamp, "%Y-%m-%dT%H:%M:%SZ")
        time.append(time_string)

    return time


def get_latitude(l2b):
    latitude = np.array(
        l2b["/HDFEOS/SWATHS/PRS_L2B_HCO/Geolocation Fields/Latitude"]
    )
    # (spatial, temporal) -> (temporal, spatial)
    latitude = np.transpose(latitude)
    return latitude


def get_longitude(l2b):
    longitude = np.array(
        l2b["/HDFEOS/SWATHS/PRS_L2B_HCO/Geolocation Fields/Longitude"]
    )
    # (spatial, temporal) -> (temporal, spatial)
    longitude = np.transpose(longitude)
    return longitude


def get_sza(l2b):
    sza = np.array(
        l2b["/HDFEOS/SWATHS/PRS_L2B_HCO/Geometric Fields/Solar_Zenith_Angle"]
    )
    # (spatial, temporal) -> (temporal, spatial)
    sza = np.transpose(sza)
    return sza


def get_vza(l2b):
    vza = np.array(
        l2b["/HDFEOS/SWATHS/PRS_L2B_HCO/Geometric Fields/Observing_Angle"]
    )
    # (spatial, temporal) -> (temporal, spatial)
    vza = np.transpose(vza)
    return vza


def get_spectrum(l1b, spectral_domain):
    # get wavelength with dimension (spectral)
    # wavelength in nm
    wavelength = np.array(
        l1b.attrs[f"List_Cw_{spectral_domain.capitalize()}"]
    )

    # get radiance
    digital_number = np.array(
        l1b["/HDFEOS/SWATHS/PRS_L1_HCO/Data Fields/"
            + f"{spectral_domain.upper()}_Cube"]
    )
    scale_factor = l1b.attrs[f"ScaleFactor_{spectral_domain.capitalize()}"]
    offset = l1b.attrs[f"Offset_{spectral_domain.capitalize()}"]

    # calculate radiance using PRISMA product specification p.107
    # radiance in W m-2 sr-1 um-1
    radiance = digital_number / scale_factor - offset

    # (spatial, spectral, temporal) -> (temporal, spatial, spectral)
    radiance = np.moveaxis(radiance, [0, 1, 2], [1, 2, 0])

    # Get radiance noise with dimension (temporal, spatial, spectral)
    # Calculate radiance noise using SNR provided in PRISMA product
    # specification p.15
    if spectral_domain == "swir":
        snr = 100/1
    elif spectral_domain == "vnir":
        snr = 200/1
    else:
        sys.exit("invalid spectral domain")
    radiance_noise = radiance / snr

    # postprocess
    # bands for which measurement failed have wavelength of zero
    # the radiance will also be zero for these bands
    invalid_spectral_indices = np.where(wavelength == 0)
    wavelength = np.delete(
        wavelength, invalid_spectral_indices)
    radiance = np.delete(
        radiance, invalid_spectral_indices, axis=2)
    radiance_noise = np.delete(
        radiance_noise, invalid_spectral_indices, axis=2)

    # wavelength dimension may be flipped.
    sort_indices = np.argsort(wavelength)
    wavelength = np.array(wavelength)[sort_indices]
    radiance = np.array(radiance)[:, :, sort_indices]
    radiance_noise = np.array(radiance_noise)[:, :, sort_indices]

    # convert units
    # radiance, radiance_noise
    # W m-2 sr-1 um-1 -> photons s-1 cm-2 sr-1 nm-1
    planck_constant = 6.62607015e-34  # J s
    light_speed = 299792458  # m s-1
    radiance = \
        radiance * 1e-7 * wavelength * 1e-9 \
        / planck_constant / light_speed
    radiance_noise = \
        radiance_noise * 1e-7 * wavelength * 1e-9 \
        / planck_constant / light_speed

    return wavelength, radiance, radiance_noise
