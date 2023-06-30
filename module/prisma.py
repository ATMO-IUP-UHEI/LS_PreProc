import h5py
import zipfile
from io import BytesIO
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import sys


def main():
    path = "/home/lscheidw/phd/RemoTeC_LS/data/tmp_preproc/l1b/prisma"
    id = "20201010072033_20201010072038"

    l1b = get_data(f"{path}/PRS_L1_STD_OFFL_{id}_0001.zip")
    l2b = get_data(f"{path}/PRS_L2B_STD_{id}_0001.zip")

    temporal = "sample"
    spatial = "line"
    spectral = "band"
    dims_1d = (spectral)
    dims_2d = (temporal, spatial)
    dims_3d = (temporal, spatial, spectral)

    prisma_data = xr.Dataset()

    prisma_data["latitude"] = (dims_2d, get_latitude(l2b))
    prisma_data.latitude.attrs["standard_name"] = "Latitude at pixel center"
    prisma_data.latitude.attrs["units"] = "degrees north"

    prisma_data["longitude"] = (dims_2d, get_longitude(l2b))
    prisma_data.longitude.attrs["standard_name"] = "Longitude at pixel center"
    prisma_data.longitude.attrs["units"] = "degrees east"

    prisma_data["sza"] = (dims_2d, get_sza(l2b))
    prisma_data.sza.attrs["standard_name"] = "Solar Zenith Angle"
    prisma_data.sza.attrs["units"] = "degrees"

    prisma_data["vza"] = (dims_2d, get_vza(l2b))
    prisma_data.vza.attrs["standard_name"] = "Viewing Zenith Angle"
    prisma_data.vza.attrs["units"] = "degrees"

    print("TODO LS: Get time.")

    print("TODO LS: Maybe get surface elevation?")
    # for surface elevation maybe check l2 product 8.2.1 HEIGHT_OFF ??
    # but maybe also page before or after, not exactly sure. This should be
    # provided starting from the L2B product.

    band1_data = xr.Dataset()

    wavelength, radiance, radiance_noise = get_spectrum(l1b, "swir")

    band1_data["wavelength"] = (dims_1d, wavelength)
    band1_data.wavelength.attrs["standard_name"] = "Wavelength"
    band1_data.wavelength.attrs["units"] = "nm"

    band1_data["radiance"] = (dims_3d, radiance)
    band1_data.radiance.attrs["standard_name"] = "At-sensor radiance"
    band1_data.radiance.attrs["units"] = "photons s-1 cm-2 sr-1 nm-1"

    band1_data["radiance_noise"] = (dims_3d, radiance_noise)
    band1_data.radiance_noise.attrs["standard_name"] = \
        "Noise of at-sensor radiance"
    band1_data.radiance_noise.attrs["units"] = "photons s-1 cm-2 sr-1 nm-1"

    band2_data = xr.Dataset()

    wavelength, radiance, radiance_noise = get_spectrum(l1b, "vnir")

    band2_data["wavelength"] = (dims_1d, wavelength)
    band2_data.wavelength.attrs["standard_name"] = "Wavelength"
    band2_data.wavelength.attrs["units"] = "nm"

    band2_data["radiance"] = (dims_3d, radiance)
    band2_data.radiance.attrs["standard_name"] = "At-sensor radiance"
    band2_data.radiance.attrs["units"] = "photons s-1 cm-2 sr-1 nm-1"

    band2_data["radiance_noise"] = (dims_3d, radiance_noise)
    band2_data.radiance_noise.attrs["standard_name"] = \
        "Noise of at-sensor radiance"
    band2_data.radiance_noise.attrs["units"] = "photons s-1 cm-2 sr-1 nm-1"

    for var in prisma_data.data_vars:
        prisma_data[var].encoding.update({"_FillValue": None})
    for var in band1_data.data_vars:
        band1_data[var].encoding.update({"_FillValue": None})
    for var in band2_data.data_vars:
        band2_data[var].encoding.update({"_FillValue": None})

    prisma_data.attrs["history"] = \
        "Created using the L1B preprocessor for RemoTeC for PRISMA data."\
        + f" Raw files used:  PRS_L1_STD_OFFL_{id}_0001.zip" \
        + f" and PRS_L2B_STD_{id}_0001.zip"

    prisma_data.to_netcdf("L1B_prisma.nc", mode="w", format="NETCDF4")
    band1_data.to_netcdf("L1B_prisma.nc", mode="a", format="NETCDF4",
                         group="BAND01")
    band2_data.to_netcdf("L1B_prisma.nc", mode="a", format="NETCDF4",
                         group="BAND02")

    return 0


def get_data(zip_file_path):
    zip_contents = zipfile.ZipFile(zip_file_path)
    # assume .zip folder contains only the he5 file and nothing else
    he5_contents = zip_contents.open(zip_contents.namelist()[0]).read()
    data = h5py.File(BytesIO(he5_contents), mode="r")
    return data


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


def get_datetime(l1b):
    # check documentation p.130 at the very bottom. This field should give
    # the UTC time in MJD2000 Decimal days format, but online converters
    # get a date in 1987 or something like that. Something seems to be wrong.
    # Instead, calculate the times for the individual frames using start and
    # stop times of the measurements and linearly interpolate between all
    # along-track pixels.
    # You could get size of (temporal) using the time field and get its size
    # even if you can't use the values.
    # Number of lines see 7.2.2, p.96 of documentation, should be in the data.
    # time = l1b["/HDFEOS/SWATHS/PRS_L1_HCO/Geolocation Fields/Time"]

    start_time = l1b.attrs["Product_StartTime"]
    stop_time = l1b.attrs["Product_StopTime"]

    print("TODO LS: Create array with size (temporal).")
    return [start_time, stop_time]


def get_spectrum(l1b, spectral_range):
    # get wavelength with dimension (spectral)
    # wavelength in nm
    wavelength = np.array(
        l1b.attrs[f"List_Cw_{spectral_range.capitalize()}"]
    )

    # get radiance
    digital_number = np.array(
        l1b["/HDFEOS/SWATHS/PRS_L1_HCO/Data Fields/"
            + f"{spectral_range.upper()}_Cube"]
    )
    scale_factor = l1b.attrs[f"ScaleFactor_{spectral_range.capitalize()}"]
    offset = l1b.attrs[f"Offset_{spectral_range.capitalize()}"]

    # calculate radiance using PRISMA product specification p.107
    # radiance in W m-2 sr-1 um-1
    radiance = digital_number / scale_factor - offset

    # (spatial, spectral, temporal) -> (temporal, spatial, spectral)
    radiance = np.moveaxis(radiance, [0, 1, 2], [1, 2, 0])

    # get radiance noise with dimension (temporal, spatial, spectral)
    # calculate radiance_noise using SNR provided in PRISMA product
    # specification p.15
    if spectral_range == "swir":
        snr = 100/1
    elif spectral_range == "vnir":
        snr = 200/1
    else:
        sys.exit("invalid spectral range")
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


if __name__ == "__main__":
    main()
