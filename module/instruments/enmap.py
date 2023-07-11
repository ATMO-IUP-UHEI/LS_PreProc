import tarfile
import zipfile
from io import BytesIO
import xml.etree.ElementTree as ET
import numpy as np
import xarray as xr
import sys
import os
from datetime import datetime
from datetime import timedelta


def import_data(config, dims):
    meta, vnir, swir = get_data(os.path.join(config["path"], config["tar_gz"]))

    enmap_data = xr.Dataset()

    enmap_data["time"] = xr.DataArray(
        data=get_time(meta),
        dims=(dims["y"]),
    )

    enmap_data["latitude"] = xr.DataArray(
        data=get_latitude(meta),
        dims=(dims["y"], dims["x"]),
    ).astype("float32")

    enmap_data["longitude"] = xr.DataArray(
        data=get_longitude(meta),
        dims=(dims["y"], dims["x"]),
    ).astype("float32")

    enmap_data["solar_zenith_angle"] = xr.DataArray(
        data=get_sza(meta),
        dims=(dims["y"], dims["x"]),
    ).astype("float32")

    enmap_data["solar_azimuth_angle"] = xr.DataArray(
        data=get_saa(meta),
        dims=(dims["y"], dims["x"]),
    ).astype("float32")

    enmap_data["viewing_zenith_angle"] = xr.DataArray(
        data=get_vza(meta),
        dims=(dims["y"], dims["x"]),
    ).astype("float32")

    enmap_data["viewing_azimuth_angle"] = xr.DataArray(
        data=get_vaa(meta),
        dims=(dims["y"], dims["x"]),
    ).astype("float32")

    print("TODO LS: Maybe get surface elevation?")

    band1_data = xr.Dataset()

    wavelength, radiance, radiance_noise = get_spectrum(meta, swir, "swir")

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

    band2_data = xr.Dataset()

    wavelength, radiance, radiance_noise = get_spectrum(meta, vnir, "vnir")

    band2_data["wavelength"] = xr.DataArray(
        data=wavelength,
        dims=(dims["z"]),
    ).astype("float32")

    band2_data["radiance"] = xr.DataArray(
        data=radiance,
        dims=(dims["y"], dims["x"], dims["z"]),
    ).astype("float32")

    band2_data["radiance_noise"] = xr.DataArray(
        data=radiance_noise,
        dims=(dims["y"], dims["x"], dims["z"]),
    ).astype("float32")

    band_list = [band1_data, band2_data]
    input_file_list = [config["tar_gz"]]

    return enmap_data, band_list, input_file_list


def get_data(tar_gz_file_path):
    tar = tarfile.open(tar_gz_file_path)

    for content in tar:
        # Find measurement. It is contained in the only zip file in the archive
        if not content.name.endswith(".ZIP"):
            continue

        # Get the zip file as a binary string and put it into memory,
        # acting as if it were a file.
        f = tar.extractfile(content)
        content = BytesIO(f.read())
        break

    # Open the .zip file which was previously inside of the .tar.gz file
    zip = zipfile.ZipFile(content)
    for member in zip.namelist():
        # Get metadata and spectral information from the desired channel(s).
        # They are contained in an .xml file and .tif file(s).
        # Get them from the .zip file and put into memory, acting as if they
        # were files.
        if "METADATA" in member and member.endswith(".XML"):
            metadata = BytesIO(zip.open(member).read())
            metadata = ET.parse(metadata).getroot()
            continue
        if "SPECTRAL_IMAGE_SWIR" in member and member.endswith(".TIF"):
            swir_band_data = BytesIO(zip.open(member).read())
            swir_band_data = xr.open_dataset(swir_band_data, engine="rasterio")
            continue
        if "SPECTRAL_IMAGE_VNIR" in member and member.endswith(".TIF"):
            vnir_band_data = BytesIO(zip.open(member).read())
            vnir_band_data = xr.open_dataset(vnir_band_data, engine="rasterio")
            continue

    return metadata, vnir_band_data, swir_band_data


def get_time(meta):
    start_time = meta.find("./base/temporalCoverage/startTime").text
    stop_time = meta.find("./base/temporalCoverage/stopTime").text
    start_time = datetime.strptime(start_time, "%Y-%m-%dT%H:%M:%S.%fZ")
    stop_time = datetime.strptime(stop_time, "%Y-%m-%dT%H:%M:%S.%fZ")
    Ntemporal = int(meta.find("./specific/heightOfScene").text)

    time_diff = timedelta(seconds=(stop_time - start_time).total_seconds())

    time = []
    for i in range(Ntemporal):
        timestamp = start_time + i/(Ntemporal - 1) * time_diff
        time_string = datetime.strftime(timestamp, "%Y-%m-%dT%H:%M:%SZ")
        time.append(time_string)

    return time


def get_latitude(meta):
    latitude_path = "./base/spatialCoverage/boundingPolygon/point/latitude"
    upper_left, lower_left, lower_right, upper_right = \
        get_pointinfo_corners(meta, latitude_path)

    latitude = interpolate_corners_to_grid(
        meta, upper_left, lower_left, lower_right, upper_right
    )

    return latitude


def get_longitude(meta):
    longitude_path = "./base/spatialCoverage/boundingPolygon/point/longitude"
    upper_left, lower_left, lower_right, upper_right = \
        get_pointinfo_corners(meta, longitude_path)

    longitude = interpolate_corners_to_grid(
        meta, upper_left, lower_left, lower_right, upper_right
    )

    return longitude


def get_sza(meta):
    sza_path = "./specific/sunElevationAngle"
    upper_left, lower_left, lower_right, upper_right = \
        get_angle_corners(meta, sza_path)

    sza = interpolate_corners_to_grid(
        meta, upper_left, lower_left, lower_right, upper_right
    )

    sza = 90 - sza

    return sza


def get_saa(meta):
    saa_path = "./specific/sunAzimuthAngle"
    upper_left, lower_left, lower_right, upper_right = \
        get_angle_corners(meta, saa_path)

    saa = interpolate_corners_to_grid(
        meta, upper_left, lower_left, lower_right, upper_right
    )

    return saa


def get_vza(meta):
    vza_path = "./specific/acrossOffNadirAngle"
    upper_left, lower_left, lower_right, upper_right = \
        get_angle_corners(meta, vza_path)

    vza = interpolate_corners_to_grid(
        meta, upper_left, lower_left, lower_right, upper_right
    )

    return vza


def get_vaa(meta):
    vaa_path = "./specific/sceneAzimuthAngle"
    upper_left, lower_left, lower_right, upper_right = \
        get_angle_corners(meta, vaa_path)

    vaa = interpolate_corners_to_grid(
        meta, upper_left, lower_left, lower_right, upper_right
    )

    return vaa


def get_pointinfo_corners(meta, pointinfo_path):
    corners = meta.findall(pointinfo_path)

    upper_left = float(corners[0].text)
    lower_left = float(corners[1].text)
    lower_right = float(corners[2].text)
    upper_right = float(corners[3].text)

    return upper_left, lower_left, lower_right, upper_right


def get_angle_corners(meta, angle_path):
    upper_left = float(
        meta.find(f"{angle_path}/upper_left").text
    )
    lower_left = float(
        meta.find(f"{angle_path}/lower_left").text
    )
    lower_right = float(
        meta.find(f"{angle_path}/lower_right").text
    )
    upper_right = float(
        meta.find(f"{angle_path}/upper_right").text
    )
    return upper_left, lower_left, lower_right, upper_right


def interpolate_corners_to_grid(
        meta, upper_left, lower_left, lower_right, upper_right):
    Ntemporal = int(meta.find("./specific/heightOfScene").text)
    Nspatial = int(meta.find("./specific/widthOfScene").text)

    da = xr.DataArray(
        data=[[lower_left, upper_left],
              [lower_right, upper_right]],
        dims=("x", "y"),
        coords={"x": [0, Nspatial - 1], "y": [0, Ntemporal - 1]},
    )

    da = da.interp(
        x=np.linspace(0, Nspatial - 1, Nspatial),
        y=np.linspace(0, Ntemporal - 1, Ntemporal)
    )

    da = da.transpose()

    return da.values


def get_spectrum(meta, l1b, spectral_domain):
    # bandIDs contain VNIR bands first, SWIR bands second.
    # Nswir can be analogously gathered from the metadata, but this script
    # just checks if the bandID is larger that Nvnir to see if the data is
    # in the SWIR range.
    Nvnir = int(meta.find("./specific/numberOfVNIRBands").text)

    # get wavelength in nm
    wavelength = np.array([])
    fwhm = np.array([])
    gain = np.array([])
    offset = np.array([])
    for band in meta.findall("./specific/bandCharacterisation/bandID"):
        # only first Nvnir bands contain information relevant for vnir spectral
        # domain
        if spectral_domain == "vnir" and float(band.attrib["number"]) > Nvnir:
            break
        # only last Nswir bands contain information relevant for swir spectral
        # domain
        if spectral_domain == "swir" and float(band.attrib["number"]) <= Nvnir:
            continue
        wavelength = np.append(
            wavelength, float(band.find("wavelengthCenterOfBand").text)
        )
        fwhm = np.append(
            fwhm, float(band.find("FWHMOfBand").text)
        )
        gain = np.append(
            gain, float(band.find("GainOfBand").text)
        )
        offset = np.append(
            offset, float(band.find("OffsetOfBand").text)
        )

    digital_number = np.array(
        l1b.band_data
    )

    # Calculate radiance using EnMAP product specification p.14
    # It doesn't say in what way the gain and offset have to be applied to the
    # digital signal, but I am following along with Luis Guanter's read
    # routine.
    # Furthermore, Section 4 mentions a background value, not sure what that
    # is for.
    # radiance in W m-2 sr-1 nm-1
    radiance = \
        digital_number[:, :, :] * gain[:, None, None] + offset[:, None, None]

    # (spectral, temporal, spatial) -> (temporal, spatial, spectral)
    radiance = np.moveaxis(radiance, [0, 1, 2], [2, 0, 1])

    # Get radiance noise with dimension (temporal, spatial, spectral)
    # Calculate radiance noise using SNR provided in EnMAP Instrument
    # Specification
    if spectral_domain == "swir":
        snr = 150/1
    elif spectral_domain == "vnir":
        snr = 500/1
    else:
        sys.exit("invalid spectral domain")
    radiance_noise = radiance / snr

    # convert units
    # radiance, radiance_noise
    # W m-2 sr-1 nm-1 -> photons s-1 cm-2 sr-1 nm-1
    planck_constant = 6.62607015e-34  # J s
    light_speed = 299792458  # m s-1
    radiance = \
        radiance * 1e-4 * wavelength * 1e-9 \
        / planck_constant / light_speed
    radiance_noise = \
        radiance_noise * 1e-4 * wavelength * 1e-9 \
        / planck_constant / light_speed

    return wavelength, radiance, radiance_noise
