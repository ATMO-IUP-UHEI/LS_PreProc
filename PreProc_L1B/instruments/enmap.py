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
from cftime import date2num


def import_data(dims):
    print(datetime.now(), "getting enmap data.")
    tile_list = []
    for file in os.listdir("tmp/spectra"):
        if file.endswith(".tar.gz"):
            tile_list.append(file)

    input_file_list = []

    for i, tile in enumerate(tile_list):
        input_file_list.append(tile)
        tar_gz_file_path = os.path.join("tmp/spectra", tile)

        tile_root, tile_band1, tile_band2 = import_tile(tar_gz_file_path, dims)

        if i == 0:
            root = tile_root
            band1 = tile_band1
            band2 = tile_band2
        else:
            root = xr.concat(
                (root, tile_root), dim="frame", data_vars="minimal")
            band1 = xr.concat(
                (band1, tile_band1), dim="frame", data_vars="minimal")
            band2 = xr.concat(
                (band2, tile_band2), dim="frame", data_vars="minimal")

    # sort by latitude
    sort_array = root.latitude.mean(dim="line")
    sort_direction = root.time[1] > root.time[0]
    root = root.sortby(sort_array, ascending=sort_direction)
    band1 = band1.sortby(sort_array, ascending=sort_direction)
    band2 = band2.sortby(sort_array, ascending=sort_direction)

    # shift swir channel
    band1["radiance"], band1["radiance_noise"] = \
        shift_swir_channel(band1.radiance, band1.radiance_noise)

    return root, [band1, band2], input_file_list


def import_tile(tar_gz_file_path, dims):
    meta, vnir, swir = get_data(tar_gz_file_path)

    enmap_data = xr.Dataset()

    time, time_units_string = get_time(meta)
    enmap_data["time"] = xr.DataArray(
        data=time,
        dims=(dims["y"]),
    ).astype("float32")
    enmap_data.time.attrs["units"] = time_units_string

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

    band1_data = xr.Dataset()

    wavelength, fwhm, radiance, radiance_noise = get_spectrum(
        meta, swir, "swir")

    band1_data["wavelength"] = xr.DataArray(
        data=wavelength,
        dims=(dims["z"]),
    ).astype("float32")

    # band1_data["fwhm"] = xr.DataArray(
    #     data=fwhm,
    #     dims=(dims["z"]),
    # ).astype("float32")

    band1_data["radiance"] = xr.DataArray(
        data=radiance,
        dims=(dims["y"], dims["x"], dims["z"]),
    ).astype("float32")

    band1_data["radiance_noise"] = xr.DataArray(
        data=radiance_noise,
        dims=(dims["y"], dims["x"], dims["z"]),
    ).astype("float32")

    band2_data = xr.Dataset()

    wavelength, fwhm, radiance, radiance_noise = get_spectrum(
        meta, vnir, "vnir")

    band2_data["wavelength"] = xr.DataArray(
        data=wavelength,
        dims=(dims["z"]),
    ).astype("float32")

    # band2_data["fwhm"] = xr.DataArray(
    #     data=fwhm,
    #     dims=(dims["z"]),
    # ).astype("float32")

    band2_data["radiance"] = xr.DataArray(
        data=radiance,
        dims=(dims["y"], dims["x"], dims["z"]),
    ).astype("float32")

    band2_data["radiance_noise"] = xr.DataArray(
        data=radiance_noise,
        dims=(dims["y"], dims["x"], dims["z"]),
    ).astype("float32")

    return enmap_data, band1_data, band2_data


def get_data(tar_gz_file_path):
    tar = tarfile.open(tar_gz_file_path)

    zip_file_number = 0
    for content in tar:
        if content.name.endswith(".ZIP"):
            zip_file = content
            zip_file_number += 1

    if zip_file_number != 1:
        sys.exit(".tar.gz must contain exactly 1 .zip file. Aborting... "
                 + "(the folder in which the .tar.gz file is saved should "
                 + "contain a script to split it into multiple tar files with "
                 + "one tile each).")

    # The measurement is contained in the only zip file in the archive
    # Get the zip file as a binary string and put it into memory,
    # acting as if it were a file.
    f = tar.extractfile(zip_file)
    content = BytesIO(f.read())

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

    reference_time = datetime.strftime(start_time, "%Y-%m-%d %H:%M:%S.%f")
    time_units_string = f"seconds since {reference_time}"

    time = []
    time_diff_total = \
        timedelta(seconds=(stop_time - start_time).total_seconds())
    for i in range(Ntemporal):
        time_measurement = start_time + i/(Ntemporal - 1) * time_diff_total
        time.append(time_measurement)

    time = date2num(time, time_units_string)

    return time, time_units_string


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

    orbit_direction = meta.find("./specific/orbitDirection").text

    if orbit_direction == "ASCENDING":
        data = [[lower_left, upper_left], [lower_right, upper_right]]
    elif orbit_direction == "DESCENDING":
        data = [[upper_left, lower_left], [upper_right, lower_right]]
    else:
        sys.exit(f"Bad orbit direction {orbit_direction}. Must be ",
                 +"ASCENDING or DESCENDING.")

    da = xr.DataArray(
        data=data,
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
    # is for, but it is zero for L1B data.
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

    return wavelength, fwhm, radiance, radiance_noise


def shift_swir_channel(radiance, radiance_noise):
    # swir data is shifted by roughly 20 spatial pixels, see
    # https://www.eoportal.org/satellite-missions/enmap#hsi-hyperspectral-imager
    # Figure 17
    shift = 20

    shifted_radiance = radiance.copy()
    shifted_radiance[:-shift, ...] = shifted_radiance[shift:, ...]
    shifted_radiance[-shift:, ...] = \
        np.zeros(shape=(
            shift,
            shifted_radiance.shape[1],
            shifted_radiance.shape[2]))

    shifted_radiance_noise = radiance_noise.copy()
    shifted_radiance_noise[:-shift, ...] = \
        shifted_radiance_noise[shift:, ...]
    shifted_radiance_noise[-shift:, ...] = \
        np.zeros(shape=(
            shift,
            shifted_radiance_noise.shape[1],
            shifted_radiance_noise.shape[2]))

    return shifted_radiance, shifted_radiance_noise
