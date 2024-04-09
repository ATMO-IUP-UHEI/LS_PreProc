import xarray as xr
import numpy as np
from datetime import datetime


def main(atm, aster_folder_path, dims):
    print(datetime.now(), "getting aster data")
    # get aster data from the folder.
    # name(s) of the required file will be constructed from the coordinates in
    # atm_data.
    aster = get_aster_data(atm, aster_folder_path)
    # drop unnecessary variables from aster_data and rename needed ones
    aster = prepare_aster_data(aster)

    # interpolation of aster_data onto atm grid
    aster = aster.interp(
        latitude=atm.latitude,
        longitude=atm.longitude
    )

    atm["surface_elevation"] = xr.DataArray(
        data=get_surface_elevation(aster),
        dims=(dims["y"], dims["x"]),
    ).astype("float32")
    atm.surface_elevation.attrs["source"] = "aster"

    return atm


def get_aster_data(atm_data, aster_folder_path):
    # get aster_data
    file_name_list = generate_file_list_from_atm_data(atm_data)

    aster_data_set = False

    for file_name in file_name_list:
        tmp_data = xr.open_dataset(
            f"{aster_folder_path}/{file_name}", engine="rasterio"
        )

        if not aster_data_set:
            aster_data = tmp_data
            aster_data_set = True
            continue

        aster_data = xr.merge([aster_data, tmp_data])

        tmp_data.close()

    return aster_data


def generate_file_list_from_atm_data(atm_data):
    min_latitude = int(np.floor(atm_data.latitude.min().values))
    max_latitude = int(np.ceil(atm_data.latitude.max().values))
    min_longitude = int(np.floor(atm_data.longitude.min().values))
    max_longitude = int(np.ceil(atm_data.longitude.max().values))

    file_name_list = []
    for latitude in range(min_latitude, max_latitude):
        for longitude in range(min_longitude, max_longitude):
            if latitude >= 0:
                lat = f"N{latitude:02}"
            else:
                lat = f"S{abs(latitude):02}"

            if longitude >= 0:
                lon = f"E{longitude:03}"
            else:
                lon = f"W{abs(longitude):03}"

            file_name_list.append(f"ASTGTMV003_{lat}{lon}_dem.tif")

    return file_name_list


def prepare_aster_data(aster_data):
    # prepare aster_data
    aster_data = aster_data.squeeze(dim="band")
    aster_data = aster_data.drop(["spatial_ref", "band"])
    aster_data = aster_data.rename(
        {"x": "longitude", "y": "latitude", "band_data": "surface_elevation"}
    )

    aster_data.latitude.attrs["standard_name"] = "latitude"
    aster_data.latitude.attrs["units"] = "degrees north"

    aster_data.longitude.attrs["standard_name"] = "longitude"
    aster_data.longitude.attrs["units"] = "degrees east"

    aster_data.surface_elevation.attrs["standard_name"] = "surface elevation"
    aster_data.surface_elevation.attrs["units"] = "m"

    # interpolation will cause problems if some variables are float32 and some
    # are float64. convert to float64
    aster_data = aster_data.assign(
        surface_elevation=aster_data.surface_elevation.astype(np.float64)
    )

    return aster_data


def get_surface_elevation(aster_data):
    return aster_data.surface_elevation.values
