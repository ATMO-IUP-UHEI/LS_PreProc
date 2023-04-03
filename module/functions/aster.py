import xarray as xr
import rioxarray
import numpy as np

def main(atm_data, aster_folder_path):
    print("getting aster data")
    # get aster data from the folder.
    # name(s) of the required file will be constructed from the coordinates in atm_data.
    aster_data = get_aster_data(atm_data, aster_folder_path)
    # drop unnecessary variables from aster_data and rename needed ones
    aster_data = prepare_aster_data(aster_data)
    # prepare variables for atm_data
    atm_data = prepare_atm_data(atm_data)

    # interpolation of aster_data and populating atm_data
    interpolated = aster_data.interp(latitude = atm_data.latitude, longitude = atm_data.longitude)
    atm_data.surface_elevation.values = interpolated.surface_elevation.values

    aster_data.close()

    return atm_data



def get_aster_data(atm_data, aster_folder_path):
    # get aster_data
    file_name_list = generate_file_list_from_atm_data(atm_data)

    aster_data_set = False

    for file_name in file_name_list:
        tmp_data = xr.open_dataset(f"{aster_folder_path}/{file_name}", engine = "rasterio")

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
    aster_data = aster_data.squeeze(dim = "band")
    aster_data = aster_data.drop(["spatial_ref", "band"])
    aster_data = aster_data.rename({"x": "longitude", "y": "latitude", "band_data": "surface_elevation"})

    aster_data.latitude.attrs["standard_name"] = "latitude"
    aster_data.latitude.attrs["units"] = "degrees north"

    aster_data.longitude.attrs["standard_name"] = "longitude"
    aster_data.longitude.attrs["units"] = "degrees east"

    aster_data.surface_elevation.attrs["standard_name"] = "surface elevation"
    aster_data.surface_elevation.attrs["units"] = "m"

    # interpolation will cause problems if some variables are float32 and some are float64. convert to float64
    aster_data = aster_data.assign(surface_elevation = aster_data.surface_elevation.astype(np.float64))

    return aster_data



def prepare_atm_data(atm_data):
    nline = len(atm_data.line)
    nsample = len(atm_data.sample)

    atm_data["surface_elevation"] = (("line", "sample"), np.empty(shape = (nline, nsample)))
    atm_data.surface_elevation.attrs["standard_name"] = "surface elevation"
    atm_data.surface_elevation.attrs["units"] = "m"

    return atm_data



if __name__ == "__main__":
    main(aster_data, aster_folder_path)
