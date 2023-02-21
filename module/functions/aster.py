import xarray as xr
import rioxarray
import numpy as np

def import_aster(atm_data, aster_folder_path):
    aster_latlon_list = ["N49E008"]
    aster_data_set = False
    for aster_latlon in aster_latlon_list:
        tmp_data = xr.open_dataset(f"{aster_folder_path}/ASTGTMV003_{aster_latlon}_dem.tif", engine = "rasterio")

        if not aster_data_set:
            aster_data = tmp_data
            aster_data_set = True
            continue

        aster_data = xr.merge([aster_data, tmp_data])

        tmp_data.close()

    aster_data = aster_data.squeeze(dim = "band")
    aster_data = aster_data.drop(["spatial_ref", "band"])
    aster_data = aster_data.rename({"x": "longitude", "y": "latitude", "band_data": "surface_elevation"})

    aster_data.latitude.attrs["standard_name"] = "latitude"
    aster_data.latitude.attrs["units"] = "degrees north"

    aster_data.longitude.attrs["standard_name"] = "longitude"
    aster_data.longitude.attrs["units"] = "degrees east"

    aster_data.surface_elevation.attrs["standard_name"] = "surface elevation"
    aster_data.surface_elevation.attrs["units"] = "m"

    # prepare atm_data
    atm_data["surface_elevation"] = (("line", "sample"), np.zeros(shape = (len(atm_data.line), len(atm_data.sample))))
    atm_data.surface_elevation.attrs["standard_name"] = "surface elevation"
    atm_data.surface_elevation.attrs["units"] = "m"

    for line in range(len(atm_data.line)):
        for sample in range(len(atm_data.sample)):
            latitude = atm_data.latitude.values[line, sample]
            longitude = atm_data.longitude.values[line, sample]

            # interpolate
            interpolated = aster_data.interp(latitude = latitude, longitude = longitude)

            # write to atm_data
            atm_data.surface_elevation.values[line, sample] = interpolated.surface_elevation.values

    aster_data.close()
################################################################################
