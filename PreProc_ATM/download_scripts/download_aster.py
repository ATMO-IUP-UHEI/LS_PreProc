import os
import sys
import xarray as xr
import numpy as np
from getpass import getpass
from requests import Session
from http import HTTPStatus
import rioxarray
from rasterio.transform import from_bounds
from pyproj import CRS


def main():
    file_name_list = generate_file_list_from_l1b_file()
    download_files(file_name_list)


def generate_file_list_from_l1b_file():
    l1b_file = "DATA_IN/L1B_DATA.nc"
    l1b_data = xr.open_dataset(l1b_file)

    min_latitude = int(np.floor(l1b_data.latitude.min().values))
    max_latitude = int(np.ceil(l1b_data.latitude.max().values))
    min_longitude = int(np.floor(l1b_data.longitude.min().values))
    max_longitude = int(np.ceil(l1b_data.longitude.max().values))

    if max_longitude - min_longitude >= 180:
        min_longitude += 360
        min_longitude, max_longitude = max_longitude, min_longitude

    file_name_list = []
    for latitude in range(min_latitude, max_latitude):
        for longitude in range(min_longitude, max_longitude):
            if latitude >= 0:
                lat = f"N{latitude:02}"
            else:
                lat = f"S{abs(latitude):02}"

            if longitude >= 180:
                # wrapping around the back
                longitude = longitude - 360
            if longitude >= 0:
                lon = f"E{longitude:03}"
            else:
                lon = f"W{abs(longitude):03}"

            file_name_list.append(f"ASTGTMV003_{lat}{lon}_dem.tif")

    return file_name_list


def download_files(file_name_list):
    with Session() as session:
        #token = getpass("Enter token. If you don't have a token, get one from "
        #                + "https://urs.earthdata.nasa.gov/. It will have an "
        #                + "expiration date.\n\ttoken (your input will be "
        #                + "hidden): ")
        print("Using token in ~/.earthdata_token")
        token = open(os.path.expanduser('~')+'/.earthdata_token').read().strip()
        session.headers = {"Authorization": f"Bearer {token}"}

        for file_name in file_name_list:
            output_file = f"tmp/meteo/aster/{file_name}"
            if os.path.exists(output_file):

               print("download_aster.py: file exists!", output_file)

            else:
                file_url = "https://data.lpdaac.earthdatacloud.nasa.gov/"\
                    + f"lp-prod-protected/ASTGTM.003/{file_name}"

                # verify login works
                response = session.get(file_url)
                match response.status_code:
                    case HTTPStatus.UNAUTHORIZED:
                        print("not authorized. maybe you entered the wrong token "
                              + "or your token expired?")
                        sys.exit()
                    case HTTPStatus.NOT_FOUND:
                        print("url not found.")

                        above_water = confirm_tile_above_water(file_name)
                        if not above_water:
                            sys.exit("Tile is not above water but data doesn't "
                                     + "exist. Some error must have occured. "
                                     + "Exiting...")

                        # If tile is fully above water, there will be no data
                        # for it in the ASTER data base. Therefore, we create
                        # our own tile object with zeros in every entry.
                        create_tile_with_zeros(
                            file_name, f"tmp/meteo/aster/{file_name}")
                        continue
                if not HTTPStatus.OK:
                    sys.exit("HTTPStatus is not OK.")

                # download file
                print(f"downloading {file_name}")
                with open(f"tmp/meteo/aster/{file_name}", "wb") as file:
                    file.write(response.content)


def confirm_tile_above_water(file_name):
    print(f"{file_name} does not exist.")
    print("This may be due to the tile being fully above water.")
    print("Please go to the following link and confirm manually.")

    lat = file_name.split("_")[1][:3]
    if lat.startswith("N"):
        lat = int(lat[1:])
    elif lat.startswith("S"):
        lat = -int(lat[1:])
    lon = file_name.split("_")[1][3:]
    if lon.startswith("E"):
        lon = int(lon[1:])
    elif lon.startswith("W"):
        lon = -int(lon[1:])

    print("https://search.earthdata.nasa.gov")
    print("Enter into the rectangle bounding box:")
    print(f"SW: {lat},{lon} NE: {lat+1},{lon+1}")

    response = input("Is this tile fully above water? (y/n) ")
    while True:
        if response.lower().startswith("y"):
            return True
        if response.lower().startswith("n"):
            return False
        response = input("Please enter y or n...")


def create_tile_with_zeros(file_name, output_file):
    print("Creating .tif file with an altitude of 0 m in all pixels...")

    lat = file_name.split("_")[1][:3]
    if lat.startswith("N"):
        lat = int(lat[1:])
    elif lat.startswith("S"):
        lat = -int(lat[1:])
    lon = file_name.split("_")[1][3:]
    if lon.startswith("E"):
        lon = int(lon[1:])
    elif lon.startswith("W"):
        lon = -int(lon[1:])

    band = [1]
    x = np.linspace(lon, lon+1, 3601)
    y = np.linspace(lat+1, lat, 3601)

    data = np.zeros(shape=(len(band), len(y), len(x))).astype("float32")

    band_data = xr.DataArray(
        data,
        dims=["band", "y", "x"],
        coords={"band": band, "x": x, "y": y},
        name="band_data"
    )

    # Assign CRS and transform
    band_data = band_data.rio.write_crs("EPSG:4326")  # WGS84 Lat/Lon

    # Compute transform from bounds (left, bottom, right, top)
    transform = from_bounds(
        west=x.min(), south=y.min(), east=x.max(), north=y.max(),
        width=len(x), height=len(y)
    )
    band_data.rio.write_transform(transform, inplace=True)

    band_data.rio.to_raster(output_file)


if __name__ == "__main__":
    main()
