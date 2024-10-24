import sys
import xarray as xr
import numpy as np
from getpass import getpass
from requests import Session
from http import HTTPStatus


def generate_file_list_from_l1b_file():
    l1b_file = "SYNTH_SPECTRA/L1B_DATA.nc"
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
        token = getpass("Enter token. If you don't have a token, get one from "
                        + "https://urs.earthdata.nasa.gov/. It will have an "
                        + "expiration date.\n\ttoken (your input will be "
                        + "hidden): ")
        session.headers = {"Authorization": f"Bearer {token}"}

        for file_name in file_name_list:
            file_url = "https://data.lpdaac.earthdatacloud.nasa.gov/"\
                     + f"lp-prod-protected/ASTGTM.003/{file_name}"

            # verify login works
            response = session.get(file_url)
            match response.status_code:
                case HTTPStatus.UNAUTHORIZED:
                    print("not authorized. maybe you entered the wrong token "
                          + "or your token expired?")
                case HTTPStatus.NOT_FOUND:
                    print("url not found.")
            if not HTTPStatus.OK:
                sys.exit("HTTPStatus is not OK.")

            # download file
            print(f"downloading {file_name}")
            with open(f"tmp/meteo/aster/{file_name}", "wb") as file:
                file.write(response.content)


file_name_list = generate_file_list_from_l1b_file()
download_files(file_name_list)
