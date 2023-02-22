from getpass import getpass
from requests import Session
from http import HTTPStatus
import sys

with Session() as session:
    token = getpass("Enter token. If you don't have a token, get one from https://urs.earthdata.nasa.gov/. It will have an expiration date.\n\ttoken: ")
    session.headers = {"Authorization": f"Bearer {token}"}

    file_url = "https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/ASTGTM.003/ASTGTMV003_N49E008_dem.tif"
    file_name = file_url.split("/")[-1]

    # verify login works
    response = session.get(file_url)
    match response.status_code:
        case HTTPStatus.UNAUTHORIZED:
            print("not authorized. maybe you entered the wrong token or your token expired?")
        case HTTPStatus.NOT_FOUND:
            print("url not found.")
    if not HTTPStatus.OK:
        sys.exit("HTTPStatus is not OK.")

    # download file
    print(f"downloading {file_name}")
    with open(file_name, "wb") as file:
        file.write(response.content)
