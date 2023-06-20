import sys
import os
import numpy as np
import xarray as xr
import cdsapi

def main():
    l1b_file = sys.argv[1]
    l1b_data = xr.open_dataset(l1b_file)

    database_name = "reanalysis-era5-complete"

    retrieve_parameters = {}

    retrieve_parameters["class"] = "ea" # do not change
    retrieve_parameters["expver"] = "1" # do not change
    retrieve_parameters["stream"] = "oper" # do not change
    retrieve_parameters["type"] = "an" # type: use "an" (analysis) unless you have a particular reason to use "fc" (forecast) https://confluence.ecmwf.int/pages/viewpage.action?pageId=85402030

    retrieve_parameters["grid"] = "0.25/0.25" # 0.25degx0.25deg grid

    yyyy_mm_dd = l1b_data.attrs["ISO 8601 datetime"].split("T")[0]
    yyyymmdd = "".join(yyyy_mm_dd.split("-"))

    retrieve_parameters["date"] = yyyy_mm_dd

    hh = int(l1b_data.attrs["ISO 8601 datetime"].split("T")[1].strip("Z").split(":")[0])
    if hh >= 23:
        sys.exit("get 23 of this day and 0 of next day. not implemented, yet.")
    time = [hh, hh+1]
    retrieve_parameters["time"] = time

    min_latitude = l1b_data.latitude.min().values - 0.25
    max_latitude = l1b_data.latitude.max().values + 0.25
    min_longitude = l1b_data.longitude.min().values - 0.25
    max_longitude = l1b_data.longitude.max().values + 0.25

    if max_latitude > 90 or min_latitude < -90:
        sys.exit("latitude out of bounds")
    if max_longitude > 180 or min_longitude < -180:
        sys.exit("longitude out of bounds")
    if max_longitude - min_longitude >= 180:
        sys.exit("measurement wraps around back of Earth?")

    retrieve_parameters["area"] = [max_latitude, min_longitude, min_latitude, max_longitude]

    os.system("cp ~/.cdsapirc_cds ~/.cdsapirc")

    c = cdsapi.Client()
    # multilayer
    retrieve_parameters["levtype"] = "ml"
    retrieve_parameters["levelist"] = "1/to/137"
    retrieve_parameters["param"] = "130/133"
    # 130 - temperature / K
    # 133 - specific humidity / kg kg-1
    output_file = f"output/era5_ml_{yyyymmdd}.grb"
    c.retrieve(
        database_name,
        retrieve_parameters,
        output_file
    )

    # surface
    retrieve_parameters["levtype"] = "sfc"
    retrieve_parameters["levelist"] = "1"
    retrieve_parameters["param"] = "129/134/165/166"
    # 129 - geopotential / m2 s-2
    # 134 - surface pressure / kg m-1 s-2
    # 165 - 10m u-component of wind / m s-1
    # 166 - 10m v-component of wind / m s-1
    output_file = f"output/era5_sfc_{yyyymmdd}.grb"
    c.retrieve(
        database_name,
        retrieve_parameters,
        output_file
    )

if __name__ == "__main__":
    main()
