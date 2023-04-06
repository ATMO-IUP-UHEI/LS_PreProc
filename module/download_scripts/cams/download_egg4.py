import sys
import os
import numpy as np
import xarray as xr
import cdsapi

def main():
    l1b_file = sys.argv[1]
    l1b_data = xr.open_dataset(l1b_file)

    database_name = "cams-global-ghg-reanalysis-egg4"

    retrieve_parameters = {}

    retrieve_parameters["variable"] = ["carbon_dioxide", "methane"]

    yyyy_mm_dd = l1b_data.attrs["ISO 8601 datetime"].split("T")[0]
    yyyymmdd = "".join(yyyy_mm_dd.split("-"))
    retrieve_parameters["date"] = f"{yyyy_mm_dd}/{yyyy_mm_dd}"

    hh = int(l1b_data.attrs["ISO 8601 datetime"].split("T")[1].strip("Z").split(":")[0])
    if hh >=  0 and hh <  3: step = ["0", "3"]
    if hh >=  3 and hh <  6: step = ["3", "6"]
    if hh >=  6 and hh <  9: step = ["6", "9"]
    if hh >=  9 and hh < 12: step = ["9", "12"]
    if hh >= 12 and hh < 15: step = ["12", "15"]
    if hh >= 15 and hh < 18: step = ["15", "18"]
    if hh >= 18 and hh < 21: step = ["18", "21"]
    if hh >= 21: sys.exit("get 21 of this day and 0 of next day. not implemented, yet.")
    retrieve_parameters["step"] = step

    min_latitude = l1b_data.latitude.min().values - 0.75
    max_latitude = l1b_data.latitude.max().values + 0.75
    min_longitude = l1b_data.longitude.min().values - 0.75
    max_longitude = l1b_data.longitude.max().values + 0.75

    if max_latitude > 90 or min_latitude < -90:
        sys.exit("latitude out of bounds")
    if max_longitude > 180 or min_longitude < -180:
        sys.exit("longitude out of bounds")
    if max_longitude - min_longitude >= 180:
        sys.exit("measurement wraps around back of Earth?")

    retrieve_parameters["area"] = [max_latitude, min_longitude, min_latitude, max_longitude]

    retrieve_parameters["format"] = "grib"

    retrieve_parameters["pressure_level"] = [
        "1", "2", "3", "5", "7", "10", "20", "30", "50", "70", "100", "150", "200", "250",
        "300", "400", "500", "600", "700", "800", "850", "900", "925", "950", "1000"
    ]

    output_file = f"output/egg4_{yyyymmdd}.grb"

    os.system("cp ~/.cdsapirc_ads ~/.cdsapirc")

    print(database_name, retrieve_parameters, output_file)
    sys.exit()

    c = cdsapi.Client()
    c.retrieve(
        database_name,
        retrieve_parameters,
        output_file
    )

if __name__ == "__main__":
    main()
