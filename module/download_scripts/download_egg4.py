import sys
import os
import xarray as xr
import numpy as np
import cdsapi


l1b_file = sys.argv[1]
l1b_data = xr.open_dataset(l1b_file)

database_name = "cams-global-ghg-reanalysis-egg4"

retrieve_parameters = {}

retrieve_parameters["variable"] = ["carbon_dioxide", "methane"]

# deal with datetime now
start = str(min(l1b_data.time).values)
stop = str(max(l1b_data.time).values)

start_yyyy_mm_dd = start[:10]
stop_yyyy_mm_dd = stop[:10]

if not start_yyyy_mm_dd == stop_yyyy_mm_dd:
    sys.exit("error: measurement wraps around midnight.")

yyyy_mm_dd = start_yyyy_mm_dd

yyyymmdd = "".join(yyyy_mm_dd.split("-"))

retrieve_parameters["date"] = f"{yyyy_mm_dd}/{yyyy_mm_dd}"

start_hh = int(start[11:13])
stop_hh = int(stop[11:13])

# get the hours that surround our measurement
available_hours = np.array([0, 3, 6, 9, 12, 15, 18, 21])
left_index = np.max(np.where(available_hours <= start_hh))
right_index = np.min(np.where(available_hours > stop_hh))

step = []
for hour in available_hours[left_index:right_index+1]:
    step.append(str(hour))

retrieve_parameters["step"] = step

min_latitude = np.round(
    l1b_data.latitude.min().values - 0.75, decimals=2
)
max_latitude = np.round(
    l1b_data.latitude.max().values + 0.75, decimals=2
)
min_longitude = np.round(
    l1b_data.longitude.min().values - 0.75, decimals=2
)
max_longitude = np.round(
    l1b_data.longitude.max().values + 0.75, decimals=2
)

if max_latitude > 90 or min_latitude < -90:
    sys.exit("latitude out of bounds")
if max_longitude > 180 or min_longitude < -180:
    sys.exit("longitude out of bounds")
if max_longitude - min_longitude >= 180:
    sys.exit("measurement wraps around back of Earth?")

retrieve_parameters["area"] = [
    max_latitude, min_longitude, min_latitude, max_longitude
]

retrieve_parameters["format"] = "grib"

retrieve_parameters["pressure_level"] = [
    "1", "2", "3", "5", "7", "10", "20", "30", "50", "70", "100", "150",
    "200", "250", "300", "400", "500", "600", "700", "800", "850", "900",
    "925", "950", "1000"
]

output_file = f"../data/egg4/download/egg4_{yyyymmdd}.grb"

os.system("cp ~/.cdsapirc_ads ~/.cdsapirc")

print(retrieve_parameters)

c = cdsapi.Client()
c.retrieve(
    database_name,
    retrieve_parameters,
    output_file
)
