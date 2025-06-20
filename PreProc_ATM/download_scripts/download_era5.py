import sys
import os
import xarray as xr
import cdsapi


l1b_file = "SYNTH_SPECTRA/L1B_DATA.nc"
l1b_data = xr.open_dataset(l1b_file)

database_name = "reanalysis-era5-complete"

retrieve_parameters = {}

# do not change
retrieve_parameters["class"] = "ea"

# do not change
retrieve_parameters["expver"] = "1"

# do not change
retrieve_parameters["stream"] = "oper"

# type: use "an" (analysis) unless you have a particular reason
# to use "fc" (forecast)
# https://confluence.ecmwf.int/pages/viewpage.action?pageId=85402030
retrieve_parameters["type"] = "an"

# 0.25degx0.25deg grid
retrieve_parameters["grid"] = "0.25/0.25"

# deal with datetime now
start = str(min(l1b_data.time).values)
stop = str(max(l1b_data.time).values)

start_yyyy_mm_dd = start[:10]
stop_yyyy_mm_dd = stop[:10]

if not start_yyyy_mm_dd == stop_yyyy_mm_dd:
    sys.exit("error: measurement wraps around midnight.")
yyyy_mm_dd = start_yyyy_mm_dd

yyyymmdd = "".join(yyyy_mm_dd.split("-"))

retrieve_parameters["date"] = yyyy_mm_dd

start_hh = int(start[11:13])
stop_hh = int(stop[11:13])

if stop_hh >= 23:
    sys.exit("get 23 of this day and 0 of next day. not implemented, yet.")
time = [start_hh, stop_hh+1]
retrieve_parameters["time"] = time

# handle latitude and longitude now
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

retrieve_parameters["area"] = [
    str(max_latitude),
    str(min_longitude),
    str(min_latitude),
    str(max_longitude)
]

os.system("cp ~/.cdsapirc_cds ~/.cdsapirc")

c = cdsapi.Client()
# multilayer
retrieve_parameters["levtype"] = "ml"
retrieve_parameters["levelist"] = "1/to/137"
retrieve_parameters["param"] = "130/133/131/132"
# 130 - temperature / K
# 133 - specific humidity / kg kg-1
# 131 - u component of wind / m s-1
# 132 - v component of wind / m s-1

output_file_ml = f"tmp/meteo/era5/era5_ml_{yyyymmdd}.grb"

c.retrieve(
    database_name,
    retrieve_parameters,
    output_file_ml
)

output_file_sfc = f"tmp/meteo/era5/era5_sfc_{yyyymmdd}.grb"

# surface
retrieve_parameters["levtype"] = "sfc"
retrieve_parameters["levelist"] = "1"
retrieve_parameters["param"] = "129/134/165/166"
# 129 - geopotential / m2 s-2
# 134 - surface pressure / kg m-1 s-2
# 165 - 10m u-component of wind / m s-1
# 166 - 10m v-component of wind / m s-1

c.retrieve(
    database_name,
    retrieve_parameters,
    output_file_sfc
)
