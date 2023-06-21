import sys
import os
import xarray as xr

path = "../../../data/scenarios/prisma/china_guanter/RAW/"
data = xr.open_dataset(os.path.join(path, "L1B_China_Guanter.nc"))
band1 = xr.open_dataset(os.path.join(path, "L1B_China_Guanter.nc"), group="BAND1")

data = data.rename({
    "nlat": "frame",
    "nlon": "line"
})
data = data.drop_vars(["time"])
data = data.transpose(
    "line", "frame"
)
data.attrs["ISO 8601 datetime"] = "2020-02-06T03:18:40Z"

band1 = band1.rename({
    "nlat": "frame",
    "nlon": "line",
    "nwave": "band",
})
band1 = band1.drop_vars([
    "radiance_error", "irradiance", "irradiance_error", "irradiance_noise"
])
band1 = band1.transpose(
    "line", "frame", "band"
)

print(data)
print(band1)
sys.exit()
data.to_netcdf(os.path.join(path, "L1B_China_Guanter_formatted.nc"))
