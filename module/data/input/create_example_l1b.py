import numpy as np

latitude = np.array([
    [49.4018, 49.4008, 49.3998, 49.3988],
    [49.4018, 49.4008, 49.3998, 49.3988],
    [49.4018, 49.4008, 49.3998, 49.3988],
])

longitude = np.array([
    [8.672434, 8.672434, 8.672434, 8.672434],
    [8.673434, 8.673434, 8.673434, 8.673434],
    [8.674434, 8.674434, 8.674434, 8.674434],
])

datetime = "2019-08-06T14:52:33Z"

#datetime = [
#    np.datetime64("2019-08-06T14:52:33"),
#]
# 
# data_vars = {
#     "latitude": (("line", "sample"), latitude),
#     "longitude": (("line", "sample"), longitude),
#     #"datetime": (("datetime"), datetime)
# }
# 
# l1b_data = xr.Dataset(data_vars = data_vars)
# 
# l1b_data.latitude.attrs["standard_name"] = "latitude"
# l1b_data.latitude.attrs["units"] = "degrees north"
# 
# l1b_data.longitude.attrs["standard_name"] = "longitude"
# l1b_data.longitude.attrs["units"] = "degrees east"
# 
# l1b_data.attrs = {"datetime": "2019-08-06T14:52:33"}
# 
# print(l1b_data)
# 
# l1b_data.to_netcdf("L1B_example.nc", mode = "w", format = "NETCDF4")

import xarray as xr
import numpy as np

l1b_data = xr.Dataset()

#l1b_data["datetime"] = datetime
l1b_data["latitude"] = (("line", "sample"), latitude)
l1b_data["longitude"] = (("line", "sample"), longitude)

l1b_data.attrs["ISO 8601 datetime"] = datetime

print(l1b_data)
l1b_data.to_netcdf("L1B_example.nc", mode = "w", format = "NETCDF4")
