import sys
import xarray as xr
import netCDF4 as nc
import matplotlib.pyplot as plt

tile_list = [*sys.argv[1:]]

# check number of groups
with nc.Dataset(tile_list[0]) as data:
    Ngroups = len(data.groups)

out_root = None
out_bands = [None]*Ngroups
out_history = ""

orbit_type = None

for tile in tile_list:
    print(f"opening tile {tile}")
    root = xr.open_dataset(tile)

    if not out_root:
        out_root = root
        out_history = root.history[:-1]
    else:
        out_root = xr.concat((out_root, root), dim="frame")
        out_history += root.history.split("file(s):")[1].strip(".")

    root.close()

    for group in range(Ngroups):
        band = xr.open_dataset(tile, group=f"BAND{group+1:02}")

        if not out_bands[group]:
            out_bands[group] = band
        else:
            out_bands[group] = xr.concat((out_bands[group], band), dim="frame")

        band.close()

print("sorting by latitude")
sort_array = out_root.latitude.mean(dim="line")
# individual tiles should determine the order direction, as the
# orbit_direction variable in the raw EnMAP data has been accounted
# for when creating the nc files for the individual tiles.
# We here copy the orbit direction from random pixels in the first tile.
sort_direction = out_root.latitude[100, 0] < out_root.latitude[200, 0]

for group in range(Ngroups):
    out_bands[group] = out_bands[group].sortby(
        sort_array, ascending=sort_direction)
out_root = out_root.sortby(
    sort_array, ascending=sort_direction)

# print("plotting")
# plt.pcolormesh(out_root.longitude, out_root.latitude, out_bands[0].radiance[:, :, 0])
# plt.show()

print("saving")
for var in out_root.data_vars:
    out_root[var].encoding.update({"_FillValue": None})
out_root.attrs["history"] = out_history
out_root.to_netcdf("L1B_DATA.nc")

for group in range(Ngroups):
    for var in out_bands[group]:
        out_bands[group][var].encoding.update({"_FillValue": None})
    out_bands[group].to_netcdf("L1B_DATA.nc", group=f"BAND{group+1:02}",
                               mode="a")
