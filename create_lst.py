import xarray as xr

data = xr.open_dataset("SYNTH_SPECTRA/L1B_DATA.nc")
print("lines:", data.sizes["line"])
print("frames:", data.sizes["frame"])

with open("LST/full.lst", "w") as file:
    for x in range(1, data.sizes["line"]+1):
        for y in range(1, data.sizes["frame"]+1):
            file.writelines(f"DATA.nc_X{x:>06}_Y{y:>06}\n")
