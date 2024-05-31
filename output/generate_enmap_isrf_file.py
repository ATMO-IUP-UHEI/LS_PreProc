import xarray as xr
import numpy as np
import matplotlib.pyplot as plt


def gauss(x, fwhm):
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    return 1/(sigma*np.sqrt(2*np.pi)) * np.exp(-0.5*x**2/sigma**2)


for group_number in range(2):  # I know that enmap has two bands
    group_number += 1

    data = xr.open_dataset("L1B_ekizak.nc", group=f"BAND{group_number:02}")

    wavelength_center = data.wavelength.values
    fwhm = data.fwhm.values

    # number of offsets must be odd for RemoTeC
    wavelength_offset = np.linspace(-3*np.max(fwhm), 3*np.max(fwhm), 999)

    response = np.empty(shape=(len(wavelength_center), len(wavelength_offset)))

    for i in range(len(wavelength_center)):
        response[i, :] = gauss(wavelength_offset, fwhm[i])

    out_data = xr.Dataset()

    out_data["wavelength_center"] = xr.DataArray(
        data=wavelength_center,
        dims="channel",
    ).astype("float32")
    out_data.wavelength_center.attrs["units"] = "nm"

    out_data["wavelength_offset"] = xr.DataArray(
        data=wavelength_offset,
        dims="d_channel",
    ).astype("float32")
    out_data.wavelength_offset.attrs["units"] = "nm"

    out_data["response"] = xr.DataArray(
        data=response,
        dims=("channel", "d_channel"),
    ).astype("float32")

    for var in out_data.data_vars:
        out_data[var].encoding.update({"_FillValue": None})

    if group_number == 1:
        mode = "w"
    else:
        mode = "a"

    out_data.to_netcdf(
        "isrf_enmap.nc",
        group=f"BAND{group_number:02}",
        mode=mode,
        format="NETCDF4"
    )

    out_data.close()
