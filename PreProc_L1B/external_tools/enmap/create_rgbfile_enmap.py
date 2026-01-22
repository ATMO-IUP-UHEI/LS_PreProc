import xarray as xr
import numpy as np
# from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import sys

import plot_functions


def main():
    path = "SYNTH_SPECTRA/L1B_DATA.nc"
    root = xr.open_dataset(path)
    # band_01 = xr.open_dataset(path, group="BAND01")
    band_02 = xr.open_dataset(path, group="BAND02")

    x = root.longitude
    y = root.latitude

    # wavelength_01 = band_01.wavelength
    # example_spectrum_01 = band_01.radiance[0, 0, :]

    wavelength_02 = band_02.wavelength
    # example_spectrum_02 = band_02.radiance[0, 0, :]

    # plt.plot(wavelength_01, example_spectrum_01)
    # plt.plot(wavelength_02, example_spectrum_02)
    # plt.xlabel("wavelength / nm")
    # plt.ylabel("radiance / photons s-1 cm-2 sr-1 nm-1")
    # plt.show()

    rgb = apply_lms_curves(wavelength_02, band_02.radiance)
    rgb = weigh(rgb)
    rgb = clip_and_normalize(rgb)
    rgb = gamma_correct(rgb)

    plot(x, y, rgb)

    out_data = xr.Dataset()

    out_data["longitude"] = x
    out_data["latitude"] = y
    out_data["rgb_values"] = xr.DataArray(
        data=rgb,
        dims=("frame", "line", "rgb")
    ).astype("float32")

    for var in out_data.data_vars:
        out_data[var].encoding.update({"_FillValue": None})

    out_data.to_netcdf("CONTRL_OUT/RGB_DATA.nc", mode="w", format="NETCDF4")


def apply_lms_curves(wavelength, radiance):
    cone_sensitivity_l = np.exp(-0.5 * ((wavelength - 600)/50)**2)
    cone_sensitivity_m = np.exp(-0.5 * ((wavelength - 550)/40)**2)
    cone_sensitivity_s = np.exp(-0.5 * ((wavelength - 450)/30)**2)

    norm_l = np.sum(cone_sensitivity_l)
    norm_m = np.sum(cone_sensitivity_m)
    norm_s = np.sum(cone_sensitivity_s)

    cone_sensitivity_l /= norm_l
    cone_sensitivity_m /= norm_m
    cone_sensitivity_s /= norm_s

    # fig, ax1 = plt.subplots()
    # ax1.plot(wavelength, radiance[0, 0, :])
    # ax2 = ax1.twinx()
    # ax2.plot(wavelength, cone_sensitivity_l, color="red", label="red")
    # ax2.plot(wavelength, cone_sensitivity_m, color="green", label="green")
    # ax2.plot(wavelength, cone_sensitivity_s, color="blue", label="blue")
    # ax1.set_xlabel("wavelength / nm")
    # ax1.set_ylabel("radiance / photons s-1 cm-2 sr-1 nm-1")
    # ax2.set_ylabel("cone sensitivity")
    # plt.legend()
    # plt.show()
    # sys.exit()

    r = np.sum(np.multiply(radiance, cone_sensitivity_l), axis=-1)
    g = np.sum(np.multiply(radiance, cone_sensitivity_m), axis=-1)
    b = np.sum(np.multiply(radiance, cone_sensitivity_s), axis=-1)

    rgb = np.dstack((r, g, b))
    return rgb


def weigh(rgb):
    return np.multiply(rgb, (1.0, 1.0, 1.0))


def clip_and_normalize(rgb):
    flat_image = rgb.flatten()
    lower_clip = np.percentile(flat_image, 0.03)
    upper_clip = np.percentile(flat_image, 99.97)

    rgb = np.clip(rgb, lower_clip, upper_clip)
    rgb = (rgb - lower_clip)/(upper_clip - lower_clip)

    return rgb


def gamma_correct(rgb):
    gamma = 2.2
    return np.power(rgb, 1/gamma)


def plot(x, y, rgb):
    fig = plt.figure()
    fig.set_size_inches(
        rgb.shape[1]/1024*5,
        rgb.shape[0]/1024*5
    )
    plt.xlabel("Longitude / deg")
    plt.ylabel("Latitude / deg")
    plt.pcolor(x, y, rgb, shading="nearest")
    data_lims = [
        np.min(x), np.max(x),
        np.min(y), np.max(y)
    ]
    aspect = plot_functions.get_aspect(*data_lims)
    plt.gca().set_aspect(aspect)
    plt.show()


if __name__ == "__main__":
    main()
