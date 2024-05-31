import xarray as xr
import matplotlib.pyplot as plt

data_path = "L1B_riyadh.nc"

data1 = xr.open_dataset(data_path, group="BAND01")
data2 = xr.open_dataset(data_path, group="BAND02")

wavelength = []
for w in data1.wavelength.values:
    wavelength.append(w)
for w in data2.wavelength.values:
    wavelength.append(w)

fwhm = []
for f in data1.fwhm.values:
    fwhm.append(f)
for f in data2.fwhm.values:
    fwhm.append(f)

plt.title("FWHM Riyadh PP09 Case1")
plt.scatter(wavelength, fwhm, marker=".", color="black")
plt.xlabel("wavelength / nm")
plt.ylabel("fwhm / nm")
plt.gca().axvspan(
    1982, 2092, color="blue", alpha=.25, label="Scheidweiler, CO2")
plt.gca().axvspan(
    2110, 2420, color="red", alpha=.25, label="Scheidweiler, CH4")
plt.axvline(
    2298.42, color="black", alpha=.25, label="Roger et al. (2023), CH4")
plt.legend()
plt.savefig("fwhm_graph.png", dpi=200)
plt.show()
