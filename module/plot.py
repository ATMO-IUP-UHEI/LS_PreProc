import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

data = xr.open_dataset("data/output/ATM_example.nc")
print(data)

plt.plot(data.h2o.values[:, 0, 0], data.pressure.values[:, 0, 0]/100, label = "h2o")
plt.plot(data.co2.values[:, 0, 0], data.pressure.values[:, 0, 0]/100, label = "co2")
plt.plot(data.ch4.values[:, 0, 0], data.pressure.values[:, 0, 0]/100, label = "ch4")
plt.xlabel(f'mole fraction / {data.co2.attrs["units"]}')
plt.ylabel(f'{data.pressure.attrs["standard_name"]} / hPa')
plt.gca().invert_yaxis()
plt.legend()
plt.show()

plt.plot(data.temperature.values[:, 0, 0], data.pressure.values[:, 0, 0]/100)
plt.xlabel(f"temperature / K")
plt.ylabel(f"pressure / hPa")
plt.yscale("log")
plt.gca().invert_yaxis()
plt.show()
