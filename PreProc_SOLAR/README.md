# PreProc_SOLAR

Solar irradiance spectrum is a standardized file with name solar_irradiance.nc
It has the dimensions:
- channel
It has the variables:
- wavelength / nm
- solar_irradiance / W m-2 nm-1
- solar_irradiance_err / W m-2 nm-1

Historically, the solar irradiance was either provided as a line list and then computed internally in the read routine, or provided in different formats with a flag defining which routine is used.
Now, there is just one standardized solar_irradiance.nc file, any new formats need to respect this file and the file needs to be created by the preprocessor.
For reference, the old read routines are provided in a folder together with this preprocessor.
