import xarray as xr

input_path = "/home/hd/hd_hd/hd_gb423/sds/data/input/solar"
input_file = "hybrid_reference_spectrum_c2021-03-04_with_unc.nc"
input_format = "TSIS-1 HSRS"

output_path = "/home/hd/hd_hd/hd_gb423/sds/data/input/solar"
output_file = "solar_irradiance_reference_spectrum.nc"


match input_format:
    case "TSIS-1 HSRS":
        print("right units already in file, only renaming variables")
        input_data = xr.open_dataset(f"{input_path}/{input_file}")

        output_data = xr.Dataset()
        output_data["wavelength"] = (("channel"), input_data["Vacuum Wavelength"].values)
        output_data["wavelength"].attrs["units"] = "nm"
        output_data["irradiance"] = (("channel"), input_data["SSI"].values)
        output_data["irradiance"].attrs["units"] = "W m-2 nm-1"
        output_data["irradiance_err"] = (("channel"), input_data["SSI_UNC"].values)
        output_data["irradiance_err"].attrs["units"] = "W m-2 nm-1"

        output_data.to_netcdf(f"{output_path}/{output_file}", mode="w")
