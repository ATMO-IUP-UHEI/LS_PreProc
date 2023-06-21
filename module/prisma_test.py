import xarray as xr
import sys

def main():
    input_file = sys.argv[1]
    input_name = input_file.split("/")[-1]

    l1b_data = xr.open_dataset(input_file)
    band1_data = xr.open_dataset(input_file, group="BAND1")

    l1b_data = set_datetime(l1b_data)
    l1b_data = rename_dims(l1b_data)
    band1_data = rename_dims(band1_data)

    l1b_data.to_netcdf(f"SYNTH_SPECTRA/{input_name}", mode="w", format="NETCDF4")
    band1_data.to_netcdf(f"SYNTH_SPECTRA/{input_name}", group="BAND1", mode="a", format="NETCDF4")



def rename_dims(group):
    group = group.rename_dims({"nlat": "line", "nlon": "sample"})

    try:
        group = group.rename_dims({"nwave": "wave"})
    except:
        pass

    return group



def set_datetime(l1b_data):
    year = f"{l1b_data.time[0]:04.0f}"
    month = f"{l1b_data.time[1]:02.0f}"
    day = f"{l1b_data.time[2]:02.0f}"
    hour = f"{l1b_data.time[3]:02.0f}"
    minute = f"{l1b_data.time[4]:02.0f}"
    second = f"{l1b_data.time[5]:02.0f}"
    datetime = f"{year}-{month}-{day}T{hour}:{minute}:{second}Z"

    l1b_data.attrs["ISO 8601 datetime"] = datetime
    l1b_data = l1b_data.drop_vars({"time"})

    return l1b_data



if __name__ == "__main__":
    main()
