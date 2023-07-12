import sys
import configparser

from instruments import prisma
from instruments import enmap
from instruments import dlr_hyspex


def main():
    try:
        config_file = sys.argv[1]
        config = configparser.ConfigParser()
        config.read(config_file)
    except IndexError:
        sys.exit("Provide settings file as command line argument.")

    dims = get_dims()
    root, band_list, input_file_list = \
        get_data(config["config"], dims)
    set_attributes(root, band_list, dims)
    write_data(config["config"], root, band_list, input_file_list)

    return


def get_dims():
    dims = {
        # across-track spatial dimension
        # These spatial pixels are measured simultaneously for a push-broom
        # imaging approach.
        "x": "line",
        # along-track spatial dimension
        # These spatial pixels are measured in sequence by the same detector
        # row for a push-broom imaging approach.
        "y": "frame",
        # spectral dimension
        # This is the spectral dimension of the detectors. This can be e.g.
        # different wavelengths or wavenumbers, depending on the instrument.
        "z": "channel",
    }
    return dims


def get_data(config, dims):
    match config["instrument"]:
        case "PRISMA":
            root, band_list, input_file_list = \
                prisma.import_data(config, dims)
        case "EnMAP":
            root, band_list, input_file_list = \
                enmap.import_data(config, dims)
        case "DLR HySpex":
            root, band_list, input_file_list = \
                dlr_hyspex.import_data(config, dims)
        case _:
            sys.exit(f"invalid instrument {config['instrument']}")

    return root, band_list, input_file_list


def set_attributes(root, band_list, dims):
    assert "time" in root.data_vars
    assert root.time.dims == (dims["y"],)
    assert root.time.dtype == "float32"
    assert root.time.attrs["units"].startswith("seconds since ")
    root.time.attrs["long_name"] = "date and time in UTC"

    assert "latitude" in root.data_vars
    assert root.latitude.dims == (dims["y"], dims["x"])
    assert root.latitude.dtype == "float32"
    root.latitude.attrs["long_name"] = \
        "latitude at pixel center"
    root.latitude.attrs["units"] = \
        "degrees north"

    assert "longitude" in root.data_vars
    assert root.longitude.dims == (dims["y"], dims["x"])
    assert root.longitude.dtype == "float32"
    root.longitude.attrs["long_name"] = \
        "longitude at pixel center"
    root.longitude.attrs["units"] = \
        "degrees east"

    assert "solar_zenith_angle" in root.data_vars
    assert root.solar_zenith_angle.dims == (dims["y"], dims["x"])
    assert root.solar_zenith_angle.dtype == "float32"
    root.solar_zenith_angle.attrs["long_name"] = \
        "solar zenith angle"
    root.solar_zenith_angle.attrs["units"] = \
        "degrees"

    if "solar_azimuth_angle" in root.data_vars:
        assert root.solar_azimuth_angle.dims == (dims["y"], dims["x"])
        assert root.solar_azimuth_angle.dtype == "float32"
        root.solar_azimuth_angle.attrs["long_name"] = \
            "solar azimuth angle"
        root.solar_azimuth_angle.attrs["units"] = \
            "degrees"

    assert "viewing_zenith_angle" in root.data_vars
    assert root.viewing_zenith_angle.dims == (dims["y"], dims["x"])
    assert root.viewing_zenith_angle.dtype == "float32"
    root.viewing_zenith_angle.attrs["long_name"] = \
        "viewing zenith angle"
    root.viewing_zenith_angle.attrs["units"] = \
        "degrees"

    if "viewing_azimuth_angle" in root.data_vars:
        assert root.viewing_azimuth_angle.dims == (dims["y"], dims["x"])
        assert root.viewing_azimuth_angle.dtype == "float32"
        root.viewing_azimuth_angle.attrs["long_name"] = \
            "viewing azimuth angle"
        root.viewing_azimuth_angle.attrs["units"] = \
            "degrees"

    if "observer_altitude" in root.data_vars:
        assert root.observer_altitude.dims == (dims["y"], dims["x"])
        assert root.observer_altitude.dtype == "float32"
        root.observer_altitude.attrs["long_name"] = \
            "observer altitude above ground"
        root.observer_altitude.attrs["units"] = \
            "m"

    for band in band_list:
        assert len(band[dims["x"]]) == len(root[dims["x"]])
        assert len(band[dims["y"]]) == len(root[dims["y"]])

        assert "wavelength" in band.data_vars
        assert band.wavelength.dims == (dims["z"],)
        assert band.wavelength.dtype == "float32"
        band.wavelength.attrs["long_name"] = "wavelength"
        band.wavelength.attrs["units"] = "nm"

        assert "radiance" in band.data_vars
        assert band.radiance.dims == (dims["y"], dims["x"], dims["z"])
        assert band.radiance.dtype == "float32"
        band.radiance.attrs["long_name"] = \
            "at-sensor radiance"
        band.radiance.attrs["units"] = \
            "photons s-1 cm-2 sr-1 nm-1"

        assert "radiance_noise" in band.data_vars
        assert band.radiance_noise.dims == (dims["y"], dims["x"], dims["z"])
        assert band.radiance_noise.dtype == "float32"
        band.radiance_noise.attrs["long_name"] = \
            "noise of at-sensor radiance"
        band.radiance_noise.attrs["units"] = \
            "photons s-1 cm-2 sr-1 nm-1"

        if "radiance_error" in band.data_vars:
            assert band.radiance_error.dims == \
                (dims["y"], dims["x"], dims["z"])
            assert band.radiance_error.dtype == "float32"
            band.radiance_error.attrs["long_name"] = \
                "realization of Gaussian noise added onto synthetic radiance"
            band.radiance_error.attrs["units"] = \
                "photons s-1 cm-2 sr-1 nm-1"


def write_data(config, root, band_list, input_file_list):
    output_instrument_string = \
        config["instrument"].lower().replace(" ", "_")

    history_string = \
        "Created using the L1B preprocessor for RemoTeC for "\
        + f"{config['instrument']} data. " \
        + f"Used input file(s): {', '.join(input_file_list)}."

    root.attrs["history"] = history_string

    output_file_name = f"output/L1B_{output_instrument_string}.nc"

    for group_num, output_data in enumerate([root, *band_list]):
        for var in output_data.data_vars:
            output_data[var].encoding.update({"_FillValue": None})

        if group_num == 0:
            output_data.to_netcdf(
                output_file_name,
                mode="w",
                format="NETCDF4"
            )
        else:
            output_data.to_netcdf(
                output_file_name,
                mode="a",
                format="NETCDF4",
                group=f"BAND{group_num:02}"
            )


if __name__ == "__main__":
    main()
