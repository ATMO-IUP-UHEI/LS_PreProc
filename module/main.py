import sys

from instruments import prisma
from instruments import enmap
from instruments import dlr_hyspex


def main():
    config_prisma = {
        "instrument":
            "PRISMA",
        "path":
            "/home/lscheidw/phd/RemoTeC_LS/data/tmp_preproc/l1b/prisma",
        "l1b":
            "PRS_L1_STD_OFFL_20201010072033_20201010072038_0001.zip",
        "l2b":
            "PRS_L2B_STD_20201010072033_20201010072038_0001.zip",
    }

    config_enmap = {
        "instrument":
            "EnMAP",
        "path":
            "/home/lscheidw/phd/RemoTeC_LS/data/tmp_preproc/l1b/enmap",
        "tar_gz":
            "dims_op_oc_oc-en_700632974_1.tar.gz",

    }

    config_dlr_hyspex = {
        "instrument":
            "DLR HySpex",
        "path":
            "/home/lscheidw/phd/RemoTeC_LS/data/tmp_preproc/l1b/dlr_hyspex",
        "l1b":
            "20180607_Pawlowice_01_SWIR_320me_SN3510_FOVx2_raw.nc",
    }

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

    config_list = [
        config_prisma,
        config_enmap,
        config_dlr_hyspex,
    ]

    for config in config_list:
        print(config["instrument"])
        root, band_list, input_file_list = get_data(config, dims)
        set_attributes(root, band_list, dims)
        write_data(config, root, band_list, input_file_list)

    return


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
    root.time.attrs["long_name"] = \
        "UTC date and time of measurement in ISO 8601 standard"
    root.time.attrs["units"] = \
        "YYYY-MM-DDThh:mm:ssZ"

    assert "latitude" in root.data_vars
    assert root.latitude.dims == (dims["y"], dims["x"])
    root.latitude.attrs["long_name"] = \
        "latitude at pixel center"
    root.latitude.attrs["units"] = \
        "degrees north"

    assert "longitude" in root.data_vars
    assert root.longitude.dims == (dims["y"], dims["x"])
    root.longitude.attrs["long_name"] = \
        "longitude at pixel center"
    root.longitude.attrs["units"] = \
        "degrees east"

    assert "solar_zenith_angle" in root.data_vars
    assert root.solar_zenith_angle.dims == (dims["y"], dims["x"])
    root.solar_zenith_angle.attrs["long_name"] = \
        "solar zenith angle"
    root.solar_zenith_angle.attrs["units"] = \
        "degrees"

    if "solar_azimuth_angle" in root.data_vars:
        assert root.solar_azimuth_angle.dims == (dims["y"], dims["x"])
        root.solar_azimuth_angle.attrs["long_name"] = \
            "solar azimuth angle"
        root.solar_azimuth_angle.attrs["units"] = \
            "degrees"

    assert "viewing_zenith_angle" in root.data_vars
    assert root.viewing_zenith_angle.dims == (dims["y"], dims["x"])
    root.viewing_zenith_angle.attrs["long_name"] = \
        "viewing zenith angle"
    root.viewing_zenith_angle.attrs["units"] = \
        "degrees"

    if "viewing_azimuth_angle" in root.data_vars:
        assert root.viewing_azimuth_angle.dims == (dims["y"], dims["x"])
        root.viewing_azimuth_angle.attrs["long_name"] = \
            "viewing azimuth angle"
        root.viewing_azimuth_angle.attrs["units"] = \
            "degrees"

    if "observer_altitude" in root.data_vars:
        assert root.observer_altitude.dims == (dims["y"], dims["x"])
        root.observer_altitude.attrs["long_name"] = \
            "observer altitude above ground"
        root.observer_altitude.attrs["units"] = \
            "m"

    for band in band_list:
        assert len(band[dims["x"]]) == len(root[dims["x"]])
        assert len(band[dims["y"]]) == len(root[dims["y"]])

        assert "wavelength" in band.data_vars
        assert band.wavelength.dims == (dims["z"],)
        band.wavelength.attrs["long_name"] = "wavelength"
        band.wavelength.attrs["units"] = "nm"

        assert "radiance" in band.data_vars
        assert band.radiance.dims == (dims["y"], dims["x"], dims["z"])
        band.radiance.attrs["long_name"] = \
            "at-sensor radiance"
        band.radiance.attrs["units"] = \
            "photons s-1 cm-2 sr-1 nm-1"

        assert "radiance_noise" in band.data_vars
        assert band.radiance_noise.dims == (dims["y"], dims["x"], dims["z"])
        band.radiance_noise.attrs["long_name"] = \
            "noise of at-sensor radiance"
        band.radiance_noise.attrs["units"] = \
            "photons s-1 cm-2 sr-1 nm-1"

        if "radiance_error" in band.data_vars:
            assert band.radiance_error.dims == \
                (dims["y"], dims["x"], dims["z"])
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
