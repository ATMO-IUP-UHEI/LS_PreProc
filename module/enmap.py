import os
import sys
import tarfile
import zipfile
from io import BytesIO
import xml.etree.ElementTree as ET
import numpy as np
import xarray as xr
import rasterio

import matplotlib.pyplot as plt

def main():
    tar_gz_file = os.path.join("/home/lscheidw/phd/RemoTeC_LS/data/tmp_preproc/l1b/enmap/dims_op_oc_oc-en_700632974_1.tar.gz")
    metadata, swir_band_data, vnir_band_data = get_data(tar_gz_file)

    read_spectrum(metadata, vnir_band_data, "vnir")
    read_spectrum(metadata, swir_band_data, "swir")



def get_data(tar_gz_file):
    tar = tarfile.open(tar_gz_file)

    for content in tar:
        # Find measurement. It is contained in the only zip file in the archive
        if not content.name.endswith(".ZIP"):
            continue

        # Get the zip file as a binary string and put it into memory, acting as if it were a file.
        f = tar.extractfile(content)
        content = BytesIO(f.read())
        break

    # open the zip file which was previously inside of the tar gz file
    zip = zipfile.ZipFile(content)
    for member in zip.namelist():
        # get metadata and spectral information from the desired channel(s)
        # they are contained in an xml file and tif file(s).
        # get them from the zip file and put into memory, acting as if they were files.
        if "METADATA" in member and member.endswith(".XML"):
            metadata = BytesIO(zip.open(member).read())
            metadata = ET.parse(metadata)
            continue
        if "SPECTRAL_IMAGE_SWIR" in member and member.endswith(".TIF"):
            swir_band_data = BytesIO(zip.open(member).read())
            swir_band_data = xr.open_dataset(swir_band_data, engine="rasterio")
            continue
        if "SPECTRAL_IMAGE_VNIR" in member and member.endswith(".TIF"):
            vnir_band_data = BytesIO(zip.open(member).read())
            vnir_band_data = xr.open_dataset(vnir_band_data, engine="rasterio")
            continue

    #xmlstr = ET.tostring(metadata.getroot(), encoding="utf8", method="xml").decode("utf8")

    return metadata, swir_band_data, vnir_band_data



def read_spectrum(metadata, band_data, spectral_domain):
    wavelength, fwhm, gain, offset = read_metadata(metadata, spectral_domain)
    plt.plot(wavelength, gain)
    plt.show()
    sys.exit()

    #plt.title(f"wavelength pixels {spectral_domain}")
    #plt.xlabel(f"pixel id {spectral_domain}")
    #plt.ylabel("wavelength / nm")
    #plt.scatter(range(len(wavelength)), wavelength, marker=".")
    #plt.show()

    print("WARN: not sure if lines and samples are swapped.")
    band_data = band_data.rename({"x": "lines", "y": "samples", "band": "bands"})
    band_data = band_data.transpose("samples", "lines", "bands")

    band_data["band_data"] = band_data["band_data"][...] * gain[None, None, :] + offset[None, None, :] * 1e+3
    band_data = band_data.rename({"band_data": "radiance"})

    print(band_data)

    return

    if spectral_domain == "swir":
        fig, ax = plt.subplots()
        ax.set_title(f"example SWIR spectrum")
        ax.set_xlabel("wavelength / nm")
        ax.set_ylabel("radiance")
        ax.scatter(wavelength, band_data["radiance"][200, 200, :], marker=".", color="black")
        ax.axvspan(1982, 2092, color="blue", alpha=0.25, label="CO$_2$ fit range")
        ax.axvspan(2110, 2450, color="red", alpha=0.25, label="CH$_4$ fit range")
        plt.legend(loc="upper right")

        zoomed = True
        if zoomed:
            ax.set_xlim(1970, 2462)
            plt.savefig("example_swir_spectrum_zoomed.png", dpi=600)
        else:
            plt.savefig("example_swir_spectrum.png", dpi=600)

        plt.show()

    # print("################################################################################")
    # print(spectral_domain)
    # print(wavelength)
    # print(f"xml file: {len(wavelength)}")
    # print(f"tif file: {band_data.dims['bands']}")
    # print("################################################################################")

    return



def read_metadata(metadata, spectral_domain):
    # Get valid spectral bands for the desired spectral domain
    for band in metadata.iter(f"{spectral_domain}ProductQuality"):
        try:
            nbands = int(band.find("numChannelsExpected").text)
        except:
            pass

    wavelength = []
    fwhm = []
    gain = []
    offset = []

    # Get calibration data for all bands
    for band in metadata.iter("bandID"):
        try:
            wavelength.append(float(band.find("wavelengthCenterOfBand").text))
            fwhm.append(float(band.find("FWHMOfBand").text))
            gain.append(float(band.find("GainOfBand").text))
            offset.append(float(band.find("OffsetOfBand").text))
        except:
            pass

    # Filter spectral domain
    if spectral_domain == "vnir":
        wavelength = np.asarray(wavelength[:nbands])
        fwhm = np.asarray(fwhm[:nbands])
        gain = np.asarray(gain[:nbands])
        offset = np.asarray(offset[:nbands])
    elif spectral_domain == "swir":
        wavelength = np.asarray(wavelength[-nbands:])
        fwhm = np.asarray(fwhm[-nbands:])
        gain = np.asarray(gain[-nbands:])
        offset = np.asarray(offset[-nbands:])

    return wavelength, fwhm, gain, offset



if __name__ == "__main__":
    main()
