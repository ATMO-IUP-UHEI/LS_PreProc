# RemoTeC_PreProc_ATM

These scripts aim to enhance your experience with the retrieval software RemoTeC.

The version of RemoTeC these scripts are written for expects measurements and some initial guesses for the atmosphere.

The measurements are provided in the form of a L1B file, generated using data from an air- or space-borne observer that looks at the surface of Earth and measures reflected sunlight.
This L1B file contains multiple spatial pixels with an across-track dimension "line" and an along-track dimension "sample".
Each of these spatial pixels contains one measured spectrum.
They further contain information about the latitude and longitude, at which they are located.
Together with the UTC time of the measurement, initial guesses for the atmosphere's parameters can be automatically downloaded from online models.

The atmospheric information needs to be provided to RemoTeC in the form of an ATM file.
These scripts will generate this ATM file using the information provided in the L1B file.
These scripts require the data to be downloaded from the different databases of the correct model and just convert this downloaded data into a format that is usable by RemoTeC.
This includes calculating variables of interest from the downloaded data or interpolating onto the spatial grid defined by the L1B file.

Contained in this repository are also the download scripts for a few different databases.
These need to be run by hand once before the PreProc scripts can then use this downloaded raw data.

Right now the scripts are hardcoded to download elevation data from ASTER, download some information about the atmosphere from ERA5, and download CO2 and CH4 information from CAMS.
