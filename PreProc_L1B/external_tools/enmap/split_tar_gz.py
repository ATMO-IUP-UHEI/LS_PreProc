import sys
import os
import tarfile
from io import BytesIO
import xml.etree.ElementTree as ET

inp_tar = sys.argv[1]
inp_folder = inp_tar.strip(".tar.gz")

if not os.path.isfile(inp_tar):
    sys.exit(f"{inp_tar} does not exist.")

if not inp_tar.endswith(".tar.gz"):
    sys.exit(f"{inp_tar} is not a .tar.gz file.")

if not os.path.isdir(inp_folder):
    os.system(f"tar -xf {inp_tar}")

iif_files = os.listdir(os.path.join(inp_folder, "iif"))

for index, iif_file in enumerate(iif_files):
    index += 1
    out_folder = f"{inp_folder}_{index:02}"

    iif_data = ET.parse(os.path.join(inp_folder, "iif", iif_file)).getroot()

    for path in iif_data.iter("path"):
        path = path.text

    zip_path = path.split("/")[1:]

    match len(zip_path):
        # Sometimes, all .zip files are in the same folder. Sometimes, this
        # folder contains subdirectories, which contain one .zip file each.
        case 1:
            zip_file = os.listdir(os.path.join(inp_folder, *zip_path))[index-1]
        case 2:
            zip_file = os.listdir(os.path.join(inp_folder, *zip_path))[0]
        case _:
            sys.exit("Unknown archive structure.")

    os.system(f"mkdir {out_folder}")
    os.system(f"cp {inp_folder}/readme.html {out_folder}")
    os.system(f"cp -r {inp_folder}/tools/ {out_folder}")
    os.system(f"mkdir {out_folder}/iif/")
    os.system(f"cp {inp_folder}/iif/{iif_file} {out_folder}/iif/")
    os.system(f"mkdir {out_folder}/{zip_path[0]}")
    match len(zip_path):
        case 1:
            os.system(f"cp {inp_folder}/{zip_path[0]}/{zip_file}"
                      + f" {out_folder}/{zip_path[0]}/")
        case 2:
            os.system(f"mkdir {out_folder}/{zip_path[0]}/{zip_path[1]}")
            os.system(f"cp {inp_folder}/{zip_path[0]}/{zip_path[1]}/{zip_file}"
                      + f" {out_folder}/{zip_path[0]}/{zip_path[1]}/")
        case _:
            sys.exit("Unknown archive structure.")
    os.system(f"tar -cf {out_folder}.tar.gz {out_folder}")
    os.system(f"rm -r {out_folder}")

os.system(f"rm -r {inp_folder}")
