#!/bin/bash
#SBATCH --partition=cpu-single
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=1gb
#SBATCH --export=NONE
#SBATCH --output=slurm/preproc/download_egg4_%j.out
#SBATCH --error=slurm/preproc/download_egg4_%j.out
python3 ~/sds/software/PreProc/PreProc_ATM/download_scripts/download_egg4.py
