#!/bin/bash
#SBATCH --partition=cpu-single         # Partition (change if needed)
#SBATCH --ntasks=1                     # Number of tasks (1 process)
#SBATCH --time=01:00:00                # Time limit (HH:MM:SS)
#SBATCH --mem=200G                     # Memory required
#SBATCH --job-name=PreProc_ATM         # Job name
#SBATCH --output=slurm/preproc/atm.log # Log file (stores stdout/stderr)
#SBATCH --error=slurm/preproc/atm.log  # Error file (if errors occur)
python3 ~/sds/software/PreProc/PreProc_ATM/main.py
