#!/bin/bash
#SBATCH --partition=cpu-single         # Partition (change if needed)
#SBATCH --ntasks=1                     # Number of tasks (1 process)
#SBATCH --time=01:00:00                # Time limit (HH:MM:SS)
#SBATCH --mem=200G                      # Memory required
#SBATCH --job-name=PreProc_L1B         # Job name
#SBATCH --output=slurm/preproc/l1b.log # Log file (stores stdout/stderr)
#SBATCH --error=slurm/preproc/l1b.log  # Error file (if errors occur)
python /home/hd/hd_hd/hd_oc152/software/RemoTeC/LS_PreProc/PreProc_L1B/main.py $1
