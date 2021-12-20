#!/bin/bash

#SBATCH --account=def-matt
#SBATCH --job-name=DVR
#SBATCH --output=DVR_%j.out
#SBATCH --error=DVR_%j.err

#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=3G     # memory per CPU core
#SBATCH --time=0-01:00:00    # time (DD-HH:MM:SS)

# Execute the program
srun ./DVR-example_advanced.sh
