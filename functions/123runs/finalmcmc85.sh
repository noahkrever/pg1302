#!/bin/sh
#SBATCH --account=astro
#SBATCH --job-name= finalmcmc85
#SBATCH --mail-type ALL
#SBATCH --mail-user ndk2115@columbia.edu
#SBATCH -c 1
#SBATCH --time=0-04:00
#SBATCH --mem-per-cpu=1gb
    
module load anaconda
source activate /rigel/astro/users/ndk2115/python-envs/pg
python3 finalmcmc85.py