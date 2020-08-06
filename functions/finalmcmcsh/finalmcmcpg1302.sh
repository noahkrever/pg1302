#!/bin/sh
#SBATCH --account=astro
#SBATCH --job-name= finalmcmcpg1302
#SBATCH --mail-type ALL
#SBATCH --mail-user ndk2115@columbia.edu
#SBATCH -c 1
#SBATCH --time=0-02:00
#SBATCH --mem-per-cpu=1gb
    
module load anaconda
source activate /rigel/astro/users/ndk2115/python-envs/pg
python3 finalmcmcpg1302.py