#!/bin/bash
#SBATCH --job-name DESY3_Covariance
#SBATCH --partition=kicp
#SBATCH --account=kicp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=36:00:00
#SBATCH --output=/home/dhayaa/DECADE/CosmicShearCosmosis/covariance/analytic/June25th2024_DESY3/%x.log
#SBATCH --mail-user=dhayaa@uchicago.edu
#SBATCH --mail-type=BEGIN,END

cd /home/dhayaa/DECADE/CosmicShearCosmosis/covariance/analytic/June25th2024_DESY3/

python -u Setup.py
