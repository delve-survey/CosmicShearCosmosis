#!/bin/bash
#SBATCH --job-name DELVE_Covariance
#SBATCH --partition=kicp
#SBATCH --account=kicp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=36:00:00
#SBATCH --output=/home/dhayaa/DECADE/CosmicShearCosmosis/covariance/analytic/20240806/%x.log
#SBATCH --mail-user=dhayaa@uchicago.edu
#SBATCH --mail-type=BEGIN,END

cd /home/dhayaa/DECADE/CosmicShearCosmosis/covariance/analytic/20240806

python -u Setup.py
