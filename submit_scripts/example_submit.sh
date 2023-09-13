#!/bin/bash
#SBATCH --job-name COSMOSIS_example
#SBATCH --output=/home/dhayaa/Cosmosis/%x.log
#SBATCH --partition=chihway
#SBATCH --account=chihway
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=30:00:00
#SBATCH --mail-user=dhayaa@uchicago.edu
#SBATCH --mail-type=BEGIN,END


#You can submit this script as sbatch <script>.sh
#then check on it  as squeue -u <user>
#to cancel, check the job id using squeue and then do
#scancel <job_id>
#You can see who else is using the nodes with squeue -p chihway

conda activate /project/chihway/envs/cosmosis3

#Run cosmosis with 40 cores (the max cores on a chihway node)
mpirun -n 40 cosmosis --mpi params.ini
