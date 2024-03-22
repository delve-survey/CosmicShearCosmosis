#!/bin/bash
#SBATCH --job-name COSMOSIS_example
#SBATCH --output=/home/dhayaa/Cosmosis/%x.log
#SBATCH --partition=caslake
#SBATCH --account=pi-chihway
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=30:00:00
#SBATCH --mail-user=dhayaa@uchicago.edu
#SBATCH --mail-type=BEGIN,END


#You can submit this script as sbatch <script>.sh
#then check on it  as squeue -u <user>
#to cancel, check the job id using squeue and then do
#scancel <job_id>
#You can see who else is using the nodes with squeue -p chihway

conda activate /project/chihway/envs/cosmosis3

#If you are running a yml file then somehow this step is needed
#to expand env variables in yml file

YML=$REPO_DIR/yml_files/delve-campaign.yml
envsubst < $YML > $TMPDIR/testing.yml


#INI file works as usual
INI=$REPO_DIR/ini_files/des-y3-shear.ini


#Campaign file runs
mpirun -n 48 cosmosis-campaign --mpi $TMPDIR/testing.yml --run fiducial


#If ini file
mpirun -n 48 cosmosis --mpi $INI