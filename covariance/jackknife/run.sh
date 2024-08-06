echo "RUNNING"
echo "slurm procid = " $SLURM_PROCID
echo "slurm ntasks = " $SLURM_NTASKS
echo "slurm cpus per task = " $SLURM_CPUS_PER_TASK

#Compute 3pt correlation measurement
cmd="python ./measure_tomo_JK.py"
date
echo $cmd
$cmd
date
