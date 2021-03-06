#!/bin/bash
 
#SBATCH --time=00:30:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1
#SBATCH --ntasks=2                 # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --tasks-per-node=40
#SBATCH --ntasks-per-socket=20
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --sockets-per-node=2
#SBATCH --cores-per-socket=20
#SBATCH --threads-per-core=1
#SBATCH -m block:cyclic          # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=MaxPerNode                   # memory required per node - amount of memory (in bytes)
#SBATCH --job-name parJob      # you can give your job a name for easier identification (same as -J)
 
 
########## Command Lines to Run ##########
  
cd $SLURM_SUBMIT_DIR 

numNodes=2

export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=close              

ompNum="OMP_NUM_THREADS=${OMP_NUM_THREADS}"
 
cmd="srun -n ${numNodes}  poissonSolvePar"

echo $cmd
echo $ompNum

$cmd             ### call your executable
