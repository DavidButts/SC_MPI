#!/bin/bash

#SBATCH --time=01:30:00             # limit of wall clock time - how long the job will run (same as -t)
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

if [ -z "$SLURM_SUBMIT_DIR" ] ; then
  SLURM_SUBMIT_DIR=$(dirname $0)
fi

echo $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR

numNodes=2
home="../../.."

loadModules="source ${home}/results/utilScripts/moduleLoad.sh"
echo $loadModules
$loadModules

######## Python  runs ##########

cmd="srun -n ${numNodes} python  ${home}/pingpong_python/ping_pong_blocking_double.py"
echo $cmd
$cmd

cmd="srun -n ${numNodes} python  ${home}/pingpong_python/ping_pong_non_blocking_double.py"
echo $cmd
$cmd

######## C++ runs ##########

buildCmd="mpic++ -O3 ${home}/pingpong_c++/ping_pong_blocking_double.cpp -o pingpongblockingdouble"
runCmd="srun -n ${numNodes}  pingpongblockingdouble"
echo $buildCmd
echo $runCmd
$buildCmd
$runCmd

buildCmd="mpic++ -O3 ${home}/pingpong_c++/ping_pong_non_blocking_double.cpp -o pingpongnoblockingdouble"
runCmd="srun -n ${numNodes}  pingpongnoblockingdouble"
echo $buildCmd
echo $runCmd
$buildCmd
$runCmd

