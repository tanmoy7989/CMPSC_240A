#!/bin/bash
#SBATCH --job-name="test_cilk"  
#SBATCH --output="testjob.out"  
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --export=ALL  
#SBATCH -t 00:05:00  
#SBATCH -A TG-ASC160059

cd /home/$USER/hw3

# unload old intel module (courtesy: Bruno Jacob)
module unload intel
module load intel/2016.3.210

# load cilk
module load cilk
export CILK_NWORKERS=10

# check races
cilkscreen ./innerproduct 10000 100

printf "\n--------------------------------\n"

# check empirical parallelism
cilkview ./innerproduct 10000 100



