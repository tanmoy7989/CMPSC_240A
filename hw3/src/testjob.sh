#!/bin/bash
#SBATCH --job-name="int_erprod"  
#SBATCH --output="testjob.out"  
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --export=ALL  
#SBATCH -t 00:05:00  
#SBATCH -A TG-ASC160059

cd /home/$USER/hw3
export CILK_NWORKERS=10
./innerproduct 1000 20



