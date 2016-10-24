#!/bin/bash
#SBATCH --job-name="coarse_sensitivity"  
#SBATCH --output="coarseness_sensitivity.out"  
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --export=ALL  
#SBATCH -t 00:30:00  
#SBATCH -A TG-ASC160059

cd /home/$USER/hw3
export CILK_NWORKERS=24

printf "Coarseness = 1\n"
./innerproduct 100000000 1
printf "\n-------------------\n\n"

printf "Coarseness = 10\n"
./innerproduct 100000000 10
printf "\n-------------------\n\n"

printf "Coarseness = 100\n"
./innerproduct 100000000 100
printf "\n-------------------\n\n"

printf "Coarseness = 1000\n"
./innerproduct 100000000 1000
printf "\n-------------------\n\n"

printf "Coarseness = 10000\n"
./innerproduct 100000000 10000
printf "\n-------------------\n\n"

printf "Coarseness = 100000\n"
./innerproduct 100000000 100000
printf "\n-------------------\n\n"

printf "Coarseness = 1000000\n"
./innerproduct 100000000 1000000
printf "\n-------------------\n\n"

printf "Coarseness = 10000000\n"
./innerproduct 100000000 10000000
printf "\n-------------------\n\n"

printf "Coarseness = 100000000\n"
./innerproduct 100000000 100000000
printf "\n-------------------\n\n"
