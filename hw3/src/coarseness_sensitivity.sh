#!/bin/bash
#SBATCH --job-name="coarse_sensitivity"  
#SBATCH --output="coarseness_sensitivity.out"  
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --export=ALL  
#SBATCH -t 00:30:00  
#SBATCH -A TG-ASC160059

cd /home/$USER/hw3

printf "Coarseness = 1\n"
printf "Single thread run--\n"
#SBATCH --ntasks-per-node=1
export CILK_NWORKERS=1
./innerproduct 1000000 1
printf "\n\n"
printf "Full 24 core run--\n"
#SBATCH --ntasks-per-node=24
export CILK_NWORKERS=24
./innerproduct 1000000 1
printf "\n-------------------\n\n"

printf "Coarseness = 10\n"
printf "Single thread run--\n"
#SBATCH --ntasks-per-node=1
export CILK_NWORKERS=1
./innerproduct 1000000 10
printf "\n\n"
printf "Full 24 core run--\n"
#SBATCH --ntasks-per-node=24
export CILK_NWORKERS=24
./innerproduct 1000000 10
printf "\n-------------------\n\n"

printf "Coarseness = 100\n"
printf "Single thread run--\n"
#SBATCH --ntasks-per-node=1
export CILK_NWORKERS=1
./innerproduct 1000000 100
printf "\n\n"
printf "Full 24 core run--\n"
#SBATCH --ntasks-per-node=24
export CILK_NWORKERS=24
./innerproduct 1000000 100
printf "\n-------------------\n\n"

printf "Coarseness = 1000\n"
printf "Single thread run--\n"
#SBATCH --ntasks-per-node=1
export CILK_NWORKERS=1
./innerproduct 1000000 1000
printf "\n\n"
printf "Full 24 core run--\n"
#SBATCH --ntasks-per-node=24
export CILK_NWORKERS=24
./innerproduct 1000000 1000
printf "\n-------------------\n\n"

printf "Coarseness = 10000\n"
printf "Single thread run--\n"
#SBATCH --ntasks-per-node=1
export CILK_NWORKERS=1
./innerproduct 1000000 10000
printf "\n\n"
printf "Full 24 core run--\n"
#SBATCH --ntasks-per-node=24
export CILK_NWORKERS=24
./innerproduct 1000000 10000
printf "\n-------------------\n\n"

printf "Coarseness = 100000\n"
printf "Single thread run--\n"
#SBATCH --ntasks-per-node=1
export CILK_NWORKERS=1
./innerproduct 1000000 100000
printf "\n\n"
printf "Full 24 core run--\n"
#SBATCH --ntasks-per-node=24
export CILK_NWORKERS=24
./innerproduct 1000000 100000
printf "\n-------------------\n\n"

printf "Coarseness = 1000000\n"
printf "Single thread run--\n"
#SBATCH --ntasks-per-node=1
export CILK_NWORKERS=1
./innerproduct 1000000 1000000
printf "\n\n"
printf "Full 24 core run--\n"
#SBATCH --ntasks-per-node=24
export CILK_NWORKERS=24
./innerproduct 1000000 1000000
printf "\n-------------------\n\n"

