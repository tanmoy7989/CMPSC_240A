#!/bin/bash
#SBATCH --job-name="arrsize_scale"  
#SBATCH --output="arrsize_scaling.out"  
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --export=ALL  
#SBATCH -t 00:30:00  
#SBATCH -A TG-ASC160059

cd /home/$USER/hw3

printf "size = 10^4\n"
printf "Single thread run:\n"
#SBATCH --ntasks-per-node=1
export CILK_NWORKERS=1
./innerproduct 10000
printf "\n\n"
printf "Full 24 core run:\n"
#SBATCH --ntasks-per-node=24
export CILK_NWORKERS=24
./innerproduct 10000
printf "\n-----------------------\n\n"

printf "size = 10^5\n"
printf "Single thread run:\n"
#SBATCH --ntasks-per-node=1
export CILK_NWORKERS=1
./innerproduct 100000
printf "\n\n"
printf "Full 24 core run:\n"
#SBATCH --ntasks-per-node=24
export CILK_NWORKERS=24
./innerproduct 100000
printf "\n-----------------------\n\n"

printf "size = 10^6\n"
printf "Single thread run:\n"
#SBATCH --ntasks-per-node=1
export CILK_NWORKERS=1
./innerproduct 1000000
printf "\n\n"
printf "Full 24 core run:\n"
#SBATCH --ntasks-per-node=24
export CILK_NWORKERS=24
./innerproduct 1000000
printf "\n-----------------------\n\n"

printf "size = 10^7\n"
printf "Single thread run:\n"
#SBATCH --ntasks-per-node=1
export CILK_NWORKERS=1
./innerproduct 10000000
printf "\n\n"
printf "Full 24 core run:\n"
#SBATCH --ntasks-per-node=24
export CILK_NWORKERS=24
./innerproduct 10000000
printf "\n-----------------------\n\n"

printf "size = 10^8\n"
printf "Single thread run:\n"
#SBATCH --ntasks-per-node=1
export CILK_NWORKERS=1
./innerproduct 100000000
printf "\n\n"
printf "Full 24 core run:\n"
#SBATCH --ntasks-per-node=24
export CILK_NWORKERS=24
./innerproduct 100000000
printf "\n-----------------------\n\n"

printf "size = 10^9\n"
printf "Single thread run:\n"
#SBATCH --ntasks-per-node=1
export CILK_NWORKERS=1
./innerproduct 1000000000
printf "\n\n"
printf "Full 24 core run:\n"
#SBATCH --ntasks-per-node=24
export CILK_NWORKERS=24
./innerproduct 1000000000
printf "\n-----------------------\n\n"

printf "size = 10^10\n"
printf "Single thread run:\n"
#SBATCH --ntasks-per-node=1
export CILK_NWORKERS=1
./innerproduct 10000000000
printf "\n\n"
printf "Full 24 core run:\n"
#SBATCH --ntasks-per-node=24
export CILK_NWORKERS=24
./innerproduct 10000000000
printf "\n-----------------------\n\n"

