#!/bin/bash
#SBATCH --job-name="ncores_scale"  
#SBATCH --output="ncores_scaling.out"  
#SBATCH --partition=compute
#SBATCH --export=ALL  
#SBATCH -t 00:30:00  
#SBATCH -A TG-ASC160059

cd /home/$USER/hw3

printf "NCores = 1\n"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
export CILK_NWORKERS=1
./innerproduct 1000000
printf "\n--------------------------------\n\n"

printf "NCores = 2\n"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
export CILK_NWORKERS=2
./innerproduct 1000000
printf "\n--------------------------------\n\n"

printf "NCores = 4\n"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
export CILK_NWORKERS=4
./innerproduct 1000000
printf "\n--------------------------------\n\n"

printf "NCores = 8\n"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
export CILK_NWORKERS=8
./innerproduct 1000000
printf "\n--------------------------------\n\n"

printf "NCores = 16\n"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
export CILK_NWORKERS=16
./innerproduct 1000000
printf "\n--------------------------------\n\n"

printf "NCores = 32\n"
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=16
export CILK_NWORKERS=32
./innerproduct 1000000
printf "\n--------------------------------\n\n"

printf "NCores = 64\n"
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
export CILK_NWORKERS=64
./innerproduct 1000000
printf "\n--------------------------------\n\n"



