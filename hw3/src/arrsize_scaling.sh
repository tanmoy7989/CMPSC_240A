#!/bin/bash
#SBATCH --job-name="arrsize_scale"  
#SBATCH --output="arrsize_scaling.out"  
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --export=ALL  
#SBATCH -t 00:30:00  
#SBATCH -A TG-ASC160059

cd /home/$USER/hw3
export CILK_WORKERS=24

printf "array size = 10^4\n"
./innerproduct 10000
printf "\n-----------------------\n\n"

printf "array size = 10^5\n"
./innerproduct 100000
printf "\n-----------------------\n\n"

printf "array size = 10^6\n"
./innerproduct 1000000
printf "\n-----------------------\n\n"

printf "array size = 10^7\n"
./innerproduct 10000000
printf "\n-----------------------\n\n"

printf "array size = 10^8\n"
./innerproduct 100000000
printf "\n-----------------------\n\n"

printf "array size = 10^9\n"
./innerproduct 1000000000
printf "\n-----------------------\n\n"

printf "array size = 10^10\n"
./innerproduct 10000000000
printf "\n-----------------------\n\n"

