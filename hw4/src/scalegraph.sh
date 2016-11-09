#!/bin/bash
#SBATCH --job-name="scalegraph"  
#SBATCH --output="scalegraph.out"  
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --export=ALL  
#SBATCH -t 00:30:00  
#SBATCH -A TG-ASC160059

cd /home/$USER/hw4

export OMP_NUM_THREADS=48

./bc -p -t 2

./bc -p -t 4

./bc -p -t 8

./bc -p -t 16

./bc -p -t 32

./bc -p -t 64

./bc -p -t 128

./bc -p -t 256

./bc -p -t 512

./bc -p -t 1024
