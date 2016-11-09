#!/bin/bash
#SBATCH --job-name="scalecores"  
#SBATCH --output="scalecores.out"  
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --export=ALL  
#SBATCH -t 00:30:00  
#SBATCH -A TG-ASC160059

cd /home/$USER/hw4

GRAPHSIZE=64

export OMP_NUM_THREADS=2
./bc -p -t $GRAPHSIZE

export OMP_NUM_THREADS=4
./bc -p -t $GRAPHSIZE

export OMP_NUM_THREADS=6
./bc -p -t $GRAPHSIZE

export OMP_NUM_THREADS=8
./bc -p -t $GRAPHSIZE

export OMP_NUM_THREADS=12
./bc -p -t $GRAPHSIZE

export OMP_NUM_THREADS=16
./bc -p -t $GRAPHSIZE

export OMP_NUM_THREADS=18
./bc -p -t $GRAPHSIZE

export OMP_NUM_THREADS=24
./bc -p -t $GRAPHSIZE

export OMP_NUM_THREADS=32
./bc -p -t $GRAPHSIZE

export OMP_NUM_THREADS=36
./bc -p -t $GRAPHSIZE

export OMP_NUM_THREADS=48
./bc -p -t $GRAPHSIZE
