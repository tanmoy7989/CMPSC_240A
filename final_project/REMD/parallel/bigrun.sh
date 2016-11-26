#!/bin/bash
#
#$ -V
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N remd

# -pe ompi 4

date
mpirun -np 4 python rex.py bigrunsettings.txt
