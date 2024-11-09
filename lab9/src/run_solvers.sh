#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=00:05:00
#SBATCH --partition=atesting
#SBATCH --output=eofiles/%j.out

module purge
module load intel
module load impi

echo "Begin execution"

mpirun -n 1 ./solvers -nPEx 1 -nPEy 1 -nCellx 64 -nCelly 64 -solver cg -nl nr -c0 5. -tau 1. -r .3

echo "Done"
