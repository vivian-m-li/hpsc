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

# TODO: change omp_set_num_threads() in solvers.cpp for each run
# mpirun -n 8 ./solvers -nPEx 2 -nPEy 4 -nCellx 45 -nCelly 45 -solver cg -nl nr -c0 5. -tau 1. -r .3 # 2 threads, 8 ranks

# mpirun -n 4 ./solvers -nPEx 2 -nPEy 2 -nCellx 64 -nCelly 64 -solver cg -nl nr -c0 5. -tau 1. -r .3 # 4 threads, 4 ranks

mpirun -n 2 ./solvers -nPEx 2 -nPEy 1 -nCellx 90 -nCelly 90 -solver cg -nl nr -c0 5. -tau 1. -r .3 # 8 threads, 2 ranks

echo "Done"
