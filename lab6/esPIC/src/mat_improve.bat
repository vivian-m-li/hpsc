#!/bin/bash

# -
# |
# | This is a batch script for running a MPI parallel job on Summit
# |
# | (o) To submit this job, enter:  sbatch --export=CODE='/home/scru5660/HPSC/codes/fd_mpi/src' ex_01.bat 
# |
# | (o) To check the status of this job, enter: squeue -u <username>
# |
# -

# -
# |
# | Part 1: Directives
# |
# -

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=00:30:00
#SBATCH --partition=atesting
#SBATCH --output=mat_improve-%j.out


# -
# |
# | Part 2: Loading software
# |
# -

module purge
module load intel
module load impi 
module load valgrind/3.17.0
module use /projects/scru5660/public/software/module
module load kcachegrind/0.7.4

# -
# |
# | Part 3: User scripting
# |
# -

echo "=="
echo "||"
echo "|| Begin Execution of esPIC in slurm batch script."
echo "||"
echo "=="

mpirun -n 4 ./esPIC -nPEx 2 -nPEy 2 -nCellx 12 -nCelly 12 -flux 1000. -vx_bdy 1. -npHat 80. -tEnd 1 -dt .01 > tty.out
 
# mpirun -n 4 valgrind --tool=callgrind --log-file="mat_tty/callgrind.log" --callgrind-out-file="mat_tty/callgrind.out" ./esPIC -nPEx 2 -nPEy 2 -nCellx 6 -nCelly 6 -flux 1000. -vx_bdy 1. -npHat 80. -tEnd 1 -dt .01 > callgrind_tty.out

# mpirun -n 4 valgrind --tool=cachegrind --log-file="mat_tty/cachegrind.log" --cachegrind-out-file="mat_tty/cachegrind.out" ./esPIC -nPEx 2 -nPEy 2 -nCellx 6 -nCelly 6 -flux 1000. -vx_bdy 1. -npHat 80. -tEnd 1 -dt .01 > cachegrind_tty.out

# ./writePlotCmd_py3.py

echo "=="
echo "||"
echo "|| Execution of esPIC in slurm batch script complete."
echo "||"
echo "=="

rm -rf *.plt
rm -rf *.sed








