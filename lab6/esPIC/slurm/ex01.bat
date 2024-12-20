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
#SBATCH --time=00:01:00
#SBATCH --partition=atesting
#SBATCH --output=ex01-%j.out

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


COMPILER := mpicc
UBUNTU   := ubuntu

# Change if on ubuntu

ifeq ($(strip $(LOGNAME)),scott)
COMPILER := mpic++ 
endif

# -----------------------------------------------------
# Make esPIC
# -----------------------------------------------------


mpicpp esPIC.cpp -g -o esPIC -std=c++11 -fopenmp

mpirun -n 4 $CODE/esPIC -n 4 ./esPIC -nPEx 2 -nPEy 2 -nCellx 10 -nCelly 10 -nPtcl 50000 -flux 200 -tEnd 2 -dt 0.01 -tPlot .2 > tty.out

echo "=="
echo "||"
echo "|| Execution of esPIC in slurm batch script complete."
echo "||"
echo "=="








