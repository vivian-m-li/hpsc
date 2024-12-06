#!/bin/bash

# -
# |
# | This is a batch script for running a MPI parallel job on Alpine
# |
# | (o) To submit this job, enter:  sbatch slurm.bat 
# | (o) To check the status of this job, enter: watch squeue -u <username>
# |
# -

# -
# |
# | Part 1: Directives
# |
# -

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:03:00
#SBATCH --partition=atesting_a100
#SBATCH --output=myslurm.out
#SBATCH --gres=gpu:1

# -
# |
# | Part 2: Loading software
# |
# -


# -
# |
# | Part 3: User scripting
# |
# -

echo "=="
echo "||"
echo "|| Begin Execution of ovenWalls in slurm batch script."
echo "||"
echo "=="

# Load the module for the pgi compiler and build the code

rm ovenWalls
module load nvhpc_sdk/2023.233
make

./ovenWalls  > tty.out


echo "=="
echo "||"
echo "|| Execution of ovenWalls in slurm batch script complete."
echo "||"
echo "=="








