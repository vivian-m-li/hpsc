#!/bin/bash
#SBATCH --job-name=lab00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vili4418@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --partition=acompile
#SBATCH --output=lab00-%j.out

module purge
module load gcc/11.2.0
module load openmpi/4.1.1

mpirun -n 4 ./main.exe # run the compiled program on 4 processors