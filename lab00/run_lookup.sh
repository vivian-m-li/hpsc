#!/bin/bash
#SBATCH --job-name=lab00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vili4418@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=00:10:00
#SBATCH --partition=atesting
#SBATCH --output=lab00-%j.out

module purge
module load gcc/10.3 openmpi

export SLURM_EXPORT_ENV=ALL

mpirun -np $SLURM_NTASKS /home/vili4418/hpsc/lab00/main.exe > output.txt
