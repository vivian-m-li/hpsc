#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=00:10:00
#SBATCH --partition=atesting
#SBATCH --output=eofiles/%j.out

module purge
module load intel
module load impi
module load python
module load anaconda

# To test the restart capabilities
# mpirun -n 16 ./transientDiffusion -nPEx 4 -nPEy 4 -nCellx  4 -nCelly 4 -solver jacobi -tEnd .02   -dt .001 -tPlot .001 | tee tmp

# read -p "First 16PE done.  Press enter to continue"

# mpirun -n 9 ./transientDiffusion  -nPEx 3 -nPEy 3 -nCellx  5 -nCelly 5 -solver jacobi -tEnd .04   -dt .001 -tPlot .001 | tee tmp

# read -p "First 9PE done.  Press enter to continue"

# mpirun -n 1 ./transientDiffusion  -nPEx 1 -nPEy 1 -nCellx  13 -nCelly 13 -solver jacobi -tEnd .06   -dt .001 -tPlot .001 -restart     
  
# read -p "1PE  done.  Press enter to continue"

# mpirun -n 16 ./transientDiffusion -nPEx 4 -nPEy 4 -nCellx  4 -nCelly 4 -solver jacobi -tEnd .10   -dt .001 -tPlot .001 -restart     

# read -p "Last 16PE  done.  Press enter to continue"

rm -rf runA runB1 runB2 pc_restart pc_end
mkdir runA runB1 runB2

echo "Run A"

mpirun -n 16 ./transientDiffusion -nPEx 4 -nPEy 4 -nCellx  4 -nCelly 4 -solver jacobi -tEnd .05 -dt .001 -tPlot .001 | tee runA_out

python writePlotCmd.py
mv *.plt pc* runA_out runA/

echo "Run B1"

mpirun -n 16 ./transientDiffusion -nPEx 4 -nPEy 4 -nCellx  4 -nCelly 4 -solver jacobi -tEnd .025 -dt .001 -tPlot .001 | tee runB1_out

python writePlotCmd.py
mv *.plt pc* runB1_out runB1/

echo "Run B2"

mpirun -n 16 ./transientDiffusion -nPEx 4 -nPEy 4 -nCellx  4 -nCelly 4 -solver jacobi -tEnd .05 -dt .001 -tPlot .001 -restart | tee runB2_out

python writePlotCmd.py
mv *.plt pc* runB2_out runB2/

python writeAllPlots.py

echo "Done"
