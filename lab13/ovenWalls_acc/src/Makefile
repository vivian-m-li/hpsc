# -----------------------------------------------------
# Set compiler based on platform
# -----------------------------------------------------

# module load gcc
#COMPILER := g++ -O3 -fopenmp -fopenacc -std=c++11

# module load nvhpc_sdk/2023.233
COMPILER := pgc++ -O3 -acc -gpu=ccnative -std=c++11


#COMPILER := mpicc -O3 -acc -std=c++11

# -----------------------------------------------------
# Make solvers
# -----------------------------------------------------

ovenWalls: ovenWalls.cpp ovenWalls.h plotter.h geom.h Makefile
	$(COMPILER) ovenWalls.cpp -o ovenWalls 



