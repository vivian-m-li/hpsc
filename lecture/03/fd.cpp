// Lecture 3

#include "fd.h"

class LaplacianOnGrid {
 public:
  double x0, x1;  // x-extent of the domain on this PE
  double y0, y1;  // y-extent of the domain on this PE
  double dx;      // Size of the cell in the x-direction
  double dy;      // Size of the cell in the y-direction
  int nRealx;     // Number of real cells in the x-direction
  int nRealy;     // Number of real cells in the y-direction
  int nField;  // Total number of unknowns on this PE, including the ghost cells
  VDD A;       //  Matrix
  VD x, y;   // Stores the x-y locations of each point. x[natural node number],
             // y[natural node number]
  VD phi;    // Our temprature solution
  VD b;      // Right hand side in our linear system
  int myPE;  // myPE number

  // Constructor: Initialize values

  LaplacianOnGrid(double _x0, double _x1, double _y0, double _y1, int ncell_x,
                  int ncell_y, mpiInfo &myMPI) {
    // Record incoming parameters
    // Set up nRealx, nRealy, and nField
    // Compute grid spacings dx and dy

    x0 = _x0;
    x1 = _x1;
    y0 = _y0;
    y1 = _y1;
    nRealx = ncellx;
    nRealy = ncelly;
    nField = (nRealx + 2) * (nRealy + 2);
    dx = (x1 - x0) / ncell_x;  // Size of the cells
    dy = (y1 - y0) / ncell_y;
  }

  // Convert i-j-k address into natural cell number / point number.
  int pid(int i, int j, int k) {
    return i + (j - 1) * nRealx + (k - 1) * nRealx * nRealy;
  }
};

int main(int argc, char *argv[]) {
  // Establish an "mpiInfo" object, myMPI
  mpiInfo myMPI;  // has domain decomposition for how we split up the mesh

  // Begin parallel region
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &myMPI.numPEs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myMPI.myPEs);

  // Get user input
  int nPEx, nPEy, nCellx, nCelly;
  for (int count = 0; count < argc; ++count) {
    if (!strcmp(argv[count], "-nPEx")) nPEx = atoi(argv[count + 1]);
    if (!strcmp(argv[count], "-nPEy")) nPEy = atoi(argv[count + 1]);
    if (!strcmp(argv[count], "-nCellx")) nCellx = atoi(argv[count + 1]);
    if (!strcmp(argv[count], "-nCelly")) nCelly = atoi(argv[count + 1]);
  }

  // Assign different parts of the global mesh to the processes
  myMPI.GridDecomposition(nPEx, nPEy, nCellx, nCelly);

  // Simulate a unit square - x0, x1 and y0, y1 are hard coded
  double length = 1.;
  double x0, x1;
  double y0, y1;

  x0 = (length / myMPI.nPEx) * myMPI.iPE;
  x1 = x0 + length / myMPI.nPEx;
  y0 = (length / myMPI.nPEy) * myMPI.jPE;
  y1 = y0 + length / myMPI.nPEy;

  // Form the solver object
  LaplacianOnGrid F(x0, x1, y0, y1, nCellx, nCelly,
                    myMPI);  // Forms the system A*phi = b for one PE

  // have each PE wait until all of their matrices are formed before they can
  // communicate
  MPI_Barrier(MPI_COMM_WORLD);

  // Each PE forms their linear system

  //            Max its RHS  Sol'n  Parallel  2 --> Jacobi, not Gauss-Seidel
  F.GS_or_Jacobi(500, F.b, F.phi, myMPI, 2);

  MPI_Barrier(MPI_COMM_WORLD);

  // Each PE writes their solution to their separate plot file

  F.plot("F_phi", F.phi, myMPI);

  MPI_Finalize();
  return 0;
}