// Lecture 2

#include "fd.h"

class LaplacianOnGrid
{
public:
  double x0, x1, y0, y1, z0, z1;
  int nRealx, nRealy, nRealz; // Number of cells in the x-y-z direction
  double dx, dy, dz;          // In case we do not want each of the cells to be a cube
  int nField;                 // The total number of temperature values in the grid. The
                              // temperature "field". field variable; spatially varying
  VDD A;
  VD b;
  VD phi;

  // constructor
  // when object gets instantiated, it gets a place in memory
  LaplacianOnGrid(double _x0, double _x1, double _y0, double _y1, double _z0,
                  double _z1, int ncell_x, int ncell_y,
                  int ncell_z) // pass in vertices of the cube
  {
    x0 = _x0;
    x1 = _x1;
    y0 = _y0;
    y1 = _y1;
    z0 = _z0;
    z1 = _z1;

    dx = (x1 - x0) / ncell_x; // size of cell's extent in that direction
    dy = (y1 - y0) / ncell_y;
    dz = (z1 - z0) / ncell_z;

    // when we do parallel processing, introduce ghost cells around the outside.
    // these are the real cells
    nRealx = ncell_x;
    nRealy = ncell_y;
    nRealz = ncell_z;

    nField = nRealx * nRealy * nRealz;

    // Acquire memory
    A.resize(nField +
             1);                   // Acquire the rows. This is a 1-based matrix. Not 0-based.
    rLOOP A[r].resize(nField + 1); // Acquire the columns.

    // Acquire memory for RHS
    b.resize(nField + 1);

    // Acquire memory for temperature (phi).
    phi.resize(nField + 1);
  }

  // Form the linear system Ax = b .... A*phi = b
  // Here, "bcs" is an array, numbered 1-6, which contains the boundary
  // conditions on the six sides of the cube. user can set the boundary
  // conditions on the command line
  void FormLS(double *bcs)
  {
    rLOOP cLOOP A[r][c] = 0.;
    rLOOP b[r] = 0.;

    // Set the boundary conditions first. Loop over the 6 planes on the surface
    // of the cube. We use the i-j-k indices for that.
    kLOOP jLOOP
    {
      int r = pid(1, j, k);
      A[r][r] = 1.;
      b[r] = bcs[1];
    } // Left BC
    kLOOP jLOOP
    {
      int r = pid(nRealx, j, k);
      A[r][r] = 1.;
      b[r] = bcs[2];
    } // Right BC
    iLOOP kLOOP
    {
      int r = pid(1, 1, k);
      A[r][r] = 1.;
      b[r] = bcs[3];
    } // Bottom BC
    iLOOP kLOOP
    {
      int r = pid(1, nRealy, k);
      A[r][r] = 1.;
      b[r] = bcs[4];
    } // Top BC

    // Front
    // Back

    // ==================
    // Interior Cells
    // ==================
    for (int i = 2; i <= nRealx - 1; ++i)
      for (int j = 2; j <= nRealy - 1; ++j)
        for (int k = 2; k <= nRealz - 1; ++k)
        {
          int p = pid(i, j, k);                                  // point in the center
          A[p][p] = -2. / dx / dx - 2. / dy / dy - 2. / dz / dz; // Cell in the center
          A[p][p - 1] = 1. / dx / dx;                            // Cell to the left
          A[p][p + 1] = 1. / dx / dx;                            // Cell to the right
          A[p][pid(i, j - 1, k)] = 1. / dy / dy;                 // Cell below (j direction)
          A[p][pid(i, j + 1, k)] = 1. / dy / dy;                 // Cell above (j direction)
        }
  }

  // Convert i-j-k address into natural cell number / point number.
  int pid(int i, int j, int k)
  {
    return i + (j - 1) * nRealx + (k - 1) * nRealx * nRealy;
  }
};

int main(int argc, char *argv[])
{

  LaplacianOnGrid F(0., 1., 0., 1., 0., 1., 10, 10, 10);

  return 0;
}