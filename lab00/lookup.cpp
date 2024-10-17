#include <mpi.h>
#include <iostream>
#include <string>
using std ::cout;
using std ::endl;
using std ::string;

void printArray(string arrayName, int n, int *array, int myPE)
{
  for (int i = 0; i < n; ++i)
    cout << "myPE: " << myPE << " " << arrayName << "[" << i << "] = " << array[i] << endl;
}

double lookupVal(int n, double *x, double *y, double xval)
{
  for (int i = 0; i < n; ++i)
    if (xval >= x[i] && xval <= x[i + 1])
      return y[i] + (xval - x[i]) * (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
  return 0.;
}

int main(int argc, char *argv[])
{

  MPI_Init(&argc, &argv);
  int numPE, myPE;

  MPI_Comm_size(MPI_COMM_WORLD, &numPE);
  MPI_Comm_rank(MPI_COMM_WORLD, &myPE);

  int n = 100;
  double x[n], y[n];

  for (int i = 0; i < n; ++i)
  {
    x[i] = i;
    y[i] = i * i;
  }

  int m = 20;
  int s = m / numPE; // this only works if m % numPE == 0; otherwise, we don't fill out all the values
  int xval[m];
  int yval[s];   // each PE only needs to handle array of length s, then they get gathered on the last step
  int output[m]; // output array

  // populate xval with values
  for (int i = 0; i < m; ++i)
  {
    xval[i] = 2. * i;
  }

  // equally split up the work amongst PEs - each PE only needs to populate yval with s values
  for (int i = 0; i < numPE; ++i)
    if (i == myPE)
      for (int j = 0; j < s; ++j)
      {
        yval[j] = lookupVal(n, x, y, xval[i * s + j]);
      }

  // gather the results from each yval array of length s into an output array of length m
  MPI_Gather(yval, s, MPI_INT,
             output, s, MPI_INT, 0, MPI_COMM_WORLD);

  printArray("yval", s, yval, myPE);
  printArray("output", m, output, myPE);

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}
