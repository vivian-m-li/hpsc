#include "main.h"

//--------------------------
// Double 1D Array
//--------------------------
double *arrayDouble(int n) {
  double *myArray;
  myArray = new double[n];
  for (int i = 0; i < n; ++i) myArray[i] = 1.;
  return myArray;
}

//--------------------------
// Double 2D Array
//--------------------------
double **arrayDouble(int nRows, int nCols) {
  double *myArray;
  myArray = new double[nRows * nCols];

  // Create a pointer that points to the beginning of each new row
  double **myArray_ptr;
  myArray_ptr = new double *[nRows];

  // Populate myArray_ptr
  int count = 0;
  for (int row = 0; row < nRows; ++row) {
    // myArray_ptr points to the memory location of myArray
    myArray_ptr[row] = &myArray[count * nCols];
    ++count;
  }

  // Initialize
  for (int i = 0; i < nRows; ++i)
    for (int j = 0; j < nCols; ++j) myArray_ptr[i][j] = 10.;

  // Return the pointer
  return myArray_ptr;
}

//--------------------------
// Integer 2D Array
//--------------------------
int **arrayInt(int nRows, int nCols) {
  int *myArray;
  myArray = new int[nRows * nCols];

  // Create a pointer that points to the beginning of each new row
  int **myArray_ptr;
  myArray_ptr = new int *[nRows];

  // Populate myArray_ptr
  int count = 0;
  for (int row = 0; row < nRows; ++row) {
    // myArray_ptr points to the memory location of myArray
    myArray_ptr[row] = &myArray[count * nCols];
    ++count;
  }

  // Initialize
  for (int i = 0; i < nRows; ++i)
    for (int j = 0; j < nCols; ++j) myArray_ptr[i][j] = 10.;

  // Return the pointer
  return myArray_ptr;
}

//--------------------------
// Standard Ax=b
//--------------------------
void matVec(double **A, double *x, double *b, int n) {
  for (int i = 1; i <= n; ++i) {
    b[i] = 0.;
    // standard, not compressed format
    for (int j = 1; j <= n; ++j) b[i] += A[i][j] * x[j];  // going down j and going across i
  }
}

// convert matrix to pointer matrix

//--------------------------
// Sparse Ax=b
//--------------------------
void matVec(double **A, int **J, double *x, double *b, int n, int bandwidth) {
  for (int i = 1; i <= n; ++i) {
    b[i] = 0.;
    // standard, not compressed format
    for (int j = 1; j <= bandwidth; ++j)
      b[i] += A[i][j] * x[J[i][j]];  // J[i][j] stores actual column number of what's in A[i][j]
  }
}

//--------------------------
// Main
//--------------------------

int main(int argc, char *argv[]) {
  // Write a finite solver that doesn't solve anything but just exercises the code
  int n, r;

  // Command line inputs:
  for (int count = 0; count < argc; ++count) {
    if (!strcmp(argv[count], "-n")) n = atoi(argv[count + 1]);  // how many nodes in the x-direction
    if (!strcmp(argv[count], "-r")) n = atoi(argv[count + 1]);  // how many repeats of the experiment
  }

  int N = n * n;  // square mesh
  //--------------------------
  // Full-matrix format
  //--------------------------

  // Get memory
  double *x = arrayDouble(N + 1);  // +1 because we are 1-based
  double *b = arrayDouble(N + 1);
  double **A = arrayDouble(N + 1, N + 1);  // pointer to a pointer

  // Matrix-vector product
  for (int count = 0; count < r; ++count) {
    matVec(A, x, b, N);
  }

  //--------------------------
  // Sparse-matrix format
  //--------------------------

  // Get memory
  int bandwidth = 5;
  double **ACoef = arrayDouble(N + 1, bandwidth + 1);
  int **JCoef = arrayInt(N + 1, bandwidth + 1);

  // Matrix-vector product
  for (int count = 0; count < r; ++count) {
    matVec(ACoef, JCoef, x, b, N, bandwidth);
  }
}