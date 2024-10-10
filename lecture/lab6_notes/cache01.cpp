#include "main.h"

// Double 1D Array
double *arrayDouble(int n) {
  double *myArray;
  myArray = new double[n];
  for (int i = 0; i < n; ++i) myArray[i] = 1.;
  return myArray;
}