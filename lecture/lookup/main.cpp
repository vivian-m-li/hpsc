#include <iostream>
using std ::cout; // print
using std ::endl; // line feed

// star means it's going to pass a address of the variable
//                      table of x-y values   interpolant x-value
//                          |          |            |
double lookupVal(int n, double *x, double *y, double xval) // Performs linear interpolation using a table of values
{
  for (int i = 0; i < n - 1; ++i)
    if (xval >= x[i] && xval <= x[i + 1])
      return y[i] + (xval - x[i]) * (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
}

int main() // returns an int
{
  cout << "Hello world!" << endl; // print out

  int n = 100;
  double x[n], y[n]; // acquire memory, 2 double precision arrays

  for (int i = 0; i < n; ++i)
  {
    x[i] = i;
    y[i] = i * 2;
  }

  double xval = 2.5;
  double yval = lookupVal(n, x, y, xval);

  cout << "For xval = " << xval << " yval = " << yval << endl;
}