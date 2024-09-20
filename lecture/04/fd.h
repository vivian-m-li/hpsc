#include <mpi.h>

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "math.h"

using std ::cout;
using std ::endl;
using std ::string;
using std ::stringstream;
using std ::vector;

typedef vector<double> VD;            // 1D array
typedef vector<vector<double> > VDD;  // 2D array
typedef vector<int> VI;               // vector of integers

// Square matrix:
#define rLOOP for (int r = 1; r <= nField; ++r)
#define cLOOP for (int c = 1; c <= nField; ++c)
#define iLOOP for (int i = 1; i <= nRealx; ++i)
#define jLOOP for (int j = 1; j <= nRealy; ++j)
#define sLOOP for (int s = 0; s <= nRealx + 1; ++s)
#define tLOOP for (int t = 0; t <= nRealy + 1; ++t)