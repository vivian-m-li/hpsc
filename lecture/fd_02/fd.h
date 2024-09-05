#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>
#include "math.h"

using std ::cout;
using std ::endl;
using std ::string;
using std ::stringstream;
using std ::vector;

typedef vector<double> VD;          // 1D array
typedef vector<vector<double> > VDD; // 2D array

// Square matrix:
#define rLOOP for (int r = 1; r <= nField; ++r) // row loop
#define cLOOP for (int c = 1; c <= nField; ++c) // column loop

// For the mesh
#define iLOOP for (int i = 1; i <= nRealx; ++i) // Loop over the i's in the mesh
#define jLOOP for (int j = 1; j <= nRealy; ++j) //               j's
#define kLOOP for (int k = 1; k <= nRealz; ++k) //               k's