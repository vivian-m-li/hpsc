#include "main.h"

void passValue(double val) {
  cout << "passValue val =" << val << endl;
  val += 1;
}

void passAddress(double *val) {
  cout << "passAddress val =" << val << endl;
  cout << "passAddress will increment the value by one." << endl;
  // val += 1;   // increments the address by 1, not the value
  *val += 1;  // increments the value by 1
}

void passValueGetAddress(double &val) {
  // I know you passed the value but I'm going to get the address
  cout << "passValueGetAddress val =" << val << endl;
  cout << "passValueGetAddress will increment it by one." << endl;
  val += 1;  // increments the value by 1
}

void passArray(double *val, int n) {
  cout << "passRray: ";  // no linefeed
  for (int i = 0; i < n; i++)
    cout << val[i];  // when you receive a value like an address, you can put brackets on it and have it as an index
  cout << "\n";
}

double **create2DArray() {  // returns a pointer to a pointer
  double *myArray;
  myArray = new double[9];

  myArray[0] = 11;
  myArray[1] = 21;
  myArray[2] = 31;

  myArray[3] = 21;
  myArray[4] = 22;
  myArray[5] = 23;

  myArray[6] = 31;
  myArray[7] = 32;
  myArray[8] = 33;

  double **myArray_ptr;
  myArray_ptr = new double *[3];  // one for each row of our 2D array
  myArray_ptr[0] = &myArray[0];
  myArray_ptr[1] = &myArray[3];
  myArray_ptr[2] = &myArray[6];

  return myArray_ptr;
}

int main(int argc, char *argv[]) {
  double val;

  cout << "\n(1) ---------- Pass a value -----------\n\n";
  val = 3.;
  passValue(val);
  cout << "main, val = " << val << endl;

  cout << "\n(2) ---------- Pass an address -----------\n\n";
  val = 3.;
  passAddress(&val);
  cout << "main, val = " << val << endl;

  cout << "\n(3) ---------- Pass a value routine grabs address -----------\n\n";
  val = 3.;
  passValueGetAddress(val);
  cout << "main, val = " << val << endl;

  cout << "\n(4) ---------- Pass array is the same as pass address[0] -----------\n\n";
  double a[2] = {3., 333.};
  passAddress(a);
  cout << "main, a[0] = " << a[0] << " a[1] = " << a[1] << endl;

  cout << "\n(5) ---------- Intro to pointer variables -----------\n\n";
  double *aptr = &a[0];
  cout << "main, aptr (address of a[0]) = " << aptr << endl;
  cout << "main, aptr content of that address = " << *aptr << endl;

  cout << "\n(6) ---------- Incrementing a pointer -----------\n\n";
  ++aptr;
  passValue(*aptr);

  cout << "\n(7) ---------- Pointing to another pointer -----------\n\n";
  double **aptr_ptr;
  --aptr;
  aptr_ptr = &aptr;
  cout << "main, aptr_ptr (address of aptr) = " << aptr_ptr << endl;
  cout << "main, aptr_ptr (content of that address) = " << *aptr_ptr << endl;

  // Incrementing another pointer
  ++aptr_ptr;
  aptr_ptr = &aptr;
  cout << "main, aptr_ptr (address of aptr) = " << aptr_ptr << endl;
  cout << "main, aptr_ptr (content of that address) = " << *aptr_ptr
       << endl;  // 8 bytes after the original one since we incremented the pointer

  cout << "\n(8) ---------- Methods of array passing -----------\n\n";
  double b[5] = {0, 100, 200, 300, 400};
  passArray(&b[0], 4);  // only pass in 4 elements
  passArray(&b[1], 4);  // starts at element 1 and only prints 4 elements

  cout << "\n(9) ---------- Pointers to row beginnings -----------\n\n";
  double c[9] = {11, 21, 31, 21, 22, 23, 31, 32, 33};  // this is a 1-D array in memory
  double *cptr[3];
  cptr[0] = &c[0];  // cptr[0] stores the address of the first row of our 3x3 matrix
  cptr[1] = &c[3];  // cptr[1]    "    "     "     "  "  second "  "   "   "    "
  cptr[2] = &c[6];  // cptr[2]    "    "     "     "  "  third  "  "   "   "    "

  cout << "print as 3 1-D arrays: " << endl;
  passArray(cptr[0], 3);
  passArray(cptr[1], 3);
  passArray(cptr[2], 3);

  cout << "print as a single 2D array: " << endl;
  for (int row = 0; row < 3; ++row) {
    for (int col = 0; col < 3; ++col) {
      cout << cptr[row][col] << " ";
    }
    cout << endl;
  }

  cout << "\n(10) ---------- Ask routine to create a 3x3 array -----------\n\n";
  double **Array2D;
  Array2D = create2DArray();

  cout << "print as a single 2D array: " << endl;
  for (int row = 0; row < 3; ++row) {
    for (int col = 0; col < 3; ++col) {
      cout << Array2D[row][col] << " ";
    }
    cout << endl;
  }

  return 0;
}