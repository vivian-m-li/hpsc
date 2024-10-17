// #include "main.h" // copied from class repo

#include <mpi.h>
#include <iostream>
#include <string>
using std ::cout; // print
using std ::endl; // line feed
using std ::string;

void printArray(string arrayName, int n, int *array, int myPE)
{
  for (int i = 0; i < n; ++i)
    cout << "myPE: " << myPE << " " << arrayName << "[" << i << "] = " << array[i] << endl;
}

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv); // & = pass the locations in memory, not the value
  int numPE = 4;          // this needs to match the number of PEs on which the code is running
  int myPE;               // IDK? not declaring this

  MPI_Comm_size(MPI_COMM_WORLD, &numPE); // PE stands for processes
  MPI_Comm_rank(MPI_COMM_WORLD, &myPE);  // rank = process

  cout << "myPE: " << myPE << " Hello world!" << endl;

  // Broadcast
  int n = 10;
  int a[n]; // initialize array 10 elements long

  for (int i = 0; i < 10; ++i)
    a[i] = 0;

  if (myPE == 0)
    for (int i = 0; i < 10; ++i)
      a[i] = 333;

  MPI_Bcast(a, n, MPI_INT, 0, MPI_COMM_WORLD); // Broadcast a to all PE, from PE 0

  printArray("a", n, a, myPE);

  // Scatter
  int m = 40;
  int b[m]; // b is an array that's length 40, and when we scattered it, we scatter it in chunks of 10 to each PE in a

  if (myPE == 0)
    for (int i = 0; i < 40; ++i)
      b[i] = i;

  MPI_Scatter(b, n, MPI_INT,
              a, n, MPI_INT, 0, MPI_COMM_WORLD); // Scatter b from process 0 to everyone

  printArray("a", n, a, myPE);

  // Reduction
  int myVal = myPE;
  int myReducedVal = -1;

  //        from each PE, the sum, how long, type, operation, receive vals here
  MPI_Reduce(&myVal, &myReducedVal, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); // sum across all processes

  // sum = 6 but only process 0 will get the answer
  cout << "myPE: " << myPE << " sum of the PE numbers = " << myReducedVal << endl;

  MPI_Allreduce(&myVal, &myReducedVal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); // sum across all processes
  cout << "myPE: " << myPE << " sum of the PE numbers = " << myReducedVal << endl;

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
}