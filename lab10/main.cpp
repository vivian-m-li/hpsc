
#include<mpi.h>  
#include <cstdlib>

// ==
// ||
// ||  Memory Allocation
// ||
// ==

float  * * Array2D_float(int nRows,int nCols)
{
  float *myArray;
  myArray = new float [ nRows * nCols ];
    
  // Create a pointer that points to the beginning of each new row

  float * * myArray_ptr;
  myArray_ptr = new float * [nRows];

  int count = 0;

  for ( int row = 0 ; row < nRows ; ++ row )
    {
      myArray_ptr[row] = &myArray[ count*nCols ];
      ++count;
    }

  // Return that pointer
  
  return myArray_ptr;
}




//  -
//  |
//  |   O V E R V I E W   O F   T H E
//  |
//  |   G R I D   A N D   P A R A L L E L   S T R U C T U R E
//  |
//  |   
//  |   In this code, the real nodes overlap at the PE boundary.
//  |   So, when we write values to the file we do not always want to write
//  |   all of the real nodes' values.  The global grid is shown below, as
//  |   it is split on each of four (2PE x 2PE) processes.  In the figure
//  |   below, the global i-j indices are shown.
//  |
//  |	   [1,5]---[2,5]---[3,5]   [3,5]---[4,5]---[5,5]
//  |	     |       |       |       |       |       |  
//  |	     |       |       |       |       |       |  
//  |	     |       |       |       |       |       |  
//  |	   [1,4]---[2,4]---[3,4]   [3,4]---[4,4]---[5,4]
//  |	     |       |       |       |       |       |  
//  |	     |       |       |       |       |       |  
//  |	     |       |       |       |       |       |  
//  |      [1,3]---[2,3]---[3,3]   [3,3]---[4,3]---[5,3]
//  |
//  |
//  |      [1,3]---[2,3]---[3,3]   [3,3]---[4,3]---[5,3]
//  |	     |       |       |       |       |       |  
//  |	     |       |       |       |       |       |  
//  |	     |       |       |       |       |       |  
//  |      [1,2]---[2,2]---[3,2]   [3,2]---[4,2]---[5,2]
//  |	     |       |       |       |       |       |  
//  |	     |       |       |       |       |       |  
//  |	     |       |       |       |       |       |  
//  |	   [1,1]---[2,1]---[3,1]   [3,1]---[4,1]---[5,1]
//  |
//  |
//  |   On each PE, the local numbers for the same grid are identical,
//  |   i.e., they look like this:
//  |
//  |
//  |      [1,3]---[2,3]---[3,3]   [1,3]---[2,3]---[3,3]   
//  |	     |       |       |       |       |       |     
//  |	     |       |       |       |       |       |     
//  |	     |       |       |       |       |       |     
//  |      [1,2]---[2,2]---[3,2]   [1,2]---[2,2]---[3,2]   
//  |	     |       |       |       |       |       |     
//  |	     |       |       |       |       |       |     
//  |	     |       |       |       |       |       |     
//  |	   [1,1]---[2,1]---[3,1]   [1,1]---[2,1]---[3,1]   
//  |
//  |
//  |      [1,3]---[2,3]---[3,3]   [1,3]---[2,3]---[3,3]   
//  |	     |       |       |       |       |       |     
//  |	     |       |       |       |       |       |     
//  |	     |       |       |       |       |       |     
//  |      [1,2]---[2,2]---[3,2]   [1,2]---[2,2]---[3,2]   
//  |	     |       |       |       |       |       |     
//  |	     |       |       |       |       |       |     
//  |	     |       |       |       |       |       |     
//  |	   [1,1]---[2,1]---[3,1]   [1,1]---[2,1]---[3,1]   
//  |
//  |
//  |   The field values at nodes on the PE boundaries exist (and are identical)
//  |   on all PEs sharing that boundary.  In other words, using local
//  |   i-j numbers, the field value A[3,2] on PE 0 (lower left) is
//  |   equal to A[1,2] on PE 1.
//  |
//  |   When writing a single dump file for the field variable, A,            
//  |   we do not want both A[3,2] on PE0 and A[1,2] on PE1 to be written.
//  |   We must adjust the MPI sub-arrays accordingly.
//  |
//  -

// ==
// ||
// ||  Main
// ||
// ==

int main(  int argc, char *argv[] )
{

  int nRealx  = 3;    // Number of real points in the x- and y-directions
  int nRealy  = 3;    //
  int myPE, numPE;    // This process' ID and the number of processes total
  int nPEx,  nPEy;    // Number of processes in the x- and y-directions, in the 2D grid of processes
  int iPE ,   jPE;    // This process' logical ID in the 2D array of processors

  // (1) Initialize MPI

  MPI_Init      ( &argc          , &argv  );
  MPI_Comm_size ( MPI_COMM_WORLD , &numPE );
  MPI_Comm_rank ( MPI_COMM_WORLD , &myPE  );

  // (2) In this example, we are not letting the user choose the number of PEs.  It is a 2x2 logical
  //     grid of PEs.  The formulas for iPE and jPE are from the previous mpiInfo.h

  nPEx = 2;
  nPEy = 2;
  
  jPE = int(myPE/nPEx);
  iPE = myPE - jPE*nPEx;

  // (4) Allocate memory for and populate the real nodes' entries in the 2D array we wish to write
  
  float **A = Array2D_float(nRealx+1, nRealy+1);
  for (int j = 1; j <= nRealy; j++)
    for (int i = 1; i <= nRealx; i++) 
      A[i][j] = 10*i + j  + myPE/10.;
  
  // (5) Create MPI derived data type that will be used to represent the real nodes in the grid.

  // (5.1) Set up specifications for the datatype, which we will call "myRealNodes"
  
  MPI_Datatype myRealNodes;
  int idxStartThisPE  [2] = { 1        , 1        };  // Index coordinates of the sub-array inside this PE's array, A
  int AsizeThisPE     [2] = { /* TO-DO in Lab */  };  // Size of the A array on this PE    
  int sub_AsizeThisPE [2] = { /* TO-DO in Lab */  };  // Size of the A-sub-array on this PE

  // Adjust sub_AsizeThisPE and idxStartThisPE to avoid writing values on PE boundaries twice

  /* TO-DO in Lab */
  /* TO-DO in Lab */

  // (5.2) Ask MPI to create that new datatype
  //
  // ------------------------------------------------------------------------------------------------------
  //
  // Number of elements in A in each dimension                   Starting indices of the subarray inside
  //                          of the sub array                   the this PE's A array
  //                                         |                   |
  //  Number of elements in each             |                   |           Order of the A array in memory, row
  // dimension of this PE's array            |                   |           major or column major
  //                            |            |                   |           |
  //    Number of dimensions    |            |                   |           |
  //                       |    |            |                   |           |     
  //                       |    |            |                   |           |     
  MPI_Type_create_subarray(2, AsizeThisPE, sub_AsizeThisPE, idxStartThisPE, MPI_ORDER_FORTRAN, MPI_FLOAT, &myRealNodes);
  //
  // ------------------------------------------------------------------------------------------------------
  
  MPI_Type_commit(&myRealNodes);
  
  // (6) Create a second derived type that will be used to describe the whole array, as it exists
  //     across PEs.  This next datatype describes this PE's relationship between its part of the
  //     array, A, and the global A array across PE.  It is used to create the "view" of the file
  //     from this PE's perspective, i.e., where in that file this PE should write its data.
  
  MPI_Datatype myPartOfGlobal;
  int idxStartInGlobal [2] = { /* TO-DO in Lab */    };  // Index cordinates of the sub-arrayinside the global array
  int AsizeGlobal      [2] = { /* TO-DO in Lab */    };  // Size of the global array

  // Adjust idxStartInGlbal to ensure proper values are written.

  /* TO-DO in Lab */
  /* TO-DO in Lab */
  
  // ------------------------------------------------------------------------------------------------------
  //
  //      Number of elements in each dimension              Starting indices of the subarray inside
  //                          of the sub array              the this PE's array
  //                                         |              |
  //  Number of elements in each             |              |                    Order of the array in memory, row
  //  dimension of this PE's array           |              |                    major or column major
  //                            |            |              |                    |
  //    Number of dimensions    |            |              |                    |
  //                       |    |            |              |                    |     
  //                       |    |            |              |                    |     
  MPI_Type_create_subarray(2, AsizeGlobal, sub_AsizeThisPE, idxStartInGlobal, MPI_ORDER_FORTRAN, MPI_FLOAT, &myPartOfGlobal);
  //
  // ------------------------------------------------------------------------------------------------------
  
  MPI_Type_commit(&myPartOfGlobal);
  
  // (7) Set the "view" of the file from this PE's perspective, i.e., where and how this PE should write data
  
  MPI_File fh;
  
  MPI_File_open(MPI_COMM_WORLD, "output.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
  
  // ------------------------------------------------------------------------------------------------------
  //
  //
  //     Displacement in file to     Contains the description
  //     accommodate other data      of how this PE's "A" array
  //     written before this data    fits into the global "A".
  //                     |                 |
  //                     |                 |
  MPI_File_set_view (fh, 0, MPI_FLOAT, myPartOfGlobal, "native", MPI_INFO_NULL);
  //                           |                                     |
  //                           |                                     |
  //                       The data type of each element     Could contain optimization hints,
  //                                                         e.g., striping information
  // ------------------------------------------------------------------------------------------------------

  //  (8) Perform the collective write operation

  
  // ------------------------------------------------------------------------------------------------------
  //                            The sub_array describing this
  //        File pointer        PE's real-nodes part of A
  //                 |                   |
  //                 |                   |
  MPI_File_write_all(fh, &A[0][0], 1, myRealNodes, MPI_STATUS_IGNORE);
  //                       |
  //                       |
  //                 Memory location of original
  //                 data set for which the
  //                 sub arrays were created.
  // ------------------------------------------------------------------------------------------------------

  MPI_File_close(&fh);
  
  // (9) Cleanup 
  
  free(A[0]);
  free(A);
  
  MPI_Type_free(&myPartOfGlobal);
  MPI_Type_free(&myRealNodes);
  MPI_Finalize();
  
  return 0;

}




