#include "omp.h"

// ==
// ||
// ||  OMP Laplacian Demo Code
// ||
// ==


int main(  int argc, char *argv[] ) 
{

  // Start the timer for the entire execution

  struct timespec t0, t1;
  StartTimer(t0);


  // Default values for user input
    
  int nCell = 190;  // nCell is for the x- and y-directions
  int numTH = 1;   // Number of threads

  // Gather user input from the command line

  for (int count =  0 ; count < argc; ++count)
    {
      if ( !strcmp(argv[count],"-nCell"    ) ) nCell     = atoi(argv[count+1]);  
      if ( !strcmp(argv[count],"-nTH"      ) ) numTH     = atoi(argv[count+1]);
    }

  // Set the number of threads for this execution

  omp_set_num_threads(numTH);

  // Set up grid and linear system information

  double dh          = 1./nCell;        // Size of each cell in each direction (dx and dy)
  int nPtsx          = nCell+1;         // Number of points in the x-direction
  int nPtsy          = nCell+1;         //   "    "    "    "   "  y-    "
  int nField         = nPtsx * nPtsy;   // Total number of points (and number of rows in the linear system)
  int bandwidth      = 5;               // Number of columns in the sparse-matrix format
  int max_iter       = nField*10;       // Number of allowed Jacobi iterations
  int numTHconverged = 0;               // Number of threads converged
    
  double   * phi;                       // The unknown for which we are solving A*phi = b
  
  phi     = Array1D_double(nField);     // Acquire memory for this shared variable.

  // This shared array will store either a 0 or 1 for each thread, 0 meaning not congverged
  // 1 meaning converged.
  
  int *THconverged = new int [numTH];   // This array records which threads have converged. done
  

// =================================================================================== //
// =======================         BEGIN PARALLEL REGION          ==================== //
// =================================================================================== //
  
#pragma omp parallel shared(numTHconverged)
  {
    int myTH = omp_get_thread_num();  

    // ------------------------------------------------------------------------------------
    // (1) Compute this thread's bounds in the global system:
    //
    //     numPerTH = number of rows per thread (the last thread may be slightly different)
    //     Lower    = the first row in the global matrix to be handled by this thread
    //     Upper    = the last   "  "   "    "      "    "  "    "     "   "      "
    // ------------------------------------------------------------------------------------

    int numPerTH = nField/numTH + 1; // assign threads, +1 to handle remainders in case it's not equally divisible, done
    int Lower    = numPerTH * myTH; // done
    int Upper    = Lower + numPerTH; // done

    //  (1.1) Adjust Upper for the last (highest numbered) thread to ensure all the rows of the global system are handled, done
    if (myTH == numPerTH - 1) { Upper = nField; }

    // ------------------------------------------------------------------------------------
    // (2) Acquire memory for this thread
    //
    //     nField_TH = the number of field variables (phi) to be handled by this thread
    // ------------------------------------------------------------------------------------

    int nField_TH = Upper - Lower; // done

    double * * Acoef;
    int    * * Jcoef;
    double   * b;
    double   * phiNew;

    //     Acoef, Jcoef, b, and phiNew are not the global matrices and arrays; here they contain only the
    //     rows being handled by this thread.
    
    Acoef    = Array2D_double(nField_TH, bandwidth); // only a local matrix, so we use nField_TH instead of nField, done
    Jcoef    = Array2D_int   (nField_TH, bandwidth); // done
    b        = Array1D_double(nField_TH); // done
    phiNew   = Array1D_double(nField_TH); // done
    
    // ------------------------------------------------------------------------------------
    // (3) Initialize the linear system and form the initial guess for phi
    // ------------------------------------------------------------------------------------

    // Each thread loops only over its local matrix.  Use the correct upper limit below.
    for ( int row = 0 ; row < nField_TH  ; ++row ) // loop through only the local matrix, done   
      {
	for ( int col = 0 ; col < bandwidth ; ++col )
	  {
	    Acoef[row][col] = 0.;
	    Jcoef[row][col] = 0 ;
	  }
	b[row] = 0.;
      }

    // Each thread is responsible for initializing only its part of the global solution array.
    
    for ( int row = Lower ; row <= Upper ; ++row )  // loop from lower to upper, done
      {
	phi[row] = 0.;
      }
    
    // ------------------------------------------------------------------------------------
    // (4) Form the linear system.  Here "pt" represents "point" number in the mesh.  It is
    //     equal to the row number in the linear system, too.
    // ------------------------------------------------------------------------------------
    
    for ( int pt = Lower ; pt <= Upper ; ++pt ) // Each thread loops only over its rows in the global matrix (row = pt)
      {

	// Using the same logic as for converting myPE to iPE and jPE,
	// compute the i,j logical coordinates of "pt" in the mesh:

	/* Compute j */
  int j = pt / nPtsx;  // describes which row we're on, done

	/* Compute i */
  int i = pt % nPtsx; // describes which column we're on, done
	
	// Compute the row number local to this thread, i.e., the row in its Acoef/Jcoef arrays
	
	int pt_TH = pt - Lower;   

	// Populate the linear system for all interior points using that local row number

	if ( i > 0 && i < nPtsx-1 )
	  if ( j > 0 && j < nPtsy-1 )
	    {
	      Acoef[pt_TH /* done */][0] = -4./dh/dh;  Jcoef[pt_TH /*  done */][0] = pt; // pt_TH is local to this thread          
	      Acoef[pt_TH /* done */][1] =  1./dh/dh;  Jcoef[pt_TH /* done */][1] = pt - 1;       
	      Acoef[pt_TH /* done */][2] =  1./dh/dh;  Jcoef[pt_TH /* done */][2] = pt + 1;       
	      Acoef[pt_TH /* done */][3] =  1./dh/dh;  Jcoef[pt_TH /* done */][3] = pt + nPtsx;   
	      Acoef[pt_TH /* done */][4] =  1./dh/dh;  Jcoef[pt_TH /* done */][4] = pt - nPtsx;   
	    }

	// Apply boundary conditions
	
	if ( i == 0       ) { Acoef[pt_TH /* done */][0] = 1. ; Jcoef[pt_TH /* done */][0] = pt; b[pt_TH /* done */] =  1.; }  // setting the column BCs, indexing is the same as above because we're modifying the local matrices
	if ( i == nPtsx-1 ) { Acoef[pt_TH /*  done  */][0] = 1. ; Jcoef[pt_TH /* done */][0] = pt; b[pt_TH /* done */] = -1.; }  
	if ( j == 0       ) { Acoef[pt_TH /* done */][0] = 1. ; Jcoef[pt_TH /* done */][0] = pt; b[pt_TH /* done */] = -1.; }  // setting the row BCs
	if ( j == nPtsy-1 ) { Acoef[pt_TH /* done */][0] = 1. ; Jcoef[pt_TH /* done */][0] = pt; b[pt_TH /* done */] =  1.; }  
      }

    // ------------------------------------------------------------------------------------
    // (5) Perform Jacobi iterations
    // ------------------------------------------------------------------------------------

    int iter       = 0;
    int thisThreadConverged;

    #pragma omp barrier

    // numTHconverged counts the number of threads converged.  We do not care which thread initializes it.
    
    #pragma omp single /* Add clause for efficiency */ // VIVIAN: check this later
    {
      numTHconverged = 0;
    }

    // Iterate until all of the threads have converged or we have exceeded the number of allowed iterations
  
    while ( numTHconverged < numTH && ++iter <= max_iter ) // check if all of our threads have converged, done
      {
	thisThreadConverged = 1;

	for ( int row = Lower ; row <= Upper ; ++row ) 
	  {
	    // Just like above, with pt_TH, compute row_TH
	    
	    int row_TH = row - Lower; // done
	    
	    phiNew[row_TH /* done */] = b[row_TH /* done */];   
	    
	    for ( int col = 1 ; col < bandwidth ; ++col ) phiNew[row_TH /* done */] -= Acoef[row_TH /*  done */][col] * phi[ Jcoef[row_TH /* done */][col] ];
	    
	    phiNew[row_TH /* done */] /= Acoef[row_TH /* done */][0];                                                                            
	    
	    if ( fabs(phiNew[row_TH /* done */] - phi[row] ) > 1.e-10 ) thisThreadConverged = 0;                                     

	  }

	// (5.1) Record in shared array if this thread converged or not
	
	THconverged[myTH /* done */ ] = thisThreadConverged;   

	// (5.2) Count the number of threads that have converged.  We do not care which thread does the counting.

#pragma omp single /* Add a clause for efficiency */ // VIVIAN 
	{
	  numTHconverged = 0;
	  for ( int i = 0 ; i < numTH ; ++i ) numTHconverged += THconverged[i /* done */ ];
	  if ( numTHconverged == numTH ) printf("Jacobi converged in %d iterations.",iter);
	}

	// (5.3) Update the shared/global array phi with this thread's values
	
	for ( int row = Lower ; row <= Upper ; ++row )  
	  {
	    int row_TH = row - Lower; /* done (same as above) */
            phi[row] = phiNew[row_TH  /* done */];  
	  }
      }

  }

// =================================================================================== //
// =======================           END PARALLEL REGION          ==================== //
// =================================================================================== //
  
  if ( numTHconverged != numTH ) printf("WARNING: Jacobi did not converge.\n");
  
  EndTimer("main", t0,t1);
  
  plot("phi",phi,nPtsx,nPtsy,dh);
  return 0;
}

