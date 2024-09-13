//  ================================================================================
//  ||                                                                            ||
//  ||              fd                                                            ||
//  ||              -------------------------------------------                   ||
//  ||              F I N I T E   D I F F E R E N C E                             ||
//  ||                                                                            ||
//  ||              D E M O N S T R A T I O N   C O D E                           ||
//  ||              -------------------------------------------                   ||
//  ||                                                                            ||
//  ||       Developed by: Scott R. Runnels, Ph.D.                                ||
//  ||                     University of Colorado Boulder                         ||
//  ||                                                                            ||
//  ||                For: CU Boulder CSCI 4576/5576 and associated labs          ||
//  ||                                                                            ||
//  ||           Copyright 2020 Scott Runnels                                     ||
//  ||                                                                            ||
//  ||                     Not for distribution or use outside of the             ||
//  ||                     this course.                                           ||
//  ||                                                                            ||
//  ================================================================================

  //  ==
  //  ||
  //  ||  Perform Jacobiiterations on nField X nField linear
  //  ||  system:  A*Solution = RHS
  //  ||
  //  ==

void SolveLinearSystem(int max_iter , VD RHS, VD &Solution , mpiInfo &myMPI )
  {
    int converged, it_converged;
    int iter = 0;
    double newval;
    double cur_delta = 0.;
    double max_delta = 0.;
    double tol;

    MPI_Status status;
    int        tag = 0;
    int        err;
    int        global_converged;
    int        zero = 0;
    int        one  = 1;
    

    VD SolutionNew; SolutionNew.resize(Solution.size());

    rLOOP Solution[r] = 0.;
    tol       = 1.e-04;
    
    // ====================================================
    // Begin Iterations
    // ====================================================
    
    converged = 0;

    
    while ( converged == 0 )
      {
	// ----------------------------------------------
	// (1) Parallel communication on PE Boundaries
	// ----------------------------------------------
	
	// TO-DO in Lab:  Activate --> myMPI.ExchangeBoundaryInfo(  /* TO-DO in Lab */  );

	// ----------------------------------------------
	// (2) Convergence
	// ----------------------------------------------
	
	if ( ++iter > max_iter )
	  {
	    cout << "Jacobi: max_iter of " << max_iter << " reached.\n";
	    return;
	  }

	max_delta    = 0.;  it_converged = 1;

	// ----------------------------------------------
	// (3) One Jacobi Iteration
	// ----------------------------------------------
	
	rLOOP
	  {

	    
	    // (3.1) Compute new guess for row r

	    newval = b[r];
	    cLOOP  if ( r != c ) newval -=  A[r][c] * Solution[c];
	    newval /= A[r][r];

	    // (3.2) Convergence check

	    cur_delta  = fabs(Solution[r] - newval);

	    if ( cur_delta > tol ) it_converged = 0;

	    // (3.3) Record new value in solution

	    Solution[r] = newval;

	  }

	// ----------------------------------------------
	// (4) Make note of the convergence state
	// ----------------------------------------------
	
	converged = it_converged;

	// ----------------------------------------------
	// (5) Gather convergence information from PEs
	// ----------------------------------------------

	int root = 0;
	err = MPI_Reduce( &converged          , &global_converged, one , MPI_INT , MPI_MIN , zero , MPI_COMM_WORLD );
	err = MPI_Bcast ( /* TO-DO in Lab */  ,                    one , MPI_INT,            zero , MPI_COMM_WORLD );
	converged = global_converged;

	if ( converged == 1 )
	  {
	    if ( myPE == 0 ) cout << "Jacobi/GS converged in " << iter << " iterations" << endl;
	    return;
	  }
	
      }

  }

