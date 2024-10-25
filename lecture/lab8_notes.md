10/24

- populate a field variable (Diag_sumPE)
  - matrix entries incomplete for half staggered cells because they only get the value from one side (PE ont he other side)
  - Diag_sumPE takes care of this
- change in (4) Begin Iterations:

  - rLOOP {
    rowSUM_sumPE[r] = b[r];
    for (int c = 2; c <= bandwidth ; ++c) RowSum_sumPE[r] -=... // not new (using sparse matrix), just using a new variable name

    myMPI.PEsum(RowSum_sumPE); // new, need to sum the boundaries
    }

  - rLOOP {
    newval = RowSum_sumPE[r] / Diag_sumPE[r]; // need to sum the diagonals
    cur_delta = fabs...
    }

- RHS not complete - not necessarily 0, use the method of manufactured solutions
  - "b is going to be doing stuff so we get a nice linear solution"
- double Dot(...)
  - dot product: r dot r
  - every single r for the global system
    - all the PEs need to know what the global r dot r is
    - but every PE needs to know its part in computing it
    - the danger is that PEs share nodes on the boundary and you're going to double count them if you're not careful
- matrix vector product: void MatVecProd(...)
  - follows PE sum/row sum business
  - but no PE has the whole matrix on the boundaries
