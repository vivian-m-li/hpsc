LaplacianOnGrid.h

- TO-DO in class: convert the following to sparse matrix format
  - diagonal = pid(i,j)
  - east = pid(i-1,j)
  - west
  - north
  - south
- for (int i = 1; i <= 4; ++i) Qval[p[i]] += w[i] \* PTCL.Qp;

  - this Qval is accumulated
  - taking each particle and dividing its weight up and giving it to the 4 surrounding nodes
  - most important line of code
  - Qval will not be complete on the processor boundary because that cell is only a half cell
    - particles on other processes contributing and sum them on the boundary
    - point on the intersection of 4 processes, sum all 4 processes together and have all 4 processes get the right answer

- looping over all point/particles
  - find out if active
    - left/right/top/bottom
      iL = int()... - takes x location, subtracts left boundary, divides by dx and adds i --> gives us the point to the left
    - assign array
    - computing weights
    - assigning particle, summing to to Qval entry for that point (line above)
