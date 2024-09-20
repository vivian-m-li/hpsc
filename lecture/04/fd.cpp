#include "fd.h"

#include "particles.h"
// #include "mpiInfo.h"

class Mesh  // Mesh = LaplacianOnGrid - Solving the Laplace's Equation +
            // Handling Particles
{           // No FormLS... No SolveLinearSystem... No LinearSolver
 public:
  double x0, x1, y0, y1;  // x-y physical extent of the mesh
  VD x, y;                // Locations of each cell center in the mesh
  int nRealx, nRealy, nField, nReal;  // Number of real cells in the mesh
  double lengthx, lengthy, dx,
      dy;  // dx/dy size of each cell in the x-y direction
  int myPE;

  // Constructor: Initialize values

  Mesh(double _x0, double _x1, double _y0, double _y, int ncell_x,
       int ncell_y) {
    // ...
  }

  // ParticlesOnMesh
  // Quantify the particles' interaction with the mesh
  // References: [1] https://www.particleincell.come/2010/es-pic-method

  // hand it a particles object
  // routine will get info on all the particles on this mesh
  void ParticlesOnMesh(particles &PTCL, mpiInfo &myMPI) {
    // Determine which particles are still on this mesh and which have left
    VI ptcl_send_list,
        ptcl_send_PE;  // These collect information about particles that have
                       // left this PE and are heading onto another PE

    for (int k = 1; k <= PTCL.n;
         ++k) {  // Loop over all spaces in the particle data structure
      if (PTCL.active[k] == 1) {  // Only consider the active ones

        // Test to see if particle k is still on this PE's mesh
        // active[k] = -1 means particle k is leaving this PE and needs to be
        // moved to another PE
        if (PTCL.x[k] < x0) {  // Leaving the left boundary
          PTCL.active[k] = -1;
          iPEnew = myMPI.iPE - 1;  // Where is it heading... to which PE?
        }

        // handle the other 3 sides in lab
        // if (PTCL.x[k] > x1) {} // Right boundary
      }

      // This particle, k, is not in our mesh anymore. Collect this particle
      // into a holding array that will be sent to the neighboring processor.

      if (PTCL.active[k] == -1) {
        // Let's make sure that k hasn't left the entire problem...
        if (iPEnew >= 0 && iPEnew < myMPI.nPEx)
          if (jPEnew >= 0 && jPEnew < myMPI.nPEy) {
            ptcl_send_list.push_back(k);  // Add k to this array.
            ptcl_send_PE.push_back(
                /* in lab */);  // This contains the PE number to which particle
                                // k will go.
                                // need to do an MPI send and receive over all
                                // of this information
          }
      }

      // Give and receive particles to/with other processors
      myMPI.ParticleExchange(ptcl_send_list, ptcl_send_PE, PTCL);

      // At this point, there are no more active[k]  = -1. They are all 0 or 1
      // (inactive or active). All particles are now on the correct PE.

      // Next week: Accumulate particles to the nodes (to be completed in the
      // esPIC code)

      // Compute forces on particles: Eventually, this class will solve the
      // Poisson equation from which we will determine the forces
      for (int k; k <= PTCL.n; ++k) PTCL.xf[k] = PTCL.yf[k] = 0.;

      // comment this loop out for debugging (so forces are only going sideways)
      for (int k; k <= PTCL.n; ++k) {
        if (PTCL.active[k] == 1) {
          PTCL.xf[k] = 0.;  // Gravity does not affect horizontal motion.
          PTCL.yf[k] = -4;  // Think of gravitation force pointing down.
        }
      }
    }
  }

  // Set up the particles
  particles PTCL(500);  // can change number of particles here

  // Time Marching Loop
  for (double t = 0.; t <= tEnd; t += dt) {  // tEnd is a user input
    if (myMPI.iPE == 0)
      PTCL.addFlux(/* ..... */);  // Particle injection on the far left
    PTCL.move(dt);                // Move the particles
    MESH.ParticlesOnMesh(PTCL, myMPI);
    PTCL.plot(...);
    MESH.plot(...);
  }