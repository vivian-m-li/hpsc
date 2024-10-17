#include "fd.h"

// gravitational field pointing down the in the mesh pulling all the particles
// down

// Lab notes:
// have a batch of particles
// collect in an array
// send to process red
// call the add function so process red can add them

class particles {
 public:
  int n;  // Maximum number of particles on this PE
  // not growing and shrinking particles on a process - create a storage
  // location static in memory, particles can move in and out of it
  VD x;       // The x-locations of the particles on this PE
  VD y;       // The y-locations of the particles on this PE
  VD vx;      // The x-velocities of the particles on this PE
  VD vy;      // The y-velocities of the particles on this PE
  VD xf;      // The forces in the x-direction on each particle
  VD yf;      // The forces in the y-direction on each particle
  VI active;  // Indicate if the particle is active on this PE
              // active[i] = 0 means particle i is not active.
  double mass;
  // memory-inefficient to be shrinking/growing particle array every time it
  // changes

  // ---------------------------------------------------------------------------
  // Constructor. The main routine will instantiate the particles class,
  // telling it how many particles, maximum, to accommodate
  // ---------------------------------------------------------------------------

  particles(int _n) {
    n = _n;           // maximum number of particles on the process
    x.resize(n + 1);  // 1-based indexing of particles.
    y.resize(n + 1);
    vx.resize(n + 1);
    vy.resize(n + 1);
    xf.resize(n + 1);
    yf.resize(n + 1);
    active.resize(n + 1);

    for (int i = 1; i <= active.size(); ++i) {
      active[i] = 0.;  // Initially, no active particles.
    }
  }

  // ---------------------------------------------------------------------------
  // Add particles to this PE
  // ---------------------------------------------------------------------------

  void add(VD &xAdd, VD &yAdd, VD &vxAdd, VD &vyAdd) {
    // Add new particles in spaces that are currently inactive
    // Search "active" array for zeros. When we find a zero,
    // add a particle in that location.
    // particle ID doesn't matter, just x-y location of the particle
    int count = 0;
    for (int j = 1; j < active.size(); ++j)
      if (active[j] == 0)  // Let's add a particle in the j-location
      {
        ++count;
        if (count > xAdd.size() - 1)
          return;  // We have added all the new particles

        x[j] = xAdd[count];
        y[j] = yAdd[count];
        vx[j] = vxAdd[count];
        vy[j] = vyAdd[count];
        active[j] = 1;
      }
  }

  // If we run out of space in x[], y[], active[], etc., we need to resize them.
  // FatalError("Out of room.");

  // -----------------------------------------------------------------------------
  // Move the particles on this PE
  // -----------------------------------------------------------------------------
  void move(double dt) {
    // dt is the time-step in seconds. It represents a small period of time
    // over which the forces are applied to the particles. Time marching

    // Apply the force to each particle to update its acceleration.
    // F = ma                   - Newton's Second Law
    // a = F/m                  - Algebra
    // a = change in velocity over time = (vnew - vold)/dt
    // vnew = vold + a*dt       - Algebra

    for (int i = 1; i < x.size(); ++i) {
      if (active[i] == 1)  // only update active particles
      {
        // vx and vy are completely independent of each other
        vx[i] += dt * xf[i] / mass;  // vnew = vold + a*dt
        vy[i] += dt * yf[i] / mass;

        // Update the locations
        x[i] += dt * vx[i];
        y[i] += dt * vy[i];
      }
    }
  }

  // Adding particles on the far-left boundary
  // Need to consider flux density, coulombs/second
  // periodically injects particles
  void addFlux(double time, double y0, double y1, double density, double vx_bdy,
               double vy_bdy) {}
}
