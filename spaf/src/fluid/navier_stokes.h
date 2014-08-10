#ifndef FLUID_SOLVE_H
#define FLUID_SOLVE_H

#include "fluid.h"
#include "mesh.h"
#include "bc.h"
#include "particle.h"

void
FluidNavierStokes(
const int     order,
const mesh_t* mesh,
const bc_t*   bc,
#ifdef VERSION_Z
      mesh_t*     mesh_o,
      double*     Vel_o[4],
#endif          
const fluid_t*    fluid,
const particle_t  particle[],
const double      tau[3],
const double      dt,
      FluidVar_t* FluidVar );

#endif
