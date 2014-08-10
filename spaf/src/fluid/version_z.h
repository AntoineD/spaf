#ifndef version_z_H
#define version_z_H

#include "includes.h"
#include "mesh.h"
#include "bc.h"
#include "operators.h"
#include "fluid.h"
#include "particle.h"

void
ParticleInitialConditions(
const mesh_t*     mesh_o,
const FluidVar_t* FluidVar_o,
      particle_t* particle );

#endif
