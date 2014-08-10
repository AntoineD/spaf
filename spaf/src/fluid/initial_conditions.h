#ifndef FLUID_INITIAL_CONDITION_H
#define FLUID_INITIAL_CONDITION_H

#include "includes.h"
#include "mesh.h"
#include "bc.h"

void
ReadFluidInitialConditions(
const char*   FileName,
const mesh_t* mesh,
      double**  Vel,
      double*   Pre );

void
GetFluidInitialConditions_z(
const double position[4],
const mesh_t*     mesh,
const mesh_t*     mesh_o,
const FluidVar_t* FluidVar_o,
      FluidVar_t* FluidVar );

#endif
