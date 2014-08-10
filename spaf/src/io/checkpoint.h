#ifndef CHECKPOINT_H
#define CHECKPOINT_H

#include "particle.h"
#include "fluid.h"

void
WriteCheckPoint(
                const bool        WriteInBinary,  ///< write in binary or ascii
                const int         TimeStep,       ///< timestep
                const FluidVar_t* FluidVar,       ///< fluid variables
                const particle_t* particles );     ///< particles variables

bool
ReadCheckPoint(
               int*        TimeStep,  
               FluidVar_t* FluidVar,  
               particle_t* particles );

#endif
