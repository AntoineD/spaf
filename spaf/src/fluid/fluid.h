#ifndef FLUID_H
#define FLUID_H

#include "includes.h"
#include "operators.h"

/// Fluid physical parameters
typedef struct
{
  double
    Re,  ///< Reynolds number
    Fr,  ///< Froude number
    Gravity[4]; ///< gravity vector
    
  bool
    FrameVelDir[4]; ///< if true then the frame moves in this direction
}
fluid_t;

/// Fluid variables structure
typedef struct
{
  int
    NbOfPreNodes, ///< number of pressure nodes
    NbOfVelNodes; ///< number of velocity nodes
  
  double
    *Pre,        ///< pressure
    *VelRHS[4],  ///< velocity rhs
    *Vel[4],     ///< velocity at current time step
    *VelOld1[4], ///< velocity at time step n
    *VelOld2[4], ///< velocity at time step n-1
    *Acc[4];     ///< velocity acceleration
}
FluidVar_t;

FluidVar_t*
FluidVarAlloc( void );

FluidVar_t*
FluidVarCreate(
const int PreNodesNb,
const int VelNodesNb );

fluid_t*
ReadFluidParameters(
const char* FileName ); 

void
FluidShiftVar(
FluidVar_t* FluidVar );

#endif
