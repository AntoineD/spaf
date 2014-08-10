#ifndef PARTICLE_RBI_H
#define PARTICLE_RBI_H

#include "fluid.h"
#include "mesh.h"
#include "particle.h"

void
ReadRigidBodyMotionParameters(
FILE* FileId );

void RigidBody(
const mesh_t*       mesh,
const double        tau[3],
      FluidVar_t*   FluidVar,
      particle_t    particle[] );

void
SetRigidBodyMotion(
const mesh_t*       mesh,     
const particle_t*   particle, 
      FluidVar_t*   FluidVar );

#endif
