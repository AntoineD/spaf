#ifndef PARTICLE_COLLISION_H
#define PARTICLE_COLLISION_H

#include "includes.h"
#include "particle.h"
#include "fluid.h"

void
ReadCollisionParameters(
FILE* FIleId );

void
GetCollisionPlanes(
const char* FileName );

void
DetectCollision(
particle_t  particle[] );

void
ParColForcePred(
const double        dt,
const fluid_t*    fluid,
      particle_t  particle[] );

void
ParColForceCor(
const fluid_t*    fluid,
      particle_t  particle[] );

#endif
