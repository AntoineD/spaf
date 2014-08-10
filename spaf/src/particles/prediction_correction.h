#ifndef PARTICLE_PREDICTION_CORRECTION_H
#define PARTICLE_PREDICTION_CORRECTION_H

#include "particle.h"
#include "fluid.h"

void
PredictParticlePosition(
const int      order,        ///< order of time integration
const mesh_t*     mesh,
const double        dt,
const fluid_t*    fluid,
      particle_t  particle[] );

void
CorrectParticlePosition(
const int      order,        ///< order of time integration
const double        dt,
const double      tau[3],       ///< time discretization coefficients
const fluid_t*    fluid,
      particle_t  particle[] );

#endif
