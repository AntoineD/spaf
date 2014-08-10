#ifndef PARTICLE_DOMAIN_MESH_H
#define PARTICLE_DOMAIN_MESH_H

#include "mesh.h"
#include "particle.h"

void
GetParticleMesh(
const mesh_t*     mesh,
      particle_t  particle[] );

bool
GetParElementJacobian(
const mesh_t*     mesh,
const Parmesh_t*  Parmesh,
const int      element,
      double*       J_times_wt );

#endif
