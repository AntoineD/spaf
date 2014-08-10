#ifndef PARTICLE_IO_H
#define PARTICLE_IO_H

#include "particle.h"
#include "mesh.h"

void
WriteParticleResults(
const bool        FullPrecision,
const double        time,
const int      TimeStep,
const particle_t  particle[] );

void
WriteParMeshVTK(
const int      iPar,
const mesh_t*     mesh,
const Parmesh_t*  ParMesh );

void
WriteParticleCheckpoint(
const bool        WriteInBinary,
const particle_t* particles );

void
ReadParticleCheckpoint(
particle_t* particles );

#endif
