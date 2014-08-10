/*=============================================================================
File: particle_mesh_subdivide.h
Description: Routines to compute the particle's equations
By: Carolina 2000/01
Last Update: 03/May/02
=============================================================================*/

#ifndef PARTICLE_DOMAIN_ELEMENT_H
#define PARTICLE_DOMAIN_ELEMENT_H

#include "particle_mesh.h"
#include "particle.h"
#include "mesh.h"

void FillInOutSurfConTab(
  const mesh_t* mesh,
  const particle_t* particle,
  const int*   partESp,
        double**  geoTSp,
        int*** surfN,
        int*   surfNIS );

void ParConnGeo(
  const mesh_t* mesh,
  const particle_t* particle,
      double** GeoTab,
      int** conTSp,
      int** fullESp,
      int* ElemPartly,
        int*** surfN,
        int*   surfNIS );

#endif
