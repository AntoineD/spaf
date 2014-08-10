/*=============================================================================
File: particle_mesh_find_elem.h

Description: Routines to compute the particle's equations

By: Carolina 2000/01
Last Update: 19/Sep/06 Antoine
=============================================================================*/

#ifndef _PARTICLE_DOMAIN_FILL_H_
#define _PARTICLE_DOMAIN_FILL_H_

#include "mesh.h"
#include "particle.h"

bool
IsPointInPar(
const particle_t* particle,
const double        nodPos[ 4 ] );

bool
IsNodeInPar(
const mesh_t*     mesh,
const particle_t* particle,
const int      node );

void
ParElemNodeTouch(
const mesh_t*     mesh,
const particle_t* particle,
const int      element,
      int**    elemT,
      int**    elemP,
      int**    NodesIn );

#endif
