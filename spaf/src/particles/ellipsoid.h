/*=============================================================================
File: ellipsoid.h

Description: Routines for ellipsoid particles

By: Antoine
Last Update: 19/Sep/06
=============================================================================*/
#ifndef PARTICLE_KIND_ELLIPSOID_H
#define PARTICLE_KIND_ELLIPSOID_H

#include "planes.h"

typedef struct
{
  double SemiAxes[4]; // semi-axes
}
ellipsoid_t;

ellipsoid_t*
EllipsoidCreate( void );

void
EllipsoidPrint(
const ellipsoid_t* ellipsoid );

ellipsoid_t*
EllipsoidRead(
FILE* FileId );

bool
PointIsInEllipsoid(
const ellipsoid_t* ellipsoid,
const double       position[ 4 ] );

double
GetEdgeIntersectParamEllipsoid(
const ellipsoid_t* ellipsoid,
      double       Pin[4],
      double       Pout[4] );

double
EllipsoidVolume(
const ellipsoid_t* ellipsoid );

void
EllipsoidInertiaTensorInv(
const ellipsoid_t*  ellipsoid,
      double        density,
      double**      InertiaTensorInv );

double
EllipsoidPlaneDist(
const ellipsoid_t* ellipsoid,
const double       ParPos[4],
      double**     RotMat,
const equation_t*  PlaneEquation );

#endif
