#ifndef FLUID_CONVECTION_H
#define FLUID_CONVECTION_H

#include "mesh.h"

void
ReadConvectionParameters(
FILE* FileId );

void
Convection(
const int       order,
const mesh_t*   mesh,
const double    dt,
const double    FrameVel[4],
      double**  VelOld1,
      double**  VelOld2,
      double**  VelTilde1,
      double**  VelTilde2 );

void
ConvectionMicroGrid(
const int       order,
const mesh_t*   mesh,
const double    dt,
const double    FrameVel[4],
const double    position[4],
const mesh_t*   mesh_o,
      double**  Vel_o,
      double**  VelOld1,
      double**  VelOld2,
      double**  VelTilde1,
      double**  VelTilde2 );

void
ConvectionMacroGrid(
const int       order,
const mesh_t*   mesh,
const double    dt,
const double    FrameVel[4],
const double    position[4],
const mesh_t*   mesh_o,
      double**  Vel_o,
      double**  VelOld1,
      double**  VelOld2,
      double**  VelTilde1,
      double**  VelTilde2 );

#endif
