#ifndef FLUID_OUTPUT_VTK_H
#define FLUID_OUTPUT_VTK_H

#include "mesh.h"

void
FluidWriteVTK(
const char*   FileName,   ///< file name
const bool    binary,     ///< file format : ascii or binary
const bool    WritePre,   ///< save or not pressure field
const mesh_t* mesh,       ///< mesh structure
const double**  Vel,        ///< fluid's velocity field
const double*   Pre );      ///< fluid's pressure field

#endif
