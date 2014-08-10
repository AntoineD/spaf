#ifndef FLUID_IO_H
#define FLUID_IO_H

#include "includes.h"
#include "mesh.h"
#include "particle.h"
#include "fluid.h"

void
WriteField(
           const bool      WriteInBinary,  ///< write in binary or ascii
           const bool      WriteAsDouble,  ///< write as double or float
           const char*     PathToDir,      ///< path to directory 
           const int       NbOfPreNodes,   ///< number of pressure nodes
           const double*   Pre,            ///< pressure field
           const int       NbOfVelNodes,   ///< number of velocity nodes
           double**  Vel );           ///< velocity field

void
ReadField(
          const char*     PathToDir,    ///< directory to read to 
          const int       NbOfPreNodes, ///< number of pressure nodes
          double*   Pre,          ///< fluid variables
          const int       NbOfVelNodes, ///< number of velocity nodes
          double**  Vel );         ///< fluid variables

void
WriteVizData(
             const bool      WriteInBinary,  ///< write in binary or ascii
             const bool      WriteAsDouble,  ///< write as double or float
             const double    time,           ///< time
             const int       TimeStep,       ///< timestep
             const FluidVar_t* FluidVar);     ///< velocity field

#endif
