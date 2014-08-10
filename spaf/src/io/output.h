#ifndef OUTPUT_H
#define OUTPUT_H

#include "includes.h"
#include "fluid.h"
#include "particle.h"

void
SetOutputDirectory(
char* PathToOutput );

char*
GoToOutputDirectory( void );

void
ReadOutputParameters(
FILE* FileId );

void
WriteOutput(
const FluidVar_t* FluidVar,
const particle_t* particle,
const double      time,
const int         TimeStep );

void
PrintTitle(
const char title[] );

void
PrintTimeStep(
const double t,
const int timestep );

#endif
