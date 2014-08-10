#ifndef OPERATORS_H
#define OPERATORS_H

#include "includes.h"
#include "mesh.h"
#include "mass_conservation.h"

void UseMacroGrid(void);
void UseMicroGrid(void);
bool UsingMicroGrid(void);

void
ReadPreconditionerParameters(
FILE* FileId );

void
GetOperators(
const double  tau0,
const double  Re,
const mesh_t* mesh );

void
UpdateOperators(
const double DeltaTau );

void
ApplyOperator(
const char*   what,
      double* x[],
      double* rhs[] );

void
SolveOperator(
const char*   what,
const bool*   OutsideNodes,
      double* rhs[],
      double* x[] );

#endif
