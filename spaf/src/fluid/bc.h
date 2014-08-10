#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H

typedef struct bc_struct bc_t;

bc_t*
ReadBC(
const char* FileName );

void
SetVelBC(
const int       NbOfNodes,
const bc_t*     bc,
      double**  Vel );

#endif
