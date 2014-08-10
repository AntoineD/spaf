#ifndef MASS_CONSERVATION_H
#define MASS_CONSERVATION_H

#include "includes.h"
#include "mesh.h"

void
GetMassConservingOperator(
                          const char*         PathToCase,
                          const mesh_t*       mesh );

void
MassConserveBC(
            const double    position[4],
            const mesh_t*   mesh_o,
            double**  Vel_o,
            const mesh_t*   mesh,
            double**  Vel);

#endif
