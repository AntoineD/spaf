#ifndef PLANES_H
#define PLANES_H

#include "includes.h"

typedef struct
{
	double coef[4];
}
equation_t;

typedef struct
{
	int n;
	equation_t *equation;
}
planes_t;

planes_t*
PlanesRead(
const char* FileName );

double
PlaneToPointDistance(
const double      point[4],
const equation_t* equation );

#endif
