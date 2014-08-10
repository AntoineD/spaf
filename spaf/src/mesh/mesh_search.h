#ifndef mesh_search_H
#define mesh_search_H

#include "mesh.h"

void
ReadMeshSearchParameters(
FILE* FileId );

bool
MeshGetLocalPosition(
const mesh_t* mesh,
const int  el,
const double    Lag_node[4],
      double    LocalPos[4] );

int
SearchElemContainingPoint(
const mesh_t* mesh,
      bool    UseNNS,
const double  point[4],
      int*   ClosestNode_ptr,
      int*   ElemSearchedNb );

#endif
