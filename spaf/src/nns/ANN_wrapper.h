#ifndef ANN_WRAPPER_H
#define ANN_WRAPPER_H

#ifdef __cplusplus
extern "C" {
#endif

void*
PrepareNNS_ANN(
const int       NbOfPoints,
      double**  points );

bool
IsInsideBoundingBox_ANN(
const void*   tree,
const double  point[4] );
  
unsigned int
DoNNS_ANN(
const void*   tree,
const double  point[4] );

#ifdef __cplusplus
}
#endif

#endif
