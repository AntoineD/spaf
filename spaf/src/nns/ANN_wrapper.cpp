/** \file
 Wrapper for the ANN C++ library.

 ANN has been modified to allow choice between double or float real values, the preprocessor flag ANN_USE_FLOAT will force ANN to use float instead of the default double. Choosing float is a little faster (5 to 10 %) and has no drawback on the result as long as the minimum point to point distance across all points is larger than FLT_EPSILON.
 A function that return whether or not a point is inside the bounding box of the set of points has been added as well.
 Finally the preprocessor macro ANN_REPORT_POINT_OUTSIDE_BOUNDING_BOX is used to force ANN to stop a search when a point is outside the bounding box, though it is not intended to be passed at compilation, search for it in ANN files in order to delete it. (it is on by default)
 
 ANN offers different tree structures and search algorithms, I choosed the fastest ones by testing with a typical mesh. Other choices have been commented out, see the wrappers below.
 
 WARNING : this library is not thread safe !!!!
 */

#include "ANN.h"
#include "ANN_wrapper.h"

using namespace std;

// filewise static variables
// As we search for only one nearest neighbor, we use those arrays for all calls
static ANNidxArray _nnIdx = NULL; // near neighbor indices
static ANNdistArray _dists = NULL;  // near neighbor distances

//==============================================================================
/** Prepare the nns structure.
 
 Return a pointer to the structure.
 */
//==============================================================================
void*
PrepareNNS_ANN(
const int       NbOfPoints, ///< number of points
      double**  points )    ///< points set
{
  _nnIdx = new ANNidx[1];
  _dists = new ANNdist[1];

  ANNpointArray dataPts = annAllocPts(NbOfPoints, 3);  // allocate data points

  for (int i = 0; i < NbOfPoints; i++)
  {
    dataPts[i][0] = points[i+1][1];
    dataPts[i][1] = points[i+1][2];
    dataPts[i][2] = points[i+1][3];
  }
  
  return (void*) new ANNkd_tree(dataPts, NbOfPoints, 3);
//  return (void*) new ANNbd_tree(dataPts, NbOfPoints, 3);
}
//==============================================================================
/** Test if a point is inside the bounding box of a mesh.
 */
//==============================================================================
bool
IsInsideBoundingBox_ANN(
const void*   tree,       ///< pointer to nns structure
const double  point[4] )  ///< point
{
#ifdef ANN_USE_FLOAT
  float _point[3] = {(float)point[1],(float)point[2],(float)point[3]};
#else
  double* _point = (double*) &point[1];
#endif
  return ((ANNkd_tree*)tree)->annIsInsideBoundingBox(_point);
}
//==============================================================================
/** Search for the nearest point.
 
 Return the node index.
 */
//==============================================================================
unsigned int
DoNNS_ANN(
const void*   tree,       ///< pointer to nns structure
const double  point[4] )  ///< point
{
#ifdef ANN_USE_FLOAT
  float _point[3] = {(float)point[1],(float)point[2],(float)point[3]};
#else
  double* _point = (double*) &point[1];
#endif
  
  // initialize in such a way we know if the point is outside the bounding box
  // if so annkSearch will return _nnIdx unchanged and DoNNS_ANN will return -1
  _nnIdx[0] = -2;
  
  ((ANNkd_tree*)tree)->annkSearch(_point, 1, _nnIdx, _dists, 0.);
//  ((ANNkd_tree*)tree)->annkPriSearch(_point, 1, _nnIdx, _dists, 0.);
  
  return _nnIdx[0]+1;
}
