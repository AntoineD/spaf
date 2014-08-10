/** \file
 Entry points to the nearesr neighbor search routines.
 */

#include "nns.h"
#include "logging.h"
#include "ANN_wrapper.h"
#include "legacy.h"
#include "parse.h"
#include "strings.h"

// static local file-wise variable
// switch between Legacy and ANN NNS
static bool _use_ann = false;
// to make sure PrepareNNS has been called before calling other routines
static bool _PrepareNNS_done = false;

// parameters for legacy nns :
// coefficients used for computing the voxel size defined as
// _coef_min * MinNode2NodeDist + _coef_max * MaxMinNodeDist + _coef_round_off
static double
  _coef_max = -DBL_MAX,
  _coef_min = -DBL_MAX,
  _coef_round_off = -DBL_MAX;

//==============================================================================
/// Read bounding box parameters
//==============================================================================
void
ReadNNSParameters(
FILE* FileId )  ///< file id
{
  // read nns type
  TokenizeFileLine( FileId );
  assert_error( TokenNameIs("type"),
               "nearest_neighbor_search parameters : type not found");
  if ( StringCompare(Token(3), "ann") )
  {
    _use_ann = true;
    
    // check we are running serial
#ifdef _OPENMP
#pragma omp parallel
    {
#pragma omp master
      {
        assert_error( omp_get_num_threads() == 1,
                     "ANN is not thread safe, fix it or run serial");
      }
    }
#endif
  }
  else if ( StringCompare(Token(3), "legacy") )
  {
    _use_ann = false;
    _coef_min = Token2double(4,'[',0.,DBL_MAX,']');
    _coef_max = Token2double(5,'[',_coef_min,DBL_MAX,']');
    _coef_round_off = Token2double(6,'[',0.,DBL_MAX,']');
  }
  else
    error("nearest_neighbor_search : bad type : %s", Token(3));

  // check end of section
  TokenizeFileLine(FileId);
  assert_error( SectionIsEnd(),
               "nearest_neighbor_search parameters : missing section = end");

  // print infos
  info("\nNearest neighbor search parameters :\n" ); 
  
  if ( _use_ann )
    info("\tusing ANN nns""\n");
  else
  {
    info("\tusing legacy nns""\n");
    info("\tcoef_min = "       "%g" "\n", _coef_min);
    info("\tcoef_max = "       "%g" "\n", _coef_max);
    info("\tcoef_round_off = " "%g" "\n", _coef_round_off);    
  }
}
//==============================================================================
/** Prepare the nns structure.
 
 Return a pointer to the structure.
 */
//==============================================================================
nns_t*
PrepareNNS(
const int       NbOfPoints,       ///< number of points
      double**  points,           ///< points set
      int**     connectivity,     ///< connectivity
const SparseMatrix_t* InvConTab ) ///< inverse connectivity
{
  _PrepareNNS_done = true;

  if ( _use_ann )
    return PrepareNNS_ANN(NbOfPoints, points);
  else
    return PrepareNNS_Legacy(_coef_min, _coef_max, _coef_round_off,
                             NbOfPoints, points, connectivity, InvConTab);
}
//==============================================================================
/** Test if a point is inside the bounding box of a mesh.
 */
//==============================================================================
bool
IsInsideBoundingBox(
const nns_t*  tree,       ///< pointer to nns structure 
const double  point[4] )  ///< point
{
  assert_error(_PrepareNNS_done, "PrepareNNS must be called first");
  
  if ( _use_ann )
    return IsInsideBoundingBox_ANN(tree, point);
  else
    return IsInsideBoundingBox_Legacy(tree, point);
}
//==============================================================================
/** Search for the nearest point.
 
 Return the node index.
 */
//==============================================================================
int
DoNNS(
const nns_t*  tree,       ///< pointer to nns structure 
const double  point[4] )  ///< point
{
  assert_error(_PrepareNNS_done, "PrepareNNS must be called first");
  
  if ( _use_ann )
    return DoNNS_ANN(tree, point);
  else
    return DoNNS_Legacy(tree, point);
}
