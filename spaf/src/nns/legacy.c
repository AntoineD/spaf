/** \file
 Legacy or original nearest neighbor search routines.
 
 Voxel tables for a fast search of the location of a point in the nns is based on a voxel edge distance = the maximum (over the mesh) of the minimum internodal distance, i.e. max i (min j (d_ij)), where d_ij is the distance between nodes i & j, j loops over all nodes touching node i, and i loops over all nodes in nns.
 
 PROBLEM : when the mesh does not have uniform sized elements, this code runs VERY slowly when searching for points in the highly refined part of the mesh.
*/

#include "output.h"
#include "memory.h"
#include "parse.h"
#include "logging.h"
#include "file_system.h"
#include "legacy.h"

/// flags for result from search
enum
{
  PointIsOutsideBBox = -1,
  VoxelIsEmpty = 0
};

struct nns_struct
  {
    double
    **points,  ///< pointer to a set of points
    min_x, ///< minimum coordinate in the mesh in x-direction
    max_x, ///< maximum coordinate in the mesh in x-direction
    min_y, ///< minimum coordinate in the mesh in y-direction
    max_y, ///< maximum coordinate in the mesh in y-direction
    min_z, ///< minimum coordinate in the mesh in z-direction
    max_z, ///< maximum coordinate in the mesh in z-direction
    step;  ///< edgelength of the bounding box voxels
    
    int
    NbOfPoints, ///< number of points in the set to search in
    size_x,  ///< size of the grid in x-direction (in units of "step")
    size_y,  ///< size of the grid in y-direction (in units of "step")
    size_z,  ///< size of the grid in z-direction (in units of "step")
    n_slices, /* size of the arrays slice_info1 & slice_info2 in the "slice direction".
     Here the idea is that the bounding structured mesh box is "sliced"
     perpendicular to the axis aligned with the longest edge of the mesh box.
     After slices, the next direction is called "rows", then "columns".*/
    d_slices,  ///< direction of the slices. Legend: 1 means x, 2 means y, 3 means z
    d_rows,    ///< direction of the rows. See above for legend
    d_columns, ///< direction of the columns. See above for legend
    n_bb,           ///< the number of non-empty voxels == the size of the array column_info
    **slice_info1,  /* Works with slice_info2.  For row number q* being searched for,
     [k][0]th entry stores an relative offset into column_info for row k.
     [k][p]th entry stores an relative offset, so that the total offset
     into column_info is entry[k][0] + entry[k][p].
     This offset is then used to find the column corresponding to the
     voxel of interest.
     Size: [n_slices][number of non-empty rows in a slice+1]
     Note that slice_info1 has exactly the same *structure* as slice_info2.
     However. entries differ between these matrices. */
    **slice_info2,  /* [k][0]th entry stores the # entries in row k of the matrix.
     The remainder of the kth row stores an ordered list of non-empty
     row numbers from the structured mesh.  Therefore, if the [k][p]th 
     entry equals q, then row q is the pth non-empty row in slice k.
     This array is binary searched for a specified q*; if found, p is then
     used as an index into slice_info1.  Size:
     [n_slices][number of non-empty rows in a slice+1] */
    **points_in_bb, /* row n_bb stores the numbers of points within a certain voxel;
     the offsets for the slices and rows are stored in slice_info1;
     the number of column to which the entry corresponds is stored in
     column_info. Size [n_bb][?] */
    *column_info;   /* If the rth entry = s, then row s in points_in_bb contains a list of
     all the points in voxel s.  The value of r is obtained as an
     offset from array slice_info1.  Size: [n_bb] */
  };

/************************************************************************
FUNCTIONS: bin_search
PURPOSE: performs a binary search in an ordered array of integers 
PARAMETERS:
   Input: array_ptr: the pointer to the 0-th element of the array
          element: the integer to be searched for
RETURNS: the position of element if it is found in the array;
          the position of the right closest entry of the array if element
          is not found in the array(with - sign); if it's to the right of all entries
          of the array it returns the size of the array + 1.
************************************************************************/
static long
bin_search(
const int* array_ptr,
const int  element,
const int  l_count,
const int  r_count )
{
  int left, right, m_count1, m_count2, middle1, middle2, shift;

  left = array_ptr[ l_count ];
  right = array_ptr[ r_count ];

  if ( element > right )
    return - ( ( long int ) r_count + 1 );
  else if ( element < left )
    return - ( long int ) l_count;

  shift = (int) ( ( ( float ) ( r_count - l_count ) ) / 2e0 + 1e-1 );
  m_count1 = (int) ( l_count + shift );
  middle1 = array_ptr[ m_count1 ];
  m_count2 = (int) ( r_count - shift );
  middle2 = array_ptr[ m_count2 ];

  if ( element == left )
    return ( long int ) l_count;
  
  else if ( element == right )
    return ( long int ) r_count;
  
  else if ( element < middle2 && element > middle1 )
    return - ( long int ) m_count2;
  
  else
  {
    while ( shift != 0 )
    {
      if ( element < middle1 )
        return bin_search( array_ptr, element, l_count, m_count1 );

      else if ( middle2 < element )
        return bin_search( array_ptr, element, m_count2, r_count );

      else if ( middle1 == element )
        return ( long int ) m_count1;
      
      else if ( middle2 == element )
        return ( long int ) m_count2;
    }

    return - r_count;
  }
}
/************************************************************************
FUNCTION: shell
 PURPOSE: sorts an array into ascending numerical order while making the corresponding rearrangment
          of another, integer array
PARAMETERS:
   Input:
     - n: the size of the arrays
     - arr: int array to be sorted
     - brr: int array to be rearranged
   Output:
     - appropriately filled arr and brr
 RETURNS: void
************************************************************************/
static void
shell(
const int     n,
      int* arr,
      int* brr )
{
  int j, i;
  int t, u;
  double ALN2I = 1.442695022,
         TINY = 1.0e-5;

  int lognb2 = ( int ) ( log( ( double ) n ) * ALN2I + TINY ),
      m = n;
  
  for ( int nn = 1 ; nn <= (int) lognb2 ; nn++ )
  {
    m >>= 1;
    for ( j = m + 1;j <= n;j++ )
    {
      i = j - m;
      t = arr[ j ];
      u = brr[ j ];
      while ( i >= 1 && arr[ i ] > t )
      {
        arr[ i + m ] = arr[ i ];
        brr[ i + m ] = brr[ i ];
        i -= m;
      }
      arr[ i + m ] = t;
      brr[ i + m ] = u;
    }
  }
}
//==============================================================================
/** Search in the bb for a voxel that contains a point.
 
 Return PointIsOutsideBBox if the point is outside the bb, the column position if it exist, otherwise returns VoxelIsEmpty.
 */
//==============================================================================
static long
GetNNSVoxel(
const nns_type* nns,        ///< bounding box
const double    point[4] )  ///< the point of interest
{
  // Is point inside bb ?
  if ( point[1] < nns->min_x || point[1] > nns->max_x ||
      point[2] < nns->min_y || point[2] > nns->max_y ||
      point[3] < nns->min_z || point[3] > nns->max_z )
    return PointIsOutsideBBox;
  
  // Compute its position in the bounding boxes grid
  int pos[ 4 ] = {0,0,0,0};
  pos[ 1 ] = (int) ( ( point[1] - nns->min_x ) / nns->step + 1 );
  pos[ 2 ] = (int) ( ( point[2] - nns->min_y ) / nns->step + 1 );
  pos[ 3 ] = (int) ( ( point[3] - nns->min_z ) / nns->step + 1 );
  
  // position of the point in the voxel grid
  int slice    = pos[ nns->d_slices ],
  row    = pos[ nns->d_rows ],
  column = pos[ nns->d_columns ];
  
  // row in the bounding box
  long pos_row = bin_search(
    nns->slice_info2[ slice ], row, 1, nns->slice_info2[ slice ][ 0 ] );
  
  if ( pos_row < 0 ) return VoxelIsEmpty;
  
  // offsets of the row in the array bb->points_in_bb
  int
  offset2 = 0,
  offset1 = nns->slice_info1[ slice ][ 0 ]
  + nns->slice_info1[ slice ][ pos_row ];
  
  if ( pos_row < ( long ) nns->slice_info2[ slice ][ 0 ] )
    offset2 = nns->slice_info1[ slice ][ 0 ]
    + nns->slice_info1[ slice ][ pos_row + 1 ] - 1;
  
  else if ( slice < nns->n_slices )
    offset2 = nns->slice_info1[ slice + 1 ][ 0 ] - 1;
  
  else
    offset2 = nns->n_bb;
  
  assert( offset1 <= offset2 );
  
  // search column in the bounding box
  long pos_column = bin_search( nns->column_info, column, offset1, offset2 );
  
  if ( pos_column < 0 ) return VoxelIsEmpty;
  
  return pos_column;
}
/************************************************************************
FUNCTION: BBFillInfo
PURPOSE: fills the information about a grid of bounding boxes
PARAMETERS:
   Input:
     - mesh->points: the array of the coordinates of the nodal points
     - tot_nodes: the number of nodes in the mesh
     - bb: the structure for the bounding box grid which is partially filled
   Output:
     - the remaining entries of the structure bb are appropriately
       filled
 RETURNS: void
************************************************************************/
/// Fills the information about a grid of bounding boxes
static void
BBFillInfo(
nns_type*  nns )  ///< nns structure
{
  int i, iloop, j, k, l, /* loop variables and counters */
        *pos,  /* store the coordinates of the point in the grid of bounding boxes */
        *sb, *sa,  /* temporary rearrangement of the node numbers in order to speed up filling the struct  */

        slice, row, column, /* store the slice, row and column number of the current point */
        offset1, offset2;  /* store the offsets in column_info */
  long int pos_row, pos_column, pos_point;


  /* decide which is slice direction, which is row direction & which is column direction */
  nns->n_slices  = nns->size_x;
  nns->d_slices  = 1;
  nns->d_rows    = 2;
  nns->d_columns = 3;
  
  if ( nns->size_y > nns->n_slices )
  {
    nns->n_slices = nns->size_y;
    nns->d_slices = 2;
    nns->d_rows = 1;
    nns->d_columns = 3;
  }
  
  if ( nns->size_z > nns->n_slices )
  {
    nns->n_slices = nns->size_z;
    nns->d_slices = 3;
    nns->d_rows   = 1;
    nns->d_columns = 2;
  }
  
  nns->n_bb = 0;

  /***  Sort the points of the mesh with increasing slice-coordinate; the convention
        is that the three principle directions are called slice, row and column
        directions;  the slice direction is chosen to be the one containing the 
        largest number of voxels (bounding boxes); it is previously given in
        bb->d_slices.  sa has the "voxel co-ordinates" of each node in the mesh.
        sb has the corresponding node #.  After the shell sort is called, sa is
        ordered in increasing order.  sb is permuted in the same way as sa.
        Why do we do this?  Having the nodes sorted  in increase slice co-ordinate order means
        that we can step through all nodes in an orderly way.  This avoids the use of a linked
        list data structure, and means that the arrays in the bb structure can be "grown"
        always from the end, without having to insert elements into the middle of an existing
        array.  ***/

  sa = (int*) calloc( nns->NbOfPoints + 1, sizeof(int) );
  assert_error( sa != NULL, "Memory allocation error" );
    
  sb = ( int * ) calloc( nns->NbOfPoints + 1, sizeof(int) );
  assert_error( sb != NULL, "Memory allocation error" );

  for ( i = 1; i <= nns->NbOfPoints; i++ )
  {
    if ( nns->d_slices == 1 )
      sa[ i ] = (int) ( ( nns->points[ i ][ 1 ] - nns->min_x ) / nns->step + 1 );
    
    else if ( nns->d_slices == 2 )
      sa[ i ] = (int) ( ( nns->points[ i ][ 2 ] - nns->min_y ) / nns->step + 1 );
    
    else if ( nns->d_slices == 3 )
      sa[ i ] = (int) ( ( nns->points[ i ][ 3 ] - nns->min_z ) / nns->step + 1 );
    
    sb[ i ] = i;
  }
  
  shell( nns->NbOfPoints, sa, sb ); /* Call sorting routine */

  /*** Check the sorting  ***/
  for ( iloop = 1; iloop < nns->NbOfPoints; iloop++ )
  {
    assert_error( sa[ iloop ] <= sa[ iloop + 1 ],
             "Wrong sorting : points %d %d.", sb[ iloop ], sb[ iloop + 1 ] );
  }

  free( sa );

  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

  pos = (int*) calloc( 4, sizeof(int) );
  assert_error( pos != NULL, "Memory allocation error" );

  /***  Alocate memory for the arrays of the structure bb ***/
  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  nns->slice_info1 = (int**) calloc( nns->n_slices + 1, sizeof(int*) );
  assert_error( nns->slice_info1 != NULL, "Memory allocation error" );

  for ( i = 1; i <= nns->n_slices; i++ )
  {
    nns->slice_info1[ i ] = (int*) calloc( 2, sizeof(int) );
    assert_error( nns->slice_info1[i] != NULL, "Memory allocation error" );
  }

  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  nns->slice_info2 = (int**) calloc( nns->n_slices + 1, sizeof(int*) );
  assert_error( nns->slice_info2 != NULL, "Memory allocation error" );

  for ( i = 1; i <= nns->n_slices; i++ )
  {
    nns->slice_info2[ i ] = (int*) calloc( 2, sizeof(int) );
    assert_error( nns->slice_info2[ i ] != NULL, "Memory allocation error" );

    nns->slice_info2[ i ][ 0 ] = 0;
  }

  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  nns->points_in_bb = (int**) calloc( 2, sizeof(int*) );
  assert_error( nns->points_in_bb != NULL, "Memory allocation error" );

  nns->points_in_bb[ 1 ] = (int*) calloc( 1, sizeof(int) );
  assert_error( nns->points_in_bb[ 1 ] != NULL, "Memory allocation error" );

  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  nns->column_info = (int*) calloc( 1, sizeof(int) );
  assert_error( nns->column_info != NULL, "Memory allocation error" );

  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  /*   Initialize the offsets for the slices  */
  for ( i = 1; i <= nns->n_slices; i++ )
    nns->slice_info1[ i ][ 0 ] = 1;

  /*** Now do the hard job - filling the arrays of the structure ***/

  for ( iloop = 1; iloop <= nns->NbOfPoints; iloop++ )
  {
    i = sb[ iloop ];

    pos[ 1 ] = (int) ( ( nns->points[ i ][ 1 ] - nns->min_x ) / nns->step + 1 );
    pos[ 2 ] = (int) ( ( nns->points[ i ][ 2 ] - nns->min_y ) / nns->step + 1 );
    pos[ 3 ] = (int) ( ( nns->points[ i ][ 3 ] - nns->min_z ) / nns->step + 1 );

    slice = pos[ nns->d_slices ];
    row = pos[ nns->d_rows ];
    column = pos[ nns->d_columns ];
    assert_error( slice <= nns->n_slices,
             "slice # larger than the total number of slices" );


    /**  Find the position of row in slice_info2 **/
    if ( nns->slice_info2[ slice ][ 0 ] < 1 )  /* No entries yet in this row */
      pos_row = -1;
    else
      pos_row = bin_search( nns->slice_info2[ slice ], row, 1, nns->slice_info2[ slice ][ 0 ] );

    if ( pos_row < 0 )
    {
      /*  Point not found; we must realloc all the arrays of the structure;
          it is supposed that the points are sorted
          into ascending slice coordinate thus only the current row offsets must be updated.
          See documentation for bin_search to understand possible returned values  */

      pos_row = -pos_row;
      nns->slice_info2[ slice ][ 0 ] += 1;
      nns->n_bb += 1;

      nns->slice_info2[ slice ] = (int*) realloc(
        nns->slice_info2[ slice ], (nns->slice_info2[slice][0]+1) * sizeof(int) );
      assert_error( nns->slice_info2[ slice ] != NULL, "Memory allocation error" );

      nns->slice_info1[ slice ] = ( int * ) realloc( nns->slice_info1[ slice ],
                                                       ( nns->slice_info2[ slice ][ 0 ] + 1 ) * sizeof(int) );
      assert_error( nns->slice_info1[ slice ] != NULL, "Memory allocation error" );

      nns->column_info = ( int * ) realloc( nns->column_info,
                                              ( nns->n_bb + 1 ) * sizeof(int) );
      assert_error( nns->column_info != NULL, "Memory allocation error" );

      nns->points_in_bb = ( int ** ) realloc( nns->points_in_bb,
                                               ( nns->n_bb + 1 ) * sizeof( int * ) );
      assert_error( nns->points_in_bb != NULL, "Memory allocation error" );

      nns->points_in_bb[ nns->n_bb ] = ( int * ) calloc( 1, sizeof(int) );
      assert_error( nns->points_in_bb[ nns->n_bb ] != NULL, "Memory allocation error" );

      /*  Shift the entries in points_in_bb and column_info so that we can fit the new entry into the
          middle somewhere */

      if ( pos_row == ( ( long int ) nns->slice_info2[ slice ][ 0 ] ) )
      {
        /*  The new row is the rightmost row in the current slice; we do not shift it; update the offset for it */
        offset1 = (int) nns->n_bb;
        nns->slice_info1[ slice ][ pos_row ] = (int) ( nns->n_bb - nns->slice_info1[ slice ][ 0 ] );
      }
      else
        offset1 = (int) ( nns->slice_info1[ slice ][ 0 ] + nns->slice_info1[ slice ][ pos_row ] );

      for ( j = nns->n_bb; j >= offset1 + 1; j-- )
        nns->column_info[ j ] = nns->column_info[ j - 1 ];

      for ( j = nns->n_bb; j >= offset1 + 1; j-- )
      {
        nns->points_in_bb[ j ] = ( int * ) realloc( nns->points_in_bb[ j ],
                                                     ( nns->points_in_bb[ j - 1 ][ 0 ] + 1 ) * sizeof(int) );
        assert_error( nns->points_in_bb[ j ] != NULL, "Memory allocation error" );

        for ( k = 0; k <= nns->points_in_bb[ j - 1 ][ 0 ]; k++ )
          nns->points_in_bb[ j ][ k ] = nns->points_in_bb[ j - 1 ][ k ];
      }
      /*  Put the number of the current point into corresponding entry of points_in_bb  */
      nns->points_in_bb[ offset1 ] = ( int * ) realloc( nns->points_in_bb[ offset1 ],
                                                         2 * sizeof(int) );
      assert_error( nns->points_in_bb[ offset1 ] != NULL, "Memory allocation error" );

      nns->points_in_bb[ offset1 ][ 0 ] = 1;
      nns->points_in_bb[ offset1 ][ 1 ] = i;
      nns->column_info[ offset1 ] = column;

      /* Shift entries in the arrays info1 and info2      */

      for ( j = nns->slice_info2[ slice ][ 0 ]; j >= (int) pos_row + 1; j-- )
      {
        nns->slice_info2[ slice ][ j ] = nns->slice_info2[ slice ][ j - 1 ];
        nns->slice_info1[ slice ][ j ] = (int) ( nns->slice_info1[ slice ][ j - 1 ] + 1 );
      }
      nns->slice_info2[ slice ][ pos_row ] = row;

      /*  Update the offsets of the slices below the current one  */

      for ( j = (int) ( slice + 1 ); j <= nns->n_slices; j++ )
        nns->slice_info1[ j ][ 0 ] = (int) ( nns->n_bb + 1 );
    } /* END if; case point is not found in the search through slices and rows */
    else
    {
      /**  The slice and row are not empty; check if the column is empty  **/
      offset1 = (int) ( nns->slice_info1[ slice ][ 0 ] + nns->slice_info1[ slice ][ pos_row ] );
      if ( pos_row < ( ( long int ) nns->slice_info2[ slice ][ 0 ] ) )
        offset2 = (int) ( nns->slice_info1[ slice ][ 0 ] + nns->slice_info1[ slice ][ pos_row + 1 ] - 1 );
      else
        offset2 = nns->n_bb;

      assert( offset1 <= offset2 );

      pos_column = bin_search( nns->column_info, column, offset1, offset2 );

      if ( pos_column < 0 )
      {
        /*  slice and row  already exist but the column not */
        pos_column = -pos_column;
        /**  Realloc column_info and points_in_bb  **/

        nns->n_bb += 1;

        nns->column_info = (int*) realloc ( nns->column_info,
                                                ( nns->n_bb + 1 ) * sizeof(int) );
        assert_error( nns->column_info != NULL, "Memory allocation error" );

        nns->points_in_bb = (int**) realloc( nns->points_in_bb,
                                                 ( nns->n_bb + 1 ) * sizeof(int*) );
        assert_error( nns->points_in_bb != NULL, "Memory allocation error" );
        
        nns->points_in_bb[ nns->n_bb ] = (int*) calloc( 1, sizeof(int) );
        assert_error( nns->points_in_bb[ nns->n_bb ] != NULL, "Memory allocation error" );

        /*  Shift the entries in points_in_bb and column_info */

        offset1 = pos_column;
        for ( j = nns->n_bb; j >= offset1 + 1; j-- )
        {
          nns->column_info[ j ] = nns->column_info[ j - 1 ];
        }

        for ( j = nns->n_bb; j >= offset1 + 1; j-- )
        {
          nns->points_in_bb[ j ] = ( int * ) realloc( nns->points_in_bb[ j ],
                                                       ( nns->points_in_bb[ j - 1 ][ 0 ] + 1 ) * sizeof(int) );
          assert_error( nns->points_in_bb[ j ] != NULL, "Memory allocation error" );

          for ( k = 0; k <= nns->points_in_bb[ j - 1 ][ 0 ]; k++ )
            nns->points_in_bb[ j ][ k ] = nns->points_in_bb[ j - 1 ][ k ];
        }
        
        /*  Put the number of the current point into corresponding entry of points_in_bb  */
        nns->points_in_bb[ offset1 ][ 0 ] = 1;
        nns->column_info[ offset1 ] = column;
        nns->points_in_bb[ offset1 ] = ( int * ) realloc( nns->points_in_bb[ offset1 ],
                                                           2 * sizeof(int) );
        assert_error( nns->points_in_bb[ offset1 ] != NULL, "Memory allocation error" );
        
        nns->points_in_bb[ offset1 ][ 1 ] = i;
        nns->column_info[ offset1 ] = column;

        /*  Update the offsets of the rows to the right of the current one */

        for ( j = nns->slice_info2[ slice ][ 0 ]; j >= (int) pos_row + 1; j-- )
          nns->slice_info1[ slice ][ j ] += 1;

        /*  Update the offsets of the slices below the current one  */

        for ( j = (int) ( slice + 1 ); j <= nns->n_slices; j++ )
          nns->slice_info1[ j ][ 0 ] = (int) ( nns->n_bb + 1 );

      } /* END if; the bounding box was empty */
      else
      {
        /**  The bounding box contains an entry already  **/
        /*  Put the number of the current point into corresponding entry of points_in_bb  */
        pos_point = nns->points_in_bb[ pos_column ][ 0 ] += 1;
        nns->points_in_bb[ pos_column ] = ( int * ) realloc( nns->points_in_bb[ pos_column ],
                                                              ( nns->points_in_bb[ pos_column ][ 0 ] + 1 ) * sizeof(int) );
        assert_error( nns->points_in_bb[ pos_column ] != NULL, "Memory allocation error" );
        nns->points_in_bb[ pos_column ][ pos_point ] = i;
      } /* END else; the bounding box was not empty */
    } /* END else; the slice  and row were not empty */
    /*
    info("\n Point : %d", i);
    info("\n Number of bounding boxes: %d", bb->n_bb);
    for( j=1; j<=bb->n_bb; j++ )
    print_ivector( bb->points_in_bb[j], bb->points_in_bb[j][0], "\n bounding box %d \n",j);
     
    for( j=1; j<=bb->n_slices; j++ )
  {
    info("\n Offset for slice: %d - %d \n", j, bb->slice_info1[j][0]);
    print_ivector( bb->slice_info1[j], bb->slice_info2[j][0], "\n offsets for slice %d \n",j);
  }
    */
  } /* END loop over the nodes of the mesh; index iloop */
  
  /*for( i=1; i<=bb->n_bb; i++ )
  print_ivector( bb->points_in_bb[i], bb->points_in_bb[i][0], "\n bounding box %d \n",i);

  for( i=1; i<=bb->n_slices; i++ )
  print_ivector( bb->slice_info1[i], bb->slice_info2[i][0], "\n offsets for slice %d \n",i);  
  */

  /***   Now check if the info is correct.  We do this by computing the voxel co-ordinates of
   each mesh point, and by then searching for the point in the bb structure.  If the voxel
   co-ordinates returned by the bb structure search match the ones known beforehand, then
   everything is OK  ***/

  int NbOfEmptySlices = 0;
  for ( i = 1; i <= nns->n_slices; i++ )
  {
    /* count is # empty slices, should be zero */
    if ( nns->slice_info2[ i ][ 0 ] == 0 ) NbOfEmptySlices++; 

    for ( j = 1; j <= nns->slice_info2[ i ][ 0 ]; j++ )
    {
      offset1 = (int) ( nns->slice_info1[ i ][ 0 ] + nns->slice_info1[ i ][ j ] );
      
      if ( j == nns->slice_info2[ i ][ 0 ] )
      { /* if true, last row in this slice */
        if ( i == nns->n_slices )
          offset2 = nns->n_bb;
        else
          offset2 = (int) ( nns->slice_info1[ i + 1 ][ 0 ] - 1 );
      }
      else
        offset2 = (int) ( nns->slice_info1[ i ][ 0 ] + nns->slice_info1[ i ][ j + 1 ] - 1 );

      for ( k = offset1; k <= offset2; k++ )
      for ( iloop = 1; iloop <= nns->points_in_bb[ k ][ 0 ]; iloop++ )
      {
        l = nns->points_in_bb[ k ][ iloop ]; /* Nodal point # */
        
        pos[ 1 ] = (int) ( ( nns->points[ l ][ 1 ] - nns->min_x ) / nns->step + 1 );
        pos[ 2 ] = (int) ( ( nns->points[ l ][ 2 ] - nns->min_y ) / nns->step + 1 );
        pos[ 3 ] = (int) ( ( nns->points[ l ][ 3 ] - nns->min_z ) / nns->step + 1 );
        
        slice = pos[ nns->d_slices ];
        row = pos[ nns->d_rows ];
        column = pos[ nns->d_columns ];
        
        if ( slice != i ||
             row != nns->slice_info2[ i ][ j ] ||
             column != nns->column_info[ k ] )
        {

          warning( "\n Directions: slices %d, rows %d, columns %d; step: %lf\n",
                   nns->d_slices, nns->d_rows, nns->d_columns, nns->step );

          error( "\n Point %d (%d, %d, %d) is in a wrong box; bb-coord: %d, %d, %d",
                   nns->points_in_bb[ k ][ iloop ], slice, row, column,
                   i, nns->slice_info2[ i ][ j ], nns->column_info[ k ] );
        }
      } /* END loop; index iloop */
    } /* END loop; index j */
  } /* END loop; index i */
  
  
  // check that all nodes are in a voxel
  for ( int node = 1 ; node <= nns->NbOfPoints ; node++ )
  {
    assert_error( GetNNSVoxel( nns, nns->points[node] ) != 0,
            "Node "INT_FMT" is not in a voxel !",node);
  }
  
  // printf some infos
  info("\t" INT_FMT " empty slices" "\n", NbOfEmptySlices );
  
  WriteNNS(nns);
  
  free( sb );
  free( pos );
}
/************************************************************************
 PURPOSE: allocate the bounding box structure and set some default values
 PARAMETERS:
 
 RETURNS: allocated ans initialized bounding box structure
 ************************************************************************/
nns_type*
PrepareNNS_Legacy(
const double    coef_min,
const double    coef_max,
const double    coef_round_off,
const int       NbOfNodes,
      double**  points,
      int**     ConTab,
const SparseMatrix_t* InvConTab )
{
  PrintTitle( "Creating bounding box" );

  nns_type* nns = (nns_type*) malloc( sizeof(nns_type) );
  
  nns->points = points;
  nns->NbOfPoints = NbOfNodes;
  
  nns->min_x = nns->min_y = nns->min_z = + DBL_MAX;
  nns->max_x = nns->max_y = nns->max_z = - DBL_MAX;
  nns->step  = 0.;
  
  nns->size_x   = 0;
  nns->size_y   = 0;
  nns->size_z   = 0;
  nns->n_slices = 0;                        
  
  nns->d_slices  = 0;
  nns->d_rows    = 0;
  nns->d_columns = 0;
  
  nns->n_bb = 0;
  nns->slice_info1  = NULL;
  nns->slice_info2  = NULL;
  nns->points_in_bb = NULL;
  nns->column_info  = NULL;
  
  double MinNode2NodeDist = DBL_MAX,  // min over all nodes of min distance between a node and its neighbours, i.e. min node to node distance
  MaxMinNodeDist = 0.,  // max over all nodes of min distance between a node and its neighbours
  MaxNode2NodeDist = 0.;  // max over all nodes of min distance between a node and its neighbours
  
  /*
   Now find the minimum distance within the points of the mesh.
   Note from CRE:  shold rewrite this at some point to use node-node connectivity table stored in MATRIX struct
   */
  
  for ( int node = 1 ; node <= nns->NbOfPoints ; node++ )
  {
    double *point = nns->points[ node ];
    
    if ( point[ 1 ] < nns->min_x ) nns->min_x = point[ 1 ];
    if ( point[ 2 ] < nns->min_y ) nns->min_y = point[ 2 ];
    if ( point[ 3 ] < nns->min_z ) nns->min_z = point[ 3 ];
    
    if ( point[ 1 ] > nns->max_x ) nns->max_x = point[ 1 ];
    if ( point[ 2 ] > nns->max_y ) nns->max_y = point[ 2 ];
    if ( point[ 3 ] > nns->max_z ) nns->max_z = point[ 3 ];
    
    double dist_min = DBL_MAX;
    
    // go over elements which contain the node
    for ( int j = 1 ; j <= NbOfEntriesInRow(InvConTab, node-1); j++ )
    {
      int element = EntryNode(InvConTab, node-1, j-1);
      
      // go over nodes of this element
      for ( int local_node = 1 ; local_node <= NODES_PER_EL ; local_node++ )
      {
        int neighbor = ConTab[ element ][ local_node ];
        
        // compute node to node distance
        if ( neighbor != node )
        {
          double* neighbor_point = nns->points[ neighbor ];
          
          double dist = sqrt( pow( neighbor_point[ 1 ] - point[ 1 ] ,2) +
                             pow( neighbor_point[ 2 ] - point[ 2 ] ,2) +
                             pow( neighbor_point[ 3 ] - point[ 3 ] ,2) );
          
          if ( dist < dist_min ) dist_min = dist;
          if ( dist > MaxNode2NodeDist ) MaxNode2NodeDist = dist;
        }
      }
    }
    
    /*
     find the minimum and maximum of all local minimums; they are used as voxel sizes
     of the "large" and "small" structured grids;  in the current version only
     the "large" structured grid is used because it seems to be enough for a fast
     search but if necessary the bb2 structure can also be used
     */
    if ( dist_min > MaxMinNodeDist ) MaxMinNodeDist = dist_min;
    if ( dist_min < MinNode2NodeDist ) MinNode2NodeDist = dist_min;
  }
  
  /*
   nns->step is set so that we can guarantee
   that all voxels lying inside the original mesh have at least one point in them.
   In other words, if a voxel is empty, then we can be sure that it is outside the mesh.
   _coef_round_off is a safety factor for numerical roundoff reasons. We may want to later
   multiply this by some number > 1 if we decide that we want larger voxels.
   */
  nns->step = coef_max * MaxMinNodeDist
  + coef_min * MinNode2NodeDist
  + coef_round_off;
  
  nns->size_x = (int) ( ( nns->max_x - nns->min_x ) / nns->step + 1 );
  nns->size_y = (int) ( ( nns->max_y - nns->min_y ) / nns->step + 1 );
  nns->size_z = (int) ( ( nns->max_z - nns->min_z ) / nns->step + 1 );
  
  // print some infos
  info("\tmin node to node distance = "DBL_FMT"\n",         MinNode2NodeDist);
  info("\tmax ( min node to node distance ) = "DBL_FMT"\n", MaxMinNodeDist);
  info("\tmax node to node distance = "DBL_FMT"\n",         MaxNode2NodeDist);
  
  
  BBFillInfo( nns );

  info("\t" INT_FMT " non-empty voxels\n", nns->n_bb );
  info("\t" INT_FMT " empty voxels\n",
       nns->size_x*nns->size_y*nns->size_z-nns->n_bb );
  
  info( "\t""Directions :""\n");
  info( "\t\t""slices  : "INT_FMT"\n",nns->d_slices);
  info( "\t\t""rows    : "INT_FMT"\n",nns->d_rows);
  info( "\t\t""columns : "INT_FMT"\n",nns->d_columns);
  
  info("\tmin point = ( "DBL_FMT" , "DBL_FMT" , "DBL_FMT" )\n",
       nns->min_x, nns->min_y, nns->min_z);
  info("\tmax point = ( "DBL_FMT" , "DBL_FMT" , "DBL_FMT" )\n",
       nns->max_x, nns->max_y, nns->max_z);
  info("\tgrid step = "DBL_FMT"\n", nns->step);
  info("\tgrid size = " INT_FMT " x " INT_FMT " x " INT_FMT "\n",
       nns->size_x, nns->size_y, nns->size_z );
  
  return nns;
}
//==============================================================================
/** Search initial guess for projecting the point.

 Output element where point (x,y,z) lies in
*/
//==============================================================================
static int
NNS(
const long      voxel,      ///< voxel index
const nns_type* nns,        ///< bounding box
const double    point[4] )  ///< the point of interest
{
  int ClosestNode = 0;
  double min_dist = DBL_MAX;
  int *vox = nns->points_in_bb[ voxel ];
  
  for ( int l = 1 ; l <= vox[ 0 ] ; l++ )
  {
    int k = vox[ l ];
    double *mesh_point = nns->points[k];
    
    double dist = pow( mesh_point[ 1 ] - point[1] ,2)
              + pow( mesh_point[ 2 ] - point[2] ,2)
              + pow( mesh_point[ 3 ] - point[3] ,2);
    
    if ( dist < min_dist )
    {
      min_dist = dist;
      ClosestNode = k;
    }
  }
  
  return ClosestNode;
}
//==============================================================================
/** Test if a point is inside the bounding box of a mesh.
 */
//==============================================================================
bool
IsInsideBoundingBox_Legacy(
const nns_type* nns,        ///< pointer to nns structure 
const double    point[4] )  ///< point
{
  return
  point[1] < nns->min_x || point[1] > nns->max_x ||
  point[2] < nns->min_y || point[2] > nns->max_y ||
  point[3] < nns->min_z || point[3] > nns->max_z;
}
//==============================================================================
/** Search for the nearest point.
 
 Return the node index.
 */
//==============================================================================
int
DoNNS_Legacy(
const nns_type* nns,        ///< pointer to nns structure 
const double    point[4] )  ///< point
{
  // search for a voxel containing point
  long voxel = GetNNSVoxel( nns, point );
  return NNS(voxel, nns, point);
}
//==============================================================================
/** Write a bounding box in vtk for visualization.
 
 When bb is NULL, switch writing on, off by default.
 */
//==============================================================================
void
WriteNNS(
const nns_type* nns )  ///< pointer to a bounding box
{
  // function wise static variable for writing
  static bool _write = false;
  
  // set writing switch on ?
  if ( nns == NULL )
  {
    _write = true;
    return;
  }
  
  if ( _write == false ) return;
  
  info("Writing bounding box to disc.""\n");
  
  // Create a vtk file for visuzalization of the bbox
  // go to output directory
  // move to output directory ans save current working directory
  char* WorkingDirectory = GoToOutputDirectory();
  assert_Directory( WorkingDirectory, "output");
  
  // open file
  FILE* FileId = fopen( "bounding_box.vtk", "w" );
  assert_error(FileId, "Failed to open bounding_box.vtk""\n");
  
  // Write header
  fprintf( FileId, 
          "# vtk DataFile Version 3.0" "\n"
          "Bounding box" "\n"
          "ASCII" "\n"
          "DATASET STRUCTURED_POINTS" "\n"
          "DIMENSIONS "INT_FMT" "INT_FMT" "INT_FMT"\n"
          "ORIGIN "DBL_FMT" "DBL_FMT" "DBL_FMT"\n"
          "SPACING "DBL_FMT" "DBL_FMT" "DBL_FMT"\n"
          "CELL_DATA "INT_FMT "\n"
          "SCALARS voxel unsigned_int 1" "\n"
          "LOOKUP_TABLE default" "\n",
          nns->size_x+1,nns->size_y+1,nns->size_z+1,
          nns->min_x   ,nns->min_y   ,nns->min_z,
          nns->step    ,nns->step    ,nns->step,
          nns->size_x*nns->size_y*nns->size_z );  
  
  // for each voxel, write the voxel index of the voxel's center point
  for ( int i = 0 ; i < nns->size_z ; i++ )
  {
    double pos[ 4 ] = {0,0,0,0};
    pos[ 3 ] = nns->min_z + ( (double)i + 0.5 ) * nns->step;
    
    for ( int j = 0 ; j < nns->size_y ; j++ )
    {
      pos[ 2 ] = nns->min_y + ( (double)j + 0.5 ) * nns->step;
      
      for ( int k = 0 ; k < nns->size_x ; k++ )
      {
        pos[ 1 ] = nns->min_x + ( (double)k + 0.5 ) * nns->step;
        
        int voxel = GetNNSVoxel( nns, pos );
        if ( voxel == PointIsOutsideBBox ) voxel = VoxelIsEmpty;
        fprintf( FileId, INT_FMT"\n", (int) voxel);
      }
    }
  }
  
  fclose(FileId);
  
  // move back to original directory
  assert_GoBackToDirectory(WorkingDirectory);
  free(WorkingDirectory);
}
