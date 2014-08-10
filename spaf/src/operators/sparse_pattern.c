/** \file
 Compute the sparsity pattern of a matrix.
 */

#include "memory.h"
#include "logging.h"
#include "sparse_pattern.h"

/// Max number of neighbors for one node
#define MAX_NB_OF_NEIGHBORS 200

//=============================================================================
/** Find the neighbors of a node.

 Find all nodes that belong to an element that includes the current node. The list of neighbors is reordered in such a way node index are increasing.
*/
//=============================================================================
static void
FindNeighbors(
const mesh_t* mesh,         ///< mesh structure
const int     theNode,      ///< node index
      bool*   NodeChecked,  ///< list of nodes checked 
      int*    neighbors )   ///< list of neighbors
{
  // Initialize the number of neighbors
  neighbors[0] = 0;

  // stores neighbors nodes before re-ordering, offset to increment list
  int NodeList[ MAX_NB_OF_NEIGHBORS ];

  // Go over the elements this node is in
  for ( int elem = 1 ; elem <= NbOfEntriesInRow(mesh->InvConTab, theNode-1); elem++ )
  {
    int theElem = EntryNode(mesh->InvConTab, theNode-1, elem-1);

    // go over local nodes of the element
    for ( int node = 1 ; node <= NODES_PER_EL ; node++ )
    {
      int otherNode = mesh->ConTab[ theElem ][ node ];

      // has this node been proccessed already ?
      if ( NodeChecked[ otherNode ] == false )
      {
        NodeChecked[ otherNode ] = true;
        
        neighbors[0] ++;
        
        assert_error( neighbors[0] <= MAX_NB_OF_NEIGHBORS,
          "Increase MAX_NB_OF_NEIGHBORS");
        
        NodeList[ neighbors[0] ] = otherNode;
      }
    }
  }

  // Order the NodeList of neighbors and initialize so that list goes increasing
  for ( int i = 1 ; i <= neighbors[0] ; i++ )
  {
    int minVal = NodeList[ i ],
        minPos = i;

    for ( int j = i + 1 ; j <= neighbors[0] ; j++ )
    {
      if ( NodeList[ j ] < minVal )
      {
        minVal = NodeList[ j ];
        minPos = j;
      }
    }

    // store min
    neighbors[ i ] = minVal;

    // swap
    NodeList[ minPos ] = NodeList[ i ];

    // initialize array for next usage
    NodeChecked[ neighbors[ i ] ] = false;
  }
}
//=============================================================================
/** Process sparse storage.
 
 Only a sub-matrix of the complete operator may be processed, which can be selected with the first and last column, and the last row. If the operator is symmetric, only the upper half part can be processed. Two maps are used to choose between processing velocity or pressure nodes. A NULL map means processing the velocity nodes. The fill switch is used to select between getting the storage size and filling it.
 When sizing (fill = false), on return the number of rows and entries have been computed and the NeighborsTable is filled.
 When filling (fill = true), on return the number of rows and entries have been computed and the offset, columns and rows arrays have been filled.
 */
//=============================================================================
static void
ProcessStorage(
const mesh_t* mesh,           ///< mesh structure
const int*    row_map,        ///< row mapping
const int     last_row,       ///< last row to be processed
const int*    col_map,        ///< column mapping
const int     first_col,      ///< last column in first submatrix
const int     last_col,       ///< last column in first submatrix
const bool    upper_half,     ///< process upper part only ?
      int*    NbOfEntries,    ///< nb of entries per sub-matrix for Grad
      int*    NbOfRows,       ///< nb of rows per sub-matrix for Grad
      int*    offset,         ///< pointer to offset array
      int*    columns,        ///< pointer to columns array
      int*    rows,           ///< pointer to rows array
      int**   NeighborsTable, ///< neighbors table
const bool    fill )          ///< switch for getting size or filling
{
  // initialize counters
  *NbOfEntries = 0,
  *NbOfRows = 0;

  // Record whether a node has already been checked
  // this is common to all FindNeighbors calls, only used there
  bool *NodeChecked = NULL;
  AllocVbool( mesh->NbOfNodes, NodeChecked );
  
  for ( int i = 1 ; i <= mesh->NbOfNodes ; i++ )
  {
    // copy loop counter as it may be overwritten
    int node = i;
    
    // skip if node is mapped to nothing
    if ( row_map && ! row_map[ node ] ) continue;
    
    // pointer to current node neighbors
    int *neighbors = NeighborsTable[ node ];
    
    // find neighbors
    if ( ! fill ) FindNeighbors( mesh, node, NodeChecked, neighbors );
    
    // get the mapped node index
    if ( row_map ) node = row_map[ node ];

    // skip if we are out of the sub-matrix
    if ( node > last_row ) continue;

    // used to know if we process a new row
    bool RowIsNew = false;
    
    // Go over all neighbor nodes
    for ( int i_neigh = 1 ; i_neigh <= neighbors[0] ; i_neigh++ )
    {
      // get neighbor index
      int neighbor = neighbors[ i_neigh ];
      
      // skip if node is mapped to nothing, otherwise get mapped node index
      if ( col_map )
      {
        neighbor = col_map[neighbor];
        if ( ! neighbor ) continue;        
      }

      // skip if we are out of the sub-matrix
      if ( neighbor > last_col || neighbor < first_col ) continue;
      if ( upper_half && neighbor < node ) continue;
      
      // update entries count
      (*NbOfEntries)++;
      
      // get entry's column
      if ( fill ) columns[ *NbOfEntries-1 ] = neighbor-1;
      
      // if it is a new row, update rows count and get entry's row
      if ( ! RowIsNew )
      {
        RowIsNew = true;
        (*NbOfRows)++;
        
        if ( fill ) rows[ *NbOfRows-1 ] = node-1;
      }
    }
    
    // get row's offset
    if ( fill ) offset[ *NbOfRows ] = *NbOfEntries;
  }
  
  free( NodeChecked );
}
//=============================================================================
/** Build the sparse pattern of an operator.
 
 Only a sub-matrix of the complete operator may be processed, which can be selected with the first and last column, and the last row. If the operator is symmetric, only the upper half part can be processed. Two maps are used to choose between processing velocity or pressure nodes. A NULL map means processing the velocity nodes.

 First the neighbor table is allocated, it holds the neighbors nodes of each nodes. Then the space required by the sparse storage is computed and the neighbor table is filled. Finally, the arrays that store the sparsity are allocated and filled.
 
 On return, the number of rows and entries, the offset, columns and rows arrays have been created and filled.
 */
//=============================================================================
void
GetSparsePattern(
const mesh_t* mesh,           ///< mesh structure
const int*    row_map,        ///< row mapping
const int     last_row,       ///< last row to be processed
const int*    col_map,        ///< column mapping
const int     first_col,      ///< last column in first submatrix
const int     last_col,       ///< last column in first submatrix
const bool    upper_half,     ///< process upper part only ?
      int*    NbOfRows,       ///< nb of rows per sub-matrix for Grad
      int*    NbOfEntries,    ///< nb of entries per sub-matrix for Grad
      int**   offset,         ///< pointer to offset array
      int**   columns,        ///< pointer to columns array
      int**   rows )          ///< pointer to rows array
                 
{
  info("\t""computing sparse pattern : ");

  // compute the maximum number of rows in order to allocate the neighbours table.
  int MaxNbOfRows = 0,
  *NbOfElemPerNode = NULL;
  AllocVint( mesh->NbOfNodes, NbOfElemPerNode );
  
  // compute the number of elements a node is in
  for ( int element = 1; element <= mesh->NbOfElements; element++ )
  for ( int LocalNodeId = 1; LocalNodeId <= NODES_PER_EL; LocalNodeId++ )
  {
    int node = mesh->ConTab[element][LocalNodeId];
    NbOfElemPerNode[ node ] ++;
  }
  
  // find the max
  for ( int NodeId = 1 ; NodeId <= mesh->NbOfNodes ; NodeId++ )
    if ( NbOfElemPerNode[ NodeId ] > MaxNbOfRows )
      MaxNbOfRows = NbOfElemPerNode[ NodeId ];
  
  free( NbOfElemPerNode );
  
  // neighbors table, this array is required by the the 2 ProcessStorage calls, this is why we allocate it here, stores neighbors nodes of each node, 10 * max row because there are 10 nodes per element and this is the maximum possible number of neighbors
  int **NeighborsTable = NULL;
  AllocMint( mesh->NbOfNodes, 10 * MaxNbOfRows, NeighborsTable );

  // compute required storage
  info("sizing, ");
  ProcessStorage(
    mesh, row_map, last_row, col_map, first_col, last_col, upper_half,
    NbOfEntries, NbOfRows, NULL, NULL, NULL, NeighborsTable, false );
  
  // allocate storage
  AllocVint( *NbOfRows, *offset );
  AllocVint( *NbOfRows-1, *rows );
  AllocVint( *NbOfEntries-1, *columns );    
  
  // fill storage
  info("filling""\n");
  ProcessStorage(
    mesh, row_map, last_row, col_map, first_col, last_col, upper_half,
    NbOfEntries, NbOfRows, *offset, *columns, *rows, NeighborsTable, true );

  FreeM( NeighborsTable );
}
