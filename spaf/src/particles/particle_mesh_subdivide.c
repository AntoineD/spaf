#include "includes.h"
#include "memory.h"
#include "linalg.h"
#include "particle_mesh_find_elem.h"
#include "logging.h"
#include "particle_mesh_subdivide.h"

// Nodes per face of element
#define NODES_PER_SURFACE 6
#define FALSE 0
#define  TRUE 1


static void
setZerosMi(
           int**  mat,
           int    m,
           int    n )
{
  for ( int i = 1 ; i <= m ; i++ )
    for ( int j = 1 ; j <= n ; j++ )
      mat[ i ][ j ] = 0;
}
/*--------------------------------------------------------------------
PURPOSE:
  get the coordinates of a node from a geometric table:
  from the fluid one or, if not in there, from the particle one.
  
INPUT:
  d: pointer to datastructure containing the fluid geometric table.
  GeoTab: geometric table for the new nodes
  node: the current node

OUTPUT:
  pos: node coordinates
  
By: Antoine 22/09/06
--------------------------------------------------------------------*/
static void
GetNodePos(
  const mesh_t* mesh,
        double**   GeoTab,
  const int     node,
        double     pos[4] )
{
  int dir;

  if ( node <= mesh->NbOfNodes ) // look in the fluid geometric table
    for ( dir = 1 ; dir <= 3 ; dir++ )
      pos[ dir ] = mesh->points[ node ][ dir ];
  else // look in the particle geometric table
    for ( dir = 1 ; dir <= 3 ; dir++ )
      pos[ dir ] = GeoTab[ dir ][ node - mesh->NbOfNodes ];
}
/*--------------------------------------------------------------------
PURPOSE:
  To solve for the intersection point between particle and edge.
  Let P = (x,y,z), Pin and Pout be the nodeIn and nodeOut coordinates.
  Equation of edge: vect(PinP) = k * vect(PinPout)
  then xi = ai + bi*t, i = 1,2,3 ; with a = Pin and b = Pout - Pin.
  Take t = 0 for nodeIn, t = 1 for nodeOut ==> get ai and bi
  From t ==> coordinates of new node (intersection point)

INPUT:
  d: pointer to datastructure
  particle: the particle
  nodeIn, nodeOut: The nodes at the ends of the edge
  GeoTab: geometric table for the new nodes

OUTPUT:
  t value

By: Veeramani
Update: 02/10/06 Antoine
--------------------------------------------------------------------*/
static double
GetEdgeIntersectParam(
const mesh_t*     mesh,
const particle_t* particle,
const int      nodeIn,
const int      nodeOut,
      double**      GeoTab )
{
  double t,
        Pin[ 4 ]  = {0.,0.,0.,0.},
        Pout[ 4 ] = {0.,0.,0.,0.};

  // get nodeIn and nodeOut positions in the regular coordinate system
  GetNodePos( mesh, GeoTab, nodeIn , Pin );
  GetNodePos( mesh, GeoTab, nodeOut, Pout );
  
  // express points in the particle coordinate system.
  ChangeFrame( particle->RotMat, particle->Pos0, Pin );
  ChangeFrame( particle->RotMat, particle->Pos0, Pout );

  // select particle kind to find intersection parameter
  switch ( particle->kind )
  {
  case KIND_ELLIPSOID :
    t = GetEdgeIntersectParamEllipsoid( particle->ellipsoid, Pin, Pout );
    break;
  
  default :
    error( "GetEdgeIntersectParam : bad particle kind\n" );
  }

  return t;
} //end GetEdgeIntersectParam()
/*--------------------------------------------------------------------
PURPOSE: To solve for the intersection point between particle and edge
     Equation of particle: (x-Xc)^2 + (y-Yc)^2 + (z-Zc)^2 = R^2
     Equation of plane: xi = ai + bi*t, i = 1,2,3.
     Take t = 0 for nodeIn, t = 1 for nodeOut ==> get ai and bi
     Solve the final Quadratic Equation for 0<t<1.
     From t ==> coordinates of new node (intersection point)
     
INPUT:
     d: pointer to datastructure
     nodeIn, nodeOut: The nodes at the ends of the edge
     t: t value for the intersection point
     GeoTab: geometric table for the new nodes
     iSp: current particle
     
OUTPUT:
     Updated GeoTabSp tables

By:  Veeramani
Update: 15/Jun/04
--------------------------------------------------------------------*/
static int
SetEdgeIntersectNode(
const mesh_t* mesh,
const int     nodeIn,
const int     nodeOut,
const double      t,
      double**    GeoTab )
{
  double a[ 4 ], b[ 4 ], XYZ[ 4 ];
  int dir;
  int nGeoNodes;

  GetNodePos( mesh, GeoTab, nodeIn , a );
  GetNodePos( mesh, GeoTab, nodeOut, b );
  
  for ( dir = 1 ; dir <= 3 ; dir++ )
    b[ dir ] -= a[ dir ];

  for ( dir = 1 ; dir <= 3 ; dir++ )
    XYZ[ dir ] = a[ dir ] + t * b[ dir ];

  // Store the intersection point in GeoTab
  GeoTab[ 0 ][ 0 ] += 1.0;
  nGeoNodes = (int) GeoTab[ 0 ][ 0 ];
  
  if ( ( nGeoNodes % MBLOCK ) == 1 )
    for ( dir = 1 ; dir <= 3 ; dir++ )
    {
      GeoTab[ dir ] = (double*) realloc ( GeoTab[ dir ], ( nGeoNodes + MBLOCK ) * sizeof(double) );
      assert_error( GeoTab[ dir ] != NULL, "Memory allocation error");
    }
  for ( dir = 1 ; dir <= 3 ; dir++ )
    GeoTab[ dir ][ nGeoNodes ] = XYZ[ dir ];

  return nGeoNodes + mesh->NbOfNodes;
} // end SetEdgeIntersectNode()
/*---------------------------------------------------------------------------
PURPOSE:
  To find the midpoint between two points by linear interpolation

INPUT:
  d: pointer to datastructure
  node1 and node2: the two edge points
  GeoTab: node table to record the new point
  iSp: current particle

OUTPUT:
  The node number of mid point created

By Veeramani
30/09/2004
----------------------------------------------------------------------------*/
static int
SetNodesMidpoint(
const mesh_t* mesh,
const int     node1,
const int     node2,
      double**    GeoTab )
{
  double xyz1[ 4 ], xyz2[ 4 ], xyzm[ 4 ];
  int dir, nGeoNodes;

  GetNodePos( mesh, GeoTab, node1, xyz1 );
  GetNodePos( mesh, GeoTab, node2, xyz2 );

  // Find the midpoint
  for ( dir = 1 ; dir <= 3 ; dir++ )
    xyzm[ dir ] = ( xyz1[ dir ] + xyz2[ dir ] ) / 2.0;

  // Create this new mid node
  GeoTab[ 0 ][ 0 ] += 1.;
  nGeoNodes = (int) GeoTab[ 0 ][ 0 ];
  
  if ( ( nGeoNodes % MBLOCK ) == 1 )
    for ( dir = 1 ; dir <= 3 ; dir++ )
    {
      GeoTab[ dir ] = (double*) realloc ( GeoTab[ dir ],
                                          ( nGeoNodes + MBLOCK ) * sizeof(double) );
      assert_error( GeoTab[ dir ] != NULL, "Memory allocation error");
    }

  for ( dir = 1 ; dir <= 3 ; dir++ )
    GeoTab[ dir ][ nGeoNodes ] = xyzm[ dir ];

  return nGeoNodes + mesh->NbOfNodes;
} // end SetNodesMidpoint()
/*--------------------------------------------------------------------
PURPOSE:
  To find the edges of the element cut and find the edge node
  We store the edge node in NodeSurf datastructure only once
  All elements with common edge will have same common edge
  node

INPUT:
  d: pointer to data structure
  particle: the particle
  ElemPartly: Elements that are partly inside particle
  GeoTab: node table for new nodes
NOTE: We have on common NodeSurf structure for all the particles

OUTPUT:
  Update NodeSurf
  NodeSurfInPar: All the nodes having connectivity to outside
  are stored here. Right now this is used only for
  reallocating NodeSurf at each timestep. Note, the number
  of surface nodes should not exceed some defined value.

By: Veeramani
Date: 02/10/06 Antoine
--------------------------------------------------------------------*/
void FillInOutSurfConTab(
  const mesh_t* mesh,
  const particle_t* particle,
  const int*   ElemPartly,
        double**  GeoTab,
        int*** NodeSurf,
        int*   NodeSurfInPar )
{
  const double TMIN = 0.01, // tolerance
              TMAX = 0.99; // tolerance
  int iElem, iNode, jNode, kNode, theElem, nSurfNodes,
        fNode, // First node
        sNode, // Second node
        eNode, // Edge node
        mNode, // Mid node
        **ElemCon = mesh->ConTabLocal;
  short flag;
  double tvalue;

  // loop over every element that are partly inside
  for ( iElem = 1 ; iElem <= ElemPartly[ 0 ] ; iElem++ )
  {
    theElem = ElemPartly[ iElem ]; // get the current element 
    
    // loop over all its P nodes, except the last one as we need another node
    for ( iNode = 1 ; iNode <= PRE_NODES_PER_EL - 1 ; iNode++ )
    {
      fNode = mesh->ConTab[ theElem ][ iNode ]; // get global first node
            
      if ( IsNodeInPar( mesh, particle, fNode ) ) // first node in the particle ?
      {// loop over remaining P nodes
        for ( jNode = iNode + 1 ; jNode <= PRE_NODES_PER_EL ; jNode++ )
        {
          sNode = mesh->ConTab[ theElem ][ jNode ]; // get global second node
                    
          if ( !IsNodeInPar( mesh, particle, sNode ) ) // seconde node not in the particle ?
          {
            flag = FALSE;
            
            // Check whether this edge has been processed before or not
            if ( NodeSurf[ fNode ][ 0 ][ 0 ] > 0 ) // There are some entries, check them
              for ( kNode = 1 ; kNode <= NodeSurf[ fNode ][ 0 ][ 0 ] ; kNode++ )
                if ( NodeSurf[ fNode ][ 0 ][ kNode ] == sNode )
                {
                  flag = TRUE; // already done
                  break;
                }

            if ( !flag )
            {
              // The edge is cut and un-processed, find the intersection point and edgeNode
              tvalue = GetEdgeIntersectParam( mesh, particle, fNode, sNode, GeoTab );
              
              if ( tvalue < TMIN )
              {
                eNode = fNode; // fNode is eNode too, edge outside
                mNode = 0;  // No mid node
              }
              else if ( tvalue > TMAX )
              {
                eNode = sNode;  // sNode is eNode too, edge inside
                mNode = mesh->ConTab[ theElem ][ ElemCon[ iNode ][ jNode ] ];
              }
              else
              {
                eNode = SetEdgeIntersectNode( mesh, fNode, sNode, tvalue, GeoTab );
                mNode = SetNodesMidpoint( mesh, fNode, eNode, GeoTab );
              }
              // Store this edge in NodeSurf
              NodeSurf[ fNode ][ 0 ][ 0 ]++;
              nSurfNodes = (int) NodeSurf[ fNode ][ 0 ][ 0 ];
              
              for ( kNode = 0 ; kNode <= 2 ; kNode++ )
              {
                NodeSurf[ fNode ][ kNode ] = (int*) realloc (
                                                                 NodeSurf[ fNode ][ kNode ], ( nSurfNodes + 1 ) * sizeof(int) );
                assert_error( NodeSurf[ fNode ][ kNode ] != NULL, "Memory allocation error");
              }
              
              NodeSurf[ fNode ][ 0 ][ nSurfNodes ] = sNode;
              NodeSurf[ fNode ][ 1 ][ nSurfNodes ] = eNode;
              NodeSurf[ fNode ][ 2 ][ nSurfNodes ] = mNode;
              // Store in node in NodeSurfInPar only once
              // We do this when nSurfNodes == 1 only
              // nSurfNodes greater than one means this in node has multiple edges cut
              if ( nSurfNodes == 1 )
              {
                NodeSurfInPar[ 0 ]++;
                
                assert( NodeSurfInPar[ 0 ] <= NODELMT );
                
                NodeSurfInPar[ NodeSurfInPar[ 0 ] ] = fNode;
              } // end if nSurfNodes == 1
            }
          } //end second node is out
          // else it is in-node, edge not cut
        } //end for jNode
      } //endif first node is in
      else
      {
        for ( jNode = iNode + 1 ; jNode <= PRE_NODES_PER_EL ; jNode++ )
        {
          sNode = mesh->ConTab[ theElem ][ jNode ];
          //Check whether the second node is inside the particle          
          if ( IsNodeInPar( mesh, particle, sNode ) )
          {
            flag = FALSE;
            // Check whether this edge has been processed before or not
            if ( NodeSurf[ sNode ][ 0 ][ 0 ] > 0 )   // There are some entries, check them
              for ( kNode = 1 ; kNode <= NodeSurf[ sNode ][ 0 ][ 0 ] ; kNode++ )
                if ( NodeSurf[ sNode ][ 0 ][ kNode ] == fNode )
                {
                  flag = TRUE;  // already done
                  break;
                }

            if ( !flag )
            {
              // The edge is cut and un-processed, find the intersection point and edgeNode
              tvalue = GetEdgeIntersectParam( mesh, particle, sNode, fNode, GeoTab );
              
              if ( tvalue < TMIN )
              {
                eNode = sNode;  // sNode is eNode too, edge outside
                mNode = 0;  // no mid node
              }
              else if ( tvalue > TMAX )
              {
                eNode = fNode;  // fNode is eNode too, edge inside
                mNode = mesh->ConTab[ theElem ][ ElemCon[ iNode ][ jNode ] ];
              }
              else
              {
                eNode = SetEdgeIntersectNode( mesh, sNode, fNode, tvalue, GeoTab );
                mNode = SetNodesMidpoint( mesh, sNode, eNode, GeoTab );
              }
              // Store this edge in NodeSurf
              NodeSurf[ sNode ][ 0 ][ 0 ]++;
              nSurfNodes = (int) NodeSurf[ sNode ][ 0 ][ 0 ];
              
              for ( kNode = 0 ; kNode <= 2 ; kNode++ )
              {
                NodeSurf[ sNode ][ kNode ] =
                (int*) realloc ( NodeSurf[ sNode ][ kNode ], ( nSurfNodes + 1 ) * sizeof(int) );
                assert_error( NodeSurf[ sNode ][ kNode ] != NULL, "Memory allocation error");
              }
              
              NodeSurf[ sNode ][ 0 ][ nSurfNodes ] = fNode;
              NodeSurf[ sNode ][ 1 ][ nSurfNodes ] = eNode;
              NodeSurf[ sNode ][ 2 ][ nSurfNodes ] = mNode;
              // Store in node in NodeSurfInPar only once
              // We do this when nSurfNodes == 1 only
              // nSurfNodes greater than one means this in node has multiple edges cut
              if ( nSurfNodes == 1 )
              {
                NodeSurfInPar[ 0 ]++;
                
                assert( NodeSurfInPar[ 0 ] <= NODELMT );
                
                NodeSurfInPar[ NodeSurfInPar[ 0 ] ] = sNode;
              } // end if nSurfNodes == 1
            }
          } //end second node is in
          // else it is out node, edge not cut
        } //end for jNode
      } //end else first node is out
    } //end for iNode
  } //end for iElem
} //end FillInOutSurfConTab()
/*----------------------------------------------------------------------------
PURPOSE:
  This routine gives which of the edges are cut for an element
  which is known to be partly inside particle and also returns number
  of edges cut

INPUT:
  d: pointer to data structure
  XYZ_Sp: center of mass of current particle
  element: the current element
  NodeSurf: surface node connectivity structure
  iSp: current particle

OUTPUT:
  inNodes: Nodes inside are set true,
           inNodes[0] gives number of pure in nodes
  onNodes: Nodes on the boundary are set true,
           onNodes[0] gives number of on nodes
  exNodes: Nodes outside are set true,
           exNodes[0] gives number of pure out nodes
  edgeCut: [4X2] matrix stores local number of [1] inNode and
           [2] outNode of the edge that is cut
  numEdgeCut: Number of edges of the element cut

By: Veeramani
Update: 21/Sept/04
----------------------------------------------------------------------------*/
static int elemEdgeCut(
  const mesh_t* mesh,
  const particle_t* particle,
  const int    element,
        int*** NodeSurf,
        int*   inNodes,
        int*   onNodes,
        int*   exNodes,
        int**  edgeCut )
{
  int iNode, jNode, kNode, theNode, isIn;
  short flag = FALSE;
  int iEdge = 0;
  //int inNodes[5] = {0,0,0,0,0};  // defined in calling function
  //int onNodes[5] = {0,0,0,0,0};
  //int exNodes[5] = {0,0,0,0,0};

  // Classify the nodes of the element into in, on and out node

  for ( iNode = 1 ; iNode <= PRE_NODES_PER_EL ; iNode++ )
  {
    theNode = mesh->ConTab[ element ][ iNode ];
    // Check if this is in node
    isIn = IsNodeInPar( mesh, particle, theNode );
    
    if ( isIn )
    {
      inNodes[ iNode ] = TRUE;
      // Now find an edge with this node belonging to this element
      for ( kNode = 1 ; kNode <= NodeSurf[ theNode ][ 0 ][ 0 ] ; kNode++ )
      {
        // Check if this node belongs to the element too
        for ( jNode = 1 ; jNode <= PRE_NODES_PER_EL ; jNode++ )
        {
          if ( ( jNode != iNode ) &&
               ( NodeSurf[ theNode ][ 0 ][ kNode ] == mesh->ConTab[ element ][ jNode ] ) )
          {
            flag = TRUE;
            break;
          } // endif
        } // end for
        // Check if this node was found in element
        // Since kNodes are all out nodes of some edge, jNode is out node
        if ( flag )
        {
          exNodes[ jNode ] = TRUE;
          //further check if either of in or out node are on node
          
          if ( theNode == NodeSurf[ theNode ][ 1 ][ kNode ] )  //in node is also on node
            onNodes[ iNode ] = TRUE;
          
          if ( NodeSurf[ theNode ][ 0 ][ kNode ] == NodeSurf[ theNode ][ 1 ][ kNode ] )  //out node is also on node
            onNodes[ jNode ] = TRUE;
        } // endif other node of edge found
        // Unset flag
        flag = FALSE;
      } // end for kNode
    } // end if in node
    else
      exNodes[ iNode ] = TRUE;
  } //end for iNode

  // Now we apply the following criterion to fill edgeCut and numEdgeCut
  // An edge is cut only between in and out node
  // Edges between in-on, in-in, out-on, out-out and on-on are not cut
  // in, on and out nodes are counted in inNodes[0], onNodes[0] and exNodes[0] res.
  // We are strictly storing only genuinely cut edges

  for ( iNode = 1 ; iNode <= PRE_NODES_PER_EL ; iNode++ )
  {
    if ( !onNodes[ iNode ] )
    {
      if ( inNodes[ iNode ] )
      { // It is pure in node
        inNodes[ 0 ] ++;
        theNode = mesh->ConTab[ element ][ iNode ];
        // Check it's edge partner, we want only pure out node
        for ( jNode = 1 ; jNode <= PRE_NODES_PER_EL ; jNode++ )
        {
          // jNode should not be same as iNode,
          // jNode should be an out node and
          // jNode should not be an on node.
          if ( jNode != iNode && exNodes[ jNode ] && !onNodes[ jNode ] )
          {
            // Find the edge node
            for ( kNode = 1 ; kNode <= NodeSurf[ theNode ][ 0 ][ 0 ] ; kNode++ )
            {
              if ( NodeSurf[ theNode ][ 0 ][ kNode ] == mesh->ConTab[ element ][ jNode ] )
              {
                // Store edge cut
                iEdge++;
                edgeCut[ iEdge ][ 0 ] = NodeSurf[ theNode ][ 1 ][ kNode ];
                edgeCut[ iEdge ][ 1 ] = iNode;
                edgeCut[ iEdge ][ 2 ] = jNode;
              } // endif
            } //end for kNode
          } // end if pure out node
        } // end for jNode
      } // end if pure in node
      else
        exNodes[ 0 ] ++;  // pure out node
    } // end if not an on node
    else
      onNodes[ 0 ] ++;    // It is an on node
  } // end for iNode
  
  return ( iEdge );
} // end elemEdgeCut()
/*-----------------------------------------------------------------------------
PURPOSE:
  To form a new edge between node1 and node2 and also creat a midpoint for this edge
  
INPUT:
   d: pointer to datastructure
  node1 and node2: the two edge points
  GeoTab: node table to record the new point
  iSp: current particle
  NodeSurf: surface nodes inside particle domain. It carries other info
  like edge node and mid node for the edge cut.
  NodeSurfInPar: log of nodes accessed in NodeSurf, for deallocating memory

NOTE:  node1 has to be the original node
    node2 can be edge/on node created before

OUTPUT:
  The node number of mid point created

By: Veeramani
30/09/2004
----------------------------------------------------------------------------*/
static int
creatEdge(
  const mesh_t *mesh,
  const int node1,
  const int node2,
        double **GeoTab,
        int ***NodeSurf,
        int *NodeSurfInPar )
{
  short flag;
  int k, midPoint;
  int nSurfNodes;

  // Check whether the edge midpoint already exist
  flag = FALSE;
  for ( k = 1; k <= NodeSurf[ node1 ][ 0 ][ 0 ]; k++ )
  {
    if ( NodeSurf[ node1 ][ 1 ][ k ] == node2 )
    {
      flag = TRUE;
      break;
    }
  }
  if ( flag )   // midpoint known
    midPoint = NodeSurf[ node1 ][ 2 ][ k ];
  else
  {    // midpoint not known
    midPoint = SetNodesMidpoint( mesh, node1, node2, GeoTab );

    // Store this edge in NodeSurf
    NodeSurf[ node1 ][ 0 ][ 0 ]++;
    nSurfNodes = (int) NodeSurf[ node1 ][ 0 ][ 0 ];
    
    for ( k = 0; k <= 2; k++ )
    {
      NodeSurf[ node1 ][ k ] =
      (int*) realloc ( NodeSurf[ node1 ][ k ], ( nSurfNodes + 1 ) * sizeof(int) );
      assert_error( NodeSurf[ node1 ][ k ] != NULL, "Memory allocation error");
    }
    
    NodeSurf[ node1 ][ 0 ][ nSurfNodes ] = 0;
    NodeSurf[ node1 ][ 1 ][ nSurfNodes ] = node2;
    NodeSurf[ node1 ][ 2 ][ nSurfNodes ] = midPoint;

    // Store in node in NodeSurfInPar only once
    // We do this when nSurfNodes == 1 only
    // nSurfNodes greater than one means this in node has multiple edges cut
    if ( nSurfNodes == 1 )
    {
      NodeSurfInPar[ 0 ]++;
      assert( NodeSurfInPar[ 0 ] <= NODELMT );
      NodeSurfInPar[ NodeSurfInPar[ 0 ] ] = node1;
    } // end if nSurfNodes == 1
  }
  return midPoint;
} // end creatEdge()
/*------------------------------------------------------------------------------
PURPOSE:
  To creat a new element from the four given corner nodes of a tetrahedron

INPUT:
  d: pointer to datastructure
  tetra: Vector of corner nodes from which new element is to be formed
  type: Gives type of tetrahedra to be formed
  theElem: the parent element from which a new element is carved
  ElemCon: local connectivity of nodes for the whole element
  SurfCon: local connectivity of nodes on the surface
  NodeSurf: surface nodes inside particle domain. It carries other info
         like edge node and mid node for the edge cut.

OUPUT:
  Updated
  GeoTab: Node Table for the particle
  conTSp: Connectivity Table for the particle
  surfTSp: Connectivity Table for the surface
  fullESp: list of refined smaller elements created

By:  Veeramani
Update: 10/Mar/04
------------------------------------------------------------------------------*/
static void
creatTetra(
  const mesh_t *mesh,
        int *tetra,
  const int type,
  const int theElem,
        int **ElemCon,
        int **SurfCon,
        int **conTSp,
        double **GeoTab,
        int **surfTSp,
        int **fullESp,
        int ***NodeSurf )
{
  short flag;
  int i, j, k;
  int tElems = mesh->NbOfElements; // Total number of original nodes and elements
  int nConElems; // dynamic number of elements in ConTabSp
  int nSrfElems; // dynamic number of surface elements in surfTSp
  int nChildren; // dynamic number of child elements created
  int iNode;

  int surfT[ 4 ] = {0, 2, 3, 4};     // Node sequence for surface element

  // We are going to creat a new element eventually,
  // so to store it we creat one more space in ConTabSp
  conTSp[ 0 ][ 0 ]++;
  nConElems = conTSp[ 0 ][ 0 ];
  if ( ( nConElems % MBLOCK ) == 1 )
    for ( iNode = 0; iNode <= NODES_PER_EL; iNode++ )
    {
      conTSp[ iNode ] = (int*) realloc ( conTSp[ iNode ], ( nConElems + MBLOCK ) * sizeof(int) );
      assert_error( conTSp[ iNode ] != NULL, "Memory allocation error");
    }

  // One surface of the new element will be external
  surfTSp[ 0 ][ 0 ]++;
  nSrfElems = surfTSp[ 0 ][ 0 ];
  if ( ( nSrfElems % MBLOCK ) == 1 )
    for ( iNode = 0; iNode <= NODES_PER_SURFACE; iNode++ )
    {
      surfTSp[ iNode ] = (int*) realloc ( surfTSp[ iNode ], ( nSrfElems + MBLOCK ) * sizeof(int) );
      assert_error( surfTSp[ iNode ] != NULL, "Memory allocation error");
    }

  // Also creat space in fullESp
  (*fullESp)[ 0 ]++;
  nChildren = (*fullESp)[ 0 ];
  if ( ( nChildren % MBLOCK ) == 1 )
  {
    *fullESp = ( int * ) realloc ( (*fullESp), ( nChildren + MBLOCK ) * sizeof(int) );
    assert_error( *fullESp != NULL, "Memory allocation error");
  }

  // We can fill in the corner nodes of the new element
  for ( i = 1; i <= PRE_NODES_PER_EL; i++ )
    conTSp[ i ][ nConElems ] = tetra[ i ];

  // There is one inNode (tetra[1]) and others are edge Nodes/ on Nodes
  // Fill the conTSp with mid nodes using NodeSurf structure
  for ( i = 2; i <= PRE_NODES_PER_EL; i++ )
  {
    // Check NodeSurf of tetra[1] for tetra[i] amongst edge nodes
    // Remember out == edge node if it is also on node
    // Similarly in == edge node if it is also on node
    // In both these cases, midnode can be obtained directly from parent element
    if ( tetra[ i ] <= mesh->NbOfNodes )
    {  // Original node
      // Find the local node number in parent element
      j = 0;
      k = 0;
      for ( iNode = 1; iNode <= PRE_NODES_PER_EL; iNode++ )
      {
        if ( tetra[ 1 ] == mesh->ConTab[ theElem ][ iNode ] )
          j = iNode;
        if ( tetra[ i ] == mesh->ConTab[ theElem ][ iNode ] )
          k = iNode;
      } // end for iNode
      conTSp[ ElemCon[ 1 ][ i ] ][ nConElems ] = mesh->ConTab[ theElem ][ ElemCon[ j ][ k ] ];
    }
    else
    {
      flag = FALSE;
      for ( k = 1; k <= NodeSurf[ tetra[ 1 ] ][ 0 ][ 0 ]; k++ )
      {
        if ( NodeSurf[ tetra[ 1 ] ][ 1 ][ k ] == tetra[ i ] )
        {
          flag = TRUE;
          break;
        }
      }

      if ( flag )    // Mid node found, just get it
        conTSp[ ElemCon[ 1 ][ i ] ][ nConElems ] = NodeSurf[ tetra[ 1 ] ][ 2 ][ k ];
      else
      {        // Mid node not found, create it
        warning( "\t >> Creating new node for edge 1-%d in TETRA-%d\n", i, type );
        conTSp[ ElemCon[ 1 ][ i ] ][ nConElems ] = SetNodesMidpoint( mesh, tetra[ 1 ], tetra[ i ], GeoTab );
      }
    }
  } // end for i

  // Fill the mid node between edge/on nodes
  conTSp[ ElemCon[ 2 ][ 3 ] ][ nConElems ] = SetNodesMidpoint( mesh, tetra[ 2 ], tetra[ 3 ], GeoTab );
  conTSp[ ElemCon[ 3 ][ 4 ] ][ nConElems ] = SetNodesMidpoint( mesh, tetra[ 3 ], tetra[ 4 ], GeoTab);
  conTSp[ ElemCon[ 4 ][ 2 ] ][ nConElems ] = SetNodesMidpoint( mesh, tetra[ 4 ], tetra[ 2 ], GeoTab);

  // Store the original element number
  conTSp[ 0 ][ nConElems ] = theElem;
  // Store the element number in fullESp
  (*fullESp)[ nChildren ] = nConElems + tElems;    // Element numbering of this child

  // Store the surface element
  surfTSp[ 0 ][ nSrfElems ] = theElem;
  // First three nodes [1,2,3] are corner nodes
  for ( i = 1; i <= 3; i++ )
    surfTSp[ i ][ nSrfElems ] = conTSp[ surfT[ i ] ][ nConElems ];
  // Next three nodes [4,5,6] are mid nodes
  for ( i = 1; i <= 2; i++ )
    for ( j = i + 1; j <= 3; j++ )
      surfTSp[ SurfCon[ i ][ j ] ][ nSrfElems ] = conTSp[ ElemCon[ surfT[ i ] ][ surfT[ j ] ] ][ nConElems ];
} //end creatTetra ()
/*------------------------------------------------------------------------------
PURPOSE:
  To creat three new elements from the six given corner nodes of a prism
  The standard prism solved here. All different conformations are mapped to this
  standard prism

INPUT:
  d: pointer to datastructure
  prismMap: vector of corner nodes from which new elements are to be formed
  type:  2 => prism with two nodes inside,
        3 => prism with three nodes inside.
  pElem: the parent element from which new elements are carved
  ElemCon: local connectivity of nodes for the whole element
  SurfCon: local connectivity of nodes on the surface
  NodeSurf: surface nodes inside particle domain. It carries other info
         like edge node and mid node for the edge cut.

OUPUT:
  Updated
  GeoTab: Node Table for the particle
  conTSp: Connectivity Table for the particle
  surfTSp: Connectivity Table for the surface
  fullESp: list of full elements inside particle

By:  Veeramani
Update: 4/May/04
------------------------------------------------------------------------------*/
static void
creatPrism(
  const mesh_t *mesh,
  const int *prismMap,
  const int type,
  const int pElem,
        int **ElemCon,
        int **SurfCon,
        int **conTSp,
        double **GeoTab,
        int **surfTSp,
        int **fullESp,
        int ***NodeSurf,
        int *NodeSurfInPar )
{

  short templat ;
  int i, j, k; // General counters
  int iTetra; // Count number of tetrahedra
  int count, scount; // Count number of nodes
  int tElems; // Total number of original nodes
  int iNode, jNode;
  int nConElems; // dynamic number of elements in ConTabSp
  int nSrfElems; // dynamic number of surface elements in surfTSp
  int nChildren;

  int surfT[ 2 ][ 4 ] = {{0, 0, 0, 0}, {0, 0, 0, 0}}; // Temporary node sequence for surface element

  int tetra[ 11 ] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; // Temporary space for storing nodes
  short pTemps[ 6 ][ 4 ][ 5 ] = {{{0, 0, 0, 0, 0}, {0, 2, 3, 1, 5}, {0, 3, 1, 5, 6}, {0, 1, 5, 6, 4}},
                                 {{0, 0, 0, 0, 0}, {0, 1, 2, 3, 6}, {0, 1, 2, 6, 5}, {0, 1, 5, 6, 4}},
                                 {{0, 0, 0, 0, 0}, {0, 1, 2, 3, 4}, {0, 2, 3, 4, 5}, {0, 3, 4, 5, 6}},
                                 {{0, 0, 0, 0, 0}, {0, 1, 2, 3, 4}, {0, 2, 3, 4, 6}, {0, 2, 6, 4, 5}},
                                 {{0, 0, 0, 0, 0}, {0, 1, 2, 3, 5}, {0, 1, 5, 3, 4}, {0, 3, 4, 5, 6}},
                                 {{0, 0, 0, 0, 0}, {0, 1, 2, 3, 6}, {0, 2, 6, 4, 5}, {0, 2, 6, 1, 4}}};

  short ptTable[ 4 ][ 5 ]; // Nodes for creating tetrahedra from prism

  int prism[7][7] = {{0,0,0,0,0,0,0},{0,0,0,0,0,0,0},{0,0,0,0,0,0,0},{0,0,0,0,0,0,0},{0,0,0,0,0,0,0},{0,0,0,0,0,0,0},{0,0,0,0,0,0,0}}; // Matrix to store the connectivity of prism
//  int **prism; // Matrix to store the connectivity of prism
//  prism = AllocM_int( 6, 6 ); // Allocate Memory for prism

  // Get nodes of the prism from prismMap
  for ( i = 1; i <= 6; i++ )
  {
    prism[ i ][ i ] = prismMap[ i ];  // Diagonal values of prism are thus corner nodes
//    printf("%u %u\n",prismMap[ i ],prism[ i ][ i ] );
  }

  tElems = mesh->NbOfElements;

  // We are going to decompose every prism into three tetrahedra,
  // So here we make sure space for three tetrahedra is available
  conTSp[ 0 ][ 0 ] += 3;
  nConElems = conTSp[ 0 ][ 0 ];
  count = ( nConElems % MBLOCK );
  if ( ( count >= 1 ) && ( count <= 3 ) )
    for ( iNode = 0; iNode <= NODES_PER_EL; iNode++ )
    {
      conTSp[ iNode ] = (int*) realloc( conTSp[ iNode ], ( nConElems + MBLOCK - count + 1 ) * sizeof(int) );
      assert_error( conTSp[ iNode ] != NULL, "Memory allocation error");
    }   

  // also creat space in childSp
  (*fullESp)[ 0 ] += 3;
  nChildren = (*fullESp)[ 0 ];
  count = ( nChildren % MBLOCK );
  
  if ( ( count >= 1 ) && ( count <= 3 ) )
  {
    *fullESp = (int*) realloc ( (*fullESp), ( nChildren + MBLOCK - count + 1 ) * sizeof(int) );
    assert_error( *fullESp != NULL, "Memory allocation error");
  }   

  if ( type == 2 )
  {  // One of the rectangular face of this prism is curved
    // THIS PRISM GIVES TWO SURFACE TRIANGLES
    // WE ALWAYS CONNECT 3-5
    surfT[ 0 ][ 1 ] = 2;
    surfT[ 0 ][ 2 ] = 3;
    surfT[ 0 ][ 3 ] = 5;
    surfT[ 1 ][ 1 ] = 3;
    surfT[ 1 ][ 2 ] = 6;
    surfT[ 1 ][ 3 ] = 5;

    // Allocate memory for two surface elements
    surfTSp[ 0 ][ 0 ] += 2;
    nSrfElems = surfTSp[ 0 ][ 0 ];
    scount = ( nSrfElems % MBLOCK );

    
    if ( ( scount >= 1 ) && ( scount <= 2 ) )
      for ( iNode = 0; iNode <= NODES_PER_SURFACE; iNode++ )
      {
        surfTSp[ iNode ] = (int*) realloc ( surfTSp[ iNode ],
                                                    ( nSrfElems + MBLOCK - scount + 1 ) * sizeof(int) );
        assert_error( surfTSp[ iNode ] != NULL, "Memory allocation error");
      }   
    
    // GET MID NODES OF PRE_DEFINED EDGES OF THE PRISM
    // THESE EDGES ARE: 1-4, 1-2, 1-3, 4-5, 4-6, 2-3, 3-6, 6-5 and 5-2
    // 1-4
    // To get the mid node from original element we need to trace back the indices
    i = 0;
    j = 0;
    for ( iNode = 1; iNode <= PRE_NODES_PER_EL; iNode++ )
    {
      if ( prism[ 1 ][ 1 ] == mesh->ConTab[ pElem ][ iNode ] )
        i = iNode;
      if ( prism[ 4 ][ 4 ] == mesh->ConTab[ pElem ][ iNode ] )
        j = iNode;
    } // end for iNode
    prism[ 1 ][ 4 ] = prism[ 4 ][ 1 ] = mesh->ConTab[ pElem ][ ElemCon[ i ][ j ] ];

//  for ( i = 0; i <= 6; i++ )
//  for ( j = 0; j <= 6; j++ )
//    printf("%u ",prism[ i ][ j ]);
//    printf("\n");

    // Get the midpoints of 1-2, 1-3, 4-5 and 4-6
    // Search for edge nodes 2 and 3 from NodeSurf
    for ( i = 1; i <= NodeSurf[ prism[ 1 ][ 1 ] ][ 0 ][ 0 ]; i++ )
    {
      if ( prism[ 2 ][ 2 ] == NodeSurf[ prism[ 1 ][ 1 ] ][ 1 ][ i ] )
        prism[ 1 ][ 2 ] = prism[ 2 ][ 1 ] = NodeSurf[ prism[ 1 ][ 1 ] ][ 2 ][ i ];
      if ( prism[ 3 ][ 3 ] == NodeSurf[ prism[ 1 ][ 1 ] ][ 1 ][ i ] )
        prism[ 1 ][ 3 ] = prism[ 3 ][ 1 ] = NodeSurf[ prism[ 1 ][ 1 ] ][ 2 ][ i ];
    }

//  for ( i = 0; i <= 6; i++ )
//  for ( j = 0; j <= 6; j++ )
//    printf("%u ",prism[ i ][ j ]);
//    printf("\n");
//    printf("%u\n",prism[ 4 ][ 4 ]);
//    printf("%u\n",NodeSurf[ prism[ 4 ][ 4 ] ][ 0 ][ 0 ]);
    
    for ( i = 1; i <= NodeSurf[ prism[ 4 ][ 4 ] ][ 0 ][ 0 ]; i++ )
    {
      if ( prism[ 5 ][ 5 ] == NodeSurf[ prism[ 4 ][ 4 ] ][ 1 ][ i ] )
        prism[ 4 ][ 5 ] = prism[ 5 ][ 4 ] = NodeSurf[ prism[ 4 ][ 4 ] ][ 2 ][ i ];
      if ( prism[ 6 ][ 6 ] == NodeSurf[ prism[ 4 ][ 4 ] ][ 1 ][ i ] )
        prism[ 4 ][ 6 ] = prism[ 6 ][ 4 ] = NodeSurf[ prism[ 4 ][ 4 ] ][ 2 ][ i ];
    }
    // Get the mid nodes of 2-3, 3-6, 6-5 and 5-2
    // Use SetNodesMidpoint routine with adjustment option
    prism[ 2 ][ 3 ] = prism[ 3 ][ 2 ] = SetNodesMidpoint( mesh, prism[ 2 ][ 2 ], prism[ 3 ][ 3 ], GeoTab );
    prism[ 3 ][ 6 ] = prism[ 6 ][ 3 ] = SetNodesMidpoint( mesh, prism[ 3 ][ 3 ], prism[ 6 ][ 6 ], GeoTab );
    prism[ 5 ][ 6 ] = prism[ 6 ][ 5 ] = SetNodesMidpoint( mesh, prism[ 5 ][ 5 ], prism[ 6 ][ 6 ], GeoTab );
    prism[ 2 ][ 5 ] = prism[ 5 ][ 2 ] = SetNodesMidpoint( mesh, prism[ 2 ][ 2 ], prism[ 5 ][ 5 ], GeoTab );

    // TO FORM THE NEW EDGES WE NEED TO DECIDE WHICH NODES TO CONNECT
    // Find which nodes to connect
    // On the face of 1-4-2-5, decide between 1-5 and 2-4
    // On the face of 1-4-3-6, decide between 1-6 and 3-4
    // We make a connect for the lower of the inner nodes
    if ( prism[ 1 ][ 1 ] < prism[ 4 ][ 4 ] )
    {
      // We connect 1-5 and 1-6
      for ( j = 5; j <= 6; j++ )
      {
        prism[ 1 ][ j ] = prism[ j ][ 1 ] = creatEdge( mesh, prism[ 1 ][ 1 ], prism[ j ][ j ], GeoTab, NodeSurf, NodeSurfInPar );
      } //end j
      prism[ 3 ][ 5 ] = prism[ 5 ][ 3 ] = SetNodesMidpoint( mesh, prism[ 3 ][ 3 ], prism[ 5 ][ 5 ], GeoTab );
    } // end if
    else
    {
      // We connect 4-2 and 4-3
      for ( j = 2; j <= 3; j++ )
      {
        prism[ 4 ][ j ] = prism[ j ][ 4 ] = creatEdge( mesh, prism[ 4 ][ 4 ], prism[ j ][ j ], GeoTab, NodeSurf, NodeSurfInPar );
      } // end j
      prism[ 3 ][ 5 ] = prism[ 5 ][ 3 ] = SetNodesMidpoint( mesh, prism[ 3 ][ 3 ], prism[ 5 ][ 5 ], GeoTab );
    } // end else

    // CHOOSE THE TEMPLATE
    if ( prism[ 1 ][ 1 ] < prism[ 4 ][ 4 ] )
      templat = 0;
    else templat = 2;
    for ( i = 1; i <= 3; i++ )
      for ( j = 1; j <= 4; j++ )
        ptTable[ i ][ j ] = pTemps[ templat ][ i ][ j ];

    // All the edges have midpoints now, use ptTable and creat tetrahedra
    for ( iTetra = 1; iTetra <= 3;  iTetra++ )
    {
      // Get the pressure nodes
      for ( iNode = 1; iNode <= PRE_NODES_PER_EL; iNode++ )
        tetra[ iNode ] = prism[ ptTable[ iTetra ][ iNode ] ][ ptTable[ iTetra ][ iNode ] ];
      // Get the midpoints of edges
      for ( iNode = 1; iNode <= PRE_NODES_PER_EL - 1; iNode++ )
        for( jNode = iNode + 1; jNode <= PRE_NODES_PER_EL; jNode++ )
          tetra[ ElemCon[ iNode ][ jNode ] ] = prism[ ptTable[ iTetra ][ iNode ] ][ ptTable[ iTetra ][ jNode ] ];
      for ( i = 1; i <= NODES_PER_EL; i++ )
        conTSp[ i ][ nConElems - iTetra + 1 ] = tetra[ i ];
      conTSp[ 0 ][ nConElems - iTetra + 1 ] = pElem;
      (*fullESp)[ nChildren - iTetra + 1 ] = nConElems - iTetra + 1 + tElems;
    } // end for iTetra

    // Store the surface elements
    for ( iTetra = 0; iTetra <= 1; iTetra++ )
    {
      // Store parent element
      surfTSp
      [ 0 ][ nSrfElems - iTetra ] = pElem;
      // Get the corner nodes
      for ( i = 1; i <= 3; i++ )
        surfTSp[ i ][ nSrfElems - iTetra ] = prism[ surfT[ iTetra ][ i ] ][ surfT[ iTetra ][ i ] ];
      // Get the midpoints
      for ( i = 1; i <= 2; i++ )
        for ( j = i + 1; j <= 3; j++ )
          surfTSp[ SurfCon[ i ][ j ] ][ nSrfElems - iTetra ] = prism[ surfT[ iTetra ][ i ] ][ surfT[ iTetra ][ j ] ];
    }
  } // end if (type == 2)
  else
  {   // One of the triangular face of prism is curved
    // This prism will give one surface element [4,5,6]
    surfT[ 0 ][ 1 ] = 4;
    surfT[ 0 ][ 2 ] = 5;
    surfT[ 0 ][ 3 ] = 6;
    // Allocate memory for this surface element
    surfTSp[ 0 ][ 0 ]++;
    nSrfElems = surfTSp[ 0 ][ 0 ];
    
    if ( ( nSrfElems % MBLOCK ) == 1 )
      for ( iNode = 0; iNode <= NODES_PER_SURFACE; iNode++ )
      {
        surfTSp[ iNode ] = (int*) realloc ( surfTSp[ iNode ],
                                                ( nSrfElems + MBLOCK ) * sizeof(int) );
        assert_error( surfTSp[ iNode ] != NULL, "Memory allocation error");
      }   
    
    // GET MID NODES OF PRE_DEFINED EDGES OF THE PRISM
    // THESE EDGES ARE: 1-2, 2-3, 3-1, 1-4, 2-5, 3-6, 4-5, 5-6 AND 6-4
    // Get the midpoints of the edges connecting original nodes
    // 1-2, 1-3 and 2-3, get index of 1 2 and 3 as i, j, k
    i = 0;
    j = 0;
    k = 0;
    for ( iNode = 1; iNode <= PRE_NODES_PER_EL; iNode++ )
    {
      if ( prism[ 1 ][ 1 ] == mesh->ConTab[ pElem ][ iNode ] )
        i = iNode;
      if ( prism[ 2 ][ 2 ] == mesh->ConTab[ pElem ][ iNode ] )
        j = iNode;
      if ( prism[ 3 ][ 3 ] == mesh->ConTab[ pElem ][ iNode ] )
        k = iNode;
    } // end for iNode
    prism[ 1 ][ 2 ] = prism[ 2 ][ 1 ] = mesh->ConTab[ pElem ][ ElemCon[ i ][ j ] ];
    prism[ 1 ][ 3 ] = prism[ 3 ][ 1 ] = mesh->ConTab[ pElem ][ ElemCon[ i ][ k ] ];
    prism[ 2 ][ 3 ] = prism[ 3 ][ 2 ] = mesh->ConTab[ pElem ][ ElemCon[ j ][ k ] ];
    // Get the midpoints of the edges connecting original and new nodes
    // 1-4, 2-5 and 3-6, from NodeSurf
    for ( j = 1; j <= 3; j++ )
    {
      for ( i = 1; i <= NodeSurf[ prism[ j ][ j ] ][ 0 ][ 0 ]; i++ )
      {
        if ( prism[ j + 3 ][ j + 3 ] == NodeSurf[ prism[ j ][ j ] ][ 1 ][ i ] )
          prism[ j ][ j + 3 ] = prism[ j + 3 ][ j ] = NodeSurf[ prism[ j ][ j ] ][ 2 ][ i ];
      }
    }
    // Get the midpoints of new edges
    // 4-5, 5-6, 6-4, by SetNodesMidpoint routine
    prism[ 4 ][ 5 ] = prism[ 5 ][ 4 ] = SetNodesMidpoint( mesh, prism[ 4 ][ 4 ], prism[ 5 ][ 5 ], GeoTab );
    prism[ 5 ][ 6 ] = prism[ 6 ][ 5 ] = SetNodesMidpoint( mesh, prism[ 5 ][ 5 ], prism[ 6 ][ 6 ], GeoTab );
    prism[ 4 ][ 6 ] = prism[ 6 ][ 4 ] = SetNodesMidpoint( mesh, prism[ 6 ][ 6 ], prism[ 4 ][ 4 ], GeoTab );

    // TO FORM THE NEW EDGES WE NEED TO DECIDE WHICH NODES TO CONNECT
    // Find which nodes to connect
    // We connect for the lowest of the inner nodes 1/2/3
    if ( prism[ 1 ][ 1 ] < prism[ 2 ][ 2 ] && prism[ 1 ][ 1 ] < prism[ 3 ][ 3 ] )
    {
      // We connect 1-5, 1-6
      prism[ 1 ][ 5 ] = prism[ 5 ][ 1 ] = creatEdge( mesh, prism[ 1 ][ 1 ], prism[ 5 ][ 5 ], GeoTab, NodeSurf, NodeSurfInPar );
      prism[ 1 ][ 6 ] = prism[ 6 ][ 1 ] = creatEdge( mesh, prism[ 1 ][ 1 ], prism[ 6 ][ 6 ], GeoTab, NodeSurf, NodeSurfInPar );
      if ( prism[ 2 ][ 2 ] < prism[ 3 ][ 3 ] )
      {
        // We connect 2-6
        prism[ 2 ][ 6 ] = prism[ 6 ][ 2 ] = creatEdge( mesh, prism[ 2 ][ 2 ], prism[ 6 ][ 6 ], GeoTab, NodeSurf, NodeSurfInPar);
        templat = 1;
      }
      else{
        // We connect 3-5
        prism[ 3 ]
        [ 5 ] = prism[ 5 ][ 3 ] = creatEdge( mesh, prism[ 3 ][ 3 ], prism[ 5 ][ 5 ], GeoTab, NodeSurf, NodeSurfInPar);
        templat = 0;
      }
    }
    else if ( prism[ 2 ][ 2 ] < prism[ 1 ][ 1 ] && prism[ 2 ][ 2 ] < prism[ 3 ][ 3 ] )
    {
      // We connect 2-4 and 2-6
      prism[ 2 ]
      [ 4 ] = prism[ 4 ][ 2 ] = creatEdge( mesh, prism[ 2 ][ 2 ], prism[ 4 ][ 4 ], GeoTab, NodeSurf, NodeSurfInPar );
      prism[ 2 ][ 6 ] = prism[ 6 ][ 2 ] = creatEdge( mesh, prism[ 2 ][ 2 ], prism[ 6 ][ 6 ], GeoTab, NodeSurf, NodeSurfInPar );
      if ( prism[ 1 ][ 1 ] < prism[ 3 ][ 3 ] )
      {
        // We connect 1-6
        prism[ 1 ][ 6 ] = prism[ 6 ][ 1 ] = creatEdge( mesh, prism[ 1 ][ 1 ], prism[ 6 ][ 6 ], GeoTab, NodeSurf, NodeSurfInPar );
        templat = 5;
      }
      else{
        // We connect 3-4
        prism[ 3 ]
        [ 4 ] = prism[ 4 ][ 3 ] = creatEdge( mesh, prism[ 3 ][ 3 ], prism[ 4 ][ 4 ], GeoTab, NodeSurf, NodeSurfInPar );
        templat = 3;
      }
    }
    else{
      // We connect 3-4 and 3-5
      prism[ 3 ][ 4 ] = prism[ 4 ][ 3 ] = creatEdge( mesh, prism[ 3 ][ 3 ], prism[ 4 ][ 4 ], GeoTab, NodeSurf, NodeSurfInPar );
      prism[ 3 ][ 5 ] = prism[ 5 ][ 3 ] = creatEdge( mesh, prism[ 3 ][ 3 ], prism[ 5 ][ 5 ], GeoTab, NodeSurf, NodeSurfInPar );
      if ( prism[ 1 ][ 1 ] < prism[ 2 ][ 2 ] )
      {
        // We connect 1-5
        prism[ 1 ]
        [ 5 ] = prism[ 5 ][ 1 ] = creatEdge( mesh, prism[ 1 ][ 1 ], prism[ 5 ][ 5 ], GeoTab, NodeSurf, NodeSurfInPar );
        templat = 4;
      }
      else{
        // We connect 2-4
        prism[ 2 ][ 4 ] = prism[ 4 ][ 2 ] = creatEdge( mesh, prism[ 2 ][ 2 ], prism[ 4 ][ 4 ], GeoTab, NodeSurf, NodeSurfInPar );
        templat = 2;
      }
    }

    // LOAD THE TEMPLATE
    for ( i = 1; i <= 3; i++ )
      for ( j = 1; j <= 4; j++ )
        ptTable[ i ][ j ] = pTemps[ templat ][ i ][ j ];

    // All the edges have midpoints now, use ptTable and creat tetrahedra
    for ( iTetra = 1; iTetra <= 3;  iTetra++ )
    {
      // Get the pressure nodes
      for ( iNode = 1; iNode <= PRE_NODES_PER_EL; iNode++ )
        tetra[ iNode ] = prism[ ptTable[ iTetra ][ iNode ] ][ ptTable[ iTetra ][ iNode ] ];
      // Get the midpoints of edges
      for ( iNode = 1; iNode <= PRE_NODES_PER_EL - 1; iNode++ )
        for( jNode = iNode + 1; jNode <= PRE_NODES_PER_EL; jNode++ )
          tetra[ ElemCon[ iNode ][ jNode ] ] = prism[ ptTable[ iTetra ][ iNode ] ][ ptTable[ iTetra ][ jNode ] ];
      for ( i = 1; i <= NODES_PER_EL; i++ )
        conTSp[ i ][ nConElems - iTetra + 1 ] = tetra[ i ];
      conTSp[ 0 ][ nConElems - iTetra + 1 ] = pElem;
      (*fullESp)[ nChildren - iTetra + 1 ] = nConElems - iTetra + 1 + tElems;
    } // end for iTetra
    // Store the surface elements
    surfTSp[ 0 ][ nSrfElems ] = pElem;
    // Get the corner nodes
    for ( i = 1; i <= 3; i++ )
      surfTSp[ i ][ nSrfElems ] = prism[ surfT[ 0 ][ i ] ][ surfT[ 0 ][ i ] ];
    // Get the midpoints
    for ( i = 1; i <= 2; i++ )
      for ( j = i + 1; j <= 3; j++ )
        surfTSp[ SurfCon[ i ][ j ] ][ nSrfElems ] = prism[ surfT[ 0 ][ i ] ][ surfT[ 0 ][ j ] ];
  } // end else (type == 3)
//  FreeM( prism );
} // end creatPrism()
/*------------------------------------------------------------------------------
PURPOSE:
  To creat two new elements from the five given corner nodes of a pyramid
  A standard pyramid is solved here. All conformation are mapped to this standard

INPUT:
  d: pointer to datastructure
  pyramMap: Vector of corner nodes from which new elements are to be formed
  pElem: the original element from which new elements are carved
  ElemCon: local connectivity of nodes for the whole element
  NodeSurf: surface nodes inside particle domain. It carries other info
  like edge node and mid node for the edge cut.

OUPUT:
  Updated
  GeoTab: Node Table for the particle
  conTSp: Connectivity Table for the particle
  fullESp: list of full elements

By:  Veeramani
Update: 31/May/04
------------------------------------------------------------------------------*/
static void
creatPyram(
  const mesh_t *mesh,
        int *pyramMap,
  const int pElem,
        int **ElemCon,
        int **SurfCon,
        int **conTSp,
        double **GeoTab,
        int **surfTSp,
        int **fullESp,
        int ***NodeSurf,
        int *NodeSurfInPar )
{

  short templat ;
  int i, j, k;
  int iTetra, count;
  int tElems; // Total number of original nodes and elements
  int iNode, jNode;
  int nConElems; // dynamic number of elements in ConTabSp
  int nSrfElems;    // dynamic number of surface elements in surfTSp
  int nChildren;

  int surfT[ 4 ] = {0, 3, 4, 5};  // Temporary sequence of surface nodes, there will be only one external surface [3,4,5]
  int tetra[ 11 ] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; // Temporary space for storing nodes
  short pTemps[ 2 ][ 3 ][ 5 ] = {{{0, 0, 0, 0, 0}, {0, 1, 2, 3, 4}, {0, 2, 3, 4, 5}},
                                 {{0, 0, 0, 0, 0}, {0, 1, 2, 3, 5}, {0, 1, 5, 3, 4}}};

  short ptTable[ 3 ][ 5 ] = { {0, 0, 0, 0, 0},
                              {0, 1, 2, 3, 4},
                              {0, 2, 3, 4, 5} }; // Creating tetrahedra from prism

//  int **pyram; // Matrix to store the connectivity of pyramid
//  pyram = AllocM_int( 5, 5 ); // Allocate memory for pyramid
  int pyram[6][6] = {{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}}; // Matrix to store the connectivity of pyramid

  // Get the corner nodes of the pyramid
  for ( i = 1; i <= 5; i++ )
    pyram[ i ][ i ] = pyramMap[ i ];

  tElems = mesh->NbOfElements;

  // We are going to decompose every pyramid into two tetrahedra,
  // So here we allocate space for two tetrahedra in one go
  conTSp[ 0 ][ 0 ] += 2;
  nConElems = conTSp[ 0 ][ 0 ];
  count = ( nConElems % MBLOCK );
  if ( ( count >= 1 ) && ( count <= 2 ) )
    for ( iNode = 0; iNode <= NODES_PER_EL; iNode++ )
  {
    conTSp[ iNode ] = (int*) realloc( conTSp[ iNode ], ( nConElems + MBLOCK - count + 1 ) * sizeof(int) );
    assert_error( conTSp[ iNode ] != NULL, "Memory allocation error");
  }   
  
  // Also allocate space for children
  (*fullESp)[ 0 ] += 2;
  nChildren = (*fullESp)[ 0 ];
  count = ( nChildren % MBLOCK );
  if ( ( count >= 1 ) && ( count <= 2 ) )
  {
    (*fullESp) = (int*) realloc ( (*fullESp), ( nChildren + MBLOCK - count + 1 ) * sizeof(int) );
    assert_error( (*fullESp) != NULL, "Memory allocation error");
  }   
      
  // Allocate memory for this surface element
  surfTSp[ 0 ][ 0 ]++;
  nSrfElems = surfTSp[ 0 ][ 0 ];
  if ( ( nSrfElems % MBLOCK ) == 1 )
    for ( iNode = 0; iNode <= NODES_PER_SURFACE; iNode++ )
      {
        surfTSp[ iNode ] = (int*) realloc ( surfTSp[ iNode ], ( nSrfElems + MBLOCK ) * sizeof(int) );
        assert_error( surfTSp[ iNode ] != NULL, "Memory allocation error");
      }   
        
  // GET THE MID NODES OF PRE-DETERMINED EDGES OF THE PYRAMID
  // THESE EDGES ARE: 1-2, 2-3, 3-1, 1-4, 2-5, 3-4, 4-5 AND 5-3
  // Get midpoints between 1,2 and 3 original nodes
  // 1-2, 1-3 and 2-3, get index of 1 2 and 3 as i, j, k
  i = 0;
  j = 0;
  k = 0;
  for ( iNode = 1; iNode <= PRE_NODES_PER_EL; iNode++ )
  {
    if ( pyram[ 1 ][ 1 ] == mesh->ConTab[ pElem ][ iNode ] )
      i = iNode;
    if ( pyram[ 2 ][ 2 ] == mesh->ConTab[ pElem ][ iNode ] )
      j = iNode;
    if ( pyram[ 3 ][ 3 ] == mesh->ConTab[ pElem ][ iNode ] )
      k = iNode;
  } // end for iNode
  pyram[ 1 ][ 2 ] = pyram[ 2 ][ 1 ] = mesh->ConTab[ pElem ][ ElemCon[ i ][ j ] ];
  pyram[ 1 ][ 3 ] = pyram[ 3 ][ 1 ] = mesh->ConTab[ pElem ][ ElemCon[ i ][ k ] ];
  pyram[ 2 ][ 3 ] = pyram[ 3 ][ 2 ] = mesh->ConTab[ pElem ][ ElemCon[ j ][ k ] ];
  // Get midpoints between edges 2-5 and 1-4 from NodeSurf
  for ( i = 1; i <= NodeSurf[ pyram[ 1 ][ 1 ] ][ 0 ][ 0 ]; i++ )
  {
    if ( pyram[ 4 ][ 4 ] == NodeSurf[ pyram[ 1 ][ 1 ] ][ 1 ][ i ] )
    {
      pyram[ 1 ][ 4 ] = pyram[ 4 ][ 1 ] = NodeSurf[ pyram[ 1 ][ 1 ] ][ 2 ][ i ];
      break;
    }
  }
  for ( i = 1; i <= NodeSurf[ pyram[ 2 ][ 2 ] ][ 0 ][ 0 ]; i++ )
  {
    if ( pyram[ 5 ][ 5 ] == NodeSurf[ pyram[ 2 ][ 2 ] ][ 1 ][ i ] )
    {
      pyram[ 2 ][ 5 ] = pyram[ 5 ][ 2 ] = NodeSurf[ pyram[ 2 ][ 2 ] ][ 2 ][ i ];
      break;
    }
  }
  // Get midpoints of the new edges
  // 3-4, 4-5, and 5-3 by SetNodesMidpoint routine
  pyram[ 3 ][ 4 ] = pyram[ 4 ][ 3 ] = SetNodesMidpoint( mesh, pyram[ 3 ][ 3 ], pyram[ 4 ][ 4 ], GeoTab );
  pyram[ 5 ][ 4 ] = pyram[ 4 ][ 5 ] = SetNodesMidpoint( mesh, pyram[ 5 ][ 5 ], pyram[ 4 ][ 4 ], GeoTab );
  pyram[ 3 ][ 5 ] = pyram[ 5 ][ 3 ] = SetNodesMidpoint( mesh, pyram[ 3 ][ 3 ], pyram[ 5 ][ 5 ], GeoTab );

  // TO FORM THE NEW EDGES WE NEED TO DECIDE WHICH NODES TO CONNECT
  // Find which nodes to connect
  // On the face of 1-4-2-5, decide between 1-5 and 2-4
  if ( pyram[ 1 ][ 1 ] < pyram[ 2 ][ 2 ] )
  {
    // We connect 1-5
    pyram[ 1 ][ 5 ] = pyram[ 5 ][ 1 ] = creatEdge( mesh, pyram[ 1 ][ 1 ], pyram[ 5 ][ 5 ], GeoTab, NodeSurf, NodeSurfInPar );
    templat = 1;
  }
  else{
    // We connect 2-4
    pyram[ 2 ][ 4 ] = pyram[ 4 ][ 2 ] = creatEdge( mesh, pyram[ 2 ][ 2 ], pyram[ 4 ][ 4 ], GeoTab, NodeSurf, NodeSurfInPar );
    templat = 0;
  }

  // LOAD ptTable WITH THE APPROPRIATE TEMPLATE
  for ( i = 1; i <= 2; i++ )
  for ( j = 1; j <= 4; j++ )
    ptTable[ i ][ j ] = pTemps[ templat ][ i ][ j ];

  // All the edges have midpoints now, use ptTable and creat tetrahedra
  for ( iTetra = 1; iTetra <= 2;  iTetra++ )
  {
    // Get the pressure nodes
    for ( iNode = 1; iNode <= PRE_NODES_PER_EL; iNode++ )
      tetra[ iNode ] = pyram[ ptTable[ iTetra ][ iNode ] ][ ptTable[ iTetra ][ iNode ] ];
    // Get the midpoints of edges
    for ( iNode = 1; iNode <= PRE_NODES_PER_EL - 1; iNode++ )
      for ( jNode = iNode + 1; jNode <= PRE_NODES_PER_EL; jNode++ )
        tetra[ ElemCon[ iNode ][ jNode ] ] = pyram[ ptTable[ iTetra ][ iNode ] ][ ptTable[ iTetra ][ jNode ] ];
    for ( i = 1; i <= NODES_PER_EL; i++ )
      conTSp[ i ][ nConElems - iTetra + 1 ] = tetra[ i ];
    conTSp[ 0 ][ nConElems - iTetra + 1 ] = pElem;
    (*fullESp)[ nChildren - iTetra + 1 ] = nConElems - iTetra + 1 + tElems;
  } // end for iTetra

  // Store the surface element
  surfTSp[ 0 ][ nSrfElems ] = pElem;
  // Get the corner nodes
  for ( i = 1; i <= 3; i++ )
    surfTSp[ i ][ nSrfElems ] = pyram[ surfT[ i ] ][ surfT[ i ] ];
  // Get the midpoints
  for ( i = 1; i <= 2; i++ )
    for ( j = i + 1; j <= 3; j++ )
      surfTSp[ SurfCon[ i ][ j ] ][ nSrfElems ] = pyram[ surfT[ i ] ][ surfT[ j ] ];

//  FreeM( pyram );
} // end for creatPyram()
/*---------------------------------------------------------------------
PURPOSE:
  To update the connectivity and Node Table with new elements
  created from partial elements list for a given particle
  NOTE: After this update, fullESp will have all the elements:

INPUT:
  d - pointer to data structures
  XYZ_Sp - center of mass of the current particle
  ElemCon - matrix with local connectivity between nodes of a element
  SurfCon - matrix with local connectivity of surface element
  fullESp[iSp] - list of elements that are fully inside particle iSp
  ElemPartly[iSp] - list of elements that are partly inside particle iSp
  NOTE: Not all elements in fullESp are original, and all the
   elements in ElemPartly are new elements refered to conTSp

OUTPUT:
  fullESp- Updated with new surface elements
  conTSp - connectivity table of new elements created
  GeoTab - Node Table of new nodes created
  surfTSp- connectivity table of new surface elements

By: Veeramani
Update: 08/Apr/04
   23/Jun/04 (Added surface elements)
--------------------------------------------------------------------*/

void
ParConnGeo(
const mesh_t* mesh,
const particle_t* particle,
      double** GeoTab,
      int** conTSp,
      int** fullESp,
      int* ElemPartly,
      int*** NodeSurf,
      int*   NodeSurfInPar )
{
  int iElem, theElem, iEdge,
        i, j, // Temporary counters
        in, out, on, fin, sin, // Temporary placeholders
        numEdgeCut,   // Number of edges cut
        numNodesInOn, // Number of nodes inside or on surface
        nFulElems,    // number of full elements
        **ElemCon = mesh->ConTabLocal,
        **SurfCon,  // Local connectivity of surface element [4][4]

        inNodes[ 5 ], // Nodes found inside particle
        onNodes[ 5 ], // Nodes found on surface of particle
        exNodes[ 5 ], // Nodes found outside particle
        edNodes[ 5 ], // Nodes found on edges
        tetra[ 5 ],   // Map sequence for tetrahedron
        pyram[ 6 ],   // Map sequence for pyramid
        prism[ 7 ],   // Map sequence for prism

        eLog[ 9 ] = {0, 0, 0, 0, 0, 0, 0, 0, 0},
        face[ 5 ][ 3 ] = { {0, 0, 0},
                           {2, 4, 3},   // Face opposite 1
                           {1, 3, 4},   // Face opposite 2
                           {1, 4, 2},   // Face opposite 3
                           {1, 2, 3} },  // Face opposite 4

        **edgeCut = AllocM_int( 4, 2 ), // Matrix to store the edges cut
        **edge = AllocM_int( 4, 4 ), // For mapping of edge nodes for prism2
        **surfTSp = NULL; //connectivity of external surface of elements created

  
  AllocMint(4,4,SurfCon);
  
  SurfCon[ 2 ][ 1 ] = SurfCon[ 1 ][ 2 ] = 4;
  SurfCon[ 3 ][ 1 ] = SurfCon[ 1 ][ 3 ] = 6;
  SurfCon[ 2 ][ 3 ] = SurfCon[ 3 ][ 2 ] = 5;
  
  
  surfTSp = (int**) calloc( 7, sizeof(int*) );
  assert_error( surfTSp != NULL, "Memory allocation error");
  
  for(int k = 0 ; k < 7 ; k++)
  {
    surfTSp[ k ] = (int*) calloc( 1, sizeof(int) );
    assert_error( surfTSp[ k ] != NULL, "Memory allocation error");
  }   
      
  for ( iElem = 1 ; iElem <= ElemPartly[ 0 ] ; iElem++ )
  {
    theElem = ElemPartly[ iElem ];
    //printf(">> Processing: Partial Element [%d] = %d\n",iElem, theElem);
    // Reset all structures
    for ( i = 0 ; i <= 4 ; i++ )
      inNodes[ i ] = 0;
    for ( i = 0 ; i <= 4 ; i++ )
      onNodes[ i ] = 0;
    for ( i = 0 ; i <= 4 ; i++ )
      exNodes[ i ] = 0;
    for ( i = 0 ; i <= 4 ; i++ )
      edNodes[ i ] = 0;
    for ( i = 0 ; i <= 4 ; i++ )
      tetra[ i ] = 0;
    for ( i = 0 ; i <= 4 ; i++ )
      pyram[ i ] = 0;
    for ( i = 0 ; i <= 4 ; i++ )
      prism[ i ] = 0;
    
    setZerosMi( edgeCut, 4, 2 );
    setZerosMi( edge, 4, 4 );

    // Find the number of edges cut.
    numEdgeCut = elemEdgeCut( mesh, particle, theElem, NodeSurf,
                              inNodes, onNodes, exNodes, edgeCut );
    numNodesInOn = ( inNodes[ 0 ] + onNodes[ 0 ] ); // Total number of nodes in or on
    //printf(" inNodes = %d, onNodes = %d, exNodes = %d\n", inNodes[0], onNodes[0], exNodes[0]);

    assert_error( ( numNodesInOn + exNodes[ 0 ] ) <= PRE_NODES_PER_EL,
             "Number of nodes exceeds PRE_NODES_PER_EL, check elemEdgeCut()" );

    switch ( numNodesInOn )
    {  // Check the number of PNODES inside

    case 0:  // The element is not cut substantially, it is completely outside
      //printf("PNodes = %d, Num. edges cut = %d : TETRA-0 \n",numNodesInOn,numEdgeCut);
      // printf("Partial Element %d processed successfully\n",iElem);
      ElemPartly[ iElem ] = 0;     // Remove the processed element
      eLog[ 0 ] ++;
      break;

    case 1:  // Check whether it is in node or on node
      if ( inNodes[ 0 ] == 1 )
      {  // It is in node
        // Only possibility relevant is number of edges cut = 3
        if ( numEdgeCut == 3 )
        {
          //printf("PNodes = %d, Num. edges cut = %d : TETRA-1 \n",numNodesInOn,numEdgeCut);
          // ------------- TETRAHEDRA CREATED ----------------//
          // Get the inNode, (there is only one inNode)
          in = edgeCut[ 1 ][ 1 ];        // in node of first edge
          edNodes[ in ] = mesh->ConTab[ theElem ][ in ];
          
          for ( iEdge = 1 ; iEdge <= numEdgeCut ; iEdge++ )
          {
            out = edgeCut[ iEdge ][ 2 ];    // edge nodes in place of outnodes
            edNodes[ out ] = edgeCut[ iEdge ][ 0 ];
          } //end for iEdge
          
          // Mapping this tetra to standard tetra
          for ( i = 0 ; i <= 2 ; i++ )
            tetra[ face[ 1 ][ i ] ] = edNodes[ face[ in ][ i ] ];
          
          tetra[ 1 ] = edNodes[ in ];
          
          creatTetra( mesh, tetra, 1, theElem, ElemCon, SurfCon, conTSp, GeoTab,
            surfTSp, fullESp, NodeSurf );
          //printf("Element No. created = %d, Nodes Up to = %d\n",conTSp[0][0],(int)GeoTab[0][0]);
          // printf("Surface Element no. = %d\n", surfTSp[0][0]);
          // printf("Partial Element %d processed successfully\n",iElem);
          
          ElemPartly[ iElem ] = 0;  // Remove the processed element
          eLog[ 1 ] ++;
        } //endif numEdgeCut == 3
      } //endif inNodes[0] == 1
      // else it is an on node, we need at least one in node to form element
      // This element remains un-processed
      break;

    case 2:  // Check the possibilities:
      // 1) Both on nodes
      //   a) numEdgeCut = 1 or 2 possible --> unprocessed
      //   b) numEdgeCut > 2 --> not possible
      // 2) One on node and one in node
      //   a) numEdgeCut = 1 --> not possible
      //   b) numEdgeCut = 2 --> Tetrahedra created
      //   c) numEdgeCut > 2 --> not possible
      // 3) Both in nodes
      //   a) numEdgeCut < 4 --> not possible
      //   b) numEdgeCut = 4 --> Prism created
      //   c) numEdgeCut > 4 --> not possible
      // So, there are only two possibilities that are relevant
      if ( inNodes[ 0 ] == 1 && onNodes[ 0 ] == 1 )
      {
        if ( numEdgeCut == 2 )
        {
          on = 1;      // Index of on node
          //printf("PNodes = %d, Num. edges cut = %d : TETRA-2 \n",numNodesInOn,numEdgeCut);
          // ---------------- TETRAHEDRA CREATED -------------------//
          // Get in node and its two edges that are cut
          in = edgeCut[ 1 ][ 1 ];
          edNodes[ in ] = mesh->ConTab[ theElem ][ in ];
          
          for ( iEdge = 1; iEdge <= numEdgeCut; iEdge++ )
          {
            out = edgeCut[ iEdge ][ 2 ];
            edNodes[ out ] = edgeCut[ iEdge ][ 0 ];   // Fill edge node in out's place
          } // end for iEdge
      
          // Get the one on node
          while ( !onNodes[ on ] )
            on++;  // Get the index of on node
          
          edNodes[ on ] = mesh->ConTab[ theElem ][ on ];
      
          // Mapping this tetra to standard tetra
          for ( i = 0; i <= 2; i++ )
            tetra[ face[ 1 ][ i ] ] = edNodes[ face[ in ][ i ] ];
          
          tetra[ 1 ] = edNodes[ in ];
      
          creatTetra( mesh, tetra, 2, theElem, ElemCon, SurfCon, conTSp, GeoTab,
            surfTSp, fullESp, NodeSurf );
          //printf("Element No. created = %d, Nodes Up to = %d\n",conTSp[0][0],(int)GeoTab[0][0]);
          // printf("Surface Element no. = %d\n", surfTSp[0][0]);
          // printf("Partial Element %d processed successfully\n",iElem);
          ElemPartly[ iElem ] = 0;  // Remove the processed element
          eLog[ 2 ] ++;
        } //endif numEdgeCut== 2
      } //end if one in and one on node
      else if ( inNodes[ 0 ] == 2 )
      {
        if ( numEdgeCut == 4 )
        {
          //printf("PNodes = %d, Num. edges cut = %d : PRISM-2 \n",numNodesInOn,numEdgeCut);
          // ------------------ PRISM CREATED -----------------------//
          // Get the two in nodes and their edge nodes
          fin = 3;
          sin = 2;    // Variables to hold first in and second in
          
          for ( iEdge = 1 ; iEdge <= numEdgeCut ; iEdge++ )
          {
            in  = edgeCut[ iEdge ][ 1 ];
            out = edgeCut[ iEdge ][ 2 ];
            
            if ( !tetra[ in ] )
              tetra[ in ] = mesh->ConTab[ theElem ][ in ];
            
            if ( in < fin )
              fin = in;
            
            if ( in > sin )
              sin = in;
            
            edge[ in ][ out ] = edge[ out ][ in ] = edgeCut[ iEdge ][ 0 ];
          }
          // Search the face for the fin node
          for ( i = 0; i <= 2; i++ )
            if ( face[ sin ][ i ] == fin )
              break;
          
          i++;
      
          // Fill prism
          for ( j = 1 ; j <= 2 ; j++ )
          {
            i = i % 3;
            prism[ j + 1 ] = edge[ fin ][ face[ sin ][ i ] ];
            prism[ j + 4 ] = edge[ sin ][ face[ sin ][ i ] ];
            i++;
          }
          prism[ 1 ] = tetra[ fin ];
          prism[ 4 ] = tetra[ sin ];
      
          // Creat tetrahedra from prism
          creatPrism( mesh, prism, 2, theElem, ElemCon, SurfCon, conTSp,
                      GeoTab, surfTSp, fullESp, NodeSurf, NodeSurfInPar );
          //printf("Element No. created = %d, Nodes Up to = %d\n",conTSp[0][0],(int)GeoTab[0][0]);
          // printf("Surface Element no. = %d\n", surfTSp[0][0]);
          // printf("Partial Element %d processed successfully\n",iElem);
          ElemPartly[ iElem ] = 0;  // Remove the processed element
          eLog[ 3 ] ++;
        } // end if numEdgeCut == 4
      } // end else if both are in nodes
      // else the element is left un-processed
      break;

    case 3:  // Check the possibilities:
      // 1) All three nodes on nodes
      //   a) numEdgeCut = 0 --> rare and un-processed
      // 2) Two on nodes and one in node
      //   a) numEdgeCut = 1 --> Tetrahedra created
      //   b) numEdgeCut > 1 --> not possible
      // 3) One on node and two in nodes
      //   a) numEdgeCut < 2 --> not possible
      //   b) numEdgeCut = 2 --> Pyramid created
      //   c) numEdgeCut > 2 --> not possible
      // 4) All three in nodes
      //   a) numEdgeCut = 3 --> Prism created
      // There are three possibilities that are relevant and
      // numEdgeCut uniquely define each case.
      if ( numEdgeCut == 1 )
      {  // Two of the nodes on surface and one inside
        //printf("PNodes = %d, Num. edges cut = %d : TETRA-3 \n",numNodesInOn,numEdgeCut);
        // ----------------- TETRAHEDRA CREATED -----------------------//
        // Get the one in node and its edge node from iEdge = 1
        in = edgeCut[ 1 ][ 1 ];
        out = edgeCut[ 1 ][ 2 ];
        edNodes[ in ] = mesh->ConTab[ theElem ][ in ];
        edNodes[ out ] = edgeCut[ 1 ][ 0 ];
        // Get the two on nodes
        for ( on = 1 ; on <= PRE_NODES_PER_EL ; on++ )
          if ( onNodes[ on ] )
            edNodes[ on ] = mesh->ConTab[ theElem ][ on ];

        // Mapping this tetra to standard tetra
        for ( i = 0; i <= 2; i++ )
          tetra[ face[ 1 ][ i ] ] = edNodes[ face[ in ][ i ] ];
        
        tetra[ 1 ] = edNodes[ in ];

        creatTetra( mesh, tetra, 3, theElem, ElemCon, SurfCon, conTSp, GeoTab,
          surfTSp, fullESp, NodeSurf );
        //printf("Element No. created = %d, Nodes Up to = %d\n",conTSp[0][0],(int)GeoTab[0][0]);
        // printf("Surface Element no. = %d\n", surfTSp[0][0]);
        // printf("Partial Element %d processed successfully\n",iElem);
        ElemPartly[ iElem ] = 0;  // Remove the processed element
        eLog[ 4 ] ++;
      } //endif numEdgeCut == 1

      else if ( numEdgeCut == 2 )
      {  // One node on surface and other two inside
        //printf("PNodes = %d, Num. edges cut = %d : PYRAM-3 \n",numNodesInOn,numEdgeCut);
        // ----------------- PYRAMID CREATED -----------------------//
        // Get the in nodes and corresponding edge nodes
        out = edgeCut[ 1 ][ 2 ];  // There is only one out node
        for ( iEdge = 1 ; iEdge <= numEdgeCut ; iEdge++ )
        {
          in = edgeCut[ iEdge ][ 1 ];
          tetra[ in ] = mesh->ConTab[ theElem ][ in ];
          edNodes[ in ] = edgeCut[ iEdge ][ 0 ];
        }
        // Get the single on node
        on = 1;
        
        while ( !onNodes[ on ] )
          on++;

        // Search the face for the on node
        for ( i = 0; i <= 2; i++ )
          if ( face[ out ][ i ] == on )
            break;
        
        i++;
        // Fill pyram
        for ( j = 1 ; j <= 2 ; j++ )
        {
          i = i % 3;
          pyram[ j ] = tetra[ face[ out ][ i ] ];  // Fill in nodes
          pyram[ j + 3 ] = edNodes[ face[ out ][ i ] ];  // Fill edge nodes
          i++;
        }
        // Fill on node
        pyram[ 3 ] = mesh->ConTab[ theElem ][ on ];

        creatPyram( mesh, pyram, theElem, ElemCon, SurfCon, conTSp, GeoTab,
          surfTSp, fullESp, NodeSurf, NodeSurfInPar );
        //printf("Element No. created = %d, Nodes Up to = %d\n",conTSp[0][0],(int)GeoTab[0][0]);
        // printf("Surface Element no. = %d\n", surfTSp[0][0]);
        // printf("Partial Element %d processed successfully\n",iElem);
        ElemPartly[ iElem ] = 0; // Remove the processed element
        eLog[ 5 ] ++;
      } // endif numEdgeCut == 2

      else if ( numEdgeCut == 3 )
      {  // All three nodes are inside
        //printf("PNodes = %d, Num. edges cut = %d : PRISM-3 \n",numNodesInOn,numEdgeCut);
        // ------------------ PRISM CREATED ------------------------//
        // Get all the in nodes and their corresponding edge nodes
        out = edgeCut[ 1 ][ 2 ];  // Use first edge for out node
        for ( iEdge = 1 ; iEdge <= numEdgeCut ; iEdge++ )
        {
          in = edgeCut[ iEdge ][ 1 ];
          tetra[ in ] = mesh->ConTab[ theElem ][ in ];
          edNodes[ in ] = edgeCut[ iEdge ][ 0 ];
        } //end for iEdge

        // Fill in the prism
        for ( i = 1 ; i <= 3 ; i++ )
        {
          prism[ i ] = tetra[ face[ out ][ i - 1 ] ];
          prism[ i + 3 ] = edNodes[ face[ out ][ i - 1 ] ];
        }
        // Construct tetrahedra from prism
        creatPrism( mesh, prism, 3, theElem, ElemCon, SurfCon, conTSp,
                    GeoTab, surfTSp, fullESp, NodeSurf, NodeSurfInPar );
        //printf("Element No. created = %d, Nodes Upto = %d\n",conTSp[0][0],(int)GeoTab[0][0]);
        // printf("Surface Element no. = %d\n", surfTSp[0][0]);
        // printf("Partial Element %d processed successfully\n",iElem);
        ElemPartly[ iElem ] = 0; // Remove the processed element
        eLog[ 6 ] ++;
      } // endif numEdgeCut == 3
      break;

    case 4:  // Check these possibilities:
      // 1)   4 on Nodes and 0 in Nodes: Low resolution
      // 2)   3 on Nodes and 1 in Node : Add to fullElemSp
      // 3)   2 on Nodes and 2 in Nodes: Add to fullElemSp
      // 4)   1 on Node  and 3 in Nodes: Add to fullElemSp
      // 5)   0 on Nodes and 4 in Nodes: Should be in fullElemSp list
      // printf("inNodes = %d and onNodes = %d\n", inNodes[0], onNodes[0]);
      // ----------------- ADDED TO FULLELEMSP  -----------------------//
      //printf("PNodes = %d, Num. edges cut = %d : TETRA-4 \n",numNodesInOn,numEdgeCut);
      (*fullESp)[ 0 ]++;
      nFulElems = (*fullESp)[ 0 ];
      
      if ( ( nFulElems % MBLOCK ) == 1 )
      {
        (*fullESp) = ( int * ) realloc ( (*fullESp), ( nFulElems + MBLOCK ) * sizeof(int) );
        assert_error( (*fullESp) != NULL, "Memory allocation error");
      }   
          
      (*fullESp)[ nFulElems ] = theElem;
      // printf("Partial Element %d processed successfully\n",iElem);
      ElemPartly[ iElem ] = 0; // Remove the processed element
      eLog[ 7 ] ++;
      break;

    default:
      error( "Consider Case: Number of PNodes in = %d, Number of edges cut = %d \n",
          numNodesInOn, numEdgeCut );
      break;
    } // end switch
  } // for iElem
  //DISPLAY eLog
  eLog[ 8 ] = eLog[ 0 ] + eLog[ 1 ] + eLog[ 2 ] + eLog[ 3 ]
            + eLog[ 4 ] + eLog[ 5 ] + eLog[ 6 ] + eLog[ 7 ];
  //printf("\t\nTETRA-0 = %d\t\nTETRA-1 = %d\t\nTETRA-2 = %d\t\nPRISM-2 = %d\t\nTETRA-3 = %d\t\nPYRAM-3 = %d\t\nPRISM-3 = %d\t\nTETRA-4 = %d\t\nTOTAL = %d\n",eLog[0],eLog[1],eLog[2],eLog[3],eLog[4],eLog[5],eLog[6],eLog[7],eLog[8]);
  
  // Free memory of local connectivity table
  FreeM( edge );
  FreeM( edgeCut );

  for(int k = 0 ; k < 7 ; k++)
    free(surfTSp[k]);

  free(surfTSp);
  FreeM(SurfCon);
}

//#endif

