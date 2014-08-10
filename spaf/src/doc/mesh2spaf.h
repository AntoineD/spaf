/** \page create_case Creating a case

 After a mesh has been created with gmsh, the mesh needs to be converted to spaf format and reordered, initial and boundary conditions have to be created as well. This is done with mesh2spaf. Depending on the physical problem, other parameters like the definition of planes for particle-wall collision detection may be necessary. Here we describe briefly how to create the mesh, the format of a case, the usage of mesh2spaf and how to define collision planes.

 \section gmsh gmsh
 Gmsh is a 3D finite element mesh generator, for manuals and up to date details, visit http://www.geuz.org/gmsh. This project is very active and new versions come often so be sure to check it once in a while. It can import geometry from CAD softwares or create it from scratch. The geometry is composed of one volume and at least one boundary surface. Spaf uses second order 10 nodes tetrahedron, those elements could be quadratic or linear.
 
 \subsection gmsh_param Parameters
 Here is an excerpt of my .gmsh-options file (located in the user's home folder, comments are preceded by a '//') with non-default options defining the mesh generation:
 
 \code
 Mesh.Algorithm = 6; // 2D mesh algorithm (1=MeshAdapt+Delaunay, 4=MeshAdapt, 5=Delaunay, 6=Frontal)
 Mesh.Binary = 1; // Write mesh files in binary format
 Mesh.ElementOrder = 2; // Element order (1=linear elements, N (<6) = elements of higher order)
 Mesh.Optimize = 1; // Optimize the mesh to improve the quality of tetrahedral elements
 Mesh.OptimizeNetgen = 1; // Optimize the mesh using Netgen to improve the quality of tetrahedral elements
 //Mesh.SecondOrderLinear = 1; // Should second order vertices simply be created by linear interpolation?
 Mesh.Smoothing = 10; // Number of smoothing steps applied to the final mesh
\endcode
 
 Lately I used the Frontal as 2D meshing algorithm which gives good results but is considered as experimental by the gmsh developpers, so be sure to check your mesh, eventhough I had no problem myself. If you don't feel confident enough to use it, just coment out the corresponding line to use the default 2D meshing algorithm. Moreover, the best is to write the mesh in binary format as it saves a good deal of space and retains the full precision of the data.
 
 \subsection gmsh_entities Physical entities
 In order to convert a mesh with mesh2spaf, physical entities have to be defined, they are used to identify the computational domain and the boundary conditions. The only volume entity has to be defined with the name "domain", and at least one entity has to be defined as a surface with a name among "wall", "inlet", "outlet", "shear" or "z".
 
 The file zone_definitions.py in mesh2spaf sources is used to define all the physical entities, you'll have to modify it in order to add different boundary conditions.
 
 \section case_format Case format
 A case is a directory that contains a mesh, a set of initial conditions and a set of boundary conditions. The extension of a case name is simply '.case'.
 
 \subsection mesh Mesh format
 A mesh is a directory named 'mesh' that contains the following files:
  - connectivity : a data file that contains the mesh connectivity.
  - points : a data file that contains the coordinates of the mesh nodes.
  - free_nodes_pressure : a data file that contains the number of free pressure nodes.
  - free_nodes_velocity : a data file that contains the number of free velocity nodes.
  - velocity_pressure_map : a data file that contains the mapping from velocity nodes indexes to pressure nodes indexes.
  - collision_planes : a file that contains the definition of at least one plane, used for particle collision detection, see \ref collision_planes section.
  - connectivity.z : created only for the 2 domains version of spaf, it contains the surface connectivity of the boundary of the inner domain.
 
 The nodes indexes will be ordered in such a way that the nodes that are not prescribed have the lowest indexes. The reason is to have a contiguous set of indexes that refer to the nodes at which the physical problem is solved. Free nodes have indexes ranging from 1 to the number of free nodes given by free_nodes_velocity for the velocity and free_nodes_pressure for the pressure. Prescribed nodes have indexes that range from the number of free nodes + 1 to the total number of nodes in the mesh. This total number of nodes is given by the number of points in 'points'. Note that no distinction is made among prescibed nodes that belongs to different physical entities.
 
 For more details on the format of the files that define mesh, see mesh.c.
 
 Eventually, a mesh can contain an 'info' file that has informations on the creation of the mesh. This file is ignored by spaf.
 
 \subsection initial_conditions Initial conditions
 A set of initial conditions is a directory named 'initial_conditions' that contains a fluid field that will be read into the previous time step at the beginning of a simulation. For more details, see fluid_io.c and initial_conditions.c.
 
 \subsection boundary_conditions Boundary conditions
 A set of boundary conditions is a directory named 'boundary_conditions' that contains a fluid field for all the prescibed nodes. Spaf can deal with the following boundary conditions:\n
 - pressure : homogeneous Dirichlet
 - velocity : Dirichlet for all components, or homogeneous Neumann for all components (treated as non-prescribed).
 
 For more details, see fluid_io.c and bc.c.

 \section collision_planes Collision planes
 It is sometimes necessary to detect the collision of a particle with a wall. This can be done easily for a plane wall by defining the planes equations in a file named collision_planes. There the number of planes and the planes equations can be defined, an example is provided inside the parameters_files directory. This file must be created by hand and placed inside the mesh directory of a case.

 \section mesh2spaf_usage mesh2spaf usage
 Mesh2spaf is written in python and requires no compilation, to use it just invoque its source file mesh2spaf.py. To get informations about how to use it, try mesh2spaf -help.
 Mesh2spaf will first read the mesh from the .msh file produced by gmsh, then it will re-order the nodes such that free nodes have the lowest indexes.
 Then if your mesh contains an inlet and outlet, you'll have to provide, for each of them, the radius, the coordinates of the center and the components of the normal pointing outside of its circumscribed disc.
 Then you will be asked to choose the type of initial conditions, if you do not know it beforehand, choose 'fluid at rest'.
 Then boundary conditions are created. For 'inlet' and 'shear' boundary conditions, you can either create them from scratch or fill them with the fields of the initial conditions which has been specified for all nodes. Generally, you will have to specify those by hand when the initial conditions are not known.
 Finally, all data are written to disc.
 
 
 */
