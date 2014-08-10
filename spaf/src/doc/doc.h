/** \mainpage Spaf user guide

\section intro_sec Introduction

Spaf is a 3D unsteady finite element code with fictitious domain formulation. It can simulate flows with or without particles. A 2 domains version can simulate very small particles motion.
 
 - Time discretization : second order backward Euler
 - Element : Taylor-Hood tetrahedron, second order for velocity, first for pressure.
 - Convection : method of characteristics
 - Velocity-pressure uncoupling : pressure correction
 - Solver : preconditioned conjugate gradient

 Algorithm: see:C. Diaz-Goano, P. Minev and K. Nandakumar, A fictitious
 domain/finite element method for particulate flows.
 J. Comp. Phys. 192 (2003), 105-123.
 
 

\section install_sec Installation

Spaf is written in C C99 standard with some portions in C++. It has been successfully compiled and used with gcc version 4.x.x and intel icc from 8.x to 10.x versions. Some old versions of the intel compiler seem buggy, make sure to use the last update of any version.
 
It relies on the following libraries:
  - Taucs: mandatory\n
    provides sparse matrix reordering and preconditioner computation. See taucs manual for installation, the library must be placed in the taucs folder of the source tree.
  - Blas: optional\n
    provides optimized linear algebra routines, standard blas, intel mkl or apple accelerate are possible known to work.

All compilation parameters are in the makefile file. The following preprocessor flags are available:
  - USE_SPARSE_MKL:\n
    use the sparse matrix vector product from intel mkl library, this requires spaf has been compiled with intel mkl blas library. So far these routines are slower than spaf's ones, it may change on newer platform and/or newer versions of mkl.
  - VERSION_Z:\n
    enables the 2 domains version of spaf for small particles.
  - PRE_INCREMENTAL_ROTATIONAL:\n
    enables the incremental pressure correction scheme, the default being non incremental.
  - ANN_USE_FLOAT:\n
    enables the use of the type float instead of double in ANN nearest neighbor search, this should be enabled when the minimum node to node distance in the mesh is bigger than the float machine epsilon.
 
Spaf can be run on shared memory parallel computer with OpenMP. To enable it, compile with the required flag (-fopenmp for gcc, -openmp for icc, see compiler manual otherwise). See usage section for how to run spaf in parallel.

To compile, just run make in the directory that contains the makefile file and the source directory src.
 
\section usage_sec Usage
 The usage is described for the 1 domain and 2 domain versions. The order of the arguments has no influence. Spaf can be interrupted and restarted from a previous state with a checkpoint (link) .
 
 
\subsection synopsys Synopsys
 spaf input_file.case input_file.num input_file.flu [input_file.par] {input_file.o.case} [-g] [-log] [-o path_to_output_directory]
 

\subsection mandatory_arg Mandatory arguments

\subsubsection case_file Case file
 The input_file.case file is the path to a directory with .case extension. This directory contains the mesh, initial and boundary conditions data. It is created with mesh2spaf.

 See \ref create_case for more informations.

 
 \subsubsection num_file Numerical parameters file
 This is a file with .num extension which contains the numerical parameters. See example.num.
 
 \subsubsection phy_file Fluid physical parameters file
 This is a file with .flu extension which contains the physical parameters of the fluid. See example.flu.
 
 
\subsection optional_arg Optional arguments
 
 \subsubsection par_file Particle parameters file
 This is a file with .par extension which contains the particles parameters. See example.par.
 
 \subsubsection case_o_file Outer Case file, for 2 domains version only
 This is only for the 2 domains  version of spaf. The input_file.o.case file is the path to a directory with .o.case extension. This directory contains the mesh, initial and boundary conditions data for the outer flow.
  
 \subsubsection g_option -g option
  Enables more output from the iteration process, use for debugging purposes.
 
 \subsubsection log_option -log option
  By default spaf will display some informations in the terminal. This option enables spaf to output these informations into a file instead. The file's name is log.date.time, where "date" is the creation date and "time" the creation time. This does not concern the output data files. This option is usefull when running spaf unattended.
 
 \subsubsection o_option -o option
  By default all ouput files are written in the directory where the spaf command was issued. This argument enables to write all output files to another directory, specified by the path to this directory provided after -o.
 
 
 \subsection parallel Run in parallel
 To run spaf in parallel, make sure it has been compiled with OpenMP, then at least set the OpenMP environment variables in the terminal such as 'export OMP_NUM_THREADS=2' for using 2 cpus. See OpenMP documentation for further parameters.
 At the time of writing these lines, the ANN nearest neighbor search is not thread safe and cannot be used with more than 1 cpu. If you try to use it with more than 1 cpu, spaf will detect it and issue an error.
 

 \section ouput_data Output data
 This section contains the description of the output data, for statistics, particles, fluid and checkpoints. The frequency at which those data are written can be specified in the \ref num_file.
 
 \subsection stat_data Statistics data
  Spaf provides some statistics about a simulation, such as the number of iterations of the solver, the number of element in a particle, the number of element searched by the characteristics method. For each kind of data, the last value recorded, the number of samples, the min, max and average are given. These data can be used to check if something is wrong with a simulation, of used to adjust some related parameters. Those data are written in a file named stat.log.
 
 See statistics.c for more informations.
 
 \subsection particle_data Particle data
 Particles trajectory and velocities history are stored in one file per particle.
 
 See the description of the function WriteParticleResults() for more informations.
 
 \subsection fluid_data Fluid data
 The fluid flow fields are stored so they can be visualized, this requires to use spaf2viz to convert those data files to a file in vtk format.
 See the file fluid_io.c as well as the description of the function WriteVizData() for more informations on the storage.
 When running the one domain version of spaf without any particle, a file called residuals.csv will be written. It contains the values of the residuals for  the 3 components of the velocity at each time steps. Use it for checking that a time independent fluid only simulation has converged in time.

 \subsection checkpoint_data Checkpoint data
 Checkpoints are used to continue a simulation that has been stopped (finished or interrupted). Simulation data are stored so that the it can be restarted just as it has not been stopped, without any loss of information or precision.
 
 See the file checkpoint.c for more informations.

 \section howto How to : 2 domains version
 
 \subsection howto1 Step 1 : get outer flow
 
 The input_file.o.case file provided when running the 2 domains version must contain both the outer domain mesh and flow field. The flow field is read from the initial conditions. If the flow field is known beforehand, create a case with it as the initial conditions. Otherwise the flow field has to be computed  by running the 1 domain version of spaf without any particle. The flow solution of this simulation will have to be placed in the 'initial_conditions' directory of the 'input_file.o.case'. Note that you can use the solution saved in the 'viz' or in the 'checkpoint' directory, this latter being always more precise as it is saved with machine precision, which is not the case by default for visualization data. Prefer the flow solution from a checkpoint (pick all the 4 fields from the directory '1' in the latest checkpoint saved).
 
 No computation is done with this case, it is just raw input data.
 
 \subsection howto2 Step 2 : create inner case
 
 The particle centroid must be at the origin of the inner mesh. In gmsh, the physical entity of the boundaries of the inner domain must be set as 'z'. This way mesh2spaf recognize that it is an inner domain and will process it accordingly, without anything to specify (no initial conditions and no boundary conditions). Moreover, the surface of this inner domain should be a parallelepiped. This way it is way more efficient to check when the foot of a characteristic is outside the inner domain. See the cube_cylinder.geo and cube_sphere.geo files for gmsh examples.
 
 \subsection howto3 Step 3 : parameter files
 
 The 'input_file.num' for the 2 domains simulation must have the 'bb_search_courant_limit_o' parameter set in the convection section. This is for enabling searching in the outer mesh.
 The initial velocity of the particle will be computed as the value of the velocity of the outer domain at the centroid position. The initial value specified in 'input_file.par' is not used at all. This is not the case for the initial angular velocity.

*/
