/** \file
  
  \author    Carolina 2002
  \author    Peter May 2003
  \author    Veeramani May 2005
  \author    Peter November 2005
  \author    Antoine Dechaume 2006-2009

  \brief 3D Unsteady Navier-Stokes solver with ellipsoidal particles
            based on the fictitious domain method.
  
  Algorithm: see:C. Diaz-Goano, P. Minev and K. Nandakumar, A fictitious
                 domain/finite element method for particulate flows.
                 J. Comp. Phys. 192 (2003), 105-123.

  Zooming box version requires the VERSION_Z preprocessor switch to be activated.
  The global domain is represented by the suffix _o, for outer domain.

*/

#include "main.h"

#define SPAF_VERSION "r120"

static void
PrintInfo(
int   argc,
char* argv[] )
{
  info( "\nSPAF version " SPAF_VERSION "\n" );
  
  info("Compiled the "__DATE__" at "__TIME__ );
  
#ifdef __INTEL_COMPILER
  info(" with icc %d",__INTEL_COMPILER);
#elif defined __GNUC__
  info(" with gcc %d.%d.%d",__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__);
#endif
  
#ifdef __sgi__
  info(" on sgi.""\n"); 
#elif defined __APPLE__
  info(" on apple.""\n");
#elif defined WIN32
  info(" on windows.""\n");
#endif

#ifdef USE_ACCELERATE
  info( "\nUsing apple blas\n" );
#elif USE_MKL
  info( "\nUsing mkl blas\n" );
#elif USE_ACML
  info( "\nUsing acml blas\n" );
#endif
  
#ifdef VERSION_Z
#ifdef TWO_WAY
  info( "\nVERSION_Z + TWO_WAY preprocessor flags activated: 2 grids version with coupling.\n" );
#else
  info( "\nVERSION_Z preprocessor flag activated: 2 grids version without coupling.\n" );
#endif
#endif
  
#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp master
{
  info( "\nOpenMP enabled : %d cpus, %d threads, dynamic is %d\n",
         omp_get_num_procs(), omp_get_num_threads(), omp_get_dynamic() );
}
}
#endif
  
#ifdef PRE_INCREMENTAL_ROTATIONAL
  info( "\nPRE_INCREMENTAL_ROTATIONAL flag enabled, using rotational pressure correction method.\n" );
#endif
  
#ifdef ANN_USE_FLOAT
  info("\nANN_USE_FLOAT flag enabled, using float for ANN library.\n");
#endif
  
  info( "\nMachine epsilon = %g\n", DBL_EPSILON );
  info( "int max = %u\n", INT_MAX );
  
  // print byte order
  if ( LittleEndian() )
    info("\nLittle endian byte order\n\n");
  else
    info("\nBig endian byte order\n\n");
  
  time_t rawtime;
  struct tm * timeinfo;
  
  time( &rawtime );
  timeinfo = localtime( &rawtime );
  info( "Run date and time:\n\t%s""\n", asctime (timeinfo) );
  
  // print working directory
  char* WorkingDirectory = GetWorkingDirectory();
  info("Working directory:\n\t%s""\n",WorkingDirectory);
  free(WorkingDirectory);
  
  // print how we were called
  info("\n""Command line:\n");
  for ( int i = 0 ; i < argc ; i++ )
    info("\t""%s""\n",argv[i]);
  
  info("\n");
}












int
main(
int   argc,
char* argv[] )
{
//=============================================================================
// Parse command line arguments
//=============================================================================

// sort input files
char *PathToCase = NULL,
#ifdef VERSION_Z
     *PathToCase_o = NULL,
#endif
// by default, the output directory is the current working directory
     *PathToOutput = GetWorkingDirectory(),
     *PathToParticlesParameters = NULL,
     *PathToFluidParameters = NULL,
     *PathToNumericalParameters = NULL;

// default logging level, writes up to info, on stdout
int LogLevel = 3;
bool WriteLogToDisc = false;
  
for ( int iarg = 1 ; iarg < argc ; iarg++  )
{
  if ( CheckFileExtension( "flu" , argv[iarg] ) )
    PathToFluidParameters = GetAbsolutePath(argv[iarg]);
  
  else if ( CheckFileExtension( "num" , argv[iarg] ) )
    PathToNumericalParameters = GetAbsolutePath(argv[iarg]);
  
  else if ( CheckFileExtension( "par" , argv[iarg] ) )
    PathToParticlesParameters = GetAbsolutePath(argv[iarg]);
  
#ifdef VERSION_Z
  else if ( CheckFileExtension( "o.case" , argv[iarg] ) )
    PathToCase_o = GetAbsolutePath(argv[iarg]);
#endif
    
  else if ( CheckFileExtension( "case" , argv[iarg] ) )
    PathToCase = GetAbsolutePath(argv[iarg]);

  else if ( StringCompare( argv[iarg],"-g" ) )
    LogLevel = 4;
  
  else if ( StringCompare( argv[iarg],"-log" ) )
    WriteLogToDisc = true;
  
  else if ( StringCompare( argv[iarg],"-o" ) )
  {
    // increment loop variable since we read 2 items here
    iarg++;
    free(PathToOutput);
    PathToOutput = GetAbsolutePath(argv[iarg]);
  }
  
  else
    warning("The extension of %s is invalid \n",argv[iarg]);
}

// set up logging, by default writes to stdout
char* PathToLog = StringDuplicate("stdout");
  
// if we write the log to a file, it will be in the output directory and in the file named log with an extension that is managed by SetLog
if ( WriteLogToDisc )
{
  free(PathToLog);
  PathToLog = StringConcat(PathToOutput,"/log");
}

SetLog(PathToLog, LogLevel);

// Print informations about platform, compiler and preprocessor flags
PrintInfo(argc,argv);
  
SetOutputDirectory(PathToOutput);
  
//=============================================================================
// Read parameters
//=============================================================================
fluid_t* fluid = ReadFluidParameters( PathToFluidParameters );

particle_t* particle = ReadParticlesParameters( PathToParticlesParameters );

double dt = 0.;         // time step increment
int TimeStepNb = 0; // time step counter

ReadNumericalParameters( PathToNumericalParameters, &dt, &TimeStepNb );

bc_t *bc = ReadBC( PathToCase );
mesh_t *mesh = ReadMesh( PathToCase );

// Coefficient for time discretization, start with first order backward euler
double tau[3] = { + 1. / dt, - 1. / dt, 0. };  

// Prepare the operators to solve fluid flow
UseMicroGrid(); // we have to keep that one for any version
GetOperators( tau[0], fluid->Re, mesh );


#ifdef VERSION_Z
//==============================================================================
// Setup some special parameters for version z
//==============================================================================
PrintTitle( "Preparing version z stuff" );
  
mesh_t *mesh_o = ReadMesh( PathToCase_o );
  
// make the operator that projects the velocity field on the surface of the inner mesh to a weakly mass conserving field
GetMassConservingOperator(PathToCase,mesh);

// Setup fluid flow fields for outer domain
FluidVar_t *FluidVar_o = FluidVarAlloc();

// Create internal structure of FluidVar_o
// we don't need all arrays of this structure for outer domain,
// only pressure and velocity
AllocVdouble( mesh_o->NbOfPressureNodes, FluidVar_o->Pre );

FluidVar_o->Vel[0] = NULL;  
for ( int dir = 1 ; dir <= 3 ; dir++ )
  AllocVdouble( mesh_o->NbOfNodes, FluidVar_o->Vel[dir] );
  
// Read outer flow field
ReadFluidInitialConditions(
  PathToCase_o, mesh_o, FluidVar_o->Vel, FluidVar_o->Pre );

// the inner box moves with the particle in all directions
for ( int dir = 1 ; dir <= 3 ; dir++ )
  fluid->FrameVelDir[dir] = true;

ParticleInitialConditions(mesh_o, FluidVar_o, particle);

#ifdef TWO_WAY
// we will have to solve the flow at order 1,
// hence we need some more fluid variables
FluidVar_o->VelOld1[0] = NULL;  
for ( int dir = 1 ; dir <= 3 ; dir++ )
  AllocVdouble( mesh_o->NbOfNodes, FluidVar_o->VelOld1[dir] );

// VelOld1 = Vel
for ( int dir = 1 ; dir <= 3; dir++ )
  copy(mesh_o->NbOfNodes+1, FluidVar_o->Vel[dir], FluidVar_o->VelOld1[dir]);
  
FluidVar_o->VelRHS[0] = NULL;  
for ( int dir = 1 ; dir <= 3 ; dir++ )
  AllocVdouble( mesh_o->NbOfNodes, FluidVar_o->VelRHS[dir] );

// used in NavierStokes.c as temporary array
FluidVar_o->Acc[0] = NULL;  
for ( int dir = 1 ; dir <= 3 ; dir++ )
  AllocVdouble( mesh_o->NbOfNodes, FluidVar_o->Acc[dir] );

// we need bc as well
bc_t *bc_o = ReadBC(PathToCase_o);

UseMacroGrid();
info("\nComputing macro grid operators\n");
GetOperators(tau[0], fluid->Re, mesh_o);
#endif

#endif

//=============================================================================
// pre-processing is finished, we move to the output directory
//=============================================================================
PrintTitle( "End of pre-processing" );
char* CurrentDirectory = GoToOutputDirectory();
assert_Directory( CurrentDirectory, "output" );
free(CurrentDirectory);

//==============================================================================
//==============================================================================
// End allocating and preparing
//==============================================================================
//==============================================================================

FluidVar_t* FluidVar = FluidVarCreate(mesh->NbOfPressureNodes, mesh->NbOfNodes);

int TimeStep = 0;  // Counter for TimeSteps
int order = 1; // time discretization order

// if a checkpoint is found, read it and continue the computation from there on.
bool restart = ReadCheckPoint(&TimeStep, FluidVar, particle);

double time = TimeStep * dt;

if ( restart == false )
{
  //=============================================================================
  // Prepare initial conditions
  //=============================================================================  
#ifndef VERSION_Z
  ReadFluidInitialConditions(
    PathToCase, mesh, FluidVar->VelOld1, FluidVar->Pre );
#else
  // compute initial fluid state (VelOld1 and Pre) by interpolation from outer domain
  GetFluidInitialConditions_z( particle[1].Pos, mesh, mesh_o, FluidVar_o, FluidVar );
#endif 

  // Vel = VelOld1
  for ( int dir = 1 ; dir <= 3; dir++ )
    copy( mesh->NbOfNodes+1, FluidVar->VelOld1[dir], FluidVar->Vel[dir] );  

  if ( particle )
  {    
    // set rigid body motion inside particles
    GetParticleMesh(mesh, particle);
    SetRigidBodyMotion(mesh, particle, FluidVar);

    // VelOld1 = Vel, because SetRigidBodyMotion works on Vel
    for ( int dir = 1 ; dir <= 3; dir++ )
      copy( mesh->NbOfNodes+1, FluidVar->Vel[dir], FluidVar->VelOld1[dir] );    
  }
  
  // save initial data
  WriteOutput( FluidVar, particle, time, TimeStep );

  //=============================================================================
  // Proceed first time step
  //=============================================================================
  time = dt;
  TimeStep = 1;

  PrintTimeStep( time, TimeStep );

  // particles prediction
  PredictParticlePosition( order, mesh, dt, fluid, particle );

  // Navier-Stokes
#ifdef VERSION_Z
#ifdef TWO_WAY
  UseMicroGrid();
#endif

  FluidNavierStokes(order, mesh, bc, mesh_o, FluidVar_o->Vel, fluid, particle, tau,
                    dt, FluidVar);
#else
  
  FluidNavierStokes(order, mesh, bc, fluid, particle, tau, dt, FluidVar);
#endif
  
  // particles RBI
  RigidBody( mesh, tau, FluidVar, particle );

  // particle correction
  CorrectParticlePosition( order, dt, tau, fluid, particle );

  // save data
  WriteOutput( FluidVar, particle, time, TimeStep );

  // fill previous time steps fields
  FluidShiftVar( FluidVar );
  ParShiftVar( particle );
}

//==============================================================================
// NOW SECOND ORDER SCHEME
//==============================================================================
order = 2;

// update time derivative coefficients
tau[0] = 3. / 2. / dt;
tau[1] =    - 2. / dt;
tau[2] = 1. / 2. / dt;
  
// Update velocity operators to second order
#ifdef VERSION_Z && TWO_WAY
UseMicroGrid();
#endif
  
UpdateOperators( dt );

//==============================================================================
//==============================================================================
// TIME STEPPING LOOP
//==============================================================================
//==============================================================================
PrintTitle( "Starting time stepping" );

for ( TimeStep = TimeStep + 1 ; 
      TimeStep <= TimeStepNb ; 
      TimeStep++ )
{
  time += dt;
  
  PrintTimeStep( time, TimeStep );
  
  // particles : prediction step
  PredictParticlePosition( order, mesh, dt, fluid, particle );
  
  // Navier-Stokes
#ifdef VERSION_Z
#ifdef TWO_WAY
  UseMacroGrid();

  const double tau_2_way[3] = {1./dt, -1./dt, 0.};

  FluidNavierStokes(1, mesh_o, bc_o, mesh, FluidVar->Vel, fluid, particle, tau_2_way,
                    dt, FluidVar_o);

  UseMicroGrid();
#endif
  
  FluidNavierStokes(order, mesh, bc, mesh_o, FluidVar_o->Vel, fluid, particle, tau,
                    dt, FluidVar);
#else
  
  FluidNavierStokes(order, mesh, bc, fluid, particle, tau, dt, FluidVar);
#endif
  
  // particles RBI
  RigidBody( mesh, tau, FluidVar, particle );

  // particles : correction step
  CorrectParticlePosition( order, dt, tau, fluid, particle );
  
  // save data
  WriteOutput( FluidVar, particle, time, TimeStep );
  
  // fill previous time steps fields
  FluidShiftVar( FluidVar );
  ParShiftVar( particle );
}

//==============================================================================
// Post computation stuff
//==============================================================================
PrintTitle("End of simulation");
  
CloseLog();

return 1;
}
