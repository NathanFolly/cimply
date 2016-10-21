static char help[] = "simply FEM program. To be extended...";

/* The PETSc packages we need: */
#include <petscksp.h>
#include <petscdmplex.h>
#include <petscds.h>
#include <petscsnes.h>
#include <petscts.h>
/* #include <petscviewerhdf5.h> */

/* Other packages we need */
#include <stdio.h>
#include <string.h>
#include "cimply.h"
#include "cimplyDF.h"

/* user defined Application Context that helps to manage the options (or the contetext) of the program  */

typedef enum {RUN_FULL, RUN_TEST} RunType;
typedef enum {NONE, LABELS, BOUNDARIES, SOLUTION} Visualization;

typedef struct  /* The Leapfrog timestepping context */
{
  PetscReal total;  /* total time in seconds */
  PetscReal dt;     /* time setp size */
  PetscBool halfstep;  /* Is the timestep a halfstep? */
} TSctx;

typedef struct  /* Material structure */
{
  PetscReal mu;
  PetscReal lbda;
}Material_type;


typedef struct {
  PetscInt debug;
  RunType runType;
  PetscBool showInitial,showSolution;
  /* definitions for mesh and problem setup */
  PetscInt dim;
  PetscBool simplex;
  PetscBool interpolate;
  PetscReal refinementLimit;
  PetscBool testPartition;
  PetscReal mu;			   /* The shear modulus */
  /* Additional options and parameters */
  PetscBool verbose;
  Visualization visualization;
  PetscBool neumann;
  PetscBool transient;
  PetscBool ldis;  /* for large displacement */
  PetscBool planestress;  /* 2D analysis type: plane stress */
  TSctx time;
  Material_type material;
  PetscViewer viewer;
  
}AppCtx;


/* Some global variables we use across functions. mainly so that we do not
 * need to create the whol FEM context from scratch every time we call the
 * Cimply solver */
SNES snes;			/* nonlinear solver */
KSP ksp;                        /* the linear sovler context */
DM dm, distributeddm;			/* problem definition */
Vec u,r;			/* solution and residual vectors */
Mat A,J,P;			/* Jacobian Matrix */
AppCtx user;			/* user-defined work context */
PetscViewer viewer;
TS ts;                          /* in case of transient analysis */
/*Done with the global variables for the FEM context  */

/* The static inline statement will tell the compiler to inline the function in the code directly instead of making multiple function calls. This improves performance especially for small functions that are called repeatedly throughout the program. */



SimmerDataStruct SimmerData  = {0,0,0,0,0,NULL,NULL,NULL};  /* Create an
                                                               instance of
                                                               SimmerData outside of
                                                               * any function
                                                               so that it
                                                               will be globally accessible */


PetscErrorCode zero_scalar(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx)
{
  u[0]=0.0;
  return 0;
}

PetscErrorCode zero_vector(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx)
{
  const PetscInt Ncomp = dim;
  PetscInt comp;
  for (comp=0;comp<Ncomp;comp++) u[comp]=0.0;
  return 0;
}

PetscErrorCode constrict_y(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx)
{
  u[1] = 0.0;
  return 0;
}

PetscErrorCode constrict_x(PetscInt dim, PetscReal time, const PetscReal x[], const PetscReal n[], PetscInt Nf, PetscScalar *u, void *ctx)
{
  u[0] = 0.0;
  return 0;
}


PetscErrorCode constrict_z(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx)
{
  u[2] = 0.0;
  return 0;
}

PetscErrorCode coordinates(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx)
{
  const PetscInt Ncomp = dim;
  PetscInt comp;
  for (comp=0; comp<Ncomp;comp++) u[comp]=x[comp];
  return 0;
}

PetscErrorCode pull(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx)
{
  PetscInt comp;
  for(comp=0;comp<dim;comp++){
    if(comp==1){u[comp]=0.1;}
    else{u[comp]=0.0;}
  }
    return 0;
}



PetscErrorCode ProcessOptions(MPI_Comm comm, AppCtx *options)
{
  const char *runTypes[2]={"full","test"};
  const char *visutype[4]={"none","labels","boundaries","solution"};
  PetscInt vis;
  PetscInt run;
  PetscErrorCode ierr;
  options->debug = 0;
  options->runType = RUN_FULL;
  options->dim = 3;
  options->interpolate = PETSC_TRUE;
  options->simplex = PETSC_TRUE;
  options->refinementLimit = 0.0;
  options->mu = 1;
  options->testPartition = PETSC_FALSE;
  options->showInitial = PETSC_FALSE;
  options->showSolution = PETSC_FALSE;
  options->verbose = PETSC_FALSE;
  options->visualization = SOLUTION;
  options->neumann = PETSC_TRUE;
  options->transient = PETSC_FALSE;
  options->ldis = PETSC_FALSE;
  options->time.total = 0.001;
  options->time.dt = 0.000005;
  options->planestress = PETSC_TRUE;

  options->material.mu = 76.9;
  options->material.lbda = 115.4;


  ierr = PetscOptionsBegin(comm,"", "Linear elasticity problem options", "DMPLEX"); CHKERRQ(ierr);
  ierr = PetscOptionsInt("-debug","The debugging level","cimply.c",options->debug, &options->debug,NULL);CHKERRQ(ierr);
  run = options->runType;
  ierr = PetscOptionsEList("-run-type","The run type","cimply.c",runTypes,2,runTypes[options->runType],&run,NULL);CHKERRQ(ierr);
  
  options->runType = (RunType) run;

  ierr = PetscOptionsInt("-dim","The topological mesh dimension", "cimply.c",options->dim,&options->dim,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-interpolate","Generate intermediate mesh elements", "cimply.c",options->interpolate,&options->interpolate,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-simplex","use simplex or tensor product cells","cimply.c",options->simplex,&options->simplex,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-refinementLimit","Largest allowable cell volume", "cimply.c",options->refinementLimit,&options->refinementLimit,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-test_partition","Use a fixed partition for testing", "cimply.c", options->testPartition,&options->testPartition,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-shear_modulus","The shear modulus","cimply.c",options->mu,&options->mu,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-show_initial","Output the initial guess for verification","cimply.c",options->showInitial,&options->showInitial,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-show_solution","Output the solution for verification","cimply.c",options->showSolution,&options->showSolution,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-verbose","Output additional information.","cimply.c",options->verbose,&options->verbose,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEList("-visualize","Visualize a certain object for debugging and learning purposes","cimply",visutype,4,visutype[options->visualization],&vis,NULL);CHKERRQ(ierr);
  options->visualization = (Visualization) vis;
  ierr = PetscOptionsBool("-neumann", "Apply Neumann boundary conditions on the rightmost face.","cimply.c",options->neumann, &options->neumann,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-transient", "Conduct a transient analysis. Default false.","cimply.c", options->transient, &options->transient,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-large_displacement", "Use Lagrangian strain and the second piola Kirchoff stress for large displacments","cimpy.c",options->ldis, &options->ldis, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-planestress", "Use plane stress formulation in case of a 2D analysis. Default true.","cimply.c",options->planestress, &options->planestress,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);
  return(0);
    }

PetscErrorCode DMVecViewLocal(DM dm, Vec v, PetscViewer viewer)
{
  Vec lv;
  PetscInt p;
  PetscMPIInt rank,numProcs;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)dm),&rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PetscObjectComm((PetscObject)dm),&numProcs); CHKERRQ(ierr);
  ierr = DMGetLocalVector(dm, &lv);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(dm, v, INSERT_VALUES, lv); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(dm, v, INSERT_VALUES, lv); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Local function \n"); CHKERRQ(ierr);
  for (p=0; p<numProcs; p++){
    if (p==rank){ierr=VecView(lv, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);}
    ierr = PetscBarrier((PetscObject)dm);CHKERRQ(ierr);
  }
  ierr = DMRestoreLocalVector(dm, &lv); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode Monitor(TS ts, PetscInt step, PetscReal crtime, Vec u, void *ctx)
{
  AppCtx *user=(AppCtx*) ctx;
  PetscErrorCode ierr;
  PetscViewer viewer;
  char filename[50];
  /* Vec glovec; */

  /* ierr = DMGetGlobalVector(dm, &glovec);CHKERRQ(ierr); */
  /* glovec = u; */
  /* ierr = DMRestoreGlobalVector(dm,&glovec);CHKERRQ(ierr); */
  sprintf(filename,"solution%i.vtu",step);
  ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  if (user->verbose) PetscPrintf(PETSC_COMM_WORLD,"Writing the solution to the vtk file. \n");
  /* ierr = PetscObjectSetName((PetscObject) u,"deformation");CHKERRQ(ierr); */
  ierr = VecView(u,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  if (user->verbose)PetscPrintf(PETSC_COMM_WORLD,"Done creating the VTK file. \n");


  /* ierr = VecView(u,user->viewer);CHKERRQ(ierr); */
    
  return(0);
}


 PetscErrorCode CreateVelocityNullSpace(DM dm, AppCtx *user, Vec *v, MatNullSpace *nullSpace)
 {
   Vec              vec;
   PetscErrorCode (*funcs[2])(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void* ctx) = {zero_vector, zero_vector};
   PetscErrorCode   ierr;

   DMGetGlobalVector(dm, &vec);
   DMProjectFunction(dm, 0.0, funcs, NULL, INSERT_ALL_VALUES, vec);
   VecNormalize(vec, NULL);
   if (user->debug) {
     PetscPrintf(PetscObjectComm((PetscObject)dm), "Velocity Null Space\n");
     VecView(vec, PETSC_VIEWER_STDOUT_WORLD);
   }
   MatNullSpaceCreate(PetscObjectComm((PetscObject)dm), PETSC_FALSE, 1, &vec, nullSpace);
   if (v) {
     DMCreateGlobalVector(dm, v);
     VecCopy(vec, *v);
   }
   DMRestoreGlobalVector(dm, &vec);
   /* New style for field null spaces */
   {
     PetscObject  velocity;
     MatNullSpace nullSpaceVel;

     DMGetField(dm, 1, &velocity);
     MatNullSpaceCreate(PetscObjectComm(velocity), PETSC_TRUE, 0, NULL, &nullSpaceVel);
     PetscObjectCompose(velocity, "nullspace", (PetscObject) nullSpaceVel);
     MatNullSpaceDestroy(&nullSpaceVel);
   }
 return(0);
}



PetscErrorCode SetupProblem(DM dm, AppCtx *user)
{
  PetscDS prob; 		   /* a DS is a "discrete system" managing discretisations */
  PetscErrorCode ierr;
  PetscInt NumFields;


  ierr = DMGetDS(dm, &prob);CHKERRQ(ierr);


  
     /* Setting up the linear elasticity problem */
    if (user->transient){
      ierr = PetscDSSetResidual(prob,0,f0_u_transient,f1_u);CHKERRQ(ierr);
      /* ierr = PetscDSSetResidual(prob,0,f0_u,f1_u);CHKERRQ(ierr); */
      ierr = PetscDSSetJacobian(prob,0,0,NULL,NULL,NULL,g3_uu_3d);CHKERRQ(ierr);
    }
    else if(user->ldis){
      if(user->verbose) ierr = PetscPrintf(PETSC_COMM_WORLD,"Setting up large displacement problem\n");
      ierr = PetscDSSetResidual(prob,0,f0_u,f1_u_ldis);CHKERRQ(ierr);
      ierr = PetscDSSetJacobian(prob,0,0,NULL,NULL,NULL,g3_uu_ldis);CHKERRQ(ierr);
    }
    else{
      
      if(user->dim==3){
        PetscPrintf(PETSC_COMM_WORLD, "Setting up Jacobian for the small strains-3D case\n" );
        ierr = PetscDSSetResidual(prob,0,f0_u,f1_u);CHKERRQ(ierr);
        ierr = PetscDSSetJacobian(prob,0,0,NULL,NULL,NULL,g3_uu_3d);CHKERRQ(ierr);
      }
      else if(user->planestress){
        ierr = PetscDSSetResidual(prob,0,f0_u,f1_u_2d_pstress);CHKERRQ(ierr);
        ierr = PetscDSSetJacobian(prob,0,0,NULL,NULL,NULL,g3_uu_2d_pstress);CHKERRQ(ierr);
      }
      else{
        PetscPrintf(PETSC_COMM_WORLD, "Setting up Jacobian for the small strains-2D case\n" );
      ierr = PetscDSSetResidual(prob,0,f0_u,f1_u);CHKERRQ(ierr);
      ierr = PetscDSSetJacobian(prob,0,0,NULL,NULL,NULL,g3_uu_2d);CHKERRQ(ierr);
      }
    }


    /* setting up the velocity field for the transient analysis */
    if(user->transient){
      ierr = PetscDSSetResidual(prob,1,f0_vel,f1_vel);CHKERRQ(ierr);
      ierr = PetscDSSetJacobian(prob,1,1,g0_velvel,NULL,NULL,NULL);CHKERRQ(ierr);
      ierr = PetscDSSetDynamicJacobian(prob,1,0,g0_velu,NULL,NULL,NULL);CHKERRQ(ierr);
      /* ierr = PetscDSSetJacobian(prob,1,0,NULL,NULL,g2_velu,NULL);CHKERRQ(ierr); */
      ierr = PetscDSSetDynamicJacobian(prob,0,1,g0_uvel,NULL,NULL,NULL);CHKERRQ(ierr);
      ierr = PetscDSSetBdResidual(prob,1,f0_vel_bd,f1_vel_bd);CHKERRQ(ierr);
      /* ierr = PetscDSSetJacobian(prob,1,1,g0_velvel_2d,NULL,NULL,NULL);CHKERRQ(ierr); */
      /* Maybe we need a dynamic jacobian? -- apparently not. */
      /* ierr = PetscDSSetDynamicJacobian(prob,1,0,g0_dyn_velu_2d,NULL,NULL,NULL);CHKERRQ(ierr); */
    }
  
  /* Setting the Neumann Boudnary Condition */
  if (user->neumann){
    if (user->ldis){
      ierr = PetscDSSetBdResidual(prob,0,f0_u_bd_ldis,f1_u_bd);CHKERRQ(ierr);
    }
    else{
      ierr = PetscDSSetBdResidual(prob,0,f0_u_bd,f1_u_bd);CHKERRQ(ierr);
      /* ierr = PetscDSSetBdJacobian(prob,0,0,NULL,g1_uu_bd_2d,NULL,NULL);CHKERRQ(ierr); */
    }
    
  }

  
  ierr = PetscDSGetNumFields(prob,&NumFields);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Problem set up with %i fields\n",NumFields);CHKERRQ(ierr);
      
  return(0);
}

PetscErrorCode SetupDiscretization(DM dm, AppCtx *user){

  DM cdm = dm;
  const PetscInt dim = user->dim; 	/* need to adapt this when changing
				   the dimensions of hte code */
  PetscFE fe, fe_bd, fe_vel, fe_vel_bd;
  PetscDS prob;
  PetscErrorCode ierr;
  PetscBool simplex = user->simplex;
  PetscQuadrature q;
  PetscInt order;
  PetscInt BSpaceOrder;
  PetscSpace space;

  /* Creating the FE */
  ierr = PetscFECreateDefault(dm, dim, dim, simplex,"def_",PETSC_DEFAULT,&fe);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fe, "deformation");CHKERRQ(ierr);
  
  if(user->neumann){
    if (user->ldis) PetscPrintf(PETSC_COMM_WORLD,"Large displacement formulation has no option for neumann conditions as of yet");
    ierr = PetscFEGetQuadrature(fe,&q);CHKERRQ(ierr);
    /* ierr = PetscQuadratureGetOrder(q,&order);CHKERRQ(ierr); */
    /* Creating BD FE */
    ierr = PetscFECreateDefault(dm,dim-1, dim, simplex, "bd_def_",PETSC_DEFAULT,&fe_bd);CHKERRQ(ierr);
    ierr = PetscFESetQuadrature(fe_bd,q);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) fe_bd, "deformation");CHKERRQ(ierr);
  }
  if(user->transient){
    if (user->ldis) PetscPrintf(PETSC_COMM_WORLD,"Large displacement formulation has no option for transient as of yet");
    ierr = PetscFEGetQuadrature(fe,&q);CHKERRQ(ierr);
    ierr = PetscQuadratureGetOrder(q,&order);CHKERRQ(ierr);
    /* Creating the FE field for the velocity */
    ierr = PetscFECreateDefault(dm,dim,dim,simplex,"vel_",PETSC_DEFAULT,&fe_vel);CHKERRQ(ierr);
    ierr = PetscFESetQuadrature(fe_vel,q);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) fe_vel, "velocity");CHKERRQ(ierr);
    ierr = PetscFECreateDefault(dm,dim-1, dim, simplex, "bd_vel_",PETSC_DEFAULT,&fe_vel_bd);CHKERRQ(ierr);
    ierr = PetscFESetQuadrature(fe_vel_bd,q);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) fe_vel_bd, "velocity");CHKERRQ(ierr);
  }
  /* Discretization and boundary conditons: */
  while (cdm)  {
    ierr = DMGetDS(cdm, &prob);CHKERRQ(ierr);
    ierr = PetscDSSetDiscretization(prob, 0, (PetscObject) fe); CHKERRQ(ierr);
    if (user->transient){
      ierr = PetscDSSetDiscretization(prob,1, (PetscObject) fe_vel);CHKERRQ(ierr);
      ierr = PetscDSSetBdDiscretization(prob,1, (PetscObject) fe_vel_bd);CHKERRQ(ierr);
    }
    if (user->neumann){
      ierr = PetscDSSetBdDiscretization(prob, 0, (PetscObject) fe_bd); CHKERRQ(ierr);
    }
    ierr = SetupProblem(cdm, user); CHKERRQ(ierr);
    /* Define the boundaries */
    
    const PetscInt Ncomp = dim;
    PetscInt components[dim];
    const PetscInt Nfid = 2;
    PetscInt d;
    PetscInt fid[Nfid];  /* fixed faces [numer of fixed faces] */
    const PetscInt Npid = 1;  /* number of pressure loaded faces */
    PetscInt pid[Npid];       /* the ids of the pressure loaded faces */
    const PetscInt Nsxid = 1; /* number of face ids for the symmetry in x
                               * direction bc */
    PetscInt sxid[Nsxid];
    const PetscInt Nszid = 1;  /* Number of face ids for the symmetriy in z
                                * direction (x-y plane) BC */
    PetscInt szid[Nszid];
    const PetscInt Nsyid = 1;
    PetscInt syid[Nsyid];

    PetscInt restrictX[] = {0};  /* array of constricted components for the
                                  * dmaddboundary routine */
    PetscInt restrictZ[] = {2}; 
    PetscInt restrictY[] = {1};
    
    for (d=0;d<dim;d++){
      components[d]=d;
    }
    if(dim==2){
      fid[0] = 5; /* {5}; */ 	/* The fixed faces.*/
      pid[0] = 4; /* {4}; */ 	/* The pressure loaded faces */
    }
    else if(dim==3){
      fid[0] = 10;  /* The fixed face */
      fid[1] = 14;  /* the second fixed face */
      pid[0]= 2;  /* The pressure loaded faces */
      sxid[0] = 4;
      szid[0] = 1;
      syid[0] = 3;  /*  */
    }

    /* ierr =  DMAddBoundary(cdm, PETSC_TRUE, "fixed", "Face Sets",0, Ncomp, components, (void (*)()) zero_vector, Nfid, fid, user);CHKERRQ(ierr); */
    /* This part is to impose the boundary conditions related to the
     * rotational symmetry. Currently works only with a pi/2 geometry. We
     * impose 0 displacement in the face normal of the fundamental planes (x-y
     * and z-y in this particular case) */
    ierr = DMAddBoundary(cdm, PETSC_TRUE, "symmx", "Face Sets",0,1,restrictX, (void (*)()) zero_scalar, Nsxid, sxid, user);CHKERRQ(ierr);
    ierr = DMAddBoundary(cdm, PETSC_TRUE, "symmy", "Face Sets",0,1,restrictY, (void (*)()) zero_scalar, Nsyid, syid, user);CHKERRQ(ierr);
    ierr = DMAddBoundary(cdm, PETSC_TRUE, "bottom", "Face Sets",0,1,restrictZ, (void (*)()) zero_scalar, Nszid, szid, user);CHKERRQ(ierr);

    /* Do the same thing for the velocity field */
    /* ierr = DMAddBoundary(cdm, PETSC_TRUE, "symmx", "Face Sets",1,1,restrictX, (void (*)()) zero_scalar, Nsxid, sxid, user);CHKERRQ(ierr); */
    /* ierr = DMAddBoundary(cdm, PETSC_TRUE, "symmy", "Face Sets",1,1,restrictY, (void (*)()) zero_scalar, Nsyid, syid, user);CHKERRQ(ierr); */
    /* ierr = DMAddBoundary(cdm, PETSC_TRUE, "bottom", "Face Sets",1,1,restrictZ, (void (*)()) zero_scalar, Nszid, szid, user);CHKERRQ(ierr); */
    
    if(user->neumann){
      ierr = DMAddBoundary(cdm, PETSC_FALSE, "load", "Face Sets",0, Ncomp, components, NULL, Npid, pid, user);CHKERRQ(ierr);
      /* ierr = DMAddBoundary(cdm, PETSC_TRUE, "constrict", "Face Sets", 0, 1, test ,(void (*)()) zero_scalar, Nfid, testid, user);CHKERRQ(ierr); /\* constrict movement of roght end *\/ */
    }
    else{
      ierr = DMAddBoundary(cdm, PETSC_TRUE, "load", "Face Sets", 0, Ncomp, components, (void(*)()) pull, Npid, pid, user);CHKERRQ(ierr);
    }
    ierr = DMGetCoarseDM(cdm, &cdm); CHKERRQ(ierr);
  }

  ierr = PetscFEGetBasisSpace(fe,&space);CHKERRQ(ierr);
  ierr = PetscSpaceGetOrder(space,&BSpaceOrder);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"The order of the basis space is %i \n", BSpaceOrder);


  {
    /* Set up the near null space (a.k.a. rigid body modes) that will be used
       by the multigrid preconditioner */
    
     /* DM           subdm; */
     /* MatNullSpace nearNullSpace; */
     /* PetscInt     fields = 0; */
     /* PetscObject  deformation, velocity; */
     /* DMCreateSubDM(dm, 1, &fields, NULL, &subdm); */
     /* DMPlexCreateRigidBody(subdm, &nearNullSpace); */
     /* DMGetField(dm, 0, &deformation); */
     /* PetscObjectCompose(deformation, "nearnullspace", (PetscObject) nearNullSpace); */
     /* DMDestroy(&subdm); */
     /* MatNullSpaceDestroy(&nearNullSpace); */

     /*  and now for the velocity */
     /* fields = 1; */
     /* ierr = DMCreateSubDM(dm, 1, &fields, NULL, &subdm);CHKERRQ(ierr); */
     /* ierr = DMPlexCreateRigidBody(subdm, &nearNullSpace);CHKERRQ(ierr); */
     /* ierr = DMGetField(dm, 1, &velocity);CHKERRQ(ierr); */
     /* ierr = PetscObjectCompose(velocity, "nearnullspace", (PetscObject) nearNullSpace);CHKERRQ(ierr); */
     /* ierr = DMDestroy(&subdm);CHKERRQ(ierr); */
     /* ierr = MatNullSpaceDestroy(&nearNullSpace);CHKERRQ(ierr); */

   }

  ierr = PetscFEDestroy(&fe); CHKERRQ(ierr);
  if (user->neumann)  ierr = PetscFEDestroy(&fe_bd); CHKERRQ(ierr);
  if (user->transient){
    ierr = PetscFEDestroy(&fe_vel); CHKERRQ(ierr);
    ierr = PetscFEDestroy(&fe_vel_bd);CHKERRQ(ierr);
  }
  return(0);
}

PetscErrorCode registerSimmerData(PetscReal PK[])
{
  int i;
  PetscErrorCode ierr;
  AppCtx user;

  /* Simmer's unit for pressure is Pa but we use GPA so:*/
  for (i=0;i<SimmerData.MMS;i++){
    SimmerData.PK[i] = PK[i]*1E-9;
  }
  return(0);
}


PetscErrorCode cimplysetup(PetscErrorCode ierr){

   
  /* Firing up Petsc */
  /* ierr= PetscInitialize(&argc, &argv,NULL,help);CHKERRQ(ierr); */
  ierr = PetscInitializeFortran();CHKERRQ(ierr);
  ierr = ProcessOptions(PETSC_COMM_WORLD,&user);CHKERRQ(ierr);

  /* Creating the hdf5 viewer */
  /* ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD,"sol.h5",FILE_MODE_WRITE,&user.viewer);CHKERRQ(ierr); */

  
  /* Reading the sim05 file */
  ierr = sim05tocimply();CHKERRQ(ierr);
  /* importing the gmsh file. Take note that only simplices give meaningful results in 2D at the moment (For which ever reasons) */

  ierr = DMPlexCreateFromFile(PETSC_COMM_WORLD,"SHammer_thin.msh", user.interpolate,&dm);CHKERRQ(ierr);

  ierr = DMGetDimension(dm,&user.dim); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"The problem dimension is %i \n",user.dim);
  ierr = DMPlexDistribute(dm,0,NULL,&distributeddm); CHKERRQ(ierr);
  if (distributeddm) {
    ierr=DMDestroy(&dm);CHKERRQ(ierr);
    dm = distributeddm;
  }
  
  ierr = DMSetFromOptions(dm);CHKERRQ(ierr);
  ierr = DMView(dm,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  /* In case of a transient analysis */
  
  if (user.transient){

    PetscPrintf(PETSC_COMM_WORLD,"starting transient analysis. Total time: %g s \n",user.time.total);
    ierr =  TSCreate(PETSC_COMM_WORLD, &ts); CHKERRQ(ierr);
    ierr = TSSetType(ts, TSBEULER); CHKERRQ(ierr);
    ierr = TSSetDM(ts,dm); CHKERRQ(ierr);
    ierr = DMSetApplicationContext(dm, &user);CHKERRQ(ierr);
    ierr = SetupDiscretization(dm,&user);CHKERRQ(ierr);
    ierr = DMPlexCreateClosureIndex(dm,NULL);CHKERRQ(ierr);

    ierr = TSMonitorSet(ts,Monitor,&user,NULL);CHKERRQ(ierr);
    ierr = DMTSSetBoundaryLocal(dm, DMPlexTSComputeBoundary, &user); CHKERRQ(ierr);
    ierr = DMTSSetIFunctionLocal(dm, DMPlexTSComputeIFunctionFEM, &user); CHKERRQ(ierr);
    ierr = DMTSSetIJacobianLocal(dm, DMPlexTSComputeIJacobianFEM, &user); CHKERRQ(ierr);


    
    ierr = DMCreateGlobalVector(dm, &u); CHKERRQ(ierr);

    ierr = VecSet(u,(PetscReal) 0.0); CHKERRQ(ierr);


    ierr = TSSetSolution(ts,u); CHKERRQ(ierr);
    
    ierr = TSSetDuration(ts,1000,user.time.total); CHKERRQ(ierr);
    ierr = TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP); CHKERRQ(ierr);
    ierr = TSSetEquationType(ts,TS_EQ_IMPLICIT);CHKERRQ(ierr);
    ierr = TSSetInitialTimeStep(ts,0.0, user.time.dt); CHKERRQ(ierr);
    ierr = TSSetFromOptions(ts); CHKERRQ(ierr);
    ierr = TSGetSNES(ts, &snes); CHKERRQ(ierr);
    ierr = SNESGetKSP(snes, &ksp); CHKERRQ(ierr);

    /* ierr = CreateVelocityNullSpace(dm,&user,NULL,&nullspace); CHKERRQ(ierr); */
    /* /\* ierr = KSPSetNullSpace(ksp,nullspace);CHKERRQ(ierr); *\/ */
    /* ierr = KSPGetOperators(ksp,&A,&P);CHKERRQ(ierr); */
    /* ierr = MatSetType(A,MATAIJ);CHKERRQ(ierr); */
    /* ierr = MatSetNullSpace(A,nullspace); CHKERRQ(ierr); */
    /* ierr = KSPSetOperators(ksp,A,P);CHKERRQ(ierr); */

    /* ierr = MatNullSpaceDestroy(&nullspace);CHKERRQ(ierr); */

      /* ierr = TSGetSolveTime(ts,&ftime); */
    /* ierr = TSGetTimeStepNumber(ts,&nsteps); CHKERRQ(ierr); */
    /* ierr = TSGetConvergedReason(ts,&ConvergedReason); CHKERRQ(ierr); */


    /* ierr = PetscPrintf(PETSC_COMM_WORLD, "%s at time %g after %D steps \n",ConvergedReason,ftime,nsteps); CHKERRQ(ierr); */ 
  }
  else{
    ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
    ierr = SNESSetDM(snes,dm);CHKERRQ(ierr);
    ierr = DMSetApplicationContext(dm, &user);CHKERRQ(ierr);
    ierr = SetupDiscretization(dm,&user);CHKERRQ(ierr);
    ierr = DMPlexCreateClosureIndex(dm,NULL);CHKERRQ(ierr);
    
    ierr = DMCreateGlobalVector(dm,&u);CHKERRQ(ierr);
    ierr = VecDuplicate(u,&r);CHKERRQ(ierr);
    
    ierr = VecSet(u,(PetscReal) 0.0);CHKERRQ(ierr);
    
    
    ierr = DMSetMatType(dm, MATAIJ);CHKERRQ(ierr);
    ierr = DMCreateMatrix(dm, &J);CHKERRQ(ierr);
    A=J;
    
    ierr = DMPlexSetSNESLocalFEM(dm,&user,&user,&user);CHKERRQ(ierr);
    
    ierr = SNESSetJacobian(snes, A, J, NULL, NULL);CHKERRQ(ierr);
    if (user.debug){  		/* Showing the Jacobi matrix for debugging purposes */
      ierr = PetscPrintf(PETSC_COMM_WORLD,"The Jacobian for the nonlinear solver ( and also the preconditioning matrix)\n");CHKERRQ(ierr);
      ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    }
    
    ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
    
    
    if (user.showInitial){ ierr = DMVecViewLocal(dm, u, PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);}
    
    if (user.debug) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"initial guess \n");CHKERRQ(ierr);
      ierr = VecView(u, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    }
  }
  return(0);
}
PetscErrorCode cimplysolve(int MMS, PetscReal  PK[],int iter, PetscErrorCode ierr){

  PetscInt its;
  PetscInt i;
  PetscDS prob;


    /* Just a test if MMS from SIMMER and MMS read by the sim05reader are the same */
  if(MMS != SimmerData.MMS){
    ierr = PetscPrintf(PETSC_COMM_WORLD,"K-cell MMS :\t  %i \n sim05 MMS: \t %i \n",MMS,SimmerData.MMS);
    /* ierr=PetscErrorPrintf("The K-cell MMS received from SIMMER and the one read from the sim05 file do not have the same value.");CHKERRQ(ierr); */
  }



  
  /* filling the global SimmerData structure */
  if (user.verbose) PetscPrintf(PETSC_COMM_WORLD,"Received Data from Simmer. Creating cimply PK array... \n");
  ierr = registerSimmerData(PK);CHKERRQ(ierr);
  if (user.verbose) {
    PetscPrintf(PETSC_COMM_WORLD,"Creation of Simmer Data object successfull. \n");
    for (i=0;i<MMS;i++){
      PetscPrintf(PETSC_COMM_WORLD,"PK [%i] \t %f \n", i, SimmerData.PK[i]);
    }
  }
  
  /* blablablablablabla */
    /* Setting the Neumann Boudnary Condition */
  /* if (user.neumann){ */
  /*   ierr = DMGetDS(dm,&prob);CHKERRQ(ierr); */
  /*   if (user.ldis){ */
  /*     ierr = PetscDSSetBdResidual(prob,0,f0_u_bd_ldis,f1_u_bd);CHKERRQ(ierr); */
  /*   } */
  /*   else{ */
  /*     ierr = PetscDSSetBdResidual(prob,0,f0_u_bd,f1_u_bd);CHKERRQ(ierr); */
  /*     /\* ierr = PetscDSSetBdJacobian(prob,0,0,NULL,g1_uu_bd_2d,NULL,NULL);CHKERRQ(ierr); *\/ */
  /*   } */
    
  /* } */

  /* blablablablabal */


  if (user.verbose) PetscPrintf(PETSC_COMM_WORLD,"PK array as registered in SimmerData: %f \n",SimmerData.PK);

  if (user.transient){
    ierr = TSSolve(ts,NULL); CHKERRQ(ierr);
    ierr = TSView(ts,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"Transient analysis finished. \n");
  }
  else{
    ierr = SNESSolve(snes,NULL,u);CHKERRQ(ierr);
    ierr = SNESGetIterationNumber(snes, &its);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of snes iterations %i \n", its);CHKERRQ(ierr);
  }

  if (user.visualization == SOLUTION) {
    char filename[50];

    sprintf(filename,"solution%i.vtk",iter);
    ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    if (user.verbose) PetscPrintf(PETSC_COMM_WORLD,"Writing the solution to the vtk file. \n");
    ierr = PetscObjectSetName((PetscObject) u,"deformation");CHKERRQ(ierr);
    ierr = VecView(u,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    if (user.verbose)PetscPrintf(PETSC_COMM_WORLD,"Done creating the VTK file. \n");
  }
  
  return(0);
}

PetscErrorCode cimplyfinalize(PetscErrorCode ierr)
{
  ierr = VecViewFromOptions(u,NULL,"-sol_vec_view");CHKERRQ(ierr);
  
  if ( A != J){MatDestroy(&A);}
  MatDestroy(&J);
  VecDestroy(&u);
  VecDestroy(&r);
  if(user.transient) TSDestroy(&ts);
  if(!user.transient) SNESDestroy(&snes);
  DMDestroy(&dm);
  PetscViewerDestroy(&user.viewer);
  return(0);
}


/* int main(int argc, char **argv){ */

/*   SNES snes;			/\* nonlinear solver *\/ */
/*   DM dm, distributeddm;			/\* problem definition *\/ */
/*   Vec u,r;			/\* solution and residual vectors *\/ */
/*   Mat A,J;			/\* Jacobian Matrix *\/ */
/*   AppCtx user;			/\* user-defined work context *\/ */
/*   PetscErrorCode ierr; */
/*   PetscViewer viewer; */
/*   TS ts; */




  
/*   /\* Firing up Petsc *\/ */
/*   ierr= PetscInitialize(&argc, &argv,NULL,help);CHKERRQ(ierr); */

/*   ierr = ProcessOptions(PETSC_COMM_WORLD,&user);CHKERRQ(ierr); */


/*   /\* ierr = CreateMesh(PETSC_COMM_WORLD,&user,&dm);CHKERRQ(ierr); *\/ */

/*   /\* importing the gmsh file. Take note that only simplices give meaningful results in 2D at the moment (For which ever reasons) *\/ */

/*   /\* testmesh_2D_box_quad.msh *\/ */
/*   /\* Beam_coarse.msh *\/ */
/*   ierr = DMPlexCreateFromFile(PETSC_COMM_WORLD,"longbeam3D.msh", PETSC_TRUE,&dm);CHKERRQ(ierr); */
/*   ierr = DMGetDimension(dm,&user.dim); CHKERRQ(ierr); */
/*   PetscPrintf(PETSC_COMM_WORLD,"The problem dimension is %i \n",user.dim); */
/*   ierr = DMPlexDistribute(dm,0,NULL,&distributeddm); CHKERRQ(ierr); */
/*   if (distributeddm) { */
/*     ierr=DMDestroy(&dm);CHKERRQ(ierr); */
/*     dm = distributeddm; */
/*   } */
  
/*   ierr = DMSetFromOptions(dm);CHKERRQ(ierr); */
/*   ierr = DMView(dm,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr); */
 




/*     /\* Write solution to a VTK file *\/ */
/*     ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD,"solution.vtk",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr); */
/*     ierr = DMPlexVTKWriteAll((PetscObject) dm, viewer); */
/*     ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr); */

/*     ierr = TSDestroy(&ts); CHKERRQ(ierr); */
/*     ierr = VecDestroy(&u); CHKERRQ(ierr); */
    
/*   } */
/*   else{ */
/*     PetscInt its;  /\* number of iterations needed by the SNES solver *\/ */
/*     ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr); */
/*     ierr = SNESSetDM(snes,dm);CHKERRQ(ierr); */
/*     ierr = DMSetApplicationContext(dm, &user);CHKERRQ(ierr); */
/*     ierr = SetupDiscretization(dm,&user);CHKERRQ(ierr); */
/*     ierr = DMPlexCreateClosureIndex(dm,NULL);CHKERRQ(ierr); */
    
/*     ierr = DMCreateGlobalVector(dm,&u);CHKERRQ(ierr); */
/*     ierr = VecDuplicate(u,&r);CHKERRQ(ierr); */
    
/*     ierr = VecSet(u,(PetscReal) 0.0);CHKERRQ(ierr); */
    
    
/*     ierr = DMSetMatType(dm, MATAIJ);CHKERRQ(ierr); */
/*     ierr = DMCreateMatrix(dm, &J);CHKERRQ(ierr); */
/*     A=J; */
    
/*     ierr = DMPlexSetSNESLocalFEM(dm,&user,&user,&user);CHKERRQ(ierr); */
    
/*     ierr = SNESSetJacobian(snes, A, J, NULL, NULL);CHKERRQ(ierr); */

/*     if (user.debug){  		/\* Showing the Jacobi matrix for debugging purposes *\/ */
/*       ierr = PetscPrintf(PETSC_COMM_WORLD,"The Jacobian for the nonlinear solver ( and also the preconditioning matrix)\n");CHKERRQ(ierr); */
/*       ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr); */
/*     } */
    
/*     ierr = SNESSetFromOptions(snes);CHKERRQ(ierr); */
    
/*     if (user.showInitial){ ierr = DMVecViewLocal(dm, u, PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);} */
    
/*     if (user.debug) { */
/*       ierr = PetscPrintf(PETSC_COMM_WORLD,"initial guess \n");CHKERRQ(ierr); */
/*       ierr = VecView(u, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr); */
/*     } */
/*     ierr = SNESSolve(snes,NULL,u);CHKERRQ(ierr); */
/*     ierr = SNESGetIterationNumber(snes, &its);CHKERRQ(ierr); */
/*     ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of snes iterations %i \n", its);CHKERRQ(ierr); */
  
/*     if (user.showSolution){ */
/*       ierr = PetscPrintf(PETSC_COMM_WORLD,"solution: \n");CHKERRQ(ierr); */
/*       ierr = VecChop(u, 3.0e-9); CHKERRQ(ierr); /\* what does vecchop do? *\/ */
/*       ierr = VecView(u, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr); */
/*     } */
    
/*     if (user.visualization ==SOLUTION){ */
        
/*       if (user.verbose)PetscPrintf(PETSC_COMM_WORLD,"Creating the vtk output file... \n"); */
/*       ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD,"solution.vtk",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr); */
/*       ierr = PetscObjectSetName((PetscObject) u,"deformation");CHKERRQ(ierr); */
/*       ierr = VecView(u,viewer);CHKERRQ(ierr); */
/*       ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr); */
/*       if (user.verbose)PetscPrintf(PETSC_COMM_WORLD,"Done creating the VTK file. \n"); */
/*     } */

/*     ierr = VecViewFromOptions(u,NULL,"-sol_vec_view");CHKERRQ(ierr); */

/*     if ( A != J){MatDestroy(&A);} */
/*     MatDestroy(&J); */
/*     VecDestroy(&u); */
/*     VecDestroy(&r); */
/*     SNESDestroy(&snes); */
/*   } */
/*   DMDestroy(&dm); */
/*   PetscFinalize(); */
/*   return 0; */
  
/* } */
