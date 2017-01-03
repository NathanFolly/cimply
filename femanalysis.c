#include "femanalysis.h"

static void * FEMAnalysis_settimestep(void * _self, PetscReal dt);

static PetscErrorCode SetupProblem(DM dm, AppCtx *user);

static PetscErrorCode SetupDiscretization(DM dm, AppCtx *user);

static PetscErrorCode Monitor(TS ts, PetscInt step, PetscReal crtime, Vec u, void *ctx)

static void * setuptransient(void * _self);

static void * setupstationary(void * _self);

static PetscErrorCode ProcessOptions(MPI_Comm comm, AppCtx *options);

static PetscErrorCode zero_vector(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx);

static PetscErrorCode constrict_y(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx);

static PetscErrorCode constrict_x(PetscInt dim, PetscReal time, const PetscReal x[], const PetscReal n[], PetscInt Nf, PetscScalar *u, void *ctx);

static PetscErrorCode constrict_z(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx);

static PetscErrorCode coordinates(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx);

static PetscErrorCode pull(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx);



static void * FEMAnalysis_ctor(void * _self, va_list * app){
  struct FEMAnalysis * self = _self;
  char meshfilename[51];
  double endtime = 0.0;
  PetscPartitioner part;
  PetscMPIInt rank, numProcs;
  DM distributeddm;
  PetscErrorCode ierr;
  
  self->user->time->iter=0;
  self->settimestep=FEMAnalysis_settimestep;
  
  meshfilename = va_arg(*app, char *);
  endtime = va_arg(*app, double);
  self->user->time->total = endtime;

  ierr = ProcessOptions(PETSC_COMM_WORLD,&(self->user));CHKERRQ(ierr);
  /* importing the gmsh file. Take note that only simplices give meaningful results in 2D at the moment (For which ever reasons) */
  ierr = DMPlexCreateFromFile(PETSC_COMM_WORLD,meshfilename, PETSC_TRUE,&(self->dm));CHKERRQ(ierr);
  ierr = DMGetDimension(self->dm,&(self->dim)); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"The problem dimension is %i \n",self->dim);

  /* distribute the mesh */
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &numProcs);
  DMPlexGetPartitioner(self->dm, &part);
  PetscPartitionerSetType(part, PETSCPARTITIONERPARMETIS);
  /* PetscPartitionerShellSetPartition(part, numProcs, NULL, NULL); */
  ierr = DMPlexDistribute(self->dm,0,NULL,&distributeddm); CHKERRQ(ierr);
  if (distributeddm) {
    ierr=DMDestroy(&(self->dm));CHKERRQ(ierr);
    self->dm = distributeddm;
    ierr = DMDestroy(&distributeddm);
  } 
  /* ierr = DMSetFromOptions(dm);CHKERRQ(ierr); */
  ierr = DMView(dm,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  if(self->user->transient) setuptransient(self);
  else setupstationary(self);

  return self;
}

static void * FEMAnalysis_dtor(void * _self){
  struct FEMAnalysis * self = _self;
  return self;
  /* TODO: clean up the whole petsc junk to avoid memory leaks */
}

static void * FEMAnalysis_update(void * _self){
  struct FEMAnalysis * self = _self;
  PetscErrorCode ierr;

   if (self->user->transient){

    if(self->user->runType==RUN_STANDALONE){
      ierr = TSSolve(self->ts,NULL); CHKERRQ(ierr);
      ierr = TSView(self->ts,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD,"Transient analysis finished. \n");
    }
    if(self->user->runType==RUN_COUPLED){
      char filename[50];
      ierr = TSSetTimeStep(self->ts,self->user->time->dt);
      ierr = TSStep(self->ts);

      sprintf(filename,"solution%i.vtu",self->user->time->iter);
      ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&(self->viewer));CHKERRQ(ierr);
      if (self->user->verbose) PetscPrintf(PETSC_COMM_WORLD,"Writing the solution to the vtk file. \n");
      /* ierr = PetscObjectSetName((PetscObject) u,"deformation");CHKERRQ(ierr); */
      ierr = VecView(self->u,self->viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&(self->viewer));CHKERRQ(ierr);
      if (self->user->verbose)PetscPrintf(PETSC_COMM_WORLD,"Done creating the VTK file. \n");
      ++ self->user->time->iter;
    }
  }
  else{
    PetscInt its;
    ierr = SNESSolve(self->snes,NULL,self->u);CHKERRQ(ierr);
    ierr = SNESGetIterationNumber(self->snes, &its);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of snes iterations %i \n", its);CHKERRQ(ierr);
  }
  
}



static const struct Class _FEMAnalysis = {sizeof(struct FEMAnalysis), FEMAnalysis_ctor, FEMAnalysis_dtor,FEMAnalysis_update};

const void * FEMAnalysis = &_FEMAnalysis;





/* ---------------------------- class specific functions  -------------------- */


static void * FEMAnalysis_settimestep(void * _self, PetscReal dt){
  struct FEMAnalysis * self = _self;

  self->user->time->dt = dt;

  return 0;
}

void * settimestep(void * _self, PetscReal dt){
  struct Class ** cp = _self;
  if(*cp != FEMAnalysis){
    fprintf(stderr,"ERROR:: error in settimestep. First argument must be of type FEMAnalysis.\n");
    return 0;
  }
  struct FEMAnalysis * self = _self;

  self->settimestep(self,dt);
  return 0;
}


static PetscErrorCode SetupProblem(DM dm, AppCtx *user)
{
  PetscDS prob; 		   /* a DS is a "discrete system" managing discretisations */
  PetscErrorCode ierr;
  PetscInt NumFields;


  ierr = DMGetDS(dm, &prob);CHKERRQ(ierr);


  
     /* Setting up the linear elasticity problem */
    if (user->transient){
      ierr = PetscDSSetResidual(prob,0,f0_u_transient,f1_u);CHKERRQ(ierr);
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
      ierr = PetscDSSetDynamicJacobian(prob,0,1,g0_uvel,NULL,NULL,NULL);CHKERRQ(ierr);
      ierr = PetscDSSetBdResidual(prob,1,f0_vel_bd,f1_vel_bd);CHKERRQ(ierr);
    }
  
  /* Setting the Neumann Boudnary Condition */
  if (user->neumann){
    if (user->ldis){
      ierr = PetscDSSetBdResidual(prob,0,f0_u_bd_ldis,f1_u_bd);CHKERRQ(ierr);
    }
    else{
      ierr = PetscDSSetBdResidual(prob,0,f0_u_bd,f1_u_bd);CHKERRQ(ierr);

    }
    
  }

  
  ierr = PetscDSGetNumFields(prob,&NumFields);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Problem set up with %i fields\n",NumFields);CHKERRQ(ierr);
      
  return(0);
}


static PetscErrorCode SetupDiscretization(DM dm, AppCtx *user){

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
  
  if(user->transient){
    if (user->ldis) PetscPrintf(PETSC_COMM_WORLD,"Large displacement formulation has no option for transient as of yet");
    ierr = PetscFEGetQuadrature(fe,&q);CHKERRQ(ierr);
    ierr = PetscQuadratureGetOrder(q,&order);CHKERRQ(ierr);
    /* Creating the FE field for the velocity */
    ierr = PetscFECreateDefault(dm,dim,dim,simplex,"vel_",PETSC_DEFAULT,&fe_vel);CHKERRQ(ierr);
    ierr = PetscFESetQuadrature(fe_vel,q);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) fe_vel, "velocity");CHKERRQ(ierr);
  }
  /* Discretization and boundary conditons: */
  while (cdm)  {
    ierr = DMGetDS(cdm, &prob);CHKERRQ(ierr);
    ierr = PetscDSSetDiscretization(prob, 0, (PetscObject) fe); CHKERRQ(ierr);
    if (user->transient){
      ierr = PetscDSSetDiscretization(prob,1, (PetscObject) fe_vel);CHKERRQ(ierr);
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

    ierr = DMAddBoundary(cdm, PETSC_TRUE, "symmx", "Face Sets",0,1,restrictX, (void (*)()) zero_scalar, Nsxid, sxid, user);CHKERRQ(ierr);
    ierr = DMAddBoundary(cdm, PETSC_TRUE, "symmy", "Face Sets",0,1,restrictY, (void (*)()) zero_scalar, Nsyid, syid, user);CHKERRQ(ierr);
    ierr = DMAddBoundary(cdm, PETSC_TRUE, "bottom", "Face Sets",0,1,restrictZ, (void (*)()) zero_scalar, Nszid, szid, user);CHKERRQ(ierr);

    if(user->neumann){
      ierr = DMAddBoundary(cdm, PETSC_FALSE, "load", "Face Sets",0, Ncomp, components, NULL, Npid, pid, user);CHKERRQ(ierr);
     }
    else{
      ierr = DMAddBoundary(cdm, PETSC_TRUE, "load", "Face Sets", 0, Ncomp, components, (void(*)()) pull, Npid, pid, user);CHKERRQ(ierr);
    }
    ierr = DMGetCoarseDM(cdm, &cdm); CHKERRQ(ierr);
  }

  ierr = PetscFEGetBasisSpace(fe,&space);CHKERRQ(ierr);
  ierr = PetscSpaceGetOrder(space,&BSpaceOrder);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"The order of the basis space is %i \n", BSpaceOrder);

  ierr = PetscFEDestroy(&fe); CHKERRQ(ierr);
  if (user->transient){
    ierr = PetscFEDestroy(&fe_vel); CHKERRQ(ierr);
  }
  return(0);
}

static PetscErrorCode Monitor(TS ts, PetscInt step, PetscReal crtime, Vec u, void *ctx)
{
  AppCtx *user=(AppCtx*) ctx;
  PetscErrorCode ierr;
  PetscViewer viewer;
  char filename[50];
 
  sprintf(filename,"solution%i.vtu",step);
  ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  if (user->verbose) PetscPrintf(PETSC_COMM_WORLD,"Writing the solution to the vtk file. \n");
  ierr = VecView(u,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  if (user->verbose)PetscPrintf(PETSC_COMM_WORLD,"Done creating the VTK file. \n");
    
  return(0);
}

static void * setuptransient(void * _self){
  struct FEMAnalysis * self = _self;
  PetscErrorCode ierr;

  PetscPrintf(PETSC_COMM_WORLD,"starting transient analysis. Total time: %g s \n",user.time.total);
  ierr = TSCreate(PETSC_COMM_WORLD, &(self->ts)); CHKERRQ(ierr);
  ierr = TSSetType(self->ts, TSBEULER); CHKERRQ(ierr);
  ierr = TSSetDM(self->ts,self->dm); CHKERRQ(ierr);
  ierr = DMSetApplicationContext(self->dm, &(self->user));CHKERRQ(ierr);
  ierr = SetupDiscretization(self->dm,&(self->user));CHKERRQ(ierr);
  ierr = DMPlexCreateClosureIndex(self->dm,NULL);CHKERRQ(ierr);
  
  ierr = TSMonitorSet(self->ts,Monitor,&(self->user),NULL);CHKERRQ(ierr);
  ierr = DMTSSetBoundaryLocal(self->dm, DMPlexTSComputeBoundary, &(self->user)); CHKERRQ(ierr);
  ierr = DMTSSetIFunctionLocal(self->dm, DMPlexTSComputeIFunctionFEM, &(self->user)); CHKERRQ(ierr);
  ierr = DMTSSetIJacobianLocal(self->dm, DMPlexTSComputeIJacobianFEM, &(self->user)); CHKERRQ(ierr);
  

  ierr = DMSetVecType(self->dm,VECMPI);CHKERRQ(ierr);    
  ierr = DMCreateGlobalVector(self->dm, &(self->u)); CHKERRQ(ierr);
  
  ierr = TSSetSolution(self->ts,self->u); CHKERRQ(ierr);
  ierr = TSSetDuration(self->ts,1000,self->user->time->total); CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(self->ts, TS_EXACTFINALTIME_MATCHSTEP); CHKERRQ(ierr);
  ierr = TSSetEquationType(self->ts,TS_EQ_IMPLICIT);CHKERRQ(ierr);
  ierr = TSSetInitialTimeStep(self->ts,0.0, user.time.dt); CHKERRQ(ierr);
  ierr = TSSetFromOptions(self->ts); CHKERRQ(ierr);
  ierr = TSGetSNES(self->ts, &(self->snes)); CHKERRQ(ierr);
  ierr = SNESGetKSP(self->snes, &(self->ksp)); CHKERRQ(ierr);

  return 0;
}

  
static void * setupstationary(void * _self){
  struct FEMAnalysis * self = _self;
  PetscErrorCode ierr;

  ierr = SNESCreate(PETSC_COMM_WORLD,&(self->snes));CHKERRQ(ierr);
  ierr = SNESSetDM(self->snes,self->dm);CHKERRQ(ierr);
  ierr = DMSetApplicationContext(self->dm, &(self->user));CHKERRQ(ierr);
  ierr = SetupDiscretization(self->dm,&(self->user));CHKERRQ(ierr);
  ierr = DMPlexCreateClosureIndex(self->dm,NULL);CHKERRQ(ierr);
  
  
  ierr = DMSetVecType(self->dm,VECMPI);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(self->dm,&(self->u));CHKERRQ(ierr);
  ierr = VecDuplicate(self->u,&(self->r));CHKERRQ(ierr);
  
  ierr = DMSetMatType(self->dm, MATMPIAIJ);CHKERRQ(ierr);
  ierr = DMCreateMatrix(self->dm, &(self->J));CHKERRQ(ierr);
  A=J;
  
  ierr = DMPlexSetSNESLocalFEM(self->dm,&(self->user),&(self->user),&(self->user));CHKERRQ(ierr);
  
  ierr = SNESSetJacobian(self->snes, self->A, self->J, NULL, NULL);CHKERRQ(ierr);
  if (self->user->debug){  		/* Showing the Jacobi matrix for debugging purposes */
    ierr = PetscPrintf(PETSC_COMM_WORLD,"The Jacobian for the nonlinear solver ( and also the preconditioning matrix)\n");CHKERRQ(ierr);
    ierr = MatView(self->A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }
  
  ierr = SNESSetFromOptions(self->snes);CHKERRQ(ierr);
  
    
  if (self->user->showInitial){ ierr = DMVecViewLocal(self->dm, self->u, PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);}
  
  if (self->user->debug) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"initial guess \n");CHKERRQ(ierr);
    ierr = VecView(self->u, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }


  return 0;
}



static PetscErrorCode ProcessOptions(MPI_Comm comm, AppCtx *options)
{
  const char *runTypes[2]={"standalone","coupled"};
  const char *visutype[4]={"none","labels","boundaries","solution"};
  PetscInt vis;
  PetscInt run;
  PetscErrorCode ierr;
  options->debug = 0;
  options->runType = RUN_COUPLED;
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
  ierr = PetscOptionsEList("-run-type","The run type: standalone or coupled","cimply.c",runTypes,2,runTypes[options->runType],&run,NULL);CHKERRQ(ierr);
  
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


/* ---------------------------- functions purely used for setting up the problem ----------------- */


static PetscErrorCode zero_scalar(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx)
{
  u[0]=0.0;
  return 0;
}

static PetscErrorCode zero_vector(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx)
{
  const PetscInt Ncomp = dim;
  PetscInt comp;
  for (comp=0;comp<Ncomp;comp++) u[comp]=0.0;
  return 0;
}

static PetscErrorCode constrict_y(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx)
{
  u[1] = 0.0;
  return 0;
}

static PetscErrorCode constrict_x(PetscInt dim, PetscReal time, const PetscReal x[], const PetscReal n[], PetscInt Nf, PetscScalar *u, void *ctx)
{
  u[0] = 0.0;
  return 0;
}


static PetscErrorCode constrict_z(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx)
{
  u[2] = 0.0;
  return 0;
}

static PetscErrorCode coordinates(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx)
{
  const PetscInt Ncomp = dim;
  PetscInt comp;
  for (comp=0; comp<Ncomp;comp++) u[comp]=x[comp];
  return 0;
}

static PetscErrorCode pull(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx)
{
  PetscInt comp;
  for(comp=0;comp<dim;comp++){
    if(comp==1){u[comp]=0.1;}
    else{u[comp]=0.0;}
  }
    return 0;
}

