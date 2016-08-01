static char help[] = "basic, linear elasticity problem solved with a nonlinear sovler. This version provides the possibility to impose natural (Neumann) boundary conditions.";


/* The PETSc packages we need: */
#include <petscdmplex.h>
#include <petscds.h>
#include <petscsnes.h>

/* We will create a user defined Application Context that helps to manage the options (or the contetext) of the program  */

typedef enum {RUN_FULL, RUN_TEST} RunType;
typedef enum {NONE, LABELS, BOUNDARIES, SOLUTION} Visualization;

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
}AppCtx;


/* The Kronecker Delta */
static const PetscInt delta2D[2*2] = {1,0,0,1};
static const PetscInt delta3D[3*3] = {1,0,0,0,1,0,0,0,1};
/* The Levi-Civita Symbol */
static const PetscInt epsilon2D[2*2] = {0,1,-1,0};

/* The static inline statement will tell the compiler to inline the function in the code directly instead of making multiple function calls. This improves performance especially for small functions that are called repeatedly throughout the program. */


PETSC_STATIC_INLINE void TensContrR44(PetscScalar C[], PetscScalar A[], PetscScalar B[], PetscInt ndim)
{				/* Tensor contraction for two rank 4 tensors
				   C=A:B  => C_ijkl = A_ijmn*B_nmkl*/
  PetscInt i,j,k,l,m,n;
  
  for (i=0;i<ndim;i++){
    for (j=0;j<ndim;j++){
      for (k=0;k<ndim;k++){
        for (l=0;l<ndim;l++){
          /* C[((i*ndim+j)*ndim+k)*ndim+l]=0; */
          for (m=0;m<ndim;m++){
            for (n=0;n<ndim;n++){
              C[((i*ndim+k)*ndim+j)*ndim+l]+=A[((i*ndim+j)*ndim+m)*ndim+n]*B[((n*ndim+m)*ndim+k)*ndim+l];
            }
          }
        }
      }
    }
  }
  
  
}

PETSC_STATIC_INLINE void TensContrR42(PetscScalar C[], PetscScalar A[],const  PetscScalar B[], PetscInt ndim)
{
  PetscInt i,j,k,l; /*Double contraction of a rank 4 and a rank 2 tensor*/
  
  for (i=0;i<ndim;i++){
    for (j=0;j<ndim;j++){
      for (k=0;k<ndim;k++){
        for (l=0;l<ndim;l++){
          C[(i*ndim+j)]=0.0;
        }
      }
    }
  }
  for (i=0;i<ndim;i++){
    for (j=0;j<ndim;j++){
      for (k=0;k<ndim;k++){
        for (l=0;l<ndim;l++){
          C[(i*ndim+j)]+= A[((i*ndim+j)*ndim+k)*ndim+l]*B[k*ndim+l];
        }
      }
    }
  }
}




PETSC_STATIC_INLINE void Det2D(PetscReal *detJ, PetscReal J[])
{  /* The determinant of a 2x2 matrix in vector form: {a_11, a_12, a_21, a_22} */
  *detJ = J[0]*J[3] - J[1]*J[2];
}

PETSC_STATIC_INLINE void Cof2D(PetscReal C[], PetscReal A[]) /* the inverse of a 2x2 matrix in vector form ( still needs to be multiplied by it's determinant) */
{							
  C[0] =  A[3];
  C[1] = -A[2];
  C[2] = -A[1];
  C[3] =  A[0];
}

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

PetscErrorCode coordinates(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx)
{
  const PetscInt Ncomp = dim;
  PetscInt comp;
  for (comp=0; comp<Ncomp;comp++) u[comp]=x[comp];
  return 0;
}

PetscErrorCode pull(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx)
{
  u[0] = 0.00;
  u[1] = 0.01;
  return 0;
}


void f0_u(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscScalar f0[])
{
  const PetscInt Ncomp = dim;
  PetscInt comp;
  for(comp=0;comp<Ncomp;comp++) f0[comp]=0.0;
  
}

void f1_u_2d_testytesty(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscScalar f1[])
{
  const PetscInt Ncomp = dim;
  const PetscReal mu = 86, lbda=115.4;
  
  /* This is a copy of the resdiual function in the original petsc ex77.c I
     made some changes to learn how the concept works. 
     In this case f1 is the cauchy stress tensor*/
  
  /*u_x = deformation gradient*/
  /* Hence is the strain epsilon_ij=0.5(u_i,j+u_j,i) => epsilon[comp*dim+d]=0.5(u_x[comp*dim+d]+u_x[d*dim+comp]) */
  PetscInt comp, d;
  for(comp=0;comp<Ncomp;comp++){
    for(d=0;d<dim;d++){
      f1[comp*dim+d]=mu*(u_x[comp*dim+d]+u_x[d*dim+comp]);
    }
    for(d=0;d<dim;d++){
      f1[comp*dim+comp]+=lbda*u_x[d*dim+d];
    }
  }
}


void g3_uu_2d(PetscInt dim, PetscInt Nf, PetscInt NfAux,
	      const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], 
	      const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
	      const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
	      const PetscScalar a_x[], PetscReal t, PetscReal u_tShift, const PetscReal x[],
	      PetscScalar g3[]){
  
  /* const PetscInt Ncomp = dim; */
  const PetscReal mu = 86, lbda=115.4;
  PetscInt i,j,k,l;
  
  /* g3 is the elasticity tensor */
  
  for (i=0;i<dim;i++){
    for (j=0;j<dim;j++){
      for (k=0; k < dim; k++){
        for (l=0;l<dim; l++){
          g3[((i*dim+j)*dim+k)*dim+l]=lbda*delta2D[i*dim+k]*delta2D[j*dim+l]+2*mu*delta2D[i*dim+k]*delta2D[j*dim+l]*delta2D[i*dim+j]+mu*(1-delta2D[i*dim+k])*(1-delta2D[j*dim+l]);
        }
      }
    }
  }
}

void f0_u_bd_2d(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
                const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
                const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
                const PetscScalar a_x[], PetscReal t, const PetscReal x[], const PetscReal n[],
                PetscScalar f0[])
{
  const PetscInt Ncomp=dim;
  /* Imposing a Neumann boundary condition is in other words, setting  the surface traction tensor eqal to the external load acting on the boundary in question  */
  const PetscScalar traction[] = {0.0, 0.01};
  PetscInt comp;
  
  for (comp=0; comp<Ncomp; ++comp){
    f0[comp] = traction[comp];
  }
}
    
void f1_u_bd_2d(PetscInt dim, PetscInt Nf, PetscInt NfAux,
        const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
        const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
        const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
        const PetscScalar a_x[], PetscReal t, const PetscReal x[], const PetscReal n[],
        PetscScalar f1[])
    {
    const  PetscInt Ncomp = dim;
    PetscInt comp,d;
    
    for (comp=0;comp<Ncomp;++comp){
      for (d = 0; d<dim;++d){
        f1[comp*dim + d] = 0.0;
      }
    }
  }
    
void g1_uu_bd_2d(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                 const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], 
                 const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
                 const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
                 const PetscScalar a_x[], PetscReal t, PetscReal u_tShift, const PetscReal x[],
                 const PetscReal n[], PetscScalar g1[]){
    
  const PetscReal mu=86, lbda=115.4;
  PetscInt Ncomp = dim;
  PetscInt i,j,k;
  
  /* Trying if we need a Jacobian on the boundary for the sovler to converge */
  /* The Jacobian in that case should be the surface traction */
  /* g1 needs to be a rank 3 tensor */
  /* should probably rethink the entire thing when converting to 3D */
    
  for (i=0; i<Ncomp; i++){
    for (j=0; j<Ncomp; j++){
      for (k=0; k<Ncomp; k++){
        g1[(i*Ncomp+j)*Ncomp+k]=n[i]*delta2D[j*dim+k]*lbda+(1-delta2D[j*dim+k])*mu*(n[k]*delta2D[i*dim+j]+n[j]*delta2D[i*dim+k])+n[i]*delta2D[j*dim+k]*delta2D[i*dim+j]*2*mu;
      }
    }
  }   
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
  options->dim = 2;
  options->interpolate = PETSC_FALSE;
  options->simplex = PETSC_TRUE;
  options->refinementLimit = 0.0;
  options->mu = 1;
  options->testPartition = PETSC_FALSE;
  options->showInitial = PETSC_FALSE;
  options->showSolution = PETSC_FALSE;
  options->verbose = PETSC_FALSE;
  options->visualization = NONE;
  options->neumann = PETSC_FALSE;

  ierr = PetscOptionsBegin(comm,"", "Linear elasticity problem options", "DMPLEX"); CHKERRQ(ierr);
  ierr = PetscOptionsInt("-debug","The debugging level","cimpleFEM.c",options->debug, &options->debug,NULL);CHKERRQ(ierr);
  run = options->runType;
  ierr = PetscOptionsEList("-run-type","The run type","cimpleFEM.c",runTypes,2,runTypes[options->runType],&run,NULL);CHKERRQ(ierr);
  
  options->runType = (RunType) run;

  ierr = PetscOptionsInt("-dim","The topological mesh dimension", "cimpleFEM.c",options->dim,&options->dim,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-interpolate","Generate intermediate mesh elements", "cimpleFEM.c",options->interpolate,&options->interpolate,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-simplex","use simplex or tensor product cells","cimpleFEM.c",options->simplex,&options->simplex,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-refinementLimit","Largest allowable cell volume", "cimpleFEM.c",options->refinementLimit,&options->refinementLimit,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-test_partition","Use a fixed partition for testing", "cimpleFEM.c", options->testPartition,&options->testPartition,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-shear_modulus","The shear modulus","cimpleFEM.c",options->mu,&options->mu,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-show_initial","Output the initial guess for verification","cimpleFEM.c",options->showInitial,&options->showInitial,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-show_solution","Output the solution for verification","cimpleFEM.c",options->showSolution,&options->showSolution,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-verbose","Output additional information.","cimpleFEM.c",options->verbose,&options->verbose,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEList("-visualize","Visualize a certain object for debugging and learning purposes","cimpleFEM",visutype,4,visutype[options->visualization],&vis,NULL);CHKERRQ(ierr);
  options->visualization = (Visualization) vis;
  ierr = PetscOptionsBool("-neumann", "Apply Neumann boundary conditions on the rightmost face.","cimpleFEM.c",options->neumann, &options->neumann,NULL);CHKERRQ(ierr);
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


PetscErrorCode SetupProblem(DM dm, AppCtx *user)
{
  PetscDS prob; 		   /* a DS is a "discrete system" managing discretisations */
  PetscErrorCode ierr;

  /* for the beginning, we set up this problem with a single field of two
     components: displacement in x direction and displacement in y direction */

  /* PetscFunctionBeginUser; */
  ierr = DMGetDS(dm, &prob);CHKERRQ(ierr);
  /* Setting up the linear elasticity problem */
  ierr = PetscDSSetResidual(prob,0,f0_u,f1_u_2d_testytesty);CHKERRQ(ierr);
  ierr = PetscDSSetJacobian(prob,0,0,NULL,NULL,NULL,g3_uu_2d);CHKERRQ(ierr);
  /* Setting the Neumann Boudnary Condition */
  if (user->neumann){
    ierr = PetscDSSetBdResidual(prob,0,f0_u_bd_2d,f1_u_bd_2d);CHKERRQ(ierr);
    /* ierr = PetscDSSetBdJacobian(prob,0,0,NULL,g1_uu_bd_2d,NULL,NULL);CHKERRQ(ierr); */

  }
  

      
  return(0);
}

PetscErrorCode SetupDiscretization(DM dm, AppCtx *user){

  DM cdm = dm;
  const PetscInt dim = user->dim; 	/* need to adapt this when changing
				   the dimensions of hte code */
  PetscFE fe, fe_bd;
  PetscDS prob;
  PetscErrorCode ierr;
  PetscBool simplex = user->simplex;
  PetscQuadrature q;
  PetscInt order;

  /* Creating the FE */
  ierr = PetscFECreateDefault(dm, dim, dim, simplex,"def_",-1,&fe);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fe, "deformation");CHKERRQ(ierr);
  if(user->neumann){
    ierr = PetscFEGetQuadrature(fe,&q);CHKERRQ(ierr);
    ierr = PetscQuadratureGetOrder(q,&order);CHKERRQ(ierr);
    /* Creating BD FE */
    ierr = PetscFECreateDefault(dm,dim-1, dim, simplex, "bd_def_",order,&fe_bd);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) fe_bd, "deformation");CHKERRQ(ierr);
  }
  /* Discretization and boundary conditons: */
  while (cdm)  {
    ierr = DMGetDS(cdm, &prob);CHKERRQ(ierr);
    ierr = PetscDSSetDiscretization(prob, 0, (PetscObject) fe); CHKERRQ(ierr);
    if (user->neumann){
      ierr = PetscDSSetBdDiscretization(prob, 0, (PetscObject) fe_bd); CHKERRQ(ierr);
    }
    ierr = SetupProblem(cdm, user); CHKERRQ(ierr);

    const PetscInt Ncomp = dim;
    const PetscInt components[] = {0,1};
    const PetscInt Nfid = 1;
    const PetscInt fid[] = {5}; /* {5}; */ 	/* The fixed faces.*/
    const PetscInt Npid  =1;
    const PetscInt pid[] = {3}; /* {3}; */ 	/* The pressure loaded faces */
  

    ierr =  DMAddBoundary(cdm, PETSC_TRUE, "fixed", "Face Sets",0, Ncomp, components, (void (*)()) zero_vector, Nfid, fid, user);CHKERRQ(ierr);
    if(user->neumann){
      ierr = DMAddBoundary(cdm, PETSC_FALSE, "load", "Face Sets",0, Ncomp, components, NULL, Npid, pid, user);CHKERRQ(ierr);
    }
    else{
      ierr = DMAddBoundary(cdm, PETSC_TRUE, "load", "Face Sets", 0, Ncomp, components, (void(*)()) pull, Npid, pid, user);CHKERRQ(ierr);
    }
        
      ierr = DMGetCoarseDM(cdm, &cdm); CHKERRQ(ierr);
  }


  ierr = PetscFEDestroy(&fe); CHKERRQ(ierr);
  if (user->neumann)  ierr = PetscFEDestroy(&fe_bd); CHKERRQ(ierr);
  return(0);
}


int main(int argc, char **argv){

  SNES snes;			/* nonlinear solver */
  DM dm, distributeddm;			/* problem definition */
  Vec u,r;			/* solution and residual vectors */
  Mat A,J;			/* Jacobian Matrix */
  AppCtx user;			/* user-defined work context */
  PetscInt its;			/* interation integer */
  PetscErrorCode ierr;
  PetscViewer viewer;
  
  /* Firing up Petsc */
  ierr= PetscInitialize(&argc, &argv,NULL,help);CHKERRQ(ierr);

  ierr = ProcessOptions(PETSC_COMM_WORLD,&user);CHKERRQ(ierr);
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
  /* ierr = CreateMesh(PETSC_COMM_WORLD,&user,&dm);CHKERRQ(ierr); */

  /* importing the gmsh file. Take note that only simplices give meaningful results in 2D at the moment (For which ever reasons) */
  ierr = DMPlexCreateFromFile(PETSC_COMM_WORLD,"testmesh_2D_box_quad.msh", PETSC_TRUE,&dm);CHKERRQ(ierr);
  ierr = DMPlexDistribute(dm,0,NULL,&distributeddm); CHKERRQ(ierr);
  if (distributeddm) {
    ierr=DMDestroy(&dm);CHKERRQ(ierr);
    dm = distributeddm;
  }
  
  ierr = DMSetFromOptions(dm);CHKERRQ(ierr);
  
  ierr = SNESSetDM(snes,dm);CHKERRQ(ierr);
  ierr = DMSetApplicationContext(dm, &user);CHKERRQ(ierr);
  ierr = SetupDiscretization(dm,&user);CHKERRQ(ierr);
  ierr = DMPlexCreateClosureIndex(dm,NULL);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(dm,&u);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&r);CHKERRQ(ierr);

  ierr = VecSet(u,(PetscReal) 1.0);CHKERRQ(ierr);


  ierr = DMSetMatType(dm, MATAIJ);CHKERRQ(ierr);
  ierr = DMCreateMatrix(dm, &J);CHKERRQ(ierr);
  A=J;
  
  ierr = DMPlexSetSNESLocalFEM(dm,&user,&user,&user);CHKERRQ(ierr);
  if (user.debug){  		/* Showing the Jacobi matrix for debugging purposes */
    ierr = PetscPrintf(PETSC_COMM_WORLD,"The Jacobian for the nonlinear solver ( and also the preconditioning matrix)\n");CHKERRQ(ierr);
    ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }
  ierr = SNESSetJacobian(snes, A, J, NULL, NULL);CHKERRQ(ierr);
  
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  if (user.showInitial){ ierr = DMVecViewLocal(dm, u, PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);}
  
  if (user.debug) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"initial guess \n");CHKERRQ(ierr);
    ierr = VecView(u, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }
  ierr = SNESSolve(snes,NULL,u);CHKERRQ(ierr);
  ierr = SNESGetIterationNumber(snes, &its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of snes iterations %i \n", its);CHKERRQ(ierr);
  if (user.showSolution){
    ierr = PetscPrintf(PETSC_COMM_WORLD,"solution: \n");CHKERRQ(ierr);
    ierr = VecChop(u, 3.0e-9); CHKERRQ(ierr); /* what does vecchop do? */
    ierr = VecView(u, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }

  if (user.visualization ==SOLUTION){
        
    if (user.verbose)PetscPrintf(PETSC_COMM_WORLD,"Creating the vtk output file... \n");
    ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD,"solution.vtk",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) u,"deformation");CHKERRQ(ierr);
    ierr = VecView(u,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    if (user.verbose)PetscPrintf(PETSC_COMM_WORLD,"Done creating the VTK file. \n");
  }

  ierr = VecViewFromOptions(u,NULL,"-sol_vec_view");CHKERRQ(ierr);

  if ( A != J){MatDestroy(&A);}
  MatDestroy(&J);
  VecDestroy(&u);
  VecDestroy(&r);
  SNESDestroy(&snes);
  DMDestroy(&dm);
  PetscFinalize();

  return 0;
  
}
