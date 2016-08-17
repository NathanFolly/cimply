
static char help[] = "simply FEM program. To be extended...";


/* The PETSc packages we need: */
#include <petscdmplex.h>
#include <petscds.h>
#include <petscsnes.h>
#include <petscts.h>

/* user defined Application Context that helps to manage the options (or the contetext) of the program  */

typedef enum {RUN_FULL, RUN_TEST} RunType;
typedef enum {NONE, LABELS, BOUNDARIES, SOLUTION} Visualization;

typedef struct  /* The Leapfrog timestepping context */
{
  PetscReal total;  /* total time in seconds */
  PetscReal dt;     /* time setp size */
  PetscBool halfstep;  /* Is the timestep a halfstep? */
} TSctx;



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
  TSctx time;
}AppCtx;



/* The static inline statement will tell the compiler to inline the function in the code directly instead of making multiple function calls. This improves performance especially for small functions that are called repeatedly throughout the program. */

/* The Kronecker Delta */
static const PetscInt delta2D[2*2] = {1,0,0,1};
static const PetscInt delta3D[3*3] = {1,0,0,0,1,0,0,0,1};

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
  u[0] = 0.0;
  u[1] = 0.1;
  return 0;
}


void f0_u(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscScalar f0[])
{
  const PetscInt Ncomp = dim;
  PetscInt comp;
  for(comp=0;comp<Ncomp;comp++) f0[comp]=0.0;
    
}
void f0_u_transient(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscScalar f0[])
{
  const PetscInt Ncomp = dim;
  const PetscReal rho = 7850;  /* Density of steel needed for the inertia term */
  PetscInt comp;
  for(comp=0;comp<Ncomp;comp++) f0[comp]= 0.0-rho*u_t[Ncomp+comp];
}


void f1_u_2d(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscScalar f1[])
{
  const PetscInt Ncomp = dim;
  const PetscReal mu = 86, lbda=115.4;
  
  /* f1 is the cauchy stress tensor*/
  
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

void f1_u_ldis(PetscInt dim, PetscInt Nf, PetscInt NfAux,
               const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
               const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
               const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
               const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscScalar f1[])
{
  const PetscInt Ncomp = dim;
  const PetscReal mu =86, lbda=115.4;
  PetscScalar F[dim*Ncomp], E[dim*Ncomp], C[dim*Ncomp];
  PetscInt d, comp,i;
  
  /* In case of large deformations, f1 should be the second piola kirchhoff stress tensor */
  /* TODO: If we have rigid body motion, we should perform a polar decomposition of the
     deformation gradient tensor. This is not yet done in this version.*/

  /* Step 1: Create the deformation gradient tensor F=I+u_x */
  for (d=0;d<dim;d++){
    for (comp=0; comp<Ncomp; comp++){
      F[d*dim+comp] = u_x[d*dim+comp];
    }
    F[d*dim+d] += 1;
  }
  
  /* Step 2: Create the Cauchy-Green Deformation Tensor C = F^T F*/
  for (d=0;d<dim;d++){
    for (comp=0; comp<Ncomp; comp++){
      for (i=0;i<dim;i++){
        C[d*dim+comp] += F[i*dim+d]*F[i*dim+comp];
      }
    }
  }
  /* Step 3: Create the Lagrangian strain tensor E=0.5*(C-1) */
  for (d=0;d<dim;d++){
    for (comp=0; comp<Ncomp; comp++){
      E[d*dim+comp] = 0.5*C[d*dim+comp];
    }
    E[d*dim+d] -= 0.5;
  }
  /* Step 4: Create the second piola - Kirchoff stress tensor f1 = lbda*tr(E)*1 +2*mu*E */
  for (d=0;d<dim;d++){
    for (comp=0; comp<Ncomp; comp++){
      f1[d*dim+comp] = 2*mu*E[d*dim+comp];
      /* PetscPrintf(PETSC_COMM_WORLD,"f1[%i,%i] = %f\n",d,comp,f1[d*dim+comp]); */
    }
    for (comp=0; comp<Ncomp; comp++){
      f1[d*dim+d] += lbda*E[comp*dim+comp];
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

void g3_uu_ldis(PetscInt dim, PetscInt Nf, PetscInt NfAux,
	      const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], 
	      const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
	      const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
	      const PetscScalar a_x[], PetscReal t, PetscReal u_tShift, const PetscReal x[],
	      PetscScalar g3[]){

  const PetscReal mu = 86, lbda = 115.4;
  const PetscInt Ncomp = dim;
  PetscInt i, j, k, l, m, n;
  PetscScalar G[dim*dim*Ncomp*Ncomp];  /* Cauchy-Green  strain tensor */
  PetscScalar D[dim*dim*Ncomp*Ncomp];  /* constitutive tensor for isotropic material */
  /* constructing the cauchy-green strain tensor */
  for (m=0;m<dim;m++){  /* primary index of the cauchy green strain tensor */
    for (j=0;j<dim;j++){  /* component of trial function */
      for (n=0;n<Ncomp;n++){  /* secondary index of the cauchy green strain tensor */
        for (l=0;l<Ncomp;l++){  /* derivative index for trial function */
          G[((m*dim+j)*dim+n)*dim+l]=delta2D[n*dim+l]*(u_x[j*dim+m]+delta2D[j*dim+m])+delta2D[m*dim+l]*delta2D[n*dim+j];
        }
      }
    }
  }

  /* constructing the constitutive tensor for isotropic material */
  for (i=0;i<dim;i++){  /* component of test function */
    for (m=0;m<dim;m++){  /* primary index of the lagrangian strain tensor */
      for (k=0;k<dim;k++){  /* derivative index of the test function */
        for (n=0;n<dim;n++){  /* secondary index of the lagrangian strain
                               * tensor*/
          D[((i*dim+m)*dim+k)*dim+n] = lbda*delta2D[i*dim+k]*delta2D[m*dim+n]+mu*(delta2D[i*dim+m]*delta2D[k*dim+n]+delta2D[i*dim+n]*delta2D[k*dim+m]);
        }
      }
    }
  }
  /* constructing the integrand for gradient of testfunction and gradient of trialfunction */
  for (i=0;i<dim;i++){  /* component of test function */
    for (j=0;j<dim;j++){  /* component of the trial function */
      for (k=0;k<dim;k++){  /* derivative index of the test funciton */
        for (l=0;l<dim;l++){  /* derivative index of the trial function */
          for (m=0;m<dim;m++){  /* primary index of the cauchy-green strain
                                 * tensor */
            for (n=0;n<dim;n++){  /* secondary index of the cauchy-green
                                   * strain tensor */
              /* Maybe it helps if I assemble the tensor directly? -- nope */
              g3[((i*dim+j)*dim+k)*dim+l] += /* 0.5*(lbda*delta2D[i*dim+k]*delta2D[m*dim+n]+mu*(delta2D[i*dim+m]*delta2D[k*dim+n]+delta2D[i*dim+n]*delta2D[k*dim+m]))*(delta2D[n*dim+l]*(u_x[j*dim+m]+delta2D[j*dim+m])+delta2D[m*dim+l]*delta2D[n*dim+j]); */

              D[((i*dim+m)*dim+k)*dim+n]*0.5*(G[((m*dim+j)*dim+n)*dim+l]);
            }
          }
          /* PetscPrintf(PETSC_COMM_WORLD,"g3 [%i %i %i %i] = %f   u0_0 = %f  u0_1 = %f   u1_0 = %f   u1_1 = %f \n", i,j,k,l, g3[((i*dim+j)*dim+k)*dim+l], u_x[0*dim+0], u_x[0*dim+1], u_x[1*dim + 0], u_x[1*dim+1]); */
          

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
  /* Setting  the surface traction tensor eqal to the external load acting on the boundary in question  */
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
  
    
  for (i=0; i<Ncomp; i++){
    for (j=0; j<Ncomp; j++){
      for (k=0; k<Ncomp; k++){
        g1[(i*Ncomp+j)*Ncomp+k]=n[i]*delta2D[j*dim+k]*lbda+(1-delta2D[j*dim+k])*mu*(n[k]*delta2D[i*dim+j]+n[j]*delta2D[i*dim+k])+n[i]*delta2D[j*dim+k]*delta2D[i*dim+j]*2*mu;
      }
    }
  }   
}


void f0_vel_2d(PetscInt dim, PetscInt Nf, PetscInt NfAux,
               const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
               const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
               const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
               const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscScalar f0[])
{
  const PetscInt Ncomp = dim;
  PetscInt comp;
  /* We're solving the equation vel - du/dt = 0 so: */

  for (comp=0;comp<Ncomp;comp++){
    f0[comp] = u[Ncomp+comp]-u_t[comp];
  }
}


void f1_vel_2d(PetscInt dim, PetscInt Nf, PetscInt NfAux,
               const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
               const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
               const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
               const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscScalar f1[])
{
  const PetscInt Ncomp = dim;
  PetscInt comp,d;

  for (comp=0;comp<Ncomp;comp++){
    for (d=0;d<dim;d++){
      f1[d*comp+d] = 0.0;
    }
  }
}

void g0_velvel_2d(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                  const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], 
                  const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
                  const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
                  const PetscScalar a_x[], PetscReal t, PetscReal u_tShift, const PetscReal x[],
                  PetscScalar g0[]){
  PetscInt Ncomp = dim;
  PetscInt i;

  for (i=0;i<Ncomp;i++){
    g0[i]= 1.0;
  }
  
}

void g0_dyn_velu_2d(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
                    const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
                    const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
                    const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscScalar g0[]){

  PetscInt Ncomp = dim;
  PetscInt i;

  for (i=0;i<Ncomp;i++){
    g0[i]=1.0;
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
  options->transient = PETSC_FALSE;
  options->ldis = PETSC_FALSE;
  options->time.total = 3.0;
  options->time.dt = 0.1;
  

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
  PetscInt NumFields;


  ierr = DMGetDS(dm, &prob);CHKERRQ(ierr);

  /* Setting up the linear elasticity problem */
  if (user->transient){
    ierr = PetscDSSetResidual(prob,0,f0_u_transient,f1_u_2d);CHKERRQ(ierr);
    ierr = PetscDSSetJacobian(prob,0,0,NULL,NULL,NULL,g3_uu_2d);CHKERRQ(ierr);
  }
  else if(user->ldis){
    if(user->verbose) ierr = PetscPrintf(PETSC_COMM_WORLD,"Setting up large displacement problem\n");
    ierr = PetscDSSetResidual(prob,0,f0_u,f1_u_ldis);CHKERRQ(ierr);
    ierr = PetscDSSetJacobian(prob,0,0,NULL,NULL,NULL,g3_uu_ldis);CHKERRQ(ierr);
  }
  else{
    ierr = PetscDSSetResidual(prob,0,f0_u,f1_u_2d);CHKERRQ(ierr);
    ierr = PetscDSSetJacobian(prob,0,0,NULL,NULL,NULL,g3_uu_2d);CHKERRQ(ierr);
  }

  /* setting up the velocity field for the transient analysis */
  if(user->transient){
    ierr = PetscDSSetResidual(prob,1,f0_vel_2d,f1_vel_2d);CHKERRQ(ierr);
    ierr = PetscDSSetJacobian(prob,1,1,g0_velvel_2d,NULL,NULL,NULL);CHKERRQ(ierr);
    /* Maybe we need a dynamic jacobian? -- apparently not. */
    ierr = PetscDSSetDynamicJacobian(prob,1,0,g0_dyn_velu_2d,NULL,NULL,NULL);CHKERRQ(ierr);
  }
  /* Setting the Neumann Boudnary Condition */
  if (user->neumann){
    ierr = PetscDSSetBdResidual(prob,0,f0_u_bd_2d,f1_u_bd_2d);CHKERRQ(ierr);
    /* ierr = PetscDSSetBdJacobian(prob,0,0,NULL,g1_uu_bd_2d,NULL,NULL);CHKERRQ(ierr); */

  }
  
  ierr = PetscDSGetNumFields(prob,&NumFields);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Problem set up with %i fields\n",NumFields);CHKERRQ(ierr);
      
  return(0);
}

PetscErrorCode SetupDiscretization(DM dm, AppCtx *user){

  DM cdm = dm;
  const PetscInt dim = user->dim; 	/* need to adapt this when changing
				   the dimensions of hte code */
  PetscFE fe, fe_bd, fe_vel;
  PetscDS prob;
  PetscErrorCode ierr;
  PetscBool simplex = user->simplex;
  PetscQuadrature q;
  PetscInt order;

  /* Creating the FE */
  ierr = PetscFECreateDefault(dm, dim, dim, simplex,"def_",-1,&fe);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fe, "deformation");CHKERRQ(ierr);
  if(user->neumann){
    if (user->ldis) PetscPrintf(PETSC_COMM_WORLD,"Large displacement formulation has no option for neumann conditions as of yet");
    ierr = PetscFEGetQuadrature(fe,&q);CHKERRQ(ierr);
    ierr = PetscQuadratureGetOrder(q,&order);CHKERRQ(ierr);
    /* Creating BD FE */
    ierr = PetscFECreateDefault(dm,dim-1, dim, simplex, "bd_def_",order,&fe_bd);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) fe_bd, "deformation");CHKERRQ(ierr);
  }
  if(user->transient){
    if (user->ldis) PetscPrintf(PETSC_COMM_WORLD,"Large displacement formulation has no option for transient as of yet");
    ierr = PetscFEGetQuadrature(fe,&q);CHKERRQ(ierr);
    ierr = PetscQuadratureGetOrder(q,&order);CHKERRQ(ierr);
    /* Creating the FE field for the velocity */
    ierr = PetscFECreateDefault(dm,dim,dim,simplex,"vel_",order,&fe_vel);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) fe_vel, "velocity");CHKERRQ(ierr);
  }
  /* Discretization and boundary conditons: */
  while (cdm)  {
    ierr = DMGetDS(cdm, &prob);CHKERRQ(ierr);
    ierr = PetscDSSetDiscretization(prob, 0, (PetscObject) fe); CHKERRQ(ierr);
    if (user->transient){
      ierr = PetscDSSetDiscretization(prob,1, (PetscObject) fe_vel);CHKERRQ(ierr);
    }
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
  if (user->transient) ierr = PetscFEDestroy(&fe_vel); CHKERRQ(ierr);
  return(0);
}



int main(int argc, char **argv){

  SNES snes;			/* nonlinear solver */
  DM dm, distributeddm;			/* problem definition */
  Vec u,r;			/* solution and residual vectors */
  Mat A,J;			/* Jacobian Matrix */
  AppCtx user;			/* user-defined work context */
  PetscErrorCode ierr;
  PetscViewer viewer;
  TS ts;




  
  /* Firing up Petsc */
  ierr= PetscInitialize(&argc, &argv,NULL,help);CHKERRQ(ierr);

  ierr = ProcessOptions(PETSC_COMM_WORLD,&user);CHKERRQ(ierr);


  /* ierr = CreateMesh(PETSC_COMM_WORLD,&user,&dm);CHKERRQ(ierr); */

  /* importing the gmsh file. Take note that only simplices give meaningful results in 2D at the moment (For which ever reasons) */
  ierr = DMPlexCreateFromFile(PETSC_COMM_WORLD,"Beam_coarse.msh", PETSC_TRUE,&dm);CHKERRQ(ierr);
  ierr = DMPlexDistribute(dm,0,NULL,&distributeddm); CHKERRQ(ierr);
  if (distributeddm) {
    ierr=DMDestroy(&dm);CHKERRQ(ierr);
    dm = distributeddm;
  }
  
  ierr = DMSetFromOptions(dm);CHKERRQ(ierr);
  
 


  if (user.transient){
    Vec u_start;
    PetscReal ftime;
    PetscInt nsteps, its;
    TSConvergedReason ConvergedReason;
    PetscBool loadedspring = PETSC_TRUE;
    
    if (loadedspring){  /* preloading the cantilever so that I can run the TS without loading to see if the problem comes from the neumann conditions in the displacement field */
      user.transient=PETSC_FALSE;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Creating initial solution for loaded cantilever\n");CHKERRQ(ierr);
      ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
      ierr = SNESSetDM(snes,dm);CHKERRQ(ierr);
      ierr = DMSetApplicationContext(dm, &user);CHKERRQ(ierr);
      ierr = SetupDiscretization(dm,&user);CHKERRQ(ierr);
      ierr = DMPlexCreateClosureIndex(dm,NULL);CHKERRQ(ierr);
      
      ierr = DMCreateGlobalVector(dm,&u_start);CHKERRQ(ierr);
      ierr = VecDuplicate(u_start,&r);CHKERRQ(ierr);
      
      
      ierr = DMSetMatType(dm, MATAIJ);CHKERRQ(ierr);
      ierr = DMCreateMatrix(dm, &J);CHKERRQ(ierr);
      A=J;

      
      ierr = DMPlexSetSNESLocalFEM(dm,&user,&user,&user);CHKERRQ(ierr);

      ierr = SNESSetJacobian(snes, A, J, NULL, NULL);CHKERRQ(ierr);

      ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
      ierr = SNESSolve(snes,NULL,u_start);CHKERRQ(ierr);

      ierr = SNESGetIterationNumber(snes, &its);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Static loaded state determined. Number of snes iterations %i \n", its);CHKERRQ(ierr);

      ierr = MatDestroy(&J); CHKERRQ(ierr);
      ierr = MatDestroy(&A); CHKERRQ(ierr);
      ierr = VecDestroy(&r); CHKERRQ(ierr);

      user.transient=PETSC_TRUE;
      user.neumann = PETSC_FALSE;

    }
    
    
    PetscPrintf(PETSC_COMM_WORLD,"starting transient analysis. Total time: %g s \n",user.time.total);
    ierr =  TSCreate(PETSC_COMM_WORLD, &ts); CHKERRQ(ierr);
    ierr = TSSetType(ts, TSCN); CHKERRQ(ierr);
    ierr = TSSetDM(ts,dm); CHKERRQ(ierr);  
    ierr = DMSetApplicationContext(dm, &user);CHKERRQ(ierr);
    ierr = SetupDiscretization(dm,&user);CHKERRQ(ierr);
    ierr = DMPlexCreateClosureIndex(dm,NULL);CHKERRQ(ierr);

    ierr = DMTSSetBoundaryLocal(dm, DMPlexTSComputeBoundary, &user); CHKERRQ(ierr);
    ierr = DMTSSetIFunctionLocal(dm, DMPlexTSComputeIFunctionFEM, &user); CHKERRQ(ierr);
    ierr = DMTSSetIJacobianLocal(dm, DMPlexTSComputeIJacobianFEM, &user); CHKERRQ(ierr);


    
    ierr = DMCreateGlobalVector(dm, &u); CHKERRQ(ierr);

    if (loadedspring){
      /* setting the loaded cantilever static solution as initial condition for the transient analysis */
      u=u_start;
    }

    ierr = TSSetSolution(ts,u); CHKERRQ(ierr);
    
    ierr = TSSetDuration(ts,1000,user.time.total); CHKERRQ(ierr);
    ierr = TSSetExactFinalTime(ts, user.time.total); CHKERRQ(ierr);
    ierr = TSSetInitialTimeStep(ts,0.0, user.time.dt); CHKERRQ(ierr);
    ierr = TSSetFromOptions(ts); CHKERRQ(ierr);

    ierr = TSSolve(ts,NULL); CHKERRQ(ierr);

    /* ierr = TSGetSolveTime(ts,&ftime); */
    /* ierr = TSGetTimeStepNumber(ts,&nsteps); CHKERRQ(ierr); */
    /* ierr = TSGetConvergedReason(ts,&ConvergedReason); CHKERRQ(ierr); */

    ierr = TSView(ts,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"Transient analysis finished. \n");
    /* ierr = PetscPrintf(PETSC_COMM_WORLD, "%s at time %g after %D steps \n",ConvergedReason,ftime,nsteps); CHKERRQ(ierr); */

    /* Write solution to a VTK file */
    ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD,"solution.vtk",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    ierr = DMPlexVTKWriteAll((PetscObject) dm, viewer);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

    ierr = TSDestroy(&ts); CHKERRQ(ierr);
    ierr = VecDestroy(&u); CHKERRQ(ierr);
    
  }  
  else{
    PetscInt its;  /* number of iterations needed by the SNES solver */
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
  }
  DMDestroy(&dm);
  PetscFinalize();
  return 0;
  
}
