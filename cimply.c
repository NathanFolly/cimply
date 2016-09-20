
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

PetscErrorCode constrict_y(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx)
{
  u[1] = 0.0;
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


void f1_u(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscScalar f1[])
{
  /* 3-D or Two-dimensional plain strain formulation */
  const PetscInt Ncomp = dim;
  const PetscReal mu =76.923076923, lbda=115.384615385;
  
  /* f1 is the cauchy stress tensor*/
  
  /*u_x = deflection gradient*/
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

void f1_u_2d_pstress(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscScalar f1[])
{
  /* Teo-dimensional plain strain formulation */
  const PetscInt Ncomp = dim;
  const PetscReal mu =76.923076923, lbda=115.384615385;
  const PetscReal lbdaBar = 2*lbda*mu/(lbda+2*mu);
  
  /* f1 is the cauchy stress tensor*/
  
  /*u_x = deformation gradient*/
  /* Hence is the strain epsilon_ij=0.5(u_i,j+u_j,i) => epsilon[comp*dim+d]=0.5(u_x[comp*dim+d]+u_x[d*dim+comp]) */
  PetscInt comp, d;
  for(comp=0;comp<Ncomp;comp++){
    for(d=0;d<dim;d++){
      f1[comp*dim+d]=mu*(u_x[comp*dim+d]+u_x[d*dim+comp]);
    }
    for(d=0;d<dim;d++){
      f1[comp*dim+comp]+=lbdaBar*u_x[d*dim+d];
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
  const PetscReal mu =76.9, lbda=115.4;
  PetscScalar F[dim*Ncomp], E[dim*Ncomp], C[dim*Ncomp];
  PetscReal detF;
  PetscInt d, comp,i;
  
  /* In case of large deformations, f1 should be the second piola kirchhoff stress tensor */
  /* TODO: If we have rigid body motion, we should perform a polar decomposition of the
     deformation gradient tensor. This is not yet done in this version.*/
  /* This is the two-dimensional plane strain formulation */

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
  /* maybe I need to multiply with J (J=detF) ?*/
  /* detF = F[0]*F[3]-F[1]*F[2]; */
  /* for (comp=0;comp<Ncomp;comp++){ */
  /*   for (d=0;d<dim;d++){ */
  /*     f1[comp*dim+d] *= (detF); */
  /*   } */
  /* } */
 }


void g3_uu_2d(PetscInt dim, PetscInt Nf, PetscInt NfAux,
	      const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], 
	      const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
	      const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
	      const PetscScalar a_x[], PetscReal t, PetscReal u_tShift, const PetscReal x[],
	      PetscScalar g3[]){
  
  /* const PetscInt Ncomp = dim; */
  const PetscReal mu =76.923076923, lbda=115.384615385;
  PetscInt i,j,k,l;
  
  /* g3 is the elasticity tensor */
  /* This is the 2-D plane strain formulation for a saint-venant kirchhoff material */
  
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

void g3_uu_2d_pstress(PetscInt dim, PetscInt Nf, PetscInt NfAux,
	      const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], 
	      const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
	      const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
	      const PetscScalar a_x[], PetscReal t, PetscReal u_tShift, const PetscReal x[],
	      PetscScalar g3[]){
  
  /* const PetscInt Ncomp = dim; */
  const PetscReal mu =76.923076923, lbda=115.384615385;
  const PetscReal lbdaBar = 2*lbda*mu/(lbda+2*mu);
  PetscInt i,j,k,l;
  
  /* g3 is the elasticity tensor */
  /* This is the 2-D plane strain formulation for a saint-venant kirchhoff material */
  
  for (i=0;i<dim;i++){
    for (j=0;j<dim;j++){
      for (k=0; k < dim; k++){
        for (l=0;l<dim; l++){
          g3[((i*dim+j)*dim+k)*dim+l]=lbdaBar*delta2D[i*dim+k]*delta2D[j*dim+l]+2*mu*delta2D[i*dim+k]*delta2D[j*dim+l]*delta2D[i*dim+j]+mu*(1-delta2D[i*dim+k])*(1-delta2D[j*dim+l]);
        }
      }
    }
  }
}

void g1_uu_ldis(PetscInt dim, PetscInt Nf, PetscInt NfAux,
	      const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], 
	      const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
	      const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
	      const PetscScalar a_x[], PetscReal t, PetscReal u_tShift, const PetscReal x[],
	      PetscScalar g1[]){

  /* We need to linearize our energy form since it is dependent on the deformation
   We have L[a(u,û)] = int_omega0 DeltaS:Ê+S:DeltaÊ dOmega0
  Now DeltaS is (partial S)/(partial E) : DeltaE which is D:DeltaE where D is the constitutive tensor for a St.Venant Kirchhoff material.
  This function here hence presents the integrand for the test function (Ê) and the trial function gradient term. It should be S, the second Piola Kirchhoff stress*/
  const PetscInt Ncomp = dim;
  const PetscReal mu =76.9, lbda=115.4;
  PetscScalar F[dim*Ncomp], E[dim*Ncomp], C[dim*Ncomp], S[dim*Ncomp];
  PetscScalar G[dim*dim*Ncomp*Ncomp];  /* Cauchy-Green  strain tensor */
  PetscReal detF;
  PetscInt d, comp,i,j,l,m,n;
 
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
      S[d*dim+comp] = 2*mu*E[d*dim+comp];
      /* PetscPrintf(PETSC_COMM_WORLD,"f1[%i,%i] = %f\n",d,comp,f1[d*dim+comp]); */
    }
    for (comp=0; comp<Ncomp; comp++){
      S[d*dim+d] += lbda*E[comp*dim+comp];
    }
  }

    /* constructing the cauchy-green strain tensor */
  for (m=0;m<dim;m++){  /* primary index of the cauchy green strain tensor */
    for (j=0;j<dim;j++){  /* component of trial function */
      for (n=0;n<Ncomp;n++){  /* secondary index of the cauchy green strain tensor */
        for (l=0;l<Ncomp;l++){  /* derivative index for trial function */
          g1[j*dim+l]+=S[n*dim+m]*0.5*(delta2D[n*dim+l]*(u_x[j*dim+m]+delta2D[j*dim+m]+delta2D[m*dim+n]*u_x[j*dim+m])+delta2D[m*dim+l]*delta2D[n*dim+j]);
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

  const PetscReal mu = 76.9, lbda = 115.4;
  const PetscInt Ncomp = dim;
  PetscInt i, j, k, l, m, n;
  PetscInt d, comp;
  PetscScalar G[dim*dim*Ncomp*Ncomp];  /* Cauchy-Green  strain tensor */
  PetscScalar D[dim*dim*Ncomp*Ncomp];  /* constitutive tensor for isotropic material */
  PetscScalar F[Ncomp*dim];            /* deformation gradient */
  PetscReal detF;                      /* determinant of F */

  PetscInt delta[Ncomp*dim];
  for (i=0;i<dim;i++){
    for (j=0;j<dim;j++){
      if(dim==3){
        delta[i*dim+j] = delta3D[i*dim+j];
      }
      else{
        delta[i*dim+j]= delta2D[i*dim+j];
      }
    }
  }
  

  /* partial derivative of the lagrange strain ensor with respect to the displacement gradient G_mjnl = (partial E_mn)/(partial u^j_l) */
  for (m=0;m<dim;m++){  /* primary index of the cauchy green strain tensor */
    for (j=0;j<dim;j++){  /* component of trial function */
      for (n=0;n<Ncomp;n++){  /* secondary index of the cauchy green strain tensor */
        for (l=0;l<Ncomp;l++){  /* derivative index for trial function */
          G[((m*dim+j)*dim+n)*dim+l]=delta[n*dim+l]*(u_x[j*dim+m]+delta[j*dim+m])+delta[m*dim+l]*(u_x[j*dim+n]+delta[n*dim+j]);
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
          D[((i*dim+m)*dim+k)*dim+n] = lbda*delta[i*dim+k]*delta[m*dim+n]+mu*(delta[i*dim+m]*delta[k*dim+n]+delta[i*dim+n]*delta[k*dim+m]);
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
              g3[((i*dim+j)*dim+k)*dim+l] += /* 0.5*(lbda*delta[i*dim+k]*delta[m*dim+n]+mu*(delta[i*dim+m]*delta[k*dim+n]+delta[i*dim+n]*delta[k*dim+m]))*(delta[n*dim+l]*(u_x[j*dim+m]+delta[j*dim+m])+delta[m*dim+l]*delta[n*dim+j]); */

                  D[((i*dim+m)*dim+k)*dim+n]*0.5*(G[((m*dim+j)*dim+n)*dim+l]);
            }
          }
          /* PetscPrintf(PETSC_COMM_WORLD,"g3 [%i %i %i %i] = %f   u0_0 = %f  u0_1 = %f   u1_0 = %f   u1_1 = %f \n", i,j,k,l, g3[((i*dim+j)*dim+k)*dim+l], u_x[0*dim+0], u_x[0*dim+1], u_x[1*dim + 0], u_x[1*dim+1]); */
        }
      }
    }
  }
  /* PetscPrintf(PETSC_COMM_WORLD,"\n \n " ); */
}

void g3_uu_3d(PetscInt dim, PetscInt Nf, PetscInt NfAux,
	      const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], 
	      const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
	      const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
	      const PetscScalar a_x[], PetscReal t, PetscReal u_tShift, const PetscReal x[],
	      PetscScalar g3[]){
  
  const PetscInt Ncomp = dim;
  const PetscReal mu =76.923076923, lbda=115.384615385;
  PetscInt i,j,k,l,m,n;
    PetscScalar G[dim*dim*Ncomp*Ncomp];  /* Cauchy-Green  strain tensor */
  PetscScalar D[dim*dim*Ncomp*Ncomp];  /* constitutive tensor for isotropic material */
    

  /* partial derivative of the cauchy strain ensor with respect to the displacement gradient G_mjnl = (partial epsilon_mn)/(partial u^j_l) */
  for (m=0;m<dim;m++){  /* primary index of the cauchy green strain tensor */
    for (j=0;j<dim;j++){  /* component of trial function */
      for (n=0;n<Ncomp;n++){  /* secondary index of the cauchy green strain tensor */
        for (l=0;l<Ncomp;l++){  /* derivative index for trial function */
          G[((m*dim+j)*dim+n)*dim+l]=0.5*(delta3D[n*dim+l]*delta3D[j*dim+m]+delta3D[m*dim+l]*delta3D[n*dim+j]);
        }
      }
    }
  }

  /* constructing the constitutive tensor for isotropic material */
  for (i=0;i<dim;i++){  /* component of test function */
    for (m=0;m<dim;m++){  /* primary index of the cauchy strain tensor */
      for (k=0;k<dim;k++){  /* derivative index of the test function */
        for (n=0;n<dim;n++){  /* secondary index of the cauchy strain
                               * tensor*/
          D[((i*dim+m)*dim+k)*dim+n] = lbda*delta3D[i*dim+k]*delta3D[m*dim+n]+mu*(delta3D[i*dim+m]*delta3D[k*dim+n]+delta3D[i*dim+n]*delta3D[k*dim+m]);
        }
      }
    }
  }
  /* constructing the integrand for gradient of testfunction and gradient of trialfunction */
  for (i=0;i<dim;i++){  /* component of test function */
    for (j=0;j<dim;j++){  /* component of the trial function */
      for (k=0;k<dim;k++){  /* derivative index of the test funciton */
        for (l=0;l<dim;l++){  /* derivative index of the trial function */
          for (m=0;m<dim;m++){  /* primary index of the cauchy strain
                                 * tensor */
            for (n=0;n<dim;n++){  /* secondary index of the cauchy
                                   * strain tensor */
              g3[((i*dim+j)*dim+k)*dim+l] +=  D[((i*dim+m)*dim+k)*dim+n]*(G[((m*dim+j)*dim+n)*dim+l]);
            }
          }
        }
      }
    }
  }
}


void f0_u_bd(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
                const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
                const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
                const PetscScalar a_x[], PetscReal t, const PetscReal x[], const PetscReal n[],
                PetscScalar f0[])
{
  const PetscInt Ncomp=dim;
  /* Setting  the surface traction tensor eqal to the external load acting on the boundary in question  */
  const PetscReal mu =76.923076923, lbda=115.384615385;
  const PetscScalar traction[] = {0.0, 0.01, 0.0};
  PetscInt comp;
  const PetscReal pressure = 0.001;
  PetscReal nonsense;

  if(x[0]<0.8)nonsense=0.0;
  else nonsense = 1.0;
  
  for (comp=0; comp<Ncomp; ++comp){
    f0[comp] = pressure*n[comp];
  }
  
}


void f0_u_bd_ldis(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
                const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
                const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
                const PetscScalar a_x[], PetscReal t, const PetscReal x[], const PetscReal n[],
                PetscScalar f0[])
{
  const PetscInt Ncomp=dim;
  /* TODO: define the boundary residual for the large displacement formulation (should be the first piola Kirchhoff stress */
  const PetscReal mu =76.9, lbda=115.4;
  const PetscScalar traction[] = {0.0,0.0,0.0,0.01,0.0,0.0,0.0,0.0,0.0};
  const PetscScalar trac[] = {0.0,0.001,0.0};
  PetscInt comp, d, i;
  PetscScalar F[Ncomp*dim], FInv[Ncomp*dim], S[Ncomp*dim];
  PetscReal J;
  PetscReal nonsense;
  
  /* step 1 construct the transpose of the deformation gradient tensor F^T */
  for(comp=0; comp<Ncomp; comp++){
    for(d=0; d<dim; d++){
      F[comp*dim+d]= u_x[comp*dim+d];
    }
    F[comp*dim+comp]+=1;
  }

  J = (F[0*3+0]*(F[1*3+1]*F[2*3+2] - F[1*3+2]*F[2*3+1])+F[0*3+1]*(F[1*3+2]*F[2*3+0] - F[1*3+0]*F[2*3+2]) + F[0*3+2]*(F[1*3+0]*F[2*3+1] - F[1*3+1]*F[2*3+0]));

  /* Inverse of F */

  FInv[0*3+0] = F[1*3+1]*F[2*3+2] - F[1*3+2]*F[2*3+1];
  FInv[0*3+1] = F[1*3+2]*F[2*3+0] - F[1*3+0]*F[2*3+2];
  FInv[0*3+2] = F[1*3+0]*F[2*3+1] - F[1*3+1]*F[2*3+0];
  FInv[1*3+0] = F[0*3+2]*F[2*3+1] - F[0*3+1]*F[2*3+2];
  FInv[1*3+1] = F[0*3+0]*F[2*3+2] - F[0*3+2]*F[2*3+0];
  FInv[1*3+2] = F[0*3+1]*F[2*3+0] - F[0*3+0]*F[2*3+1];
  FInv[2*3+0] = F[0*3+1]*F[1*3+2] - F[0*3+2]*F[1*3+1];
  FInv[2*3+1] = F[0*3+2]*F[1*3+0] - F[0*3+0]*F[1*3+2];
  FInv[2*3+2] = F[0*3+0]*F[1*3+1] - F[0*3+1]*F[1*3+0];
  
  /* Step 2 construct the inverse of F^T: */
  /* FInv[0]=F[3]*1/(F[0]*F[3]-F[1]*F[2]); */
  /* FInv[1]=-F[1]*1/(F[0]*F[3]-F[1]*F[2]); */
  /* FInv[2]=-F[2]*1/(F[0]*F[3]-F[1]*F[2]); */
  /* FInv[3]=F[0]*1/(F[0]*F[3]-F[1]*F[2]); */

  if(x[0]<0.8){nonsense = 0.0;}
  else{nonsense=1.0;}

  /* Step 3 Create the second piola Kirchhoff stress tensor */
  for (comp=0; comp<Ncomp; ++comp){
    for (d=0;d<dim;d++){
      for (i=0;i<dim;i++){
        S[comp*dim+d] += traction[comp*dim+i]*FInv[i*dim+d]/J;
      }
    }
  }
  /* Step 4 S*n */

  for (comp=0; comp<Ncomp; ++comp){
    for (d=0;d<dim;d++){
      f0[comp]+=FInv[comp*dim+d]/J*n[comp];
      /* PetscPrintf(PETSC_COMM_WORLD,"S[%i] = %f      n[%i]  = %f \n", comp*dim+d, S[comp*dim+d],comp,n[comp]); */
    }
    /* PetscPrintf(PETSC_COMM_WORLD,"f0[%i] = %f \n", comp, f0[comp]);  */
    f0[comp]*=trac[comp];
  }
  /* f0[1] =0.06; */
}
    

    
void f1_u_bd(PetscInt dim, PetscInt Nf, PetscInt NfAux,
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
    
  const PetscReal mu=76.9, lbda=115.4;
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
  options->dim = 3;
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


PetscErrorCode SetupProblem(DM dm, AppCtx *user)
{
  PetscDS prob; 		   /* a DS is a "discrete system" managing discretisations */
  PetscErrorCode ierr;
  PetscInt NumFields;


  ierr = DMGetDS(dm, &prob);CHKERRQ(ierr);


  
     /* Setting up the linear elasticity problem */
    if (user->transient){
      ierr = PetscDSSetResidual(prob,0,f0_u_transient,f1_u);CHKERRQ(ierr);
      ierr = PetscDSSetJacobian(prob,0,0,NULL,NULL,NULL,g3_uu_2d);CHKERRQ(ierr);
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
      ierr = PetscDSSetResidual(prob,1,f0_vel_2d,f1_vel_2d);CHKERRQ(ierr);
      ierr = PetscDSSetJacobian(prob,1,1,g0_velvel_2d,NULL,NULL,NULL);CHKERRQ(ierr);
      /* Maybe we need a dynamic jacobian? -- apparently not. */
      ierr = PetscDSSetDynamicJacobian(prob,1,0,g0_dyn_velu_2d,NULL,NULL,NULL);CHKERRQ(ierr);
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
  PetscFE fe, fe_bd, fe_vel;
  PetscDS prob;
  PetscErrorCode ierr;
  PetscBool simplex = user->simplex;
  PetscQuadrature q;
  PetscInt order;
  PetscInt BSpaceOrder;
  PetscSpace space;

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
    /* Define the boundaries */
    
    const PetscInt Ncomp = dim;
    PetscInt components[dim];
    const PetscInt Nfid = 2;
    PetscInt d;
    PetscInt fid[Nfid];  /* fixed faces [numer of fixed faces] */
    const PetscInt Npid = 1;  /* number of pressure loaded faces */
    PetscInt pid[Npid];       /* the ids of the pressure loaded faces */

    PetscInt test[] = {0};
    PetscInt testid[] = {4};
    
    for (d=0;d<dim;d++){
      components[d]=d;
    }
    if(dim==2){
      fid[0] = 5; /* {5}; */ 	/* The fixed faces.*/
      pid[0] = 4; /* {4}; */ 	/* The pressure loaded faces */
    }
    else if(dim==3){
      fid[0] = 6;  /* The fixed face */
      fid[1] = 14;  /* the second fixed face */
      pid[0]= 5;  /* The pressure loaded faces */
    }

    ierr =  DMAddBoundary(cdm, PETSC_TRUE, "fixed", "Face Sets",0, Ncomp, components, (void (*)()) zero_vector, Nfid, fid, user);CHKERRQ(ierr);
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
    

  ierr = PetscFEDestroy(&fe); CHKERRQ(ierr);
  if (user->neumann)  ierr = PetscFEDestroy(&fe_bd); CHKERRQ(ierr);
  if (user->transient) ierr = PetscFEDestroy(&fe_vel); CHKERRQ(ierr);
  return(0);
}




void callthis(){
  SNES snes;			/* nonlinear solver */
  DM dm, distributeddm;			/* problem definition */
  Vec u,r;			/* solution and residual vectors */
  Mat A,J;			/* Jacobian Matrix */
  AppCtx user;			/* user-defined work context */
  PetscErrorCode ierr;
  PetscViewer viewer;
  PetscInt its;

 /* Firing up Petsc */
  /* ierr= PetscInitialize(&argc, &argv,NULL,help);CHKERRQ(ierr); */
  ierr = PetscInitializeFortran();CHKERRQ(ierr);
  ierr = ProcessOptions(PETSC_COMM_WORLD,&user);CHKERRQ(ierr);

  /* importing the gmsh file. Take note that only simplices give meaningful results in 2D at the moment (For which ever reasons) */

  ierr = DMPlexCreateFromFile(PETSC_COMM_WORLD,"longbeam3D_prepped.msh", PETSC_TRUE,&dm);CHKERRQ(ierr);
  ierr = DMGetDimension(dm,&user.dim); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"The problem dimension is %i \n",user.dim);
  ierr = DMPlexDistribute(dm,0,NULL,&distributeddm); CHKERRQ(ierr);
  if (distributeddm) {
    ierr=DMDestroy(&dm);CHKERRQ(ierr);
    dm = distributeddm;
  }
  
  ierr = DMSetFromOptions(dm);CHKERRQ(ierr);
  ierr = DMView(dm,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
 
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
  /* PetscFinalize(); */
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
 


/*   if (user.transient){ */
/*     Vec u_start; */
/*     PetscReal ftime; */
/*     PetscInt nsteps, its; */
/*     TSConvergedReason ConvergedReason; */
/*     PetscBool loadedspring = PETSC_FALSE; */
    
/*     if (loadedspring){  /\* preloading the cantilever so that I can run the TS without loading to see if the problem comes from the neumann conditions in the displacement field *\/ */
/*       user.transient=PETSC_FALSE; */
/*       ierr = PetscPrintf(PETSC_COMM_WORLD,"Creating initial solution for loaded cantilever\n");CHKERRQ(ierr); */
/*       ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr); */
/*       ierr = SNESSetDM(snes,dm);CHKERRQ(ierr); */
/*       ierr = DMSetApplicationContext(dm, &user);CHKERRQ(ierr); */
/*       ierr = SetupDiscretization(dm,&user);CHKERRQ(ierr); */
/*       ierr = DMPlexCreateClosureIndex(dm,NULL);CHKERRQ(ierr); */
      
/*       ierr = DMCreateGlobalVector(dm,&u_start);CHKERRQ(ierr); */
/*       ierr = VecDuplicate(u_start,&r);CHKERRQ(ierr); */
      
      
/*       ierr = DMSetMatType(dm, MATAIJ);CHKERRQ(ierr); */
/*       ierr = DMCreateMatrix(dm, &J);CHKERRQ(ierr); */
/*       A=J; */

      
/*       ierr = DMPlexSetSNESLocalFEM(dm,&user,&user,&user);CHKERRQ(ierr); */

/*       ierr = SNESSetJacobian(snes, A, J, NULL, NULL);CHKERRQ(ierr); */

/*       ierr = SNESSetFromOptions(snes);CHKERRQ(ierr); */
/*       ierr = SNESSolve(snes,NULL,u_start);CHKERRQ(ierr); */

/*       ierr = SNESGetIterationNumber(snes, &its);CHKERRQ(ierr); */
/*       ierr = PetscPrintf(PETSC_COMM_WORLD,"Static loaded state determined. Number of snes iterations %i \n", its);CHKERRQ(ierr); */

/*       ierr = MatDestroy(&J); CHKERRQ(ierr); */
/*       ierr = MatDestroy(&A); CHKERRQ(ierr); */
/*       ierr = VecDestroy(&r); CHKERRQ(ierr); */

/*       user.transient=PETSC_TRUE; */
/*       user.neumann = PETSC_FALSE; */

/*     } */
    
    
/*     PetscPrintf(PETSC_COMM_WORLD,"starting transient analysis. Total time: %g s \n",user.time.total); */
/*     ierr =  TSCreate(PETSC_COMM_WORLD, &ts); CHKERRQ(ierr); */
/*     ierr = TSSetType(ts, TSCN); CHKERRQ(ierr); */
/*     ierr = TSSetDM(ts,dm); CHKERRQ(ierr); */
/*     ierr = DMSetApplicationContext(dm, &user);CHKERRQ(ierr); */
/*     ierr = SetupDiscretization(dm,&user);CHKERRQ(ierr); */
/*     ierr = DMPlexCreateClosureIndex(dm,NULL);CHKERRQ(ierr); */

/*     ierr = DMTSSetBoundaryLocal(dm, DMPlexTSComputeBoundary, &user); CHKERRQ(ierr); */
/*     ierr = DMTSSetIFunctionLocal(dm, DMPlexTSComputeIFunctionFEM, &user); CHKERRQ(ierr); */
/*     ierr = DMTSSetIJacobianLocal(dm, DMPlexTSComputeIJacobianFEM, &user); CHKERRQ(ierr); */


    
/*     ierr = DMCreateGlobalVector(dm, &u); CHKERRQ(ierr); */

/*     if (loadedspring){ */
/*       /\* setting the loaded cantilever static solution as initial condition for the transient analysis *\/ */
/*       u=u_start; */
/*     } */

/*     ierr = TSSetSolution(ts,u); CHKERRQ(ierr); */
    
/*     ierr = TSSetDuration(ts,1000,user.time.total); CHKERRQ(ierr); */
/*     ierr = TSSetExactFinalTime(ts, user.time.total); CHKERRQ(ierr); */
/*     ierr = TSSetInitialTimeStep(ts,0.0, user.time.dt); CHKERRQ(ierr); */
/*     ierr = TSSetFromOptions(ts); CHKERRQ(ierr); */

/*     ierr = TSSolve(ts,NULL); CHKERRQ(ierr); */

/*     /\* ierr = TSGetSolveTime(ts,&ftime); *\/ */
/*     /\* ierr = TSGetTimeStepNumber(ts,&nsteps); CHKERRQ(ierr); *\/ */
/*     /\* ierr = TSGetConvergedReason(ts,&ConvergedReason); CHKERRQ(ierr); *\/ */

/*     ierr = TSView(ts,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr); */
/*     PetscPrintf(PETSC_COMM_WORLD,"Transient analysis finished. \n"); */
/*     /\* ierr = PetscPrintf(PETSC_COMM_WORLD, "%s at time %g after %D steps \n",ConvergedReason,ftime,nsteps); CHKERRQ(ierr); *\/ */

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
