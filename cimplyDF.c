/* functions for creating the Jacobians and residuals of the nonlinear finite element analysis */

#include "cimplyDF.h"

/* defining the kronecker deltas */
const PetscInt delta2D[2*2] = {1,0,0,1};
const PetscInt delta3D[3*3] = {1,0,0,0,1,0,0,0,1};


/* defining the material parameters */
const PetscReal mu =76.923076923, lbda=115.384615385, rho_steel = 7850E-9;



void f0_u(PetscInt dim, PetscInt Nf, PetscInt NfAux,
          const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
          const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
          const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
          const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscScalar f0[])
{
  const PetscInt Ncomp = dim;
  PetscInt comp;
  for(comp=0;comp<Ncomp;comp++) f0[comp]=0.0;
    
}




void f0_u_transient(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
                    const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
                    const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
                    const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscScalar f0[])
{
  const PetscInt Ncomp = dim;
  /* const PetscReal rho_steel = 7850E-9;  /\* Density of steel needed for the inertia term *\/ */
  PetscInt comp;
  for(comp=0;comp<Ncomp;comp++) f0[comp]= 1.0*rho_steel*u_t[uOff[1]+comp];
}


void f1_u(PetscInt dim, PetscInt Nf, PetscInt NfAux,
          const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
          const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
          const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
          const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscScalar f1[])
{
  /* 3-D or Two-dimensional plain strain formulation */
  const PetscInt Ncomp = dim;
  /* const PetscReal mu =76.923076923, lbda=115.384615385; */
  
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

void f1_u_2d_pstress(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                     const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
                     const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
                     const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
                     const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscScalar f1[])
{
  /* Teo-dimensional plain strain formulation */
  const PetscInt Ncomp = dim;
  /* const PetscReal mu =76.923076923, lbda=115.384615385; */
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
  /* const PetscReal mu =76.9, lbda=115.4; */
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
  /* const PetscReal mu =76.923076923, lbda=115.384615385; */
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
  /* const PetscReal mu =76.923076923, lbda=115.384615385; */
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
  /* const PetscReal mu =76.9, lbda=115.4; */
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

  /* const PetscReal mu = 76.9, lbda = 115.4; */
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
  /* const PetscReal mu =76.923076923, lbda=115.384615385; */
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

  const PetscScalar traction[] = {0.0, 0.00, 0.0};
  PetscInt comp;
  PetscReal pressure = 0.000;

  /* get the pressure from the Simmer data according to the location of the
   * boundary element */
  getSimmerPressure((void * ) simmeranal,x,&pressure);
    
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
  /* const PetscReal mu =76.9, lbda=115.4; */
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
    
  /* const PetscReal mu=76.9, lbda=115.4; */
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


void f0_vel(PetscInt dim, PetscInt Nf, PetscInt NfAux,
               const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
               const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
               const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
               const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscScalar f0[])
{
  const PetscInt Ncomp = dim;
  PetscInt comp, d;
  /* We're solving the equation vel - du/dt = 0 so: */
  for (comp=0;comp<Ncomp;comp++){
    f0[comp] = u[uOff[1]+comp]-u_t[comp];
    /* This is for the material derivative */
    /* for (d=0;d<dim;d++){ */
    /*   f0[comp]-=u_t[comp]*u_x[d*dim+comp]; */
    /* } */
  }
}


void f1_vel(PetscInt dim, PetscInt Nf, PetscInt NfAux,
               const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
               const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
               const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
               const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscScalar f1[])
{
  const PetscInt Ncomp = dim;
  PetscInt comp,d;

  for (comp=0;comp<Ncomp;comp++){
    for (d=0;d<dim;d++){
      f1[comp*Ncomp+d] = 0.0;
    }
  }
}

void g0_velvel(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                  const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], 
                  const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
                  const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
                  const PetscScalar a_x[], PetscReal t, PetscReal u_tShift, const PetscReal x[],
                  PetscScalar g0[]){
  PetscInt Ncomp = dim;
  PetscInt i;

  for (i=0;i<Ncomp;i++){
    g0[i*Ncomp+i]= 1.0;
  }
  
}

void g0_uvel(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                  const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], 
                  const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
                  const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
                  const PetscScalar a_x[], PetscReal t, PetscReal u_tShift, const PetscReal x[],
                  PetscScalar g0[]){
  /* const PetscReal rho_steel = 7850E-9;  /\* Density of steel needed for the inertia term *\/ */
  PetscInt Ncomp = dim;
  PetscInt i;

  for (i=0;i<Ncomp;i++){
    g0[i*Ncomp+i]= 1.0*rho_steel;
  }
  /* PetscPrintf(PETSC_COMM_WORLD,"t = %f\n",t); */
  /* u_tShift = 1.0; */
  
}

void g2_velu(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                  const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], 
                  const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
                  const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
                  const PetscScalar a_x[], PetscReal t, PetscReal u_tShift, const PetscReal x[],
                  PetscScalar g2[]){

  const PetscInt NComp=dim;
  PetscInt i,j,k;
  /* This is for the material derivative */
  for (i=0;i<dim;i++){
    for (j=0; j < dim; ++j) {
      for (k=0; k<dim; k++){
        g2[(i*dim+j)*dim+k] = -u_t[i];
      }
    }
  }

}


void g0_velu(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                  const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], 
                  const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
                  const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
                  const PetscScalar a_x[], PetscReal t, PetscReal u_tShift, const PetscReal x[],
                  PetscScalar g0[]){

  PetscInt Ncomp = dim;
  PetscInt i,j;

  for (i=0;i<Ncomp;i++){
    /* for (j=0;j<dim;j++){ */
    /*   g0[i*Ncomp+i]-=u_x[i*Ncomp+j];  /\* Lagrangian formulation -> substanital derivative D/Dt *\/ */
    /* } */
    g0[i*Ncomp+i]=-1.0;
  }
}

void f0_vel_bd(PetscInt dim, PetscInt Nf, PetscInt NfAux,
        const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
        const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
        const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
        const PetscScalar a_x[], PetscReal t, const PetscReal x[], const PetscReal n[],
        PetscScalar f0[])
{
  const  PetscInt Ncomp = dim;
  PetscInt comp,d;
  
  for (comp=0;comp<Ncomp;++comp){
    f0[comp] = 0.0;
  }
}

void f1_vel_bd(PetscInt dim, PetscInt Nf, PetscInt NfAux,
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
