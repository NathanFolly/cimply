/* Residual functions and handcoded Jacobians for the CIMPLY FEM structural code for SIMMER */
#ifndef CIMPLYDF_H
#define CIMPLYDF_H

#include "cimply.h"
#include "simmerutilities.h"
/* declaring The Kronecker Delta */
extern const PetscInt delta2D[2*2];
extern const PetscInt delta3D[3*3];

/* declaring the material parameters */
extern const PetscReal mu;
extern const PetscReal lbda;
extern const PetscReal rho_steel;

void f0_u(PetscInt dim, PetscInt Nf, PetscInt NfAux,
          const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
          const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
          const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
          const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscScalar f0[]);

    
void f0_u_transient(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
                    const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
                    const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
                    const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscScalar f0[]);



void f1_u(PetscInt dim, PetscInt Nf, PetscInt NfAux,
          const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
          const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
          const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
          const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscScalar f1[]);


void f1_u_2d_pstress(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[],
                     const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[],
                     const PetscScalar u_x[], const PetscInt aOff[], const PetscInt aOff_x[],
                     const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                     PetscReal t, const PetscReal x[], PetscScalar f1[]);


void f1_u_ldis(PetscInt dim, PetscInt Nf, PetscInt NfAux,
               const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
               const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
               const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
               const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscScalar f1[]);


void g3_uu_2d(PetscInt dim, PetscInt Nf, PetscInt NfAux,
	      const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], 
	      const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
	      const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
	      const PetscScalar a_x[], PetscReal t, PetscReal u_tShift, const PetscReal x[],
	      PetscScalar g3[]);


void g3_uu_2d_pstress(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                      const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], 
                      const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
                      const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
                      const PetscScalar a_x[], PetscReal t, PetscReal u_tShift, const PetscReal x[],
                      PetscScalar g3[]);

void g1_uu_ldis(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], 
                const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
                const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
                const PetscScalar a_x[], PetscReal t, PetscReal u_tShift, const PetscReal x[],
                PetscScalar g1[]);

void g3_uu_ldis(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], 
                const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
                const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
                const PetscScalar a_x[], PetscReal t, PetscReal u_tShift, const PetscReal x[],
                PetscScalar g3[]);


void g3_uu_3d(PetscInt dim, PetscInt Nf, PetscInt NfAux,
	      const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], 
	      const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
	      const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
	      const PetscScalar a_x[], PetscReal t, PetscReal u_tShift, const PetscReal x[],
	      PetscScalar g3[]);


void f0_u_bd(PetscInt dim, PetscInt Nf, PetscInt NfAux,
             const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
             const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
             const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
             const PetscScalar a_x[], PetscReal t, const PetscReal x[], const PetscReal n[],
             PetscScalar f0[]);


void f0_u_bd_ldis(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                  const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
                  const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
                  const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
                  const PetscScalar a_x[], PetscReal t, const PetscReal x[], const PetscReal n[],
                  PetscScalar f0[]);

    
void f1_u_bd(PetscInt dim, PetscInt Nf, PetscInt NfAux,
             const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
             const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
             const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
             const PetscScalar a_x[], PetscReal t, const PetscReal x[], const PetscReal n[],
             PetscScalar f1[]);
   
    
void g1_uu_bd_2d(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                 const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], 
                 const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
                 const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
                 const PetscScalar a_x[], PetscReal t, PetscReal u_tShift, const PetscReal x[],
                 const PetscReal n[], PetscScalar g1[]);


void f0_vel(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
            const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
            const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
            const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscScalar f0[]);


void f1_vel(PetscInt dim, PetscInt Nf, PetscInt NfAux,
            const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
            const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
            const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
            const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscScalar f1[]);


void g0_velvel(PetscInt dim, PetscInt Nf, PetscInt NfAux,
               const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], 
               const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
               const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
               const PetscScalar a_x[], PetscReal t, PetscReal u_tShift, const PetscReal x[],
               PetscScalar g0[]);


void g0_uvel(PetscInt dim, PetscInt Nf, PetscInt NfAux,
             const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], 
             const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
             const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
             const PetscScalar a_x[], PetscReal t, PetscReal u_tShift, const PetscReal x[],
             PetscScalar g0[]);

void g2_velu(PetscInt dim, PetscInt Nf, PetscInt NfAux,
             const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], 
             const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
             const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
             const PetscScalar a_x[], PetscReal t, PetscReal u_tShift, const PetscReal x[],
             PetscScalar g2[]);


void g0_velu(PetscInt dim, PetscInt Nf, PetscInt NfAux,
             const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], 
             const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
             const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
             const PetscScalar a_x[], PetscReal t, PetscReal u_tShift, const PetscReal x[],
             PetscScalar g0[]);

void f0_vel_bd(PetscInt dim, PetscInt Nf, PetscInt NfAux,
               const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
               const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
               const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
               const PetscScalar a_x[], PetscReal t, const PetscReal x[], const PetscReal n[],
               PetscScalar f0[]);

void f1_vel_bd(PetscInt dim, PetscInt Nf, PetscInt NfAux,
               const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
               const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[],
               const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[],
               const PetscScalar a_x[], PetscReal t, const PetscReal x[], const PetscReal n[],
               PetscScalar f1[]);
#endif
