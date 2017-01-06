#ifndef FEMANALYSIS_H
#define FEMANALYSIS_H

#include "cimplyobjects.h"
#include "femapplicationcontext.h"
#include "cimplyDF.h"
#include <petscdm.h>
#include <petscdmlabel.h>
#include <petscds.h>
#include <petscdmplex.h>
#include <petscksp.h>
#include <petscsnes.h>
#include <petscts.h>


struct FEMAnalysis{
  const struct Class * class;  /* first attribute of all objects (and classes) */

  SNES snes;			/* nonlinear solver */
  KSP ksp;                        /* the linear sovler context */
  DM dm, distributeddm;			/* problem definition */
  Vec u,r;			/* solution and residual vectors */
  Mat A,J,P;			/* Jacobian Matrix */
  AppCtx user;			/* user-defined work context */
  PetscViewer viewer;
  TS ts;                          /* in case of transient analysis */
  int dim;                        /* dimension of the anlysis */

  const char * immersedboundaryname;   /**/


  void * (* getPhantomFractions) (void * _self, float * PhantomFractions);
  void * (* settimestep) (void * _self, double dt);
  void * (* selectfsinterface) (void * _self, const char interfacename[], const PetscInt interfacelabelid);
  void * (* copyFSInterface) (const void * _self, void * _b);
};

extern const void * FEMAnalysis;

void * settimestep (void * _self, double dt);
void * selectfsinterface(void * _self, const char interfacename[], const PetscInt interfacelabelid);
void * copyFSInterface(const void * _self, void * _b);
  

#endif
