#ifndef FEMAPPLICATIONCONTEXT_H
#define FEMAPPLICATIONCONTEXT_H
#include <petscsys.h>

typedef enum {RUN_STANDALONE, RUN_COUPLED} RunType;
typedef enum {NONE, LABELS, BOUNDARIES, SOLUTION} Visualization;

typedef struct
{
  PetscReal mu;
  PetscReal lbda;
}Material_Type;

typedef struct  /* The time context */
{
  PetscReal total;  /* total time in seconds */
  PetscReal dt;                      /* the current time step size */
  PetscInt iter;                     /* number of iterations so far */
  /* PetscBool halfstep;  /\* Is the timestep a halfstep? *\/ */
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
  PetscBool planestress;  /* 2D analysis type: plane stress */
  TSctx time;
  Material_Type material;
  PetscViewer viewer;
  
}AppCtx;


#endif
