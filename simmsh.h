#ifndef SIMMSH_H
#define SIMMSH_H

#include "cimplyobjects.h"
#include <petscsys.h>
struct SIMMsh{
  const struct Class* class;  /* first attribute of every object */
  PetscInt IB;  /* Number of cells in radial direction */
  PetscInt JB;  /* Number of cells in vertical direction */
  PetscInt IBP2;  /* IB plus 2 shadow cells */
  PetscInt JBP2;  /* JB plus 2 shadow cells */
  PetscInt MMS;   /* The total length of a K-cell scalar data array */
  PetscReal * DRINP;  /*Array of cell heigths in radial direction  */
  PetscReal * DZINP;  /*Array of cell heights in z direction */
  PetscReal * PK;  /* The pressure in each K-cell */
  PetscReal TWFIN; /* The total time of the analysis */

  void * (* getPhantomFractions)(void * _self, float ** PhantomFractions);
  void * (* prepare)(void * _self);
  void * (* assign)(void * _self, void * _b);
};

extern const void * SIMMsh;

  
#endif
