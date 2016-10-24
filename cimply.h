#ifndef CIMPLY_H
#define CIMPLY_H
#include <petscsys.h>
/* Create a data structure for the data from SIMMER  */

typedef struct {
  PetscInt IB;  /* Number of cells in radial direction */
  PetscInt JB;  /* Number of cells in vertical direction */
  PetscInt IBP2;  /* IB plus 2 shadow cells */
  PetscInt JBP2;  /* JB plus 2 shadow cells */
  PetscInt MMS;   /* The total length of a K-cell scalar data array */
  PetscReal *DRINP;  /*Array of cell heigths in radial direction  */
  PetscReal *DZINP;  /*Array of cell heights in z direction */
  PetscReal *PK;  /* The pressure in each K-cell */
  PetscReal TWFIN; /* The total time of the analysis */
}SimmerDataStruct;

extern SimmerDataStruct SimmerData;
#endif 
