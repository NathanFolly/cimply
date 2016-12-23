#ifndef SIMMERANALYSIS_H
#define SIMMERANALYSIS_H

#include "cimplyobjects.h"
#include <petscsys.h>

struct SimmerAnalysis{
  const struct Class * class;  /* first attribute of every object */

  int IB;  /* Number of cells in radial direction */
  int JB;  /* Number of cells in vertical direction */
  int IBP2;  /* IB plus 2 shadow cells */
  int JBP2;  /* JB plus 2 shadow cells */
  int MMS;   /* The total length of a K-cell scalar data array */
  double * DRINP;  /*Array of cell heigths in radial direction  */
  double * DZINP;  /*Array of cell heights in z direction */
  double TWFIN; /* The total time of the analysis */
  /* simmer variables */
  double * PK;  /* The pressure in each K-cell */

};

extern const void * SimmerAnalysis;

#endif
