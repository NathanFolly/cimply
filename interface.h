#ifndef INTERFACE_H
#define INTERFACE_H

#include "cimplyobjects.h"
#include "femanalysis.h"
#include "simmsh.h"

/* Interface for the coupling between SIMMER and a FEM analysis (cimply) */

struct Interface{
  const struct Class * class;  /* always needs to be the first attribute of any object */
  struct FEMAnalysis * FEM;    /* the FEM analysis associated with the interface */
  struct SIMMsh * SIMMsh;  /* the SIMMER mesh associated with this analysis */

  void * (* assign)(void * _self, void * _b);
  void * (* getPhantomFractions)(void * _self, float ** PhantomFractions);
  void * (* pull)(void * _self);
  void * (* prepare)(void * _self);
  
};

extern const void * Interface;
#endif