#ifndef FEMANALYSIS_H
#define FEMANALYSIS_H

#include "cimplyobjects.h"


struct FEMAnalysis{
  const struct Class * class;  /* first attribute of all objects (and classes) */
  const char * somestring;   /* just for test purposes */
  int updatecounter;     /* just for test purposes */


  void * (* getPhantomFractions) (void * _self, float * PhantomFractions);
};

extern const void * FEMAnalysis;

#endif
