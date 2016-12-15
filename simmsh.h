#ifndef SIMMSH_H
#define SIMMSH_H

#include "cimplyobjects.h"
struct SIMMsh{
  const struct Class* class;  /* first attribute of every object */
  int updatecounter;   /* test purposes */

  void * (* getPhantomFractions)(void * _self, float ** PhantomFractions);
};

extern const void * SIMMsh;

  
#endif
