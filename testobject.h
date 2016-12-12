#ifndef TESTOBJECT_H
#define TESTOBJECT_H
#include "cimplyobjects.h"

struct testobject{
  const void * class;
  int somenumber;
  int (* givevalue)(const void * _self);
 };

extern const void * testobject;

void * getvalue(const void * _self);

#endif
