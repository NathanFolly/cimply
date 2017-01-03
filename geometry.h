#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "cimplyobjects.h"

struct Geometry{
  const struct Class * class;  /* first attribute of all objects  */
  float currentposition[3];
  
  void * (* getposition)(const void * _self, float ** position, const char * coordinatesystem );

};

extern const void * Geometry;
#endif