#ifndef BOUNDARYVERTEX_H
#define BOUNDARYVERTEX_H

#include "cimplyobjects.h"
#include "geometry.h"
#include <petscdm.h>
#include <petscdmlabel.h>

struct BoundaryVertex{
  const struct Geometry _;  /* part of geometry superclass */
  /* float currentposition[3]; */ /* current position of this vertex [x,y,z]*/
  void * (* updateposition)(void * self, const float * coords);
};


extern const void * BoundaryVertex;


#endif
