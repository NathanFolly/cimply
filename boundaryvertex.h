#ifndef BOUNDARYVERTEX_H
#define BOUNDARYVERTEX_H

#include "cimplyobjects.h"
#include "geometry.h"

struct BoundaryVertex{
  const struct Geometry _;  /* part of geometry superclass */
  /* float currentposition[3]; */ /* current position of this vertex [x,y,z]*/
};

extern const void * BoundaryVertex;

#endif
