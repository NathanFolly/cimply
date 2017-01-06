#ifndef BOUNDARYVERTEX_H
#define BOUNDARYVERTEX_H

#include "cimplyobjects.h"
#include "geometry.h"
#include <petscdm.h>
#include <petscdmlabel.h>

struct BoundaryVertex{
  const struct Class * class;  /* first attributa of any object */
  /* float currentposition[3]; */ /* current position of this vertex [x,y,z]*/
  /* void * (* updateposition)(void * self, const float * coords); */

  float currentposition[3];
  
  int vertexnumber;


  void * (* getposition)(const void * _self, float ** position, const char * coordinatesystem );
  void * (* updateposition)(void * _self, const float * position);
};


extern const void * BoundaryVertex;

int vertexnumber(const void * _self);
#endif
