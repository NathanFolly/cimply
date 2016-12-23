#include "boundaryvertex.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

static void * BoundaryVertex_getposition(const void * _self, float ** position, const char * coordinatesystem);

static void * BoundaryVertex_ctor(void * _self, va_list * app){
  struct Geometry * self = _self;
  /* temporary solution: create vertices with position as arument */
  self->currentposition[0] = va_arg(*app, double);
  self->currentposition[1] = va_arg(*app, double);
  self->currentposition[2] = va_arg(*app, double);

  self->getposition=BoundaryVertex_getposition;

  (struct BoundaryVertex *) self;
  return self;
}

static void * BoundaryVertex_dtor(void * _self){
  struct BoundaryVertex * self = _self;
  return self;
}

static void * BoundaryVertex_update(void * _self){
  struct BoundaryVertex * self = _self;
  return 0;
}

static const struct Class _BoundaryVertex = {sizeof(struct BoundaryVertex), BoundaryVertex_ctor, BoundaryVertex_dtor, BoundaryVertex_update};

const void * BoundaryVertex = &_BoundaryVertex;


static void * BoundaryVertex_getposition(const void * _self, float ** position, const char * coordinatesystem){
  /* phi = 0 is the x-z plane, phi in radians */
  const struct Geometry * self = _self;
  float x = self->currentposition[0];
  float y = self->currentposition[1];
  float z = self->currentposition[2];
  
  if (strcmp(coordinatesystem,"cylind")==0){
    *position = (float *) calloc(3,sizeof(float));
    (*position)[0] = sqrt(pow(x,2)+pow(y,2));
    (*position)[1] = atan(y/x);
    (*position)[2] = z;
  }
  else{
    fprintf(stderr,"ERROR :: error when getting the position of a vertex. coordinatesystem not recognized. try cylind or cart. \n");
  }
  return 0;
}

