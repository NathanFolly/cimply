#include "geometry.h"
#include <stdio.h>
static void * Geometry_getposition(const void * _self, float ** position, const char * coordinatesystem);
static void * Geometry_updateposition(void * self, const float * coords);


static void * Geometry_ctor(void * _self, va_list * app){
  struct Geometry * self = _self;
  self->currentposition[0] = va_arg(*app, double);
  self->currentposition[1] = va_arg(*app, double);
  self->currentposition[2] = va_arg(*app, double);

  self->getposition=Geometry_getposition;
   
  self->updateposition = Geometry_updateposition;
  
  return self;
}

static void * Geometry_dtor(void * _self){
  struct Geometry * self = _self;
  return self;
}

static void * Geometry_update(void * _self){
  struct Geometry * self = _self;
  return 0;
}


static const struct Class _Geometry = {sizeof(struct Geometry), Geometry_ctor, Geometry_dtor, Geometry_update};
const void * Geometry = &_Geometry;


/*  ------------------------ class specific functions ----------------------------*/

static void * Geometry_getposition(const void * _self, float ** position, const char * coordinatesystem){
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
    fprintf(stderr,"ERROR :: error when getting the position of a geometry. coordinatesystem not recognized. try cylind or cart. \n");
  }
  return 0;
}


static void * Geometry_updateposition(void * _self, const float * coords){
  struct Geometry * self = _self;
  int i;
  printf("Point 2. \n");
  for (i=0;i<3;i++)
  {
  self->currentposition[i] = coords[i];
  printf("Point 3. \n");
  }

  return 0;
}
