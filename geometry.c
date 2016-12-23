#include "geometry.h"

static void * Geometry_ctor(void * _self, va_list * app){
  struct Geometry * self = _self;
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
