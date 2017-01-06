#include <assert.h>
#include "cimplyobjects.h"
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include "interface.h"
#include "simmsh.h"
#include "femanalysis.h"
#include "geometry.h"
#include "phantomcell.h"
#include "boundaryvertex.h"

void * new(const void * _class, ...){
  const struct Class * class = _class;
  void * p = calloc(1, class -> size);

  assert(p);
  * (const struct Class **) p = class;  /* convert p to a Class pointerpointer and
                                  * assign the class to it. This works
                                  * because the class attribute is the first
                                  * attr of any object. &object is hence
                                  * equal to &object.class */
  if (class->ctor){
    va_list ap;
    va_start (ap, class);
    p = class->ctor(p,&ap);
    va_end(ap);
  }

  return p;
}

void delete(void * self){
  const struct Class ** cp = self;
  if (self && *cp && (*cp) ->dtor){
    self = (*cp)->dtor(self);
  }
  free(self);
}

/* int differ(void * self, void * b){ */
/*   const struct Class * const * cp = self; */
/*   assert(self && * cp && (*cp) -> differ); */
/*   return (*cp) -> differ(self, b); */
/* } */

/* void * clone(const void * self){ */
/*   const struct Class * const * classptr = self; */
/*   assert (self && *classptr && (*classptr)->clone); */
/*     return (*classptr)->clone(self); */
/* } */
  
  
void update(void * self){
  const struct Class ** cp = self;
  assert(self && *cp && (*cp)->update);
  (*cp)->update(self);
}


/* The function getPhantomFractions is defined for interfaces, simmermeshes and phantomcells */

void * getPhantomFractions(void * _self, float ** PhantomFractions){
  /* todo: still a bit ugly. might be good idea to define superclass container
   * that has a getPhantomFractions function and the subclasses interface,
   * phantomcell and simmermesh. this will do for now */
  const struct Class ** cp = _self;  /* the class pointer to find out which class we are dealing with */
  if(*cp==Interface){
    struct Interface * self= _self;
    assert(self->getPhantomFractions);
    self->getPhantomFractions(self, PhantomFractions);
  }
  else if(*cp==PhantomMesh){
    struct PhantomMesh * self= _self;
    assert(self->getPhantomFractions);
    self->getPhantomFractions(self, PhantomFractions);
  }
  /* else if(*cp== PhantomCell){ */
  /*  struct PhantomCell * self = _self; */
  /*  assert(self->getPhantomFractions);
      self->getPhantomFractions(self, PhantomFractions);
    } */
  else{
    fprintf(stderr, "ERROR:: error in getPhantomFractions. First Argument not of expected type. Expected interface, PhantomMesh or phantomCell.\n");
  }

  return 0;
      
}

/* the function prepare is defined for interfaces, simmermeshes*/

void * prepare(void * _self){
  const struct Class ** cp = _self;
  if(*cp == Interface){
    struct Interface * self = _self;
    self->prepare(self);
  }
  else if(*cp == SIMMsh){
    struct SIMMsh * self = _self;
    self->prepare(self);
  }
  else if(*cp == PhantomCell){
    struct PhantomCell * self = _self;
    self->prepare(self);
  }
  else{
    fprintf(stderr, "ERROR:: error in prepare. Argument not of expected class (Interface, Simmermesh).\n");
    return 0;
  }
  return 0;
}

/* the function getposition is defined for all objects of the superclass
 * geometry (Cells, Vertices). getposition(const void * geometryobject, float *
 * position ,const char coordinatesystem)
  geometryobject -- the object we want the current position of (vcenter of
 volume position for non-vertex objects)
  position -- the pointer to the position vector. will be allocated according
 to the number of spatial dimensions related to the object
  coordinatesystem : "cart" / "cylind" -- type of coordinates the vector will
 be given in*/



void * getposition(const void * _self, float ** position, const char * coordsystem){
  const struct BoundaryVertex * self = _self;
  self->getposition((void * )self, position, coordsystem);
  return 0;
}

void * updateposition(void * _self, const float * coords)
{
  struct BoundaryVertex * self = _self;
  assert(self->updateposition);
  self->updateposition(self, coords);
  return 0;
}
      



void * assign(void * _self,void * _b){
  const struct Class ** cp = _self;
  if(*cp == Interface){
    struct Interface * self = _self;
    self->assign(_self,_b);
  }
  else if (*cp == PhantomCell){
    struct PhantomCell * self = _self;
    self->assign(_self, _b);
  }
  else{
    fprintf(stderr,"ERROR:: error in assign, first argument not of expected type (interface or phantomcell)\n");
  }
    return 0;
}

float phantomfraction(void * _self){
  struct PhantomCell * self = _self;
  float pfrac=0;
  self->givePhantomFraction(self, &pfrac);

  return pfrac;
}

void * generateTestSphere(void * _self, float radius, int nvertices){
  const struct Class ** cp = _self;
  if(*cp!=Interface){
    fprintf(stderr,"ERROR:: error in generateTestSphere. First argument not of type Interface.\n");
    return 0;
  }
  struct Interface * self = _self;
  self->phantommesh->generateTestSphere(self->phantommesh, radius , nvertices);
  return 0;
}
