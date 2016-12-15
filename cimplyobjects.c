#include <assert.h>
#include "cimplyobjects.h"
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>

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
  else if(*cp==SIMMsh){
    struct SIMMsh * self= _self;
    assert(self->getPhantomFractions);
    self->getPhantomFractions(self, PhantomFractions);
  }
  /* else if(*cp== PhantomCell){ */
  /*  struct PhantomCell * self = _self; */
  /*  assert(self->getPhantomFractions);
      self->getPhantomFractions(self, PhantomFractions);
    } */
  else{
    fprintf(stderr, "ERROR:: error in getPhantomFractions. First Argument not of expected type. Expected interface, SIMMsh or phantomCell.");
  }

  return 0;
      
}
