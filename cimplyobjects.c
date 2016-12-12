#include <assert.h>
#include "cimplyobjects.h"
#include <stdarg.h>
#include <stdlib.h>


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

int differ(void * self, void * b){
  const struct Class * const * cp = self;
  assert(self && * cp && (*cp) -> differ);
  return (*cp) -> differ(self, b);
}

void * clone(const void * self){
  const struct Class * const * classptr = self;
  assert (self && *classptr && (*classptr)->clone);
    return (*classptr)->clone(self);
}
  
  
