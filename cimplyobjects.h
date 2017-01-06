/* cimply mimics object oriented programming. This header file contains basic
 * data types and function declarations necessary to implement this. */
#ifndef CIMPLYOBJECTS_H
#define CIMPLYOBJECTS_H
#include <stdarg.h>
#include <stdlib.h>


struct Class{
  size_t size;
  void * (* ctor) (void *self, va_list *app);  /* constructor */
  void * (* dtor) (void * self);               /* destructor */
  void * (* update) (void * self);
  /* int (* differ) (const void * self, const void * b); */
};

/* struct Set { */
/*   const void *class;  /\* must always be the first attribute in any object *\/ */
/*   unsigned count;  /\* number of objects in this set *\/ */
/* }; */


void * new(const void * class, ...);
void delete(void * self);
void update(void * self);
/* void * clone(void * self); */
/* int differ(void * self, void * b); */

float phantomfraction(void * _self);
float testfunction(void);
void * generateTestSphere(void * _self, float radius, int nvertives);
void * getposition(const void * _self, float ** position, const char * coordsystem);
void * updateposition(void * _self, const float * coords);    

#endif
