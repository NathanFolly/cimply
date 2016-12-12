#include <assert.h>
#include <stdio.h>
#include "testobject.h"

static int testobject_givevalue(const void * _self){
  const struct testobject * self = _self;
  assert(self);
  return self->somenumber;
        
}

static void * testobject_ctor(void * _self, va_list * app){
  struct testobject * self = _self;
  const int number = va_arg( * app, const int );
  self -> somenumber = number;
  self -> givevalue = testobject_givevalue;
  return self;
}

static void * testobject_dtor(void * _self){
  struct testobject * self = _self;
  return self;
}

static void * testobject_clone(const void * _self){
  const struct testobject * self = _self;
  return new(testobject,self->somenumber);
}

static int testobject_diff(const void * _self,const void * _b){
  const struct testobject * self = _self;
  const struct testobject * b = _b;
  if (self==b) {
    return 0;
  }
  if (!b ||b->class!=testobject){
    return 1;
  }
  return (self->somenumber!=b->somenumber);
}



static const struct Class _testobject = {sizeof(struct testobject),testobject_ctor,testobject_dtor,testobject_clone,testobject_diff};

const void * testobject = & _testobject;

void * getvalue(const void * _self){
  const struct testobject * self = _self;
  assert(self->givevalue);
  return self->givevalue(self);
}


