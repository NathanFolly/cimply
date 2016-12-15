#include "simmsh.h"
#include <stdio.h>

static void * SIMMsh_getPhantomFractions(void * _self, float ** PhantomFractions);

static void * SIMMsh_ctor(void * _self, va_list * app){
  struct SIMMsh * self = _self;
  self->updatecounter = 0;
  self->getPhantomFractions = SIMMsh_getPhantomFractions;
  return self;
}

static void * SIMMsh_dtor(void * _self){
  struct SIMMsh * self = _self;

  return self;
}


static void * SIMMsh_update(void * _self){
  struct SIMMsh * self = _self;
  ++self->updatecounter;
}




static void * SIMMsh_getPhantomFractions(void * _self, float ** PhantomFractions){
  struct SIMMsh * self = _self;
  int i, ncells=10;  /* test purposes */
  
  *PhantomFractions = (float *) calloc(ncells,sizeof(float));
  for(i=0;i<ncells;i++){
    /* temporary: code for test purposes */
    (*PhantomFractions)[i]= self->updatecounter;
    /* end of temporary code */
  }
  return 0;
}

static const struct Class _SIMMsh = {sizeof(struct SIMMsh), SIMMsh_ctor, SIMMsh_dtor, SIMMsh_update};

const void * SIMMsh = &_SIMMsh;
