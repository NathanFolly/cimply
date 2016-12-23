#include "interface.h"
#include "femanalysis.h"
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

static void * Interface_getPhantomFractions(void * _self, float ** PhantomFractions);
static void * Interface_prepare(void * _self);
static void * Interface_assign(void * _self, void * _b);

static void * Interface_ctor(void * _self, va_list * app){
  struct Interface * self = _self;

  self->getPhantomFractions = Interface_getPhantomFractions;
  self->prepare = Interface_prepare;
  self->assign = Interface_assign;
  
  return self;
}

static void * Interface_dtor(void * _self){
  struct Interface * self = _self;
  if(self->SIMMsh){
    delete(self->SIMMsh);
  }
  return self;
}

static void * Interface_update(void * _self){
  struct Interface * self = _self;
  assert(self->FEM);
  update(self->FEM);
  if(!self->SIMMsh){
    fprintf(stderr,"ERROR :: Interface has no SIMMER Mesh assigned to it yet\n");
    return 0;
  }
  update(self->SIMMsh);
}

static const struct Class _Interface = {sizeof(struct Interface), Interface_ctor, Interface_dtor, Interface_update};

const void * Interface = &_Interface;
/* ------------------- End of Class definition ------------- */








/* -------------------- Class specific functions: ---------- */

static void * Interface_assign(void * _self,void * _b){
  struct Interface * self = _self;
  const struct Class ** cp = _self;
  const struct Class ** b = _b;  /* class pointerpointer (first element in every class) */
  assert(*cp = Interface);
  if(*b == FEMAnalysis){         /* check if second argument has correct class */
    struct FEMAnalysis * FEM = _b;
    self->FEM = FEM;
  }
  else if(*b == SimmerAnalysis){
    struct SimmerAnalysis * simmeranalysis =_b;
    self->simmeranalysis = simmeranalysis;
  }
  else{
    fprintf(stderr,"Second input argument: unsupported format. Expected FEMAnalysis or simmeranalysis. \n");
  }
  return 0;
}


/* this is ugly. we should be able to call getPhantomFractions on interfaces, SimmerMeshes and PhantomCells. like this it is only possible for interfaces*/
static void * Interface_getPhantomFractions(void * _self, float ** PhantomFractions){
 const struct Interface * self = _self;
 if(!self->SIMMsh){
   fprintf(stderr,"ERROR :: Interface has no SIMMER Mesh assigned to it yet\n");
   return 0;
 }
 struct SIMMsh * sm = self->SIMMsh;
  assert(self->class==Interface);
  getPhantomFractions(sm,PhantomFractions);
  return 0;
}

void * pull(void * _self) {
  printf("pulling \n");
  return 0;
}

static void * Interface_prepare(void * _self) {
  struct Interface * self = _self;
  assert(self->class==Interface);
  printf("preparing the interface... \n");
  self->SIMMsh=new(SIMMsh);  /* creating a ne simmermesh */
  prepare(self->SIMMsh);     /* preparing the new simmermesh */
  return 0;
}

