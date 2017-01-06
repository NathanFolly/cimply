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

  self->phantommesh = new(PhantomMesh);

  self->getPhantomFractions = Interface_getPhantomFractions;
  self->prepare = Interface_prepare;
  self->assign = Interface_assign;
  
  
  return self;
}

static void * Interface_dtor(void * _self){
  struct Interface * self = _self;
  if(self->phantommesh){
    delete(self->phantommesh);
  }
  return self;
}

static void * Interface_update(void * _self){
  struct Interface * self = _self;
  if(self->FEM){
    assert(self->FEM);
    update(self->FEM);
  }
  /* if(!self->phantommesh){ */
  /*   fprintf(stderr,"ERROR :: Interface has no SIMMER Mesh assigned to it yet\n"); */
  /*   return 0; */
  /* } */
  /* update(self->phantommesh); */
  update(self->phantommesh);
  return 0;
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
    /* TODO: this does not confirm to the concept of information hiding. better create a simmermesh object that both other objects share or functions to assign meshes to objects */
    assert(self->phantommesh);
    self->phantommesh->DRINP = self->simmeranalysis->DRINP;
    self->phantommesh->DZINP = self->simmeranalysis->DZINP;
    self->phantommesh->MMS = self->simmeranalysis->MMS;
    self->phantommesh->JB = self->simmeranalysis->JB;
    self->phantommesh->IB = self->simmeranalysis->IB;
    self->phantommesh->phantomcell = (void **) calloc(self->phantommesh->MMS,sizeof(void *));

  }
  else{
    fprintf(stderr,"Second input argument: unsupported format. Expected FEMAnalysis or simmeranalysis. \n");
  }
  return 0;
}


/* this is ugly. we should be able to call getPhantomFractions on interfaces, SimmerMeshes and PhantomCells. like this it is only possible for interfaces*/
static void * Interface_getPhantomFractions(void * _self, float ** PhantomFractions){
 const struct Interface * self = _self;
 if(!self->phantommesh){
   fprintf(stderr,"ERROR :: Interface has no SIMMER Mesh assigned to it yet\n");
   return 0;
 }
 /* struct PhantomMesh * phantommesh = self->phantommesh; */
  assert(self->class==Interface);
  getPhantomFractions(self->phantommesh,PhantomFractions);
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
  /* self->phantommesh=new(PhantomMesh);  /\* creating a ne simmermesh *\/ */
  /* prepare(self->phantommesh);     /\* preparing the new simmermesh *\/ */

  /* TODO: make this nicer */
  FEMAnalysis_copyFSInterface(self->FEM,self->phantommesh);
  
  return 0;
}

