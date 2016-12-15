#include "femanalysis.h"

static void * FEMAnalysis_ctor(void * _self, va_list * app){
  struct FEMAnalysis * self = _self;
  self->somestring = "blablabla";
  self->updatecounter = 0;

  return self;
}

static void * FEMAnalysis_dtor(void * _self){
  struct FEMAnalysis * self = _self;
  return self;
}

static void * FEMAnalysis_update(void * _self){
  struct FEMAnalysis * self = _self;
  ++self->updatecounter;
}



static const struct Class _FEMAnalysis = {sizeof(struct FEMAnalysis), FEMAnalysis_ctor, FEMAnalysis_dtor,FEMAnalysis_update};

const void * FEMAnalysis = &_FEMAnalysis;
