/* This program tests the capabilities of the femanalysis class */
#include "femanalysis.h"
#include "simmeranalysis.h"
#include "cimply.h"
#include <petscsys.h>

  void * simmeranal;


int main(int argc, char *argv[])
{
  void * fem;
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc, &argv, NULL, NULL); CHKERRQ(ierr);

  simmeranal = new(SimmerAnalysis);
  setUniformSimmerPressure(simmeranal, 0.5);
  fem = new(FEMAnalysis, "SHammer_thin_fine.msh");
  update(fem);
  delete(simmeranal);
  delete(fem);

  PetscFinalize();
  
  return 0;
}
