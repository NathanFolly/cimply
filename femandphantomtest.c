/* This program tests the capabilities of femanalysis and phantommesh in combination */

#include "femanalysis.h"
#include "phantommesh.h"
#include "simmeranalysis.h"
#include "interface.h"
#include "cimply.h"
#include <petscsys.h>

void * simmeranal;

int main(int argc, char *argv[])
{
  void * fem;
  void * interface;
  PetscErrorCode ierr;
  double * phantomfractions=NULL;

  ierr = PetscInitialize(&argc, &argv, NULL, NULL); CHKERRQ(ierr);

  interface = new(Interface);

  simmeranal = new(SimmerAnalysis);
  assign(interface, simmeranal);
  setUniformSimmerPressure(simmeranal, 0.1);
  fem = new(FEMAnalysis, "SHammer_thin.msh");
  selectfsinterface(fem,"Face Sets",2);
  assign(interface, fem);
  prepare(interface);
  update(interface);
  printf("done updating \n");
  getPhantomFractions(interface,&phantomfractions);
  printf("got the phantom fractions \n");
  int i;
  for (i=0; i < 39; ++i)
  {
    printf("Phatomfraction %i       %f\n",i,phantomfractions[i]) ;
  }
  setUniformSimmerPressure(simmeranal, 0.7);
  update(interface);
  getPhantomFractions(interface,&phantomfractions);

  for (i=0; i < 39; ++i)
  {
    printf("Phatomfraction %i       %f\n",i,phantomfractions[i]) ;
  }
  delete(interface);
  delete(fem);
  delete(simmeranal);
  free(phantomfractions);

  PetscFinalize();
  
  return 0;
}
