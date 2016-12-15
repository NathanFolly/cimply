#include <petscsys.h>
#include <petscdmplex.h>
#include "interface.h"

 int main(int argc, char *argv[])
{
  void * interface;  /* the interface */
  void * femanalysis;  /* the FEM analysis */
  float * ALPPHK;  /* the phantom fraction vector coming from and going to SIMMER */
  
  /* prior to analysis */
  interface = new(Interface);  /* create a new interface */
  femanalysis = new(FEMAnalysis);  /* create a new FEM Analysis */
  assign(interface, femanalysis);  /* assign the FEM analysis to the SIMMER interface */
  prepare(interface);              /* prepares the interface, reads the sim05 file, sets the SIMMER mesh up tec */

  /* during analysis*/
  pull(interface);  /*get information from SIMMER */
  /* conduct(FEMAnalysis);  /\* conduct FEM Analysis with the updated information from SIMMER *\/ */
  update(interface);     /* update the phantomcells, phantom fractions etc */
  printf("Update succeeded.\n");
  getPhantomFractions(interface, &ALPPHK);
  delete(femanalysis);
  delete(interface);
  printf("number of updates %0.6f \n", ALPPHK[1]);
  
  return 0;
}
