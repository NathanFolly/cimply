#include <petscsys.h>
#include <petscdmplex.h>
#include "interface.h"

 int main(int argc, char *argv[])
{
  void * interface;  /* the interface */
  void * femanalysis;  /* the FEM analysis */
  void * simmeranalysis;  /* the simmer analysis */
  float * ALPPHK;  /* the phantom fraction vector coming from and going to SIMMER */
  
  /* prior to analysis */
  interface = new(Interface);  /* create a new interface */
  femanalysis = new(FEMAnalysis);  /* create a new FEM Analysis */
  simmeranalysis = new(SimmerAnalysis);  /* create new Simmer Analysis (just a collection of simmer variables. the actual simmer code reamains the seperate fortran code.) */
  assign(interface, femanalysis);  /* assign the FEM analysis to the SIMMER interface */
  assign(interface, simmeranalysis);  /* assign the simmeranalysis to the interface */
  prepare(interface);              /* prepares the interface, sets up the SIMMER mesh, etc. */

  /* during analysis*/
  /* run(simmeranalysis);  */ /*get information from SIMMER */
  /* conduct(FEMAnalysis);  /\* conduct FEM Analysis with the updated information from SIMMER *\/ */
  update(interface);     /* update the phantomcells, phantom fractions etc */
  printf("Update succeeded.\n");
  getPhantomFractions(interface, &ALPPHK);
  delete(femanalysis);
  delete(interface);
  printf("number of updates %0.6f \n", ALPPHK[1]);
  
  return 0;
}
