/* this program test the capabilities of the phantommesh functionality
 it creates a spherical body inside a cylindrical mesh that is read from a sim05 file.*/
#include "cimplyobjects.h"
#include "boundaryvertex.h"
#include "phantommesh.h"
#include "simmeranalysis.h"
#include "interface.h"

int main(int argc, char *argv[])
{
  void * interface;
  void * simmeranalysis;
  void * femanalysis;

  float * phantomfractions = NULL;
  int i;
  int nboffractions;

  const float radius = 0.7;
  

  interface= new(Interface);
  printf("Interface created\n");
  simmeranalysis = new(SimmerAnalysis);
  printf("Simmer Analysis created\n");
  assign(interface, simmeranalysis);
  printf("Simmer Analysis assigned\n");
  generateTestSphere(interface, radius ,100);
  printf("Test Sphere generated\n");
  update(interface);
  printf("Interface updated\n");
  getPhantomFractions(interface,&phantomfractions);
  for (i=0; i<48; i++){
    printf("Phantomfraction of cell %i =      %f\n",i,phantomfractions[i]);
  }
  delete(interface);
  delete(simmeranalysis);
      
  return 0;
}
