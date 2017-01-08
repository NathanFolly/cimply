#include "cimply.h"
#include "cimplyobjects.h"
#include "simmeranalysis.h"
#include "femanalysis.h"
#include "interface.h"
#include <petscsys.h>

void * simmeranal;
void * fem;
void * interface;
/* double * phantomfractions=NULL; */


void cimplysetup()
{
  PetscErrorCode ierr;

  ierr = PetscInitializeFortran();CHKERRQ(ierr);

  interface = new(Interface);
  simmeranal = new(SimmerAnalysis);
  assign(interface,simmeranal);
  fem=new(FEMAnalysis,"SHammer_thin.msh");
  selectfsinterface(fem,"Face Sets",2);
  assign(interface,fem); 
}

void runcimply(const int SimmerMMS, const PetscReal PK[], const int iter, const PetscReal dt, PetscReal phantomfractions[])
{
  int MMS;
  getMMS(simmeranal, &MMS);
  double pressure[MMS];
  double * phfractions = NULL;
  int i;
  
  for (i=0;i<MMS;i++)
    {
      pressure[i] = PK[i];
    }
  setpressure(simmeranal,pressure);
  setiteration(fem,iter);
  settimestep(fem,dt);
  update(interface);
  getphantomfractions(interface,&phfractions);
  for(i=0;i<MMS;i++)
    {
      phantomfractions[i] = (PetscReal) phfractions[i];
    }
  free(phfractions);

}

void cimplyfinalize()
{
  delete(interface);
  delete(fem);
  delete(simmeranal);

}
