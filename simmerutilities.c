#include "simmerutilities.h"

void findSIMMERCell(const void * _self, const PetscReal x[], PetscBool isoutside,
                    PetscInt *Ipos, PetscInt *Jpos, PetscInt *CellNr)
{  /* takes coordinate vector. gives back cell that the position is in or
    * closest to */
  const struct SimmerAnalysis * self = _self;
  PetscBool routside=PETSC_FALSE, zoutside=PETSC_FALSE;  /* do the positions fall outside the SIMMER
                                  * grid */
  PetscBool foundI=PETSC_FALSE, foundJ=PETSC_FALSE;
  PetscReal cummulR = 0, cummulZ=0;  /* total R and Z position of the IJ outer
                                     * cell face */
  const PetscReal r = sqrt(x[0]*x[0]+x[1]*x[1]);
  /* Ipos=0; */
  /* Jpos=0; */
  while (!foundI){
    if (*Ipos>=self->IB){
      routside=PETSC_TRUE;
      foundI = PETSC_TRUE;
      *Ipos=*Ipos-1;
      break;
    }
    cummulR = cummulR+self->DRINP[*Ipos];
    if (cummulR>r){
      foundI=PETSC_TRUE;
      break;
    }
    *Ipos=*Ipos+1;
  }
  if (x[2]<0)
  {
    foundJ=PETSC_TRUE;
    zoutside=PETSC_TRUE;
  }
  while (!foundJ){
    if (*Jpos>=self->JB){
      zoutside=PETSC_TRUE;
      foundJ = PETSC_TRUE;
      *Jpos=*Jpos-1;
      break;
    }
    cummulZ = cummulZ + self->DZINP[*Jpos];
    if (cummulZ>x[2]){
      foundJ = PETSC_TRUE;
      break;
    }
    *Jpos=*Jpos+1;
  }
  if(routside || zoutside){
    isoutside = PETSC_TRUE;
    /* PetscPrintf(PETSC_COMM_WORLD,"WARNING: location x=%f,y=%f,z=%f lies not within the mesh of SIMMER\n",x[0],x[1],x[2]); */
  }
  *CellNr = self->IBP2*(*Jpos+1)+*Ipos+1;
}
