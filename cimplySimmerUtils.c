
#undef __FUNCT__
#define __FUNCT__ "findSIMMERCell"
void findSIMMERCell(SimmerData SiDat, const PetscReal x[], PetscBool isoutside, PetscInt *Ipos, PetscInt *Jpos, PetscInt *CellNr)
{  /* takes coordinate vector. gives back cell that the position is in or
    * closest to */
  PetscBool routside=PETSC_FALSE, zoutside=PETSC_FALSE;  /* do the positions fall outside the SIMMER
                                  * grid */
  PetscBool foundI=PETSC_FALSE, foundJ=PETSC_FALSE;
  PetscReal cummulR = 0, cummulZ=0;  /* total R and Z position of the IJ outer
                                     * cell face */
  const PetscReal r = sqrt(x[0]*x[0]+x[1]*x[1]);
  /* Ipos=0; */
  /* Jpos=0; */
  while (!foundI){
    if (*Ipos>=SiDat.IB){
      routside=PETSC_TRUE;
      foundI = PETSC_TRUE;
      *Ipos=*Ipos-1;
      break;
    }
    cummulR = cummulR+SiDat.DRINP[*Ipos];
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
    if (*Jpos>=SiDat.JB){
      zoutside=PETSC_TRUE;
      foundJ = PETSC_TRUE;
      *Jpos=*Jpos-1;
      break;
    }
    cummulZ = cummulZ + SiDat.DZINP[*Jpos];
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
  *CellNr = SiDat.IBP2*(*Jpos+1)+*Ipos+1;
}

#undef __FUNCT__
#define __FUNCT__ "getSimmerPressure"

void getSimmerPressure(SimmerData SiDat, const PetscReal x[], PetscReal *pressure)
{
  PetscInt Ipos=0, Jpos=0, CellNr=0;  /* The i and j indices of the cell from which we should
                       * read the pressure information */
  PetscBool isoutside;

  findSIMMERCell(SiDat,x,isoutside,&Ipos,&Jpos,&CellNr);
  *pressure = SiDat.PK[CellNr];  /* +1 at the end is missing
                                                  * because c array indices
                                                  * start from 0 */
}
