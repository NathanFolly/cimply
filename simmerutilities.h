#ifndef SIMMERUTILITIES_H
#define SIMMERUTILITIES_H

#include <petscsys.h>
#include "simmeranalysis.h"

#undef __FUNCT__
#define __FUNCT__ "findSIMMERCell"
void findSIMMERCell(const void * _self, const PetscReal x[], PetscBool isoutside, PetscInt *Ipos, PetscInt *Jpos, PetscInt *CellNr);

#endif
