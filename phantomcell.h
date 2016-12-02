#ifndef PHANTOMCELL_H
#define PHANTOMCELL_H
#include <petscsys.h>

typedef struct {
  unsigned count;  /* number of objects in this set */
} Set;

typedef struct {
  /* The count attribute allows us to just reduce the count of the element in on set and delete the element only once the counter reaches zero */
  unsigned count;         /* number of times this object was added to a set */
  Set *in;            /* set this object is part of */
  PetscInt SimmerCellNr;  /* The number of the phantom cell according to
                           * SIMMER numbering */
  PetscReal RLB, RUB, ZLB, ZUB;  /* The cell boundaries (radius upper/lower;
                                  * Zpos upper/lower) */
  PetscInt IntFace;           /* The face over which we integrate */
  PetscInt quadnodes[3];     /* The quadrature node numbers. */
  PetscInt Ndaughters;       /* the number of daughter cells (see below) */
  /* PhantomCell *daughtercell;  /\* we may need to subdivide some phantom cells for the sake of integration *\/ */
  DM ParentDM;  /* The DM to which the IB nodes belong*/
  PetscReal PhantomFraction;  /* The volume fraction that is occupied by the
                               * phantom material, ergo not part of the CFD domain */

  PetscReal (*getPhantomFraction)();  /* function pointer to
                                                        * calculate the
                                                        * phantom fraction */
  PetscErrorCode (*updateQuadNodes)();  /* function to update the quadrature nodes */
  
}PhantomCell;

PhantomCell * new_PhantomCell(PetscInt SimmerCellNr,DM ParentDM, PetscInt NIBNodes);

void destroy_PhantomCell(PhantomCell *PhC);


#endif
