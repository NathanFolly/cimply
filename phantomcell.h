#ifndef PHANTOMCELL_H
#define PHANTOMCELL_H
#include <petscsys.h>
#include <petscdmplex.h>
#include "cimplyobjects.h"

struct PhantomCell{
  const void * class;  /* must be the first attribute of all objects */
  /* The count attribute allows us to just reduce the count of the element in on set and delete the element only once the counter reaches zero */
  struct PhantomCell * next;
  unsigned count;         /* number of times this object was added to a set (=
                           * number of interface nodes in this cell) */
  PetscInt SimmerCellNr;  /* The number of the phantom cell according to
                           * SIMMER numbering. Negative integers declare subcells. */
  const PetscReal RLB, RUB, ZLB, ZUB;  /* The cell boundaries (radius upper/lower;
                                  * Zpos upper/lower) */
  const PetscReal Cellvolume;  /* total volume of the cell */
  PetscReal PhantomVolume;     /* total volome of phantom material in this
                                * cell */
  PetscReal PhantomFraction;  /* The volume fraction that is occupied by the
                               * phantom material, ergo not part of the CFD domain */
  struct Node *nodes;  /* nodes inside the phantomcell */
  struct PhantomCell * subcells;  /* subcells of the phantom cell */


  PetscReal (*getPhantomFraction)(void * _self);  /* function pointer to
                                                        * calculate the
                                                        * phantom fraction */
  void (*subdivide)(void * _self, PetscReal Rpos, PetscReal Zpos);  /* to
                                                                     * subdivide
                                                                     * the
                                                                     * cell
                                                                     * at
                                                                     * specified
                                                                     * locations
                                                                     */
  void (*addnode)(void * _self, void * _node);  /* add node to this phantomcell */
};


struct Node{
  const void * class;  /* must be first attribute of all objects */
  const PetscInt DMNodeNbr;  /* number of this done in the DMPLEX numbering
                              * scheme */
  const struct Node * neighbours;  /* neighbouring nodes */
  const PetscReal CoordsUndef[3];  /* undeformed coordinates of the node */
  PetscReal CoordsDef[3];
  PetscInt inSimmerCell;

  void (* updateCoords)(void * _self);
};
  
static struct PhantomCell * Mothercells;


static void * PhantomCell_ctor(void * _self, va_list * app);


PhantomCell * new_PhantomCell(PetscInt SimmerCellNr,DM ParentDM, PetscInt NIBNodes);

void destroy_PhantomCell(PhantomCell *PhC);


#endif
