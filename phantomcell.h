#ifndef PHANTOMCELL_H
#define PHANTOMCELL_H
#include <petscsys.h>
#include <petscdmplex.h>
#include "cimplyobjects.h"

struct PhantomCell{
  const void * class;  /* must be the first attribute of all objects */
  /* The count attribute allows us to just reduce the count of the element in on set and delete the element only once the counter reaches zero */
  struct PhantomCell * next;  /* nect phantomcell in ring */
  struct PhantomCell * ring;  /* ring this phantomcell belongs to (highest
                               * level PC are = SimmerCells, then first level
                               * of subcells ...) */
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
  struct nodenbr *nodes;  /* nodes inside the phantomcell */
  struct PhantomCell * subcells;  /* subcells of the phantom cell */


  PetscReal (*getPhantomFraction)(void * _self);  /* function pointer to
                                                        * return the
                                                        * phantom fraction */
  void (* calculatePhantomFraction)(void * _self);                 /* calculate
                                                                    * the
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


struct  nodenbr{  /* for linked list of nodes that belong to the respective phantomcell */
  const PetscInt Nodenumber;
  struct nodenbr * next;
};

static struct PhantomCell * Mothercells;


static void * PhantomCell_ctor(void * _self, va_list * app);


PhantomCell * new_PhantomCell(PetscInt SimmerCellNr,DM ParentDM, PetscInt NIBNodes);

void destroy_PhantomCell(PhantomCell *PhC);


#endif
