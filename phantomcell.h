#ifndef PHANTOMCELL_H
#define PHANTOMCELL_H
#include <petscsys.h>
#include <petscdmplex.h>
#include "cimplyobjects.h"

struct PhantomCell{
  const struct Class * class; /* part of Geometry superclass */
  /* The count attribute allows us to just reduce the count of the element in on set and delete the element only once the counter reaches zero */
  unsigned vertexcount;               /* number of vertices in this phantomcell*/
  PetscReal RLB, RUB, ZLB, ZUB;  /* The cell boundaries (radius upper/lower;
                                  * Zpos upper/lower) */
  float cellvolume;  /* total volume of the cell */
  float phantomvolume;     /* total volome of phantom material in this
                                * cell */
  float phantomfraction;  /* The volume fraction that is occupied by the
                               * phantom material, ergo not part of the CFD domain */
  struct Node * vertexring;  /* nodes inside the phantomcell */
  struct PhantomCell ** subcell;  /* subcells of the phantom cell */


  void (*givePhantomFraction)(void * _self, float*  phantomfraction);  /* function pointer to
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
  void (*assign)(void * _self, void * _boundaryvertex);  /* add node to this phantomcell */

  void (*prepare)(void * _self);  /* prepare for calculation of the phantom fraction: subdivide phantomcell, order boundary vertices or whatever */
  
};


struct Node{  /* for linked list*/
  const void * content;
  struct Node * next;
};

extern const void * PhantomCell;

/* static void * PhantomCell_ctor(void * _self, va_list * app); */


/* PhantomCell * new_PhantomCell(PetscInt SimmerCellNr,DM ParentDM, PetscInt NIBNodes); */

/* void destroy_PhantomCell(PhantomCell *PhC); */


#endif
