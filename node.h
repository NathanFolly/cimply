#ifndef NODE_H
#define NODE_H
#include "cimplyobjects.h"


struct Node{                
  const void * class;  /* must be first attribute of all objects */
  const PetscInt DMNodeNbr;  /* number of this done in the DMPLEX numbering
                              * scheme */
  const struct Node * neighbours;  /* neighbouring nodes */
  const PetscReal CoordsUndef[3];  /* undeformed coordinates of the node */
  PetscReal CoordsDef[3];
  PetscReal SimmerCoords[3];  /* simmer uses polar coordinates */
  struct PhantomCell * HomeCell;/* a node has to belong to one phantomcell. and one
                                 * phantomcell alone */
  struct Node * next;  /* since the node belongs to one specific phantomcell,
                        * we can put them in a linked list. next points to
                        * the next node in the linked list */
  
  void (* updateCoords)(void * _self);
  void (* findCell)(void * _self);  /* find the cell this node is in */
};

extern const void * node;

#endif
