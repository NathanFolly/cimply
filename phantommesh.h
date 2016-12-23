#ifndef PHANTOMMESH_H
#define PHANTOMMESH_H

#include "cimplyobjects.h"

struct PhantomMesh{
  const struct Class * class;  /* first attribute of all objects */
  /* what do I need to construct a phantom mesh?
   - MMS = how big is the mesh
   - vector or array containing the phantomcells
   - routine or function to get the phantom fractions
   - routine or function to find out which cells are phantomcells
   */
  int MMS;  /* same MMS as in SIMMER */
  int IB;
  int JB;
  int nbofvertices;  /* number of vertices */
  struct PhantomCell * phantomcell;
  struct BoundaryVertex * boundary;

  void * (* getPhantomFractions)(void * _self, float ** phantomfractions);
  
  
};



#endif
