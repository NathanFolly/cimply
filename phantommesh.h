#ifndef PHANTOMMESH_H
#define PHANTOMMESH_H

#include "cimplyobjects.h"
#include "phantomcell.h"
#include "boundaryvertex.h"
#include <petscdm.h>
#include <petscdmlabel.h>

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
  double * DRINP;  /*Array of cell heigths in radial direction  */
  double * DZINP;  /*Array of cell heights in z direction */
  int nbofvertices;  /* number of vertices */
  void ** phantomcell;
  void ** boundary;

  /* information necessary to update the vertex locations */

  DM dm;
  Vec coords,u;
  PetscSection CoordSect, DispSection;

  char * immersedboundaryname;

  /* DM * dmptr;  /\* the Petsc DM this phantom mesh is associated with *\/ */
  /* void * sltnptr;  /\* pointer to the solution of the FEM analysis *\/ */
  

  void * (* getPhantomFractions)(void * _self, float ** phantomfractions);
  void * (* generateTestSphere)(void * _self, float radius, int nvertices);
  
};

extern const void * PhantomMesh;
void * PhantomMeshSetDM(void * _self, DM dm);
void * PhantomMeshSetSolution(void * _self, Vec solution);
void * PhantomMeshSetupInterface(void * _self);



#endif
