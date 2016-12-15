#include <assert.h>
#include <stdio.h>
#include <petscdm.h>
#include <petscdmlabel.h>
/* #include <petscdmplex.h> */
#include "phantomcell.h"
#include "cimplySimmerUtils.h"


static void * PhantomCell_ctor(void * _self, va_list *app ){
  struct PhantomCell * self = _self;
  /* phantom cells are in a linked list. If the SIMMER cell in question is
   * already a phantom cell, we don't create a new phantomcell but instead up the
   * count of the existing one

   call: new(PhantomCell,struct PhantomCell * partof, PetscReal RLB,
   PetscReal RUB, PetscReal ZLB, PetscReal ZUB)
   
   partof: collection this PhantomCell belongs to
   RLB: Radial lower boundary of the cell
   RUB: radial upper boundary of the cell
   ZLB: Z-direction lower boundary of the cell
   ZUB: Z-direction upper boundary of the cell
  */

  
  const struct PhantomCell * ring = va_arg(*app, struct PhantomCell *);
  if(ring){
    struct PhantomCell * b = ring;
    const PetscReal RLB = va_arg(*app, PetscReal);
    const PetscReal RUB = va_arg(*app, PetscReal);
    const PetscReal ZLB = va_arg(*app, PetscReal);
    const PetscReal ZUB = va_arg(*app, PetscReal);
    
    while(b->next!=ring){
      b=b->next;
      if(RLB==b->RLB && RUB==b->RUB && ZLB==b->ZLB && ZUB==b->ZUB){
        ++b->count;
        free(self);
        return b;
      }
    }
  }
  else{
    ring = self;
  }
  self->next = ring->next, ring->next=self;
  self->ring = ring;
  self->count = 1;
}

static void * PhantomCell_dtor(void * _self){
  struct PhantomCell *self = _self;
  struct PhantomCell *ring = self->ring;
  -- self->count;
  if (self->count > 0){
    return 0;
  }
  assert(ring);
  if(ring==self){
    ring = self->next;
  }
  if(ring==self){
    ring=NULL;
  }
  else{
    struct PhantomCell *b= ring;
    while(b->next!=self){
      b=b->next;
      assert(b!=ring);
    }
    b->next=self->next;
  }

  if(self->nodes){  /* deleting all the nodepointers from the phantomcell */
    struct Node node=self->nodes;
    struct Node nextnode = node->next;
    printf("WARNING:: deleting PhantomCell with nodes in it. Those nodes will be lost forever. Hope you know what you\'re doing \n");
    while (node){
      nextnode = node->next;
      delete(node);
      node=nextnode;
    }
  }
  if(self->subcells){  /* deleting all Phatom Subcells */
    struct PhantomCell * subcell = self->subcells;
    struct PhantomCell * nextsubcell = subcell->next;
    while (subcell){
      nextsubcell = subcell->next;
      delete(subcell);
      subcell = nextsubcell;
    }
  }
  return self;
}

static void * PhantomCell_diff(const void *_self, const void * _b){
  const struct PhantomCell * self = _self;
  const struct PhantomCell * b = _b;

  if(self==b) return 0;
  return 1;
}

static void * PhantomCell_clone(const void * _self){
  const struct PhantomCell * self = _self;
  ++ self->count;
  return self;
}

static const struct Class _PhantomCell{sizeof(struct PhantomCell), PhantomCell_ctor, PhantomCell_dtor, PhantomCell_clone, PhantomCell_diff};

const void * PhatomCell = & _PhantomCell;



/* ---------------Nodes---------------------------- */

static void * Node_ctor(void * _self, va_list * app){
  struct Node * self = _self;

  const PetscReal CoordsUndef[3];
  CoordsUndef[0] = va_arg(*app, PetscReal);
  CoordsUndef[1] = va_arg(*app, PetscReal);
  CoordsUndef[2] = va_arg(*app, PetscReal);

  self->findCell(self);
  /* could put a check if the node already exists. if yes it would be in the
   * same cell */
  return self;
}

static void * Node_dtor(void * _self){
  struct Node *self = _self;
  struct Node * ring = self->HomeCell->nodes;
  assert(ring);
  if(ring==self){
    ring=ring->next;
  }
  if(ring==self){
    ring=NULL;
  }
  else{
    struct Node *b = ring;
    while (b->next!=self){
      b=b->next;
      assert(b!=ring);
    }
    b->next = self->next;
  }
  return self;
}

static void * Node_diff(const void * _self, const void * _b){
  const struct Node * self = _self;
  const struct Node * b = _b;
  const PetscReal eps = 10E-9;
  PetscInt i;
  for(i=0;i<3;i++){
    if(fabs(self->CoordsUndef[i]-b->CoordsUndef[i])>eps){
      return 1;
    }
  }
  return 0;
}

static void * Node_clone(const void * _self){
  struct Node *self = _self;
  printf("WARNING:: You should not clone nodes. They are unique.")
  return self;
}

static const struct Class _Node{sizeof(struct Node), Node_ctor, Node_dtor,Node_clone, Node_diff};

const void * Node = &_Node;
  
    
  

/* PhantomCell * new_PhantomCell(PetscInt SimmerCellNr,DM ParentDM, PetscInt NIBNodes){ */
/*   PhantomCell *PhC= calloc(1,sizeof(PhantomCell)); */
/*   /\*This is supposed to allocate the space for a Phantom Cell type *\/ */
/*   /\* PhC->SimmerCellNr = SimmerCellNr; *\/ */
/*   /\* PhC->ParentDM = ParentDM; *\/ */
/*   /\* PhC->NIBNodes = NIBNodes; *\/ */
/*   /\* check if this creator works *\/ */
/*   assert(PhC); */
/*   return PhC; */
/* } */

/* void destroy_PhantomCell(PhantomCell *PhC){ */
/*   free (PhC); */
/* } */


/* /\*@ */
/*   Function to add an Element (Phantom Cell) to a set.  */

/*   @*\/ */
/* void * Set_add_PhantomCell(Set *set, PhantomCell *PhC){ */
/*   assert (set); */
/*   assert (PhC); */

/*   if(!PhC->in) PhC->in = set; */
/*   else { */
/*     assert(PhC->in == set); */
/*     set->count +=1; */
/*     PhC->count +=1; */
/*   } */

/*   return PhC; */
/* } */

/* void * findPhantomCell(const PhantomCell * Element, const void * set){ */
/*   /\* returns the element if it is in the set. 0 else *\/ */
/*   PhantomCell *PhC= Element; */
/*   assert(Element); */
/*   assert(set); */
/*   return Element->in == set ? (void *) PhC : 0 ; */

/* } */

/* int containsPhantomCell(const void * set, const PhantomCell * Element){ */
/*   return findPhantomCell(Element, set) != 0; */
/* } */

/* void * dropPhantomCell(const void * set, const PhantomCell * Element) { */
/*   Set * theset = set; */
/*   PhantomCell * PhC = findPhantomCell(Element,theset); */

/*   if(PhC){ */
/*     --PhC->count; */
/*     if (PhC->count==0){ */
/*       PhC->in = NULL; */
/*       -- theset->count; */
/*     } */
/*   } */
/* } */

/* unsigned countElements(const Set * set){ */
/*   assert(set); */
/*   return set->count; */
/* } */
  

typedef struct {
  PetscInt NPHC;  /* Number of phantom cells */
  PetscReal *PHCELL;  /* Array with the cell numbers of phantom cells*/  
}CouplingDataStruct;

extern CouplingDataStruct CouplingData;

PetscErrorCode generatePhantomCells(SimmerData Sidat, DM dm, const char interfacename[]){
  /* this routine scans through all boundary nodes, looks which SIMMER cells have a fluid solid interface (FSInt) in them and finds adequate quadrature points */
  PetscBool exists;  /* does the phantom cell exist already? */
  DMLabel label;
  PetscErrorCode ierr;
  PetscInt NInterfaceNodes;  /* number of interface nodes */
  Vec coordinates;           /* local coordinate vector */
  DM cdm;                    /* DM that prescribes coordinate layout and
                              * scattering between local/global coordintes */
  
  ierr = DMGetLabel(dm, interfacename, &label); CHKERRQ(ierr);
  if(!label) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "No boundary with the name %s exists in this DM", interfacename);


  return 0;        
}

PetscErrorCode getNextDepth(DM dm, IS ids){
  const PetscInt *Points;
  PetscInt NextDepthPoints;  /* these are not vertices. a point can be vertex, line,
                      * surf-elem or 3D-elem according to its depth */
  PetscInt NPoints;
  PetscErrorCode ierr;
  
  ierr = ISGetLocalSize(ids, &NPoints); CHKERRQ(ierr);
  ierr = ISGetIndices(ids, &Points); CHKERRQ(ierr);
  
  
  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "makeBoundaryLabel"
  /*@
    makeBoundaryLabel - from a set of points defined by a label and an id in that label create a
   * new label with prescribed name. The ids are the depthss of the points
   * (0=vertex, 1=line, 2=surface, 3=cell)
   @*/

PetscErrorCode makeBoundaryLabel(DM dm, const char OldLabelName[], const PetscInt OldLabelId,const char NewLabelName[]){


  PetscErrorCode ierr;
  IS OldIS;
  DMLabel NewLabel, DepthLabel;
  const PetscInt *Pntids;
  PetscInt NOldPnts;
  PetscInt depth;
  PetscBool exists, contains;
  PetscInt i;

  /* Make sure old label exists, new label does note yet exist and the dmlpex
   * is stratified */
  ierr = DMHasLabel(dm, OldLabelName, &exists); CHKERRQ(ierr);
  if (!exists) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"There is no label with the name %s in the specified DM \n",OldLabelName);
  ierr = DMHasLabel(dm, NewLabelName, &exists); CHKERRQ(ierr);
  if(exists) SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "A label with the name %s already exists \n",NewLabelName);
  ierr = DMHasLabel(dm, "depth", &exists); CHKERRQ(ierr);
  if(!exists) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Given DM has no depth label. Please use DMPlexStratify() before calling this function.\n");

  ierr = DMPlexGetDepth(dm, &depth); CHKERRQ(ierr);  /* check depth*/
  ierr = DMGetStratumIS(dm, OldLabelName, OldLabelId, &OldIS); CHKERRQ(ierr);
  ierr = ISGetLocalSize(OldIS, &NOldPnts); CHKERRQ(ierr);
  ierr = ISGetIndices(OldIS, &Pntids); CHKERRQ(ierr);
  ierr = DMCreateLabel(dm, NewLabelName); CHKERRQ(ierr);
  ierr = DMGetLabel(dm, NewLabelName, &NewLabel); CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "depth", &DepthLabel); CHKERRQ(ierr);
  
  for (i = 0; i < NOldPnts; ++i)
  {
    PetscInt *cone=NULL;
    PetscInt conesize=0;
    PetscInt j=0;
    ierr = DMLabelStratumHasPoint(DepthLabel,depth,Pntids[i], &contains); CHKERRQ(ierr);
    if(contains) SETERRQ2(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"The Label Stratum %s %i contains none-surface cells \n",OldLabelName,OldLabelId );

    ierr = DMLabelStratumHasPoint(DepthLabel,depth-1,Pntids[i], &contains); CHKERRQ(ierr);
    if(contains){  /* skip points with wrong depth */
      ierr = DMLabelSetValue(NewLabel, Pntids[i], depth-1); CHKERRQ(ierr);
      ierr = DMPlexGetCone(dm, Pntids[i], &cone); CHKERRQ(ierr);
      ierr = DMPlexGetConeSize(dm,Pntids[i], &conesize); CHKERRQ(ierr);

      for (j=0; j<conesize; j++)
      {
        PetscInt *VertexCone=NULL;
        PetscInt VertexConeSize=0;
        PetscInt k=0;
        
        ierr = DMLabelSetValue(NewLabel, cone[j], depth-2); CHKERRQ(ierr);
        if(depth==3){
          ierr = DMPlexGetCone(dm, cone[j], &VertexCone); CHKERRQ(ierr);
          ierr = DMPlexGetConeSize(dm,cone[j], &VertexConeSize); CHKERRQ(ierr);
          for (k=0;k<VertexConeSize;k++){
            ierr = DMLabelSetValue(NewLabel, VertexCone[k], depth-3); CHKERRQ(ierr);
          }
        }
      }
    }
  }
  ierr = ISRestoreIndices(OldIS, &Pntids); CHKERRQ(ierr);

  ierr = ISDestroy(&OldIS); CHKERRQ(ierr);
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "getVertexCoordinatsDeformed"
  /*@
    getVertexCoordinatesDeformed - Return the coordinates of a vertex in the
    deformed state
    
    CoordSection:     The Petsc Section for the coordinate DM
    Def Section:      Petsc subsection of the Solution Section for the
                      deformation /displacement field. Obtain with
                      PetscSectionGetField(). Solution section might have more
                      DOF than just the spatial directions.
    coords:           Vector holding the dmplex coordinates
    solution:         Vector holding the solution of the FE analysis
                      (displacements, velocities and/or more)
    VertexNo:         The vertex (in the global numbering scheme) for which
                      you want the deformed coordinates.
    VertexCoords:     Array holding the coordinates of the specified vertex
                      in the deformed state. Should be preallocated with a
                      size reflecting the spatial dimension
     @*/

PetscErrorCode getVertexCoordinatsDeformed(PetscSection CoordSect, PetscSection DispSection, Vec coords, Vec solution, PetscInt VertexNo, PetscReal VertexCoords[]){

  PetscInt i,dim,DispDof, CoordOffset, DispOffset;
  PetscInt *CoordIdx=NULL, *DispIdx=NULL;
  PetscReal *displacement=NULL;
  PetscErrorCode ierr;

  ierr = PetscSectionGetOffset(CoordSect,VertexNo,&CoordOffset); CHKERRQ(ierr);
  ierr = PetscSectionGetDof(CoordSect, VertexNo, &dim); CHKERRQ(ierr);
  ierr = PetscSectionGetDof(DispSection, VertexNo, &DispDof); CHKERRQ(ierr);
  ierr = PetscSectionGetOffset(DispSection, VertexNo, &DispOffset); CHKERRQ(ierr);

  /* due to unfortunate placement of dirichlet boundary conditions, the
   * deformation might lack a degree of freedom. We currently handle this
   * with an error since we don't know which spatial dimension is
   * restricted. This issue should be adressed in the future.*/
  if (!dim==DispDof) SETERRQ3(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Vertex No %i has only %i degrees of freedom in a %i dimensional domain.",VertexNo,DispDof,dim);

  ierr = PetscCalloc1(DispDof, &displacement); CHKERRQ(ierr);
  ierr = PetscCalloc1(DispDof, &DispIdx); CHKERRQ(ierr);
  ierr = PetscCalloc1(dim,&CoordIdx); CHKERRQ(ierr);

  for (i=0; i<dim; i++){
    DispIdx[i] = i+DispOffset;
    CoordIdx[i] = i+CoordOffset;
  }
  ierr = VecGetValues(coords, dim, CoordIdx, VertexCoords); CHKERRQ(ierr);
  ierr = VecGetValues(solution, DispDof, DispIdx, displacement); CHKERRQ(ierr);
  
  for (i=0; i<dim; i++){
    VertexCoords[i] = VertexCoords[i] + displacement[i];
  }

  ierr = PetscFree(displacement); CHKERRQ(ierr);
  ierr = PetscFree(DispIdx); CHKERRQ(ierr);
  ierr = PetscFree(CoordIdx); CHKERRQ(ierr);
  
  return 0;
  
   
}


PetscErrorCode getFaceCoordinatesDeformed(DM dm, const char LabelName[], Vec *FaceCoords){
  /* Takes a DM and the name of a label. The label IDs should represent the
   * depth of the respective strata i.e. Label[0] = Vertices, Label[1] lines etc */
  DM cdm;
  PetscSection sect, ClosSect;  /* section defining the coordinates and the
                                 * closure size of each element in the label */
  PetscInt NSurfElements,dim;
  const PetscInt *SurfElement;
  Vec coordinates, coordinatessubset;
  PetscScalar *coords=NULL;
  IS SurfVertices;
  PetscErrorCode ierr;
  DMLabel FaceLabel;
  PetscInt i;
  PetscInt *ClsrPoints= NULL;
  PetscInt NClsrPoints=0;
  PetscBool exists=PETSC_FALSE;
  PetscInt pstart, pend;

  /* check whether a label with the given name exists in the DM */
  ierr = DMHasLabel(dm, LabelName, &exists); CHKERRQ(ierr);
  if(!exists) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"The specified DM does not contain a label with the name %s.\n",LabelName);
  
  
  /* The coordinates are related to the vertices of the mesh. We hence need
   * the node numbers of the vertices from the label that the function
   * receives*/
  
  ierr = DMGetDimension(dm,&dim); CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(dm,&coordinates); CHKERRQ(ierr);
  ierr = DMGetCoordinateDM(dm, &cdm); CHKERRQ(ierr);
  ierr = DMGetDefaultSection(cdm, &sect); CHKERRQ(ierr);


  /* ierr = DMGetLabel(dm, LabelName, &FaceLabel);CHKERRQ(ierr); */
  /* ierr = DMLabelView(FaceLabel, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr); */
  /* ierr = DMPlexLabelComplete(dm, FaceLabel); CHKERRQ(ierr); */
  /* ierr = DMLabelView(FaceLabel, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr); */

  /* ierr = DMPlexCreateClosureIndex(dm,sect); CHKERRQ(ierr); */
  /* ierr = PetscSectionGetClosureIndex(sect, (PetscObject)FaceLabel, &ClosSect, &ClPnts); CHKERRQ(ierr); */
  /* ierr = ISView(ClPnts, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr); */
  
  
  ierr = PetscSectionGetChart(sect,&pstart,&pend); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "pstart %i     pend %i", pstart, pend); CHKERRQ(ierr);
  ierr = DMGetStratumIS(dm, LabelName, 0, &SurfVertices); CHKERRQ(ierr);
  ierr = VecGetSubVector(coordinates, SurfVertices, &coordinatessubset); CHKERRQ(ierr);
  ierr = VecDuplicate(coordinatessubset, FaceCoords); CHKERRQ(ierr);
  ierr = VecCopy(coordinatessubset, *FaceCoords); CHKERRQ(ierr);
  ierr = VecRestoreSubVector(coordinates, SurfVertices, &coordinatessubset); CHKERRQ(ierr);

  ierr = ISDestroy(&SurfVertices); CHKERRQ(ierr);


  /* ierr = DMPlexLabelComplete(dm,LabelName); CHKERRQ(ierr); */
  /* ierr = DMGetStratumIS(dm, LabelName,FaceId,&ids); CHKERRQ(ierr); */
  /* ierr = ISGetLocalSize(ids,&NSurfElements); CHKERRQ(ierr); */
  /* ierr = ISGetIndices(ids,&SurfElement); CHKERRQ(ierr); */
  /* ierr = VecCreateMPI(PETSC_COMM_WORLD,NSurfElements*dim,PETSC_DECIDE, FaceCoords); CHKERRQ(ierr); */

  /* PetscPrintf(PETSC_COMM_WORLD, "%i\n", NSurfElements); */
  
  /* for ( i = 0; i < NSurfElements; ++i) */
  /* { */
  /*   PetscInt ClsrSize=0; */
  /*   PetscInt j=0; */
    
  /*   ierr = DMPlexVecGetClosure(dm,sect,coordinates, SurfElement[i], &ClsrSize, &coords); CHKERRQ(ierr); */
  /*   /\* ierr = DMPlexGetTransitiveClosure(dm, SurfElement[i], PETSC_FALSE, &NClsrPoints, &ClsrPoints); CHKERRQ(ierr); *\/ */
  /*   PetscPrintf(PETSC_COMM_WORLD, "Point %i has a closure of size %i\n", SurfElement[i],ClsrSize); */
  /*   /\* ierr = DMPlexRestoreTransitiveClosure(dm, SurfElement[i], PETSC_FALSE, &NClsrPoints, &ClsrPoints); CHKERRQ(ierr); *\/ */
  /*   for (j=0; j < ClsrSize; ++j) */
  /*   { */
  /*     PetscPrintf(PETSC_COMM_WORLD, "%f\n", coords[j]); */
  /*   } */
  /*   ierr = DMPlexVecRestoreClosure(dm,sect,coordinates, SurfElement[i], &ClsrSize, &coords); CHKERRQ(ierr); */
  /* } */


  /* ierr = VecGetArrayRead(coordinates,&coords); CHKERRQ(ierr); */

  /* ierr = VecRestoreArrayRead(coordinates, &coords);CHKERRQ(ierr); */
  /* ierr = ISView(ids, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr); */
  /* ierr = VecView(coordinates,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr); */


  /* ierr = VecGetSubVector(coordinates, ids,&coordinatessubset); CHKERRQ(ierr); */
  /* ierr = VecDuplicate(coordinatessubset,FaceCoords); CHKERRQ(ierr); */
  /* ierr = VecCopy(coordinatessubset,*FaceCoords); CHKERRQ(ierr); */
  /* ierr = VecRestoreSubVector(coordinates,ids,&coordinatessubset); CHKERRQ(ierr); */

  return 0;
  
}
