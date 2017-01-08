#include <petscdm.h>
#include <petscdmlabel.h>
#include <petscdmplex.h>


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
      ierr = DMPlexGetCone(dm, Pntids[i],(const PetscInt **)  &cone); CHKERRQ(ierr);
      ierr = DMPlexGetConeSize(dm,Pntids[i], &conesize); CHKERRQ(ierr);

      for (j=0; j<conesize; j++)
      {
        PetscInt *VertexCone=NULL;
        PetscInt VertexConeSize=0;
        PetscInt k=0;
        
        ierr = DMLabelSetValue(NewLabel, cone[j], depth-2); CHKERRQ(ierr);
        if(depth==3){
          ierr = DMPlexGetCone(dm, cone[j], (const PetscInt **) &VertexCone); CHKERRQ(ierr);
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
  if (!(dim==DispDof)) SETERRQ3(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Vertex No %i has only %i degrees of freedom in a %i dimensional domain.",VertexNo,DispDof,dim);

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
  PetscSection sect;  /* section defining the coordinates and the
                                 * closure size of each element in the label */
  PetscInt dim;
  /* const PetscInt *SurfElement; */
  Vec coordinates, coordinatessubset;
  /* PetscScalar *coords=NULL; */
  IS SurfVertices;
  PetscErrorCode ierr;
  /* DMLabel FaceLabel; */
  /* PetscInt i; */
  /* PetscInt *ClsrPoints= NULL; */
  /* PetscInt NClsrPoints=0; */
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

/* #undef __FUNCT__ */
/* #define __FUNCT__ "getVertexSurfaceNormal" */
/* /\* given a dm, a surface label and a vertex number, gives the face normal vector of the surface at that point in space*\/ */

/* PetscErrorCode getVertexSurfaceNormal(DM dm, const char * labelname, const int vertexnumber, double * facenormal) */
/* { */
/*   PetscErrorCode ierr; */
  
