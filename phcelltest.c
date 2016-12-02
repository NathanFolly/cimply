#include <petscsys.h>
#include <petscdmplex.h>


int main(int argc, char **argv){
  DM dm, distributeddm;
  PetscErrorCode ierr;
  Vec u, coords;
  DMLabel interface;
  IS is;
  PetscInt i, size, dim;
  const PetscInt *indices=NULL;
  PetscReal *VertexCoords=NULL;
  PetscSection CoordSect;
  

  ierr = PetscInitialize(&argc, &argv, NULL, NULL); CHKERRQ(ierr);

  ierr = DMPlexCreateFromFile(PETSC_COMM_WORLD, "testmesh.msh", PETSC_TRUE,&dm); CHKERRQ(ierr);
  ierr = DMPlexDistribute(dm,0,NULL,&distributeddm); CHKERRQ(ierr);

  if (distributeddm) {
    ierr=DMDestroy(&dm);CHKERRQ(ierr);
    dm = distributeddm;
  }

  ierr = DMGetDimension(dm, &dim); CHKERRQ(ierr);
  ierr = PetscMalloc1(dim, &VertexCoords); CHKERRQ(ierr);

  ierr = makeBoundaryLabel(dm, "Face Sets", 1, "Interface"); CHKERRQ(ierr);
  
  ierr = DMGetLabel(dm, "Interface", &interface); CHKERRQ(ierr);
  /* ierr = DMPlexLabelComplete(dm, interface); CHKERRQ(ierr); */
  ierr = DMView(dm, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  ierr = DMGetCoordinateSection(dm, &CoordSect); CHKERRQ(ierr);
  ierr = DMSetDefaultSection(dm, CoordSect); CHKERRQ(ierr);
  ierr = DMGetCoordinates(dm, &coords); CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(dm, &u); CHKERRQ(ierr);
  ierr = VecSet(u, 0.5); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(u); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(u); CHKERRQ(ierr);
  
  ierr = DMGetStratumIS(dm, "Interface", 0, &is); CHKERRQ(ierr);
  ierr = ISGetIndices(is, &indices); CHKERRQ(ierr);
  ierr = ISGetSize(is, &size); CHKERRQ(ierr);

  for (i=0; i<size; i++){
    PetscInt j;
    ierr = getVertexCoordinatsDeformed(CoordSect, CoordSect, coords, u, indices[i], VertexCoords); CHKERRQ(ierr);
    for(j=0;j<dim;j++){
      PetscPrintf(PETSC_COMM_WORLD,"Vertex No: %i       dim: %i     value. %f\n",indices[i],j,VertexCoords[j]);
    }
  }
  
  ierr = ISRestoreIndices(is, &indices); CHKERRQ(ierr);
  
  
  /* ierr = getFaceCoordinatesDeformed(dm, "Interface",&coords); CHKERRQ(ierr); */
  
  /* ierr = VecView(coords,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr); */

  VecDestroy(&coords);
  DMDestroy(&dm);
  PetscFinalize();
  
  return 0;
}
