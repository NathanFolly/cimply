#include "phantommesh.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


static void * distributevertices(void * _self);
/* distributes the boundary vertices and assigns them to the adequate phantom cells */
static void * PhantomMesh_getPhantomFractions(void * _self, float ** phantomfractions);
/* returns an allocated array with the phantomfractions of each cell. 0 if the cell is not a phantom cell */
static void * findPhantomCell(const void * _self, float * position);
  /*takes phantommesh and 3-D position in cylindrical coordinates
  returns pointer to phantomcell in phantommesh at specified position */
static void * PhantomMesh_generatetestsphere(void * _self, float radius, int nvertices);
/* when running a test, generates a random set of boundary vertices on a sphere with the indicated radius. center of the sphere is 0/0/0 */
static void * handlePhantomCellsWaterHammerlike(void * _self);
/* updates the phantom cells and creates new phantom cells on the upper right side of the fluid solid interface. works only for waterhammer like problems */





static void * PhantomMesh_ctor(void * _self, va_list * app){
  struct PhantomMesh * self = _self;
  self->getPhantomFractions = PhantomMesh_getPhantomFractions;
  self->generateTestSphere = PhantomMesh_generatetestsphere;
  return self;
}

static void * PhantomMesh_dtor(void * _self){
  struct PhantomMesh * self = _self;
  int i;
  
  if(self->phantomcell){
    for (i=0; i<self->MMS; i++){
      if (self->phantomcell[i]){
        delete(self->phantomcell[i]);
      }
    }
    free(self->phantomcell);
  }
  if(self->boundary){
    for(i=0; i<self->nbofvertices; i++){
      if (self->boundary[i]){
        delete(self->boundary[i]);
      }
    }
    free(self->boundary);
  }
  return self;
}

static void * PhantomMesh_update(void * _self){
  struct PhantomMesh * self = _self;
  int i;
    /* TODO: write update routine for the Phantom Mesh */

  for (i =0; i< self->MMS; i++){
    if(self->phantomcell[i]){
      delete(self->phantomcell[i]);
    }
  }
  for (i = 0; i< self->nbofvertices; i++){
    update(self->boundary[i]);
  }
  
  distributevertices(self);  /* reassign vertices to phantom cells according to their updated position */

  
  /* TODO: the following code only works for waterhammer-like problems (if the upper right part of the mesh is the phantom part)*/

  handlePhantomCellsWaterHammerlike(self);
  return 0;
}

static const struct Class _PhantomMesh={sizeof(struct PhantomMesh), PhantomMesh_ctor, PhantomMesh_dtor, PhantomMesh_update};
const void * PhantomMesh = &_PhantomMesh;



/* -------------- Class specific functions -------------- */


static void * PhantomMesh_getPhantomFractions(void * _self, float ** phantomfractions){
  const struct PhantomMesh * self = _self;
  int i;
  free(*phantomfractions);
  (*phantomfractions) = (float *) calloc(self->MMS,sizeof(float));
   for (i = 0; i<self->MMS; i++){
    if(self->phantomcell[i]){
      (*phantomfractions)[i] = phantomfraction(self->phantomcell[i]);
    }
  }
  return 0;
}
    


static void * distributevertices(void * _self){
  struct PhantomMesh * self = _self;
  int i;
  for (i = 0; i<self->nbofvertices; i++){
    float * position;
    struct PhantomCell * mothercell;
    getposition(self->boundary[i], &position, "cylind");  /* get
                                                               * coordinates
                                                               * of current
                                                               * position in
                                                               * cartesian /
                                                               * cylindrical coordinates */
    mothercell = findPhantomCell(self, position);
    assign(mothercell, self->boundary[i]);
  }
  return 0;
}


static void * findPhantomCell(const void * _self, float * position){
  /*takes phantommesh and 3-D position in cylindrical coordinates
  returns pointer to phantomcell in phantommesh at specified position */
  const struct PhantomMesh * self = _self;
  struct PhantomCell * phantomcell;
  /* TODO : this is only for SIMMER 2D at the moment */
  PetscBool routside=PETSC_FALSE, zoutside=PETSC_FALSE;  /* do the positions fall outside the SIMMER
                                  * grid */
  PetscBool foundI=PETSC_FALSE, foundJ=PETSC_FALSE, isoutside=PETSC_FALSE;
  PetscReal cummulR = 0, cummulZ=0;  /* total R and Z position of the IJ outer
                                     * cell face */
  int CellNr, Ipos=0, Jpos=0;
  const PetscReal r = position[0];
  
  while (!foundI){
    if (Ipos>=self->IB){
      routside=PETSC_TRUE;
      foundI = PETSC_TRUE;
      Ipos=Ipos-1;
      break;
    }
    cummulR = cummulR+self->DRINP[Ipos];
    if (cummulR>r){
      foundI=PETSC_TRUE;
      break;
    }
    Ipos=Ipos+1;
  }
  if (position[2]<0)
  {
    foundJ=PETSC_TRUE;
    zoutside=PETSC_TRUE;
  }
  while (!foundJ){
    if (Jpos>=self->JB){
      zoutside=PETSC_TRUE;
      foundJ = PETSC_TRUE;
      Jpos=Jpos-1;
      break;
    }
    cummulZ = cummulZ + self->DZINP[Jpos];
    if (cummulZ>position[2]){
      foundJ = PETSC_TRUE;
      break;
    }
    Jpos=Jpos+1;
  }
  if(routside || zoutside){
    isoutside = PETSC_TRUE;
   }
  CellNr = (self->IB+2)*(Jpos+1)+Ipos+1;  

  if(!self->phantomcell[CellNr]){
    float RLB=0,RUB=0,ZLB=0,ZUB=0;
    int i;
    for (i=0;i<Ipos;i++){
      RLB+=self->DRINP[i];
    }
    RUB=RLB+self->DRINP[Ipos];
    for (i=0;i<Jpos;i++){
      ZLB+=self->DZINP[i];
    }
    ZUB=ZLB+self->DZINP[Jpos];
    self->phantomcell[CellNr]=new(PhantomCell,RLB,RUB,ZLB,ZUB);
  }
  phantomcell = self->phantomcell[CellNr];
  return phantomcell;
}


static void * PhantomMesh_assignFEM(void * _self, void * _fem){
  struct PhantomMesh * self = _self;
  struct FEMAnalysis * fem = _fem;

  FEMGetDMPointer(fem,self->dmptr);
  FEMGetSolutionPointer(fem,self->sltnptr);
  FEMGetImmersedBoundaryName(fem,self->immersedboundaryname);

  return 0;
}

static void * PhantomMesh_generatetestsphere(void * _self, float radius, int nvertices){
  struct PhantomMesh * self = _self;
  int i;

  self->nbofvertices = nvertices;
  self->boundary = (void **) calloc(nvertices,sizeof(void *));
  for (i=0;i<nvertices;i++){
    float x,y,z;
    x = ((float)rand()/(float)(RAND_MAX)) * radius;  /* generate random number between 0 and radius */
    y = ((float)rand()/(float)(RAND_MAX)) * sqrt(pow(radius,2)-pow(x,2));
    z = sqrt(pow(radius,2)-pow(x,2)-pow(y,2));
    self->boundary[i] = new(BoundaryVertex, x,y,z);
    /* printf("x = %f        y= %f,        z=%f \n",x,y,z); */
  }
  return 0;
};


/* TODO:: this ugly part works only as long as the problem looks like a water hammer with the upper right side of the SIMMER mesh being the phantom region */
static void * handlePhantomCellsWaterHammerlike(void * _self){
  struct PhantomMesh * self = _self;
  int i;
    int fsicellnb=-1;  /* cell number of fluid/solid interface cell last encountered. determining  the 100% phantomfraction cells */
  for (i =0; i< self->MMS; i++){
    if(self->phantomcell[i]){
      update(self->phantomcell[i]);
      if(fsicellnb!=0) fsicellnb=i;
      if(i%(self->IB+2)==1) fsicellnb=0;
    }
   else if((i+1)%(self->IB+2)==0&&fsicellnb!=0){  /* check if new line and if in past line, first cell was not a boundary cell */
      fsicellnb=-1;
    }
   else if(0<=fsicellnb&&i>fsicellnb){
      int Ipos=i%(self->IB+2);
      int Jpos=(i-Ipos-1)/(self->IB+2);
      float RLB=0,RUB=0,ZLB=0,ZUB=0;
      int j;
      if(Jpos<self->JB){
        for (j=0;j<Jpos;j++){
          RLB+=self->DRINP[j];
        }
        RUB=RLB+self->DRINP[Jpos];
        for (j=0;j<Jpos;j++){
          ZLB+=self->DZINP[j];
        }
        ZUB=ZLB+self->DZINP[Jpos];
        self->phantomcell[i]=new(PhantomCell,RLB,RUB,ZLB,ZUB);
        update(self->phantomcell[i]);
      }
    }
  }
  return 0; 
}
