#include <assert.h>
#include <stdio.h>
#include <petscdm.h>
#include <petscdmlabel.h>
/* #include <petscdmplex.h> */
#include "phantomcell.h"
#include "cimplySimmerUtils.h"
#include <math.h>

static void PhantomCell_calculatePhantomFraction_rectangleslabs(void * _self);
static void PhantomCell_addvertex(void * _self, void * _boundaryvertex);
static void PhantomCell_prepare(void * _self);
static void PhantomCell_givePhantomFraction(void * _self, float * phantomfraction);


static void * PhantomCell_ctor(void * _self, va_list *app ){
  struct PhantomCell * self = _self;
  /*
   call: new(PhantomCell, PetscReal RLB,
   PetscReal RUB, PetscReal ZLB, PetscReal ZUB)
   
   RLB: Radial lower boundary of the cell
   RUB: radial upper boundary of the cell
   ZLB: Z-direction lower boundary of the cell
   ZUB: Z-direction upper boundary of the cell
  */
  self->RLB = va_arg(*app, double);
  self->RUB = va_arg(*app, double);
  self->ZLB = va_arg(*app, double);
  self->ZUB = va_arg(*app, double);

  self->cellvolume=M_PI*(pow(self->RUB,2)-pow(self->RLB,2))*(self->ZUB-self->ZLB);
      
  self->prepare=PhantomCell_prepare;
  self->calculatePhantomFraction = PhantomCell_calculatePhantomFraction_rectangleslabs;
  self->assign=PhantomCell_addvertex;
  self->givePhantomFraction = PhantomCell_givePhantomFraction;
  
  return self;
}

static void * PhantomCell_dtor(void * _self){
  struct PhantomCell *self = _self;

  
  if(self->vertexring){  /* deleting all the nodepointers from the phantomcell */
    /* actually not deleting the vertices but the pointers to them */
    struct Node * vertex = self->vertexring;
    struct Node * nextvertex = vertex->next;
    while (vertex){
      nextvertex = vertex->next;
      free(vertex);
      vertex=nextvertex;
    }
  }
  if(self->subcell){  /* deleting all Phatom Subcells */
    struct PhantomCell ** subcell = self->subcell;
    int i;
    for (i = 0; i<sizeof(self->subcell)/sizeof(self->subcell[0]); i++) {
      delete (self->subcell[i]);
    }
  }
  return self;
}

static void * PhantomCell_update(void * _self){
  struct PhantomCell * self = _self;
  self->prepare(self);  /* find penetration surfaces, subdivide cell if necessary */
  self->calculatePhantomFraction(self);
  return 0;
}
static const struct Class _PhantomCell = {sizeof(struct PhantomCell), PhantomCell_ctor, PhantomCell_dtor, PhantomCell_update};

const void * PhantomCell = & _PhantomCell;


static void PhantomCell_givePhantomFraction(void * _self, float *  phantomfraction){
  struct PhantomCell * self = _self;
  /* TODO: maybe add uptodate attribute in PhantomCell and update here if not already up to date */
  * phantomfraction = self->phantomfraction;
  return;
}


static void PhantomCell_calculatePhantomFraction_rectangleslabs(void * _self){
  struct PhantomCell * self = _self;
  struct BoundaryVertex * boundaryvertex;
  struct Node * node=self->vertexring;
  double * position;
  double rmax, zmax;  /* the maximum radial and vertical vertex positions encountered so far */
  double reductionvolume=0;  /* volume by which the phantomvolume is reduced*/
  int i=0;

    if(!self->vertexring){
    /* fprintf(stderr,"ERROR:: error calculating the phantom fractions. Cell has no vertices appointed to it \n"); */
    reductionvolume=0.0;
  }

  while(node){
    /* TODO: this is ugly. there must be a better solution to this */
    double rclosest=self->RLB, zclosest=self->ZLB;  /* r- and z- positions of closest subsquares */
    boundaryvertex =(struct BoundaryVertex *) node->content;
    getposition(boundaryvertex,&position,"cylind");
     struct Node * historynode=self->vertexring;
    while(historynode!=node){
      double * historyposition;
      double rdist, zdist;
      getposition(historynode->content, &historyposition, "cylind");
      rdist = position[0]-historyposition[0];
      zdist = position[2]-historyposition[2];
      if(rdist>0&&rdist<(position[0]-rclosest)){
        rclosest = historyposition[0];
      }
      if(zdist>0&&zdist<(position[2]-zclosest)){
        zclosest = historyposition[2];
      }
      historynode=historynode->next;
    }
    reductionvolume += M_PI*(pow(position[0],2)-pow(rclosest,2))*(position[2]-zclosest);
    
    node = node->next;
    i++;
  }
  self->phantomvolume = self->cellvolume - reductionvolume;
  self->phantomfraction = self->phantomvolume/self->cellvolume;
  /* printf("phantomvolume = %f \n", self->phantomvolume); */
  /* printf("cellvolume = %f \n", self->cellvolume); */
  /* printf("PHANTOMFRACTION = %f \n", self->phantomfraction); */
  

  
  
  /* run  */
}




  

static void PhantomCell_prepare(void * _self){
  struct PhantomCell * self = _self;

  return;
}


static void PhantomCell_addvertex(void * _self, void * _boundaryvertex){
  struct PhantomCell * self = _self;
  struct BoundaryVertex * boundaryvertex = _boundaryvertex;
  struct Node * node;
        
  if(!self->vertexring){
    self->vertexring = (struct Node *) calloc(1,sizeof(struct Node));
    self->vertexring->content=boundaryvertex;
    self->vertexcount++;
    return;
  }
  node = self->vertexring;
  while(node->next){
    node = node->next;
  }
  node->next=(struct Node * ) calloc(1,sizeof(struct Node));
  node->next->content=boundaryvertex;
  self->vertexcount++;
  return;
}



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
