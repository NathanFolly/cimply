#include "phantommesh.h"

static void * distributevertices(void * _self);



static void * PhantomMesh_ctor(void * _self, va_list * app){
  struct PhantomMesh * self = _self;

  return self;
}

static void * PhantomMesh_dtor(void * _self){
  struct PhantomMesh * self = _self;
  
  if(self->phantomcell){
    for (int i=0; i<self->MMS; i++){
      if (self->phantomcell[i]){
        delete(self->phantomcell[i]);
      }
    }
    free self->phantomcell;
  }
  if(self->boundary){
    for(int i=0; i<self->nbofvertices; i++){
      if (self->boundary[i]){
        delete(self->boundary[i]);
      }
    }
    free self->boundary;
  }
  return self;
}

static void * PhantomMesh_update(void * _self){
  struct PhantomMesh * self = _self;
  /* TODO: write update routine for the Phantom Mesh */
  for (int i = 0; i< self->nbofvertices; i++){
    update(self->boundary[i]);
  }
  
  distributevertices(self);  /* reassign vertices to phantom cells according to their updated position */

  for (int i =0; i< self->MMS; i++){
    update(self->phantomcell[i]);
  }
  

  return 0;
}

static const struct Class _PhantomMesh={sizeof(struct PhantomMesh), PhantomMesh_ctor, PhantomMesh_dtor, PhantomMesh_update};
const void * PhantomMesh = &_PhantomMesh;



/* -------------- Class specific functions -------------- */
static void * distributevertices(void * _self){
  struct PhantomMesh * self = _self;

  for (int i = 0; i<self->nbofvertices; i++){
    position = getposition(self->boundary[i], cart / cylind);  /* get
                                                               * coordinates
                                                               * of current
                                                               * position in
                                                               * cartesian /
                                                               * polar coordinates */
    mothercell = findPhatomCell(position);
    appoint(mothercell, self->boundary[i]);
  }
  return 0;
}
