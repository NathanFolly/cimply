#include "node.h"


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
  
