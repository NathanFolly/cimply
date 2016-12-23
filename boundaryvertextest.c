#include "cimplyobjects.h"
#include "boundaryvertex.h"
#include "phantomcell.h"
#include <math.h>



 int main(int argc, char *argv[])
{
  void * phantomcell;
  float * position;
  const int nvertices = 50;
  void * boundaryvertex[nvertices];
  int i;
  double RLB=0,RUB=1,ZLB=0,ZUB=1;
  double lafraction = 0;


  phantomcell = new(PhantomCell,RLB,RUB,ZLB,ZUB);

  for (i=0;i<nvertices;i++){
    float x,y,z;
    x = ((float)rand()/(float)(RAND_MAX)) * 0.5;  /* generate random number between 0 and 0.5 */
    y = ((float)rand()/(float)(RAND_MAX)) * sqrt(pow(0.5,2)-pow(x,2));
    z = sqrt(pow(0.5,2)-pow(x,2)-pow(y,2));
    boundaryvertex[i] = new(BoundaryVertex, x,y,z);
    /* printf("x = %f        y= %f,        z=%f \n",x,y,z); */
    assign(phantomcell,boundaryvertex[i]);
  }
  getposition(boundaryvertex[0], &position, "cylind");
  for ( i = 0; i<3; i++){
    printf("position [%i] = %f \n",i,position[i]);
  }
  update(phantomcell);
  lafraction = phantomfraction(phantomcell);
  printf("phantom fraction is %f \n",lafraction);
  for (i=0;i<nvertices;i++){
    delete(boundaryvertex[i]);
  }
  delete(phantomcell);
  return 0;
}
