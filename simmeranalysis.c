#include "simmeranalysis.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>

static void * readfromsim05file(void * _self);
static void * SimmerAnalysis_givepressurefromlocation(const void * _self, const PetscReal x[], PetscReal * pressure);
static void * SimmerAnalysis_setuniformpressure(void * _self, const double pressure);  /* for test reasons  */
static void * SimmerAnalysis_setpressure(void * _self, const double * pressure);
static void * SimmerAnalysis_getMMS(const void * self, int *  MMS);


static void * SimmerAnalysis_ctor(void * _self, va_list * app){
  struct SimmerAnalysis * self = _self;
  /* intializing the variables */
  self->PK = NULL;
  self->givepressurefromlocation=SimmerAnalysis_givepressurefromlocation;
  self->setuniformpressure= SimmerAnalysis_setuniformpressure;
  self->setpressure = SimmerAnalysis_setpressure;
  self->getMMS = SimmerAnalysis_getMMS;
  
  /* reading the sim05 file */
  readfromsim05file(self);
  /* allocating the correct size of the variable vectors */
  self->PK = (double *) calloc(self->MMS, sizeof(double));
  return self;
}

static void * SimmerAnalysis_dtor(void * _self){
  struct SimmerAnalysis * self = _self;
  if(self->DRINP){
    free(self->DRINP);
  }
  if(self->DZINP){
    free(self->DZINP);
  }
  if(self->PK){
    free(self->PK);
  }
  return self;
}

static void * SimmerAnalysis_update(void * _self){
  struct SimmerAnalysis * self = _self;
  /* TODO: what shoudl this do? */
  return 0;
}

static const struct Class _SimmerAnalysis = {sizeof(struct SimmerAnalysis),SimmerAnalysis_ctor,SimmerAnalysis_dtor, SimmerAnalysis_update};

const void * SimmerAnalysis = &_SimmerAnalysis;










/*-------------- reader for the sim05 file ----------- */

void * getSimmerPressure(const void * _self, const PetscReal x[], PetscReal * pressure)
{
  const struct Class * const * cp = _self;
  if(*cp != SimmerAnalysis)
  {
    fprintf(stderr,"ERROR:: error in getSimmerPressure. First argument needs to be of type SimmerAnalysis. \n");
    return 0;
  }
  
  const struct SimmerAnalysis * self = _self;
  self->givepressurefromlocation(self, x, pressure);
  return 0;
}

void * setUniformSimmerPressure(void * _self, const double pressure){
  struct Class ** cp = _self;
  if(*cp!=SimmerAnalysis)
  {
    fprintf(stderr,"ERROR:: error in setUniformSimmerPressure. First argument needs to be of type SimmerAnalysis. \n");
    return 0;
  }
  
  struct SimmerAnalysis * self = _self;
  self->setuniformpressure(self, pressure);
  return 0;
}


static void * SimmerAnalysis_givepressurefromlocation(const void * _self,
                                                      const PetscReal x[], PetscReal * pressure)
{
  const struct SimmerAnalysis * self = _self;
  PetscBool isoutside;  /* don't do anything with this information atm */
  PetscInt Ipos=0, Jpos=0;  /* don't use this either atm */
  PetscInt CellNr=0;
  findSIMMERCell(self, x, isoutside, &Ipos, &Jpos, &CellNr);
  *pressure=self->PK[CellNr];
  return 0;
  
}

static void * SimmerAnalysis_setuniformpressure(void * _self, const double pressure)
{
  struct SimmerAnalysis * self = _self;
  int i=0;

  for(i=0;i<self->MMS;i++)
  {
    self->PK[i] = pressure;
  }
  return 0;
}

static void * SimmerAnalysis_setpressure(void * _self, const double * pressure)
{
  struct SimmerAnalysis * self = _self;
  int i;
  /* assert(sizeof(pressure)/sizeof(double)==self->MMS); */
  for (i=0;i<self->MMS;i++)
    {
      self->PK[i]=pressure[i];
    }
  return 0;
}
    
  
static void * SimmerAnalysis_getMMS(const void * self, int *  MMS)
{
  const struct SimmerAnalysis * self = _self;
  *MMS = self->MMS;
  return 0;
}


static void * readfromsim05file(void * _self){

  struct SimmerAnalysis * self = _self;

  
  /* basically copied from sim05reader */
  /* TODO: make this more robust. Likely to break if format of input varies. */
  /* TODO: introduce simmerdata class. put TWFIN and PK etc there. simmerdata should have a simmermesh class object as attribute */
  FILE *sim05;
  char rstrng[20], valstrng[20];  /* rstrng is an actual string used to
                                   * identify which variables is defined in
                                   * the sim05 file
                                   valstrng is actually a value. sometimes it
                                   has the fortran style double notation
                                   (2.50000D-2) which c cannot read so we
                                   have to manipulate. Sometimes there is a multiplier*/
  int rint;         /* The integer we read in from the file */
  PetscReal rreal;  /* The float we read in from the file */
  int startat, NSSize=1;      /* start point for the grid size definition,
                                * number of cells with the same size */
  int IBFilled=0, JBFilled=0;  /* number of cells that we have registered the
                            * heigths for */
  double hcell;     /* cell height */
  int xmshread = 0;  /* 0= Tag &XMSH not encountered yet in sim05 file
                        1= in XMSH region; next job: read number of cells in
                        each dir
                        2= in XMSH region; next job: read cell sizes
                        3= completed read
                        4= encountered &END after &XMSH before completing
                        data read*/
  int xtmeread = 0;  /* Same principle as for xmshread. Just for the time
                      * context */
  PetscInt IB;
  PetscInt JB;
  PetscInt IBP2;
  PetscInt JBP2;
  PetscInt MMS;
  PetscReal * DRINP;
  PetscReal * DZINP;
  PetscReal TWFIN;
  
  PetscErrorCode ierr;
  sim05 = fopen("sim05","r");
  if (sim05){
    while(xmshread==0)
    {
      fscanf(sim05,"%s",rstrng);
      if(strcmp(rstrng,"&XMSH")==0){
        xmshread=1;
      }
    }
    while(xmshread==1)
    {
      fscanf(sim05,"%2s %*[=] %d %*[,]",rstrng, &rint);
      if(strcmp(rstrng,"IB")==0) self->IB=rint;
      if(strcmp(rstrng,"JB")==0) self->JB=rint;
      if(self->IB*self->JB!=0) xmshread=2;  /* once both variables
                                                       * are assigned values:continue */
    }
    self->DRINP = (double *) calloc(self->IB, sizeof(double));
    self->DZINP = (double *) calloc(self->JB, sizeof(double));
    /* strcp(frmt, "%5c %u") */
    while(xmshread==2)
    {
      int i = 0;
      fscanf(sim05,"%5s %*[(] %d %*[)] %*[=] %19s",rstrng,&startat,valstrng);
      NSSize=1;
      while (valstrng[i]!='\0') 
      {
        if (valstrng[i]=='D') valstrng[i]='E';  /* replacing D with E so that C can read the
                                                 * float  */
        if (valstrng[i]=='*') NSSize=0;  /* we have more than one
                                          * cell with the same height */
        i++;
      }
      if (NSSize==0) sscanf(valstrng, "%d %*[*] %lG",&NSSize, &hcell);
      else sscanf(valstrng, "%lG", &hcell);
      for (i=startat-1;i<startat+NSSize-1;i++)
      {
        if(strcmp(rstrng,"DRINP")==0)
        {
          self->DRINP[i]=hcell;
          IBFilled++;
        }
        if(strcmp(rstrng,"DZINP")==0)
        {
          self->DZINP[i]=hcell;
          JBFilled++;
        }
      }
      if ((IBFilled==self->IB)&&(JBFilled==self->JB)) xmshread=3;
      if (strcmp(rstrng,"&END")==0) {
        xmshread=4;

      }
    }
      
    /* while(xmshread==1) */
    /* { */
    if(xmshread==4) {
      fprintf(stderr, "ERROR reading the meshdata from the sim05 file.\n");
      exit(0);
    }
    
 while(xtmeread==0)
    {
      fscanf(sim05,"%s",rstrng);
      if(strcmp(rstrng,"&XTME")==0){
        xtmeread=1;
      }
    }
    while(xtmeread==1)
    {
      /* fscanf(sim05,"%5s %*[=] %f %*s %i %*[,]",rstrng, &rreal, &power); */
      int i =0;
      fscanf(sim05,"%5s %*[=] %19s",rstrng,valstrng);
      NSSize=1;
      while (valstrng[i]!='\0') 
      {
        int EOFortDoub;
        if (valstrng[i]=='D') valstrng[i]='E';  /* replacing D with E so that C can read the
                                                 * float  */
        if (valstrng[i]==',') EOFortDoub=1;  /* we have more than one
                                              * cell with the same height */
        if (EOFortDoub==1) valstrng[i] = ' ';
        i++;
      }
      sscanf(valstrng, "%lG", &rreal);
      
      if(strcmp(rstrng,"TWFIN")==0) {
        self->TWFIN=rreal;
      }
      if(self->TWFIN!=0) xtmeread=2;
      if (strcmp(rstrng,"&END")==0) {
        xtmeread=4;
      }
    }

    /* while(xmshread==1) */
    /* { */
    if(xtmeread==4) {
      fprintf(stderr, "ERROR reading the time data from the sim05 file.\n");
      exit(0);
    }
    fclose(sim05);
  }
  self->JBP2 = self->JB+2;
  self->IBP2 = self->IB+2;
  self->MMS = self->JBP2*self->IBP2;
  free(self->PK);

  return 0;
 
}
