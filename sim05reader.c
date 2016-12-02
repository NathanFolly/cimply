/* Reader for sim05 data in c. Data to be used in the structural FE code
 * cimply that is part of the SIMMER code. Reads in strict format and may be
 * sensitive to changes in the sim05 syntax. */
#include <petscsys.h>
#include <stdio.h>
#include <string.h>
#include "cimplySimmerUtils.h"
PetscErrorCode sim05tocimply(SimmerData SiDat)
{

  /* TODO: make this more robust. Likely to break if fromat of input varies. */
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
  int xtmeread = 0;  /* Same principle as for xmshread. Just for the time context */
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
      if(strcmp(rstrng,"IB")==0) SiDat.IB=rint;
      if(strcmp(rstrng,"JB")==0) SiDat.JB=rint;
      if(SiDat.IB*SiDat.JB!=0) xmshread=2;  /* once both variables
                                                       * are assigned values:continue */
    }
    ierr = PetscMalloc1(SiDat.IB,&SiDat.DRINP);CHKERRQ(ierr);
    ierr = PetscMalloc1(SiDat.JB,&SiDat.DZINP);CHKERRQ(ierr);
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
          SiDat.DRINP[i]=hcell;
          IBFilled++;
        }
        if(strcmp(rstrng,"DZINP")==0)
        {
          SiDat.DZINP[i]=hcell;
          JBFilled++;
        }
      }
      if ((IBFilled==SiDat.IB)&&(JBFilled==SiDat.JB)) xmshread=3;
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
        SiDat.TWFIN=rreal;
      }
      if(SiDat.TWFIN!=0) xtmeread=2;
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
  SiDat.JBP2 = SiDat.JB+2;
  SiDat.IBP2 = SiDat.IB+2;
  SiDat.MMS = SiDat.JBP2*SiDat.IBP2;
  ierr = PetscFree(SiDat.PK);CHKERRQ(ierr);
  ierr = PetscMalloc1(SiDat.MMS, &SiDat.PK);CHKERRQ(ierr);

  return(0);
}
