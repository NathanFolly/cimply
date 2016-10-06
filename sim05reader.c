/* Reader for sim05 data in c. Data to be used in the structural FE code
 * cimply that is part of the SIMMER code. Reads in strict format and may be
 * sensitive to changes in the sim05 syntax. */
#include <petscsys.h>
#include <stdio.h>
#include <string.h>
#include "cimply.h"
PetscErrorCode sim05tocimply()
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
      if(strcmp(rstrng,"IB")==0) SimmerData.IB=rint;
      if(strcmp(rstrng,"JB")==0) SimmerData.JB=rint;
      if(SimmerData.IB*SimmerData.JB!=0) xmshread=2;  /* once both variables
                                                       * are assigned values:continue */
    }
    ierr = PetscMalloc1(SimmerData.IB,&SimmerData.DRINP);CHKERRQ(ierr);
    ierr = PetscMalloc1(SimmerData.JB,&SimmerData.DZINP);CHKERRQ(ierr);
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
          SimmerData.DRINP[i]=hcell;
          IBFilled++;
        }
        if(strcmp(rstrng,"DZINP")==0)
        {
          SimmerData.DZINP[i]=hcell;
          JBFilled++;
        }
      }
      if ((IBFilled==SimmerData.IB)&&(JBFilled==SimmerData.JB)) xmshread=3;
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
    fclose(sim05);
  }
  SimmerData.JBP2 = SimmerData.JB+2;
  SimmerData.IBP2 = SimmerData.IB+2;
  SimmerData.MMS = SimmerData.JBP2*SimmerData.IBP2;
  ierr = PetscFree(SimmerData.PK);CHKERRQ(ierr);
  ierr = PetscMalloc1(SimmerData.MMS, &SimmerData.PK);CHKERRQ(ierr);

  return(0);
}
