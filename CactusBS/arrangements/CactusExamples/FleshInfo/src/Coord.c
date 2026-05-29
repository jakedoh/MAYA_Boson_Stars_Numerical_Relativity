/*@@
   @file      Coord.c
   @date      Sat Dec 29 2001
   @author    Gabrielle Allen
   @desc
   Write information about coordinate systems
   @enddesc
   @version   $Id$
@@*/

#include <stdio.h>

#include "cctk.h"
#include "cctk_Arguments.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusExamples_FleshInfo_Coord_c)

void CoordInfo(CCTK_ARGUMENTS);

 /*@@
   @routine    IOInfo
   @date       Thu Dec 29 2001
   @author     Gabrielle Allen
   @desc 
   Write information about coordinate systems
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
@@*/

void CoordInfo(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS

  const char *name;
  const char *imp;
  const char *thorn;
  const char *coord;
  int number;
  int handle;
  int dir;
  int dim;
  int ierr,iperr;
  int ilower,iupper;
  CCTK_REAL lower,upper;

  number = CCTK_NumCoordSystems();

  printf("Coordinate Systems\n");
  printf("------------------\n\n");
  printf(" There are %d registered coordinate systems\n\n",number);
  
  if (number > 0)
  {
    printf(" Handle  | Name                      | Dim | Thorn           | Implementation\n");
    printf(" ---------------------------------------------------------------------------\n");

    for (handle = 0; handle < number; handle++)
      {
        name = CCTK_CoordSystemName(handle);
        if (name)
        {
          imp  = CCTK_CoordSystemImplementation(handle);
          thorn = CCTK_ImplementationThorn(imp);
          dim = CCTK_CoordSystemDim(name);
          printf(" %7d | %-25.25s | %1d   | %-15.15s | %-15.15s \n",handle,name,dim,thorn, imp);
        }
      }
  }

  if (number > 0)
  {
    for (handle = 0; handle < number; handle++)
    {
      name = CCTK_CoordSystemName(handle);
      if (name)
      {
        printf("\n");
        printf("    System: %s\n",name);
        printf("    -------\n\n");
        printf("      Direction | Name              | Computational Range | Physical Index Range\n");
        printf("      ---------------------------------------------------------------------------\n");

        dim = CCTK_CoordSystemDim(name);
          
        for (dir = 1; dir <= dim; dir++)
        {
          coord = CCTK_CoordName(dir,name);
          ierr = CCTK_CoordRange(cctkGH,&lower,&upper,dir,NULL,name);
          iperr = CCTK_CoordRangePhysIndex(cctkGH,&ilower,&iupper,dir,NULL,name);
          if (coord)
          {
            printf("      %9d | %-17.17s |",dir,coord);
            if (ierr == 0)
            {
              printf(" %-6.5f to %-6.5f |",lower,upper);
            }
            else
            {
              printf(" - N/A -             |");
            }

            if (iperr == 0)
            {
              printf(" %5d to %5d\n",ilower,iupper);
            }
            else
            {
              printf(" - N/A -             \n");
            }

          }
        }
      }
    }
  }

  printf("\n");

  return;

}
