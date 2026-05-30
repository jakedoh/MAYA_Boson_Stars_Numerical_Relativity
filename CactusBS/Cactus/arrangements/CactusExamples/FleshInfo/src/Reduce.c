/*@@
   @file      Reduce.c
   @date      Thu Dec 27 2001
   @author    Gabrielle Allen
   @desc
   Write information about reduction operators
   @enddesc
   @version   $Id$
@@*/

#include <stdio.h>

#include "cctk.h"
#include "cctk_Arguments.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusExamples_FleshInfo_Reduce_c)

void ReduceInfo(CCTK_ARGUMENTS);

 /*@@
   @routine    ReduceInfo
   @date       Thu Dec 27 2001
   @author     Gabrielle Allen
   @desc 
   Write information about reduction operators
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
@@*/

void ReduceInfo(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS

  const char *name;
  const char *imp;
  const char *thorn;
  int number;
  int handle;

  number = CCTK_NumReduceOperators();

  printf("Reduction Operators\n");
  printf("-------------------\n\n");
  printf(" There are %d registered reduction operators\n\n",number);
  
  if (number > 0)
  {
    printf(" Handle  | Name                           | Thorn           | Implementation\n");
    printf(" ---------------------------------------------------------------------------\n");

    for (handle = 0; handle < number; handle++)
      {
        name = CCTK_LocalArrayReduceOperator(handle);
        if (name)
        {
          imp  = CCTK_LocalArrayReduceOperatorImplementation(handle);
          thorn = CCTK_ImplementationThorn(imp);
            printf(" %7d | %-30.30s | %-15.15s | %-15.15s \n",handle,name,thorn, imp);
        }
      }
  }

  printf("\n");

  return;

  printf("\n");

}
