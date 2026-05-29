/*@@
   @file      IO.c
   @date      Thu Dec 27 2001
   @author    Gabrielle Allen
   @desc
   Write information about reduction operators
   @enddesc
@@*/

#include <stdio.h>

#include <cctk.h>
#include <cctk_Arguments.h>

 /*@@
   @routine    IOInfo
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

void IOInfo(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;

  int number = CCTK_NumIOMethods();

  printf("IOMethods\n");
  printf("---------\n\n");
  printf(" There are %d registered IO methods\n\n", number);
  
  if (number > 0)
  {
    printf(" Handle  | Name                           | Thorn           | Implementation\n");
    printf(" ---------------------------------------------------------------------------\n");

    for (int handle = 0; handle < number; handle++)
    {
      const struct IOMethod *method = CCTK_IOMethod(handle);
      if (method)
      {
        const char *name = CCTK_IOMethodName(handle);
        const char *imp  = CCTK_IOMethodImplementation(handle);
        const char *thorn = CCTK_ImplementationThorn(imp);
        printf(" %7d | %-30.30s | %-15.15s | %-15.15s \n",
               handle, name, thorn, imp);
      }
    }
  }

  printf("\n");
}
