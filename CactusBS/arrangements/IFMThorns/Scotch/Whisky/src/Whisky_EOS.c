 /*@@
   @file      Whisky_EOS.c
   @date      Wed Feb  6 18:25:33 2002
   @author    
   @desc 
   Sets the EOS handle number from CactusEOS for use by
   all the Whisky routines
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include <stdio.h>
#include <stdlib.h>

#include "EOS_Base.h"

 /*@@
   @routine    Whisky_EOSHandle
   @date       Wed Feb  6 18:28:01 2002
   @author     Ian Hawke
   @desc 
   Sets the EOS handle number
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void Whisky_EOSHandle(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int ierr = 0;
  char *infoline;
  
  *whisky_eos_handle = EOS_Handle(whisky_eos_table);
  
  if (*whisky_eos_handle < 0) {
    ierr = *whisky_eos_handle;
    CCTK_WARN(0, "Table handle is incorrect");
  }

  *whisky_polytrope_handle = EOS_Handle("2D_Polytrope");
  
  if (*whisky_polytrope_handle < 0) {
    ierr = *whisky_polytrope_handle;
    CCTK_WARN(0, "Whisky requires the 2D_Polytrope thorn to be active for the atmosphere");
  }

  infoline = (char *)malloc(100*sizeof(char));
  if (!infoline) {
    CCTK_WARN(0, "Failed to malloc an array of 100 chars in Whisky_EOS.");
  }
  else {
    sprintf(infoline,"Whisky will use the %s equation of state.",
            whisky_eos_table);
    CCTK_INFO(infoline);
  }
  free(infoline);
  infoline = NULL;

}
