/*@@
   @file      Boundary.c
   @date      Friday 18th July 2003
   @author    Gabrielle Allen
   @desc
              Use any provided boundary condition with WaveToy
   @enddesc
   @version   $Id$
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "util_Table.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusWave_WaveToyExtra_Boundary_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

void WaveToyExtra_Boundary(CCTK_ARGUMENTS);

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

static int handle=-1;

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    WaveToyC_Boundary
   @date       Friday 18th July 2003
   @author     Gabrielle Allen
   @desc 
               Mark wavetoy variables for custom boundary condition
   @enddesc 
   @history 
   @endhistory 
@@*/

void WaveToyExtra_Boundary(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr=0;

  if (CCTK_EQUALS(bound,"custom"))
  {
    if (!CCTK_EQUALS(custom_options,""))
    {
      if (handle == -1)
      {
        handle = Util_TableCreateFromString(custom_options);
        if (handle < 0)
        {
          CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
                     "WaveToyC_Boundaries: Error creating table for "
                     "boundary condition %s",custom_bound);
        }
      }
    }
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle, 
                                   "wavetoy::phi",custom_bound);  
  }    

  if (ierr < 0) 
  {
    CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
               "WaveToyC_Boundaries: Error selecting boundary "
               "condition %s",bound);
  }

  return;
}


 /*@@
   @routine    WaveToyC_Terminate
   @date       Friday 18th July 2003
   @author     Gabrielle Allen
   @desc 
               Tidy up wavetoy
   @enddesc 
   @history 
   @endhistory 
@@*/
void WaveToyC_Terminate(CCTK_ARGUMENTS)
{
  Util_TableDestroy(handle);
}
