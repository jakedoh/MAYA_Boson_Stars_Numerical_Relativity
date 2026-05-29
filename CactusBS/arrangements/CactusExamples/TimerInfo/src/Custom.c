/*@@
   @file      Custom.c
   @date      Wed Oct 17 2001
   @author    Gabrielle Allen
   @desc
   Customised output of timer data
   @enddesc
   @version $Header$
 @@*/

#include <stdio.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusExamples_TimerInfo_Custom_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/
 
/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

void PrintTimerData (const char *timername, const cTimerData *info);

 /*@@
   @routine    PrintTimerData
   @date       Wed Oct 17 2001
   @author     Gabrielle Allen
   @desc
   Customise output of timer data
   @enddesc
@@*/

void PrintTimerData (const char *timername, const cTimerData *info)
{
  int i;

  printf("\nResults for timer called '%s'\n", timername);
  for (i = 0; i < info->n_vals; i++)
  {
    switch (info->vals[i].type)
    {
      case val_int:
        printf("  %s: %d %s\n", 
               info->vals[i].heading,info->vals[i].val.i, 
               info->vals[i].units);
        break;
        
      case val_long:
        printf("  %s: %d %s\n", 
               info->vals[i].heading,(int) info->vals[i].val.l, 
               info->vals[i].units);
        break;
        
      case val_double:
        printf("  %s: %.3f %s\n", 
               info->vals[i].heading,info->vals[i].val.d, 
               info->vals[i].units);
        break;
        
      default:
        CCTK_WARN(1, "Unknown data type for timer info");
        break;
    }
  }
}
