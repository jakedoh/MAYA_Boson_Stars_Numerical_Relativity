/*@@
   @file      TimerInfo.c
   @date      Wed Oct 17 2001
   @author    Gabrielle Allen
   @desc
   Scheduled rountines for printing all information from timers
   @enddesc
   @version $Header$
 @@*/

#include <stddef.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Timers.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusExamples_TimerInfo_TimerInfo_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

void TimerInfo(CCTK_ARGUMENTS);
 
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
   @routine    TimerInfo
   @date       Wed Oct 17 2001
   @author     Gabrielle Allen
   @desc
   Print timing information for all clocks
   @enddesc
@@*/

void TimerInfo(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  const char *newtimer=NULL;
  const char *newclock=NULL;

  if (every > 0 && (cctk_iteration % every)== 0) 
  {
    if (!CCTK_Equals(timer,"all"))
    {
      newtimer = timer;
    }

    if (!CCTK_Equals(clock,"all"))
    {
      newclock = clock;
    }

    CCTK_TimerPrintData(newtimer,newclock);
  }
}

