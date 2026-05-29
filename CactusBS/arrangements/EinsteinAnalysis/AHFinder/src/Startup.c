 /*@@
   @file      Startup.c
   @date      Sat 10th November 2001
   @author    Gabrielle Allen
   @desc
              Startup routines for AHFinder.
   @enddesc
   @version   $Id$
 @@*/

#include "cctk.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusEinstein_AHFinder_Startup_c)

int AHFinder_TimeForOutput (const cGH *GH, int vindex);

/********************************************************************
 ********************    External Routines   ************************
 ********************************************************************/
int AHFinder_Startup (void);

 /*@@
   @routine   AHfinder_Startup
   @date      Sat 10th November 2001
   @author    Gabrielle Allen
   @desc
              The startup registration routine for AHFinder.
              Registers AHFinder as an IO Method and provide the
              only necessary method TimeForOutput
   @enddesc
   @calls     CCTK_RegisterGHExtensionSetupGH
@@*/
int AHFinder_Startup (void)
{
  int handle;

  handle = CCTK_RegisterIOMethod ("AHFinder");
  CCTK_RegisterIOMethodTimeToOutput (handle, AHFinder_TimeForOutput);

  return 0;
}




