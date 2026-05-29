 /*@@
   @file      Startup.c
   @date      Wed Sep 13 21:26:56 2000
   @author    Tom Goodale
   @desc
              Thorn HTTPDExtra's startup routines.
   @enddesc
   @version   $Id$
 @@*/

#include "cctk.h"
#include "httpextra_HostNames.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusConnect_HTTPDExtra_Startup_c)


/********************************************************************
 ********************    External Routines   ************************
 ********************************************************************/
int HTTPUTILS_Startup(void);


 /*@@
   @routine    HTTP_Startup
   @date       Wed Sep 13 21:26:56 2000
   @author     Tom Goodale
   @desc
               Startup routine for the webserver.
   @enddesc
   @calls      HTTPDExtra_CollateHostData
               HTTPUTILS_RegisterPages
               HTTP_util_RegisterIOPages
               HTTPDExtra_RegisterProcessorsPages
	       HTTPDExtra_RegisterTimerInfoPages

   @returntype int
   @returndesc
               0 for success
   @endreturndesc
@@*/
int HTTPUTILS_Startup(void)
{
  /* FIXME: should put this prototype into a header file */
  extern int HTTPUTILS_RegisterPages(void);
  extern int HTTP_util_RegisterIOPages(void);
  extern int HTTPDExtra_RegisterProcessorsPages(void);
  extern int HTTPDExtra_RegisterTimerInfoPages(void);


  HTTPDExtra_CollateHostData();

  HTTPUTILS_RegisterPages();

  HTTP_util_RegisterIOPages();

  HTTPDExtra_RegisterProcessorsPages();

  HTTPDExtra_RegisterTimerInfoPages();
  
  return 0;
}
