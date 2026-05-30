 /*@@
   @file      Startup.c
   @date      
   @author    Gabrielle Allen
   @desc 
              Register banner 
   @enddesc 
   @version $Header$
 @@*/

#include "cctk.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusWave_WaveToyC_Startup_c)

int WaveToyC_Startup(void);

 /*@@
   @routine    WaveToyC_Startup
   @date       
   @author     Gabrielle Allen
   @desc 
   
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
int WaveToyC_Startup(void)
{

   const char *banner = "WaveToyC: Evolutions of a Scalar Field";

   CCTK_RegisterBanner(banner);

   return 0;
}
