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

 /*@@
   @routine    WaveToyCXX_Startup
   @date       
   @author     Gabrielle Allen
   @desc 
   
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
extern "C" int WaveToyCXX_Startup(void)
{
   CCTK_RegisterBanner("WaveToyC++: Evolutions of a Scalar Field");
   return 0;
}
