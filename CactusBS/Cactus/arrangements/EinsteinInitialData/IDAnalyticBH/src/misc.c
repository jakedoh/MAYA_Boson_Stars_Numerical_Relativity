 /*@@
   @file      misc.c
   @date      23.Mar.2003
   @author    Jonathan Thornburg <jthorn@aei.mpg.de>
   @desc 
      Misc support routines for this thorn.
   @enddesc 
   @version $Header$
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"

#include "IDAnalyticBH.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusEinstein_IDAnalyticBH_Misc_c)


/*@@
  @routine    IDAnalyticBH_zero_CCTK_REAL_array
  @date       23.Mar.2003
  @author     Jonathan Thornburg <jthorn@aei.mpg.de>
  @desc       Sets a CCTK_REAL array to all (floating-point) zeros.
  @enddesc 
  @@*/
void IDAnalyticBH_zero_CCTK_REAL_array(int N, CCTK_REAL A[])
{
int i;
  for (i = 0 ; i < N ; ++i)
  {
    A[i] = 0.0;
  }
}
