/* $Header$ */

#include <math.h>

#include "cctk.h"

int we_isnan (const CCTK_REAL x)
{
#ifdef HAVE_ISNAN
  return isnan(x);
#else
  return 0;
#endif
}

int CCTK_FCALL CCTK_FNAME(we_isnan) (const CCTK_REAL * restrict const x)
{
  return we_isnan(*x);
}
