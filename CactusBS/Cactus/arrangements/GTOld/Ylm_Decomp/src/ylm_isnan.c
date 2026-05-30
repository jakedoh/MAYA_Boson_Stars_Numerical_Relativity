// get isnan

#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

int ylm_isnan (const CCTK_REAL x)
{
#ifdef HAVE_ISNAN
  return isnan(x);
#else
  if (x!=x) {
    return 1;
  }
  else {
    return 0;
  }
#endif
}

void CCTK_FCALL CCTK_FNAME(ylm_isnan) (CCTK_INT retval,const CCTK_REAL * restrict const x)
{
  retval=ylm_isnan(*x);
}
