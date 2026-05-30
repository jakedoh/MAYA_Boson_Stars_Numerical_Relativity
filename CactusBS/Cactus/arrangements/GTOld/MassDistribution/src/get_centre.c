#include <assert.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "MassDistribution.h"

/*******************************
 ****** External Routines ******
 *******************************/
void MassDistribution_get_centre(CCTK_ARGUMENTS, int c, CCTK_REAL *xc, CCTK_REAL *yc, CCTK_REAL *zc);

void MassDistribution_get_centre(CCTK_ARGUMENTS, int c, CCTK_REAL *xc, CCTK_REAL *yc, CCTK_REAL *zc)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  assert(c >= 0);
  assert(c < intvolumes);
  assert(xc != 0);
  assert(yc != 0);
  assert(zc != 0);

  if (!CCTK_EQUALS(centre_from[c], "spherical surface"))
  {
    if (surface_index[c] != -1)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Parameter %s::surface_index[%d] has been set, but centre %d "
                  "is not being set from a spherical surface, so this parameter will be "
                  "ignored.", CCTK_THORNSTRING, c, c);
    }
  }

  if (CCTK_EQUALS(centre_from[c], "parameter"))
  {
    *xc = centre_x[c]; *yc = centre_y[c]; *zc = centre_z[c];
  }
  else if (CCTK_EQUALS(centre_from[c], "spherical surface"))
  {
    int sn=surface_index[c];
    if (sn < 0) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "surface number sn=%d is invalid for volume %d", sn,c);
    } else if (sn>=sphericalsurfaces_nsurfaces) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "surface number sn=%d too large, increase SphericalSurface::nsurfaces from its current value %d",
                 sn,sphericalsurfaces_nsurfaces);
    }
    if (sf_valid[sn]<=0) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "didn't find valid volume centre for sn=%d, vol=%d",sn,c);
    }

    *xc = sf_centroid_x[surface_index[c]];
    *yc = sf_centroid_y[surface_index[c]];
    *zc = sf_centroid_z[surface_index[c]];
  }
  else
  {
    assert(0); /* Should never get here */
  }

  if (verbose >= 3)
  {
    CCTK_VInfo(CCTK_THORNSTRING,
               "centre %d, at (%g,%g,%g) radius %g", c, *xc, *yc, *zc, int_wi_rad[c]);
  }
}
