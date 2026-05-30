#include <assert.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

//////////////////////////////////////////////////////////
// Centres: copied from Ian Hinder's RegridBoxes thorn
//////////////////////////////////////////////////////////

namespace AHMask {

void get_centre_parameters(CCTK_ARGUMENTS, bool initial, int c, 
                                  CCTK_REAL *radius, CCTK_REAL coords[3])
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  assert(c >= 0);
  assert(c < number_of_centres);
  assert(radius != 0);
  assert(coords != 0);

  if (!CCTK_EQUALS(centre_from[c], "spherical surface"))
  {
    if (surface_index[c] != -1)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Parameter %s::surface_index[%d] has been set, but centre %d "
                  "is not being set from a spherical surface, so this parameter "
                  "will be ignored.", CCTK_THORNSTRING, c, c);
    }
  }

  if (CCTK_EQUALS(centre_from[c], "parameter") || initial)
  {
    coords[0] = centre_x[c]; coords[1] = centre_y[c]; coords[2] = centre_z[c];
    switch(c)
    {
      case 0: *radius = centre_0_radius; break;
      case 1: *radius = centre_1_radius; break;
      case 2: *radius = centre_2_radius; break;
      default:  assert(0);
    }
    if (*radius == -1.)
    {
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Parameter %s::centre_%d_radius has not been set.",
                  CCTK_THORNSTRING, c);
    }
  }
  else if (CCTK_EQUALS(centre_from[c], "spherical surface"))
  {
    assert(surface_index[c] != -1);

    assert(surface_index[c] < nsurfaces);
    assert(sf_valid[surface_index[c]]);
    
    coords[0] = sf_centroid_x[surface_index[c]];
    coords[1] = sf_centroid_y[surface_index[c]];
    coords[2] = sf_centroid_z[surface_index[c]];
    *radius   = sf_max_radius[surface_index[c]];
  }
  else
  {
    assert(false); // Should never get here
  }

  if (verbose > 0)
  {
    CCTK_VInfo(CCTK_THORNSTRING,
               "centre %d, at (%g,%g,%g) radius %g", c, coords[0], coords[1], coords[2], *radius);
  }
}

}
