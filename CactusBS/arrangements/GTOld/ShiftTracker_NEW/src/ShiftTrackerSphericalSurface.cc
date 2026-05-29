#include <cassert>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"


extern "C" void ShiftTracker_Set_SphericalSurface(CCTK_ARGUMENTS);

extern "C" void ShiftTracker_Set_SphericalSurface(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if( verbose > 2)
    CCTK_VInfo(CCTK_THORNSTRING,
        "ShiftTracker_Set_SphericalSurface: updating the sphericalsurfaces");

  // For each tracker set spherical surface information
  CCTK_INT sind=-1;
  for( int i=0; i < num_trackers; i++ )  {
    sind = surface_index[i];
    if (sind==-1) {
      continue;
    }
    else {
      assert(sind >=0 && sind < nsurfaces);
      sf_valid[sind]=1;
      sf_active[sind]=1;
      sf_centroid_x[sind]=st_x[i];
      sf_centroid_y[sind]=st_y[i];
      sf_centroid_z[sind]=st_z[i];

      if( verbose > 0 ) { 
        CCTK_VInfo(CCTK_THORNSTRING,
          "ShiftTracker_Set_SphericalSurface(): Tracker %d in Surface: %d",
          i,sind);
      }
    } 
  }
}

