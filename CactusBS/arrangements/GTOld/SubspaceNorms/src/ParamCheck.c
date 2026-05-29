#include <string.h>
#include <assert.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

//*******************************//
//****** External Routines ******//
//*******************************//
void SubspaceNorms_ParamCheck(CCTK_ARGUMENTS);

void SubspaceNorms_ParamCheck(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;

  for (int vol=0; vol<intvolumes; vol++) {
        int my_calc_norms_every;

        /* ignore any volume that ist not actually computed */
        my_calc_norms_every = calc_norms_every_vol[vol] >= 0 ? calc_norms_every_vol[vol] : calc_norms_every;
        assert(my_calc_norms_every >= 0);
        if (my_calc_norms_every == 0)
        {
                continue;
        }

	/* Attempt to get variable index of weight */
	int ierrW=CCTK_VarIndex(weight_var[vol]);
	if (ierrW < 0)
	{
		CCTK_VParamWarn(CCTK_THORNSTRING,"SubspaceNorms: Error checking weight variable %s.", weight_var[vol]);
	}

        /* Check scalar floor variables */
        if ( scalar_floor[vol] && strcmp(scalar_floor_var[vol],"")>0 ) {
           int ierrF = CCTK_VarIndex(scalar_floor_var[vol]);
           if ( ierrF < 0 )
              CCTK_VParamWarn(CCTK_THORNSTRING,"SubspaceNorms: Error checking variable %s for use as a floor in a scalar integration.", scalar_floor_var[vol] );
        }

        /* check if spherical surface index is valid */
        if (CCTK_Equals(centre_from[vol], "spherical surface"))
        { 
                int sn, sf_error;

                sf_error = 0;
                sn = surface_index[vol];
                if (sn < 0) {
                  CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                             "surface number sn=%d is invalid for volume %d", sn,vol);
                  sf_error = 1;
                } else if (sn >= sphericalsurfaces_nsurfaces) {
                  CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                             "surface number sn=%d too large, increase SphericalSurface::nsurfaces from its current value %d",
                             sn,sphericalsurfaces_nsurfaces);
                  sf_error = 1;
                }
                if (sf_error)
		  CCTK_PARAMWARN("SubspaceNorms: Error in checking spherical surfaces.");
        }
  }

}
