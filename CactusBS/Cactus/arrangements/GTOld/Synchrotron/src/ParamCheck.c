#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include <assert.h>

//*******************************//
//****** External Routines ******//
//*******************************//
void Synchro_Startup(CCTK_ARGUMENTS);

void Synchro_Startup(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;

  for (int vol=0; vol<intvolumes; vol++) {

        /* ignore any volume that ist not actually computed */
        int my_calc_syn_every;
        my_calc_syn_every = calc_syn_every_vol[vol] >= 0 ? calc_syn_every_vol[vol] : calc_syn_every;
        assert(my_calc_syn_every >= 0);
        if (my_calc_syn_every == 0)
                continue;

	/* Attempt to get variable index */
	int ierr=CCTK_VarIndex(weight_var[vol]);
	if (ierr < 0)
	{
		CCTK_PARAMWARN("Synchrotron: Error in checking weight variable.");
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
		  CCTK_PARAMWARN("Synchrotron: Error in checking spherical surfaces.");
        }
  }

}
