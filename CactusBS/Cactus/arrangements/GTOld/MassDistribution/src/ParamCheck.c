#include <string.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "MassDistribution.h"

/*******************************
 ****** External Routines ******
 *******************************/
void MassDistribution_ParamCheck(CCTK_ARGUMENTS);

void MassDistribution_ParamCheck(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  int ierr;

  for (int vol=0; vol<intvolumes; vol++) {
    /* Attempt to get variable indices */
    ierr = CCTK_VarIndex(weight_variable[vol]);
    if (ierr < 0)
    {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING, "Could not find variable '%s'", 
          weight_variable[vol]);
      CCTK_PARAMWARN("MassDistribution: Error in checking weight variable.");
    }
    if (0 != MassDistribution_validate_special_variable(cctkGH,
          bin_variable[vol], NUM_SPECIAL_VARIABLE_TYPES, special_variables))
    {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING, "Could not find variable '%s'", 
          bin_variable[vol]);
      CCTK_PARAMWARN("MassDistribution: Error in checking binning variable.");
    }
    if (0 != MassDistribution_validate_special_variable(cctkGH,
          mass_variable[vol], NUM_MASS_VARIABLE_TYPES, mass_variables))
    {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING, "Could not find variable '%s'", 
          mass_variable[vol]);
      CCTK_PARAMWARN("MassDistribution: Error in checking mass variable.");
    }
    if (CCTK_EQUALS(bin_mode[vol], "manual") && 
        bin_variable_maximum[vol] <= bin_variable_minimum[vol])
    {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING, "Min-Max range invalid for %s mode bining of volume %d", 
          bin_mode[vol], vol);
      CCTK_PARAMWARN("MassDistribution: maximum of binning variable not larger than minimum.");
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
        CCTK_PARAMWARN("MassDistribution: Error in checking spherical surfaces.");
    }
  }
}
