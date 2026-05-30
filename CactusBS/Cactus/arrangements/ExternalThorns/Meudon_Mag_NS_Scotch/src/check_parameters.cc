#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>



extern "C"
void ID_Mag_NS_Scotch_check_parameters (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (not CCTK_EQUALS (initial_data,    "ID_Mag_NS_Scotch") or
      not CCTK_EQUALS (initial_lapse,   "ID_Mag_NS_Scotch") or
      not CCTK_EQUALS (initial_shift,   "ID_Mag_NS_Scotch") or
      not (CCTK_EQUALS (initial_dtlapse, "ID_Mag_NS_Scotch") or
           CCTK_EQUALS (initial_dtlapse, "none")) or
      not (CCTK_EQUALS (initial_dtshift, "ID_Mag_NS_Scotch") or
           CCTK_EQUALS (initial_dtshift, "none")) or
      not CCTK_EQUALS (initial_hydro,   "ID_Mag_NS_Scotch"))
  {
    CCTK_PARAMWARN ("The parameters ADMBase::initial_data, ADMBase::initial_lapse, ADMBase::initial_shift, ADMBase::initial_dtlapse, ADMBase::initial_dtshift, and HydroBase::initial_hydro must all be set to the value \"ID_Mag_NS_Scotch\" or \"none\"");
  }

  if (not CCTK_EQUALS(initial_Bvec,  "ID_Mag_NS_Scotch") and
      not CCTK_EQUALS (initial_Avec,  "ID_Mag_NS_Scotch"))
   {
    CCTK_PARAMWARN ("Either HydroBase::initial_Bvec or HydroBase::initial_Avec must be set to \"ID_Mag_NS_Scotch\" or \"none\"");
   }
}
