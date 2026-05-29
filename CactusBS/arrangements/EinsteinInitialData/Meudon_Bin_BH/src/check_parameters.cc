#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>



extern "C"
void ID_Bin_BH_check_parameters (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (not CCTK_EQUALS (initial_data,    "ID_Bin_BH") or
      not CCTK_EQUALS (initial_lapse,   "ID_Bin_BH") or
      not CCTK_EQUALS (initial_shift,   "ID_Bin_BH") or
      not (CCTK_EQUALS (initial_dtlapse, "ID_Bin_BH") or
           CCTK_EQUALS (initial_dtlapse, "none")) or
      not (CCTK_EQUALS (initial_dtshift, "ID_Bin_BH") or
           CCTK_EQUALS (initial_dtshift, "none")))
  {
    CCTK_PARAMWARN ("The parameters ADMBase::initial_data, ADMBase::initial_lapse, ADMBase::initial_shift, ADMBase::initial_dtlapse, and ADMBase::initial_dtshift must all be set to the value \"ID_Bin_BH\" or \"none\"");
  }
}
