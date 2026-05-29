#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

/*
 * External functions
 */
void idscalarfield_ParamCheck(CCTK_ARGUMENTS);

/*
 * Implementations
 */

void idscalarfield_ParamCheck(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;
  
  /* check that variables given in scalarPhi and scalarPi at least exist */
  CCTK_INT scalarphi = CCTK_VarIndex(scalarPhi);
  CCTK_INT scalarpi = CCTK_VarIndex(scalarPi);

  if(scalarphi < 0 || scalarpi < 0)
  {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Could not find variables scalarPhi='%s' and scalarPi='%s', found indices: %d,%d",
               scalarPhi, scalarPi, scalarphi, scalarpi);
    CCTK_PARAMWARN("Error when checking variables.");
  }
}
