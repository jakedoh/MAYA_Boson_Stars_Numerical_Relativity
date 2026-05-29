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
  CCTK_INT scalarphia = CCTK_VarIndex(scalarPhia);
  CCTK_INT scalarpia = CCTK_VarIndex(scalarPia);
  CCTK_INT scalarphib = CCTK_VarIndex(scalarPhib);
  CCTK_INT scalarpib = CCTK_VarIndex(scalarPib);

  if(scalarphia < 0 || scalarpia < 0 || scalarphib < 0 || scalarpib < 0)
  {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Could not find variables scalarPhia='%s', scalarPia='%s', scalarPhib='%s' and scalarPib='%s', found indices: %d,%d",
               scalarPhia, scalarPia, scalarphib, scalarpib);
    CCTK_PARAMWARN("Error when checking variables.");
  }
}
