#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include <assert.h>

/*****************************************************************************
 ************** External routines ********************************************
 *****************************************************************************/
void TestMoL_InitVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  cGroupDynamicData groupdata;

  const int group = CCTK_GroupIndex(CCTK_THORNSTRING "::evolved_ga");
  const int ierr = CCTK_GroupDynamicData(cctkGH, group, &groupdata);
  if(ierr != 0)
  {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Failed to obtain information about evolved_ga group: %d",
                ierr);
  }
  assert(groupdata.dim == 1);

  for(int i = 0 ; i < groupdata.lsh[0] ; i++)
  {
    evolved_ga[i] = -1.;
    sandr_ga[i] = 1.;
    sandr_ga_p[i] = 1.5;
    sandr_ga_p_p[i] = 1.75;
    constrained_ga[i] = 2.;
    constrained_ga_p[i] = 2.5;
    constrained_ga_p_p[i] = 2.75;
  }

#pragma omp parallel
  CCTK_LOOP3_ALL(TestMoL_InitVars, cctkGH, i,j,k)
  {
    const int idx = CCTK_GFINDEX3D(cctkGH, i,j,k);
    evolved_gf[idx] = -1.;
    evolvedslow_gf[idx] = -1.;
    sandr_gf[idx] = 1.;
    sandr_gf_p[idx] = 1.5;
    sandr_gf_p_p[idx] = 1.75;
    constrained_gf[idx] = 2.;
    constrained_gf_p[idx] = 2.5;
    constrained_gf_p_p[idx] = 2.75;
  } CCTK_ENDLOOP3_ALL(TestMoL_InitVars);
}
