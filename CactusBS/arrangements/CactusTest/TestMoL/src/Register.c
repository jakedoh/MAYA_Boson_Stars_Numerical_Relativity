#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include <assert.h>

/*****************************************************************************
 ************** External routines ********************************************
 *****************************************************************************/
void TestMoL_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;

  if(evolve_grid_function)
  {
    const int evolvedslow_gf_group =
      CCTK_GroupIndex(CCTK_THORNSTRING "::evolvedslow_gf");
    const int rhsslow_gf_group =
      CCTK_GroupIndex(CCTK_THORNSTRING "::rhsslow_gf");
    const int evolved_gf_group =
      CCTK_GroupIndex(CCTK_THORNSTRING "::evolved_gf");
    const int rhs_gf_group = CCTK_GroupIndex(CCTK_THORNSTRING "::rhs_gf");
    const int constrained_gf_group =
      CCTK_GroupIndex(CCTK_THORNSTRING "::constrained_gf");
    const int sandr_gf_group = CCTK_GroupIndex(CCTK_THORNSTRING "::sandr_gf");

    const int ierr_evolved_gf_slow =
      MoLRegisterEvolvedGroupSlow(evolvedslow_gf_group, rhsslow_gf_group);
    const int ierr_evolved_gf =
      MoLRegisterEvolvedGroup(evolved_gf_group, rhs_gf_group);
    const int ierr_constrained_gf =
      MoLRegisterConstrainedGroup(constrained_gf_group);
    const int ierr_saveandrestore_gf =
      MoLRegisterSaveAndRestoreGroup(sandr_gf_group);

    assert(!ierr_evolved_gf_slow && !ierr_evolved_gf && !ierr_constrained_gf &&
           !ierr_saveandrestore_gf);
  }

  if(evolve_grid_array)
  {
    const int evolved_ga_group =
      CCTK_GroupIndex(CCTK_THORNSTRING "::evolved_ga");
    const int rhs_ga_group = CCTK_GroupIndex(CCTK_THORNSTRING "::rhs_ga");
    const int constrained_ga_group =
      CCTK_GroupIndex(CCTK_THORNSTRING "::constrained_ga");
    const int sandr_ga_group = CCTK_GroupIndex(CCTK_THORNSTRING "::sandr_ga");

    const int ierr_evolved_ga =
      MoLRegisterEvolvedGroup(evolved_ga_group, rhs_ga_group);
    const int ierr_constrained_ga =
      MoLRegisterConstrainedGroup(constrained_ga_group);
    const int ierr_saveandrestore_ga =
      MoLRegisterSaveAndRestoreGroup(sandr_ga_group);

    assert(!ierr_evolved_ga && !ierr_constrained_ga && !ierr_saveandrestore_ga);
  }

  // TODO: add code to test changing of variable types
  // TODO: add code to test upgrades from sandr -> constrained -> evolved
}
