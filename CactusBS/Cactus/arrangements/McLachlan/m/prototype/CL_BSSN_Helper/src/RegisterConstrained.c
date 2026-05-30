#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <cctk.h>
#include <cctk_Arguments.h>

static void
register_constrained (char const * restrict const gn);

void
CL_BSSN_RegisterConstrained (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  
  register_constrained ("ADMBase::metric");
  register_constrained ("ADMBase::curv");
  register_constrained ("ADMBase::lapse");
  register_constrained ("ADMBase::shift");
  register_constrained ("ADMBase::dtshift");
}

static void
register_constrained (char const * restrict const gn)
{
  assert (gn);
  
  int const gi = CCTK_GroupIndex (gn);
  int const ierr = MoLRegisterConstrainedGroup (gi);
  assert (! ierr);
}
