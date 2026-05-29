/* $Header$ */

#include <assert.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void CCTK_FCALL
CCTK_FNAME(NoExcision_Reduce) (int const * cctk_iteration,
                               int const * cctk_lsh,
                               CCTK_REAL * rhs,
                               CCTK_REAL const * x,
                               CCTK_REAL const * y,
                               CCTK_REAL const * z);

void
NoExcision_Reduce (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  int var, rhs;
  void * ptr;
  
  for (var=0; var<CCTK_NumVars(); ++var) {
    rhs = MoLQueryEvolvedRHS (var);
    if (rhs >= 0) {
      ptr = CCTK_VarDataPtrI (cctkGH, 0, rhs);
      assert (ptr);
      CCTK_FNAME(NoExcision_Reduce) (& cctk_iteration, cctk_lsh, ptr, x, y, z);
    }
  }
}
