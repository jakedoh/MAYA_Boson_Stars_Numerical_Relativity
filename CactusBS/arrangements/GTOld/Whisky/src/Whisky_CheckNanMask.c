#include <assert.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void Whisky_CheckNanMask(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL output_value;
  
  const int varindex = CCTK_VarIndex("Whisky::whisky_nankind"); 
  assert(varindex >= 0);

  const int max_op = CCTK_ReductionHandle ("maximum");
  if (max_op < 0)  CCTK_WARN (CCTK_WARN_ABORT, "error");

  const int ierr = CCTK_Reduce (cctkGH, -1, max_op, 1, CCTK_VARIABLE_REAL, &output_value, 1, varindex);
  assert(ierr == 0);

  if(output_value > 0.0)
  {
    CCTK_OutputVarAsByMethod (cctkGH,
                              "Whisky::whisky_nankind",
                              "IOHDF5", "Whisky_NaNMask");
  }
}

