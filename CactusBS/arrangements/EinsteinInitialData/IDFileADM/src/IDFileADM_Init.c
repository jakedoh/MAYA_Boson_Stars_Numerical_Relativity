/* $Header$ */

#include <assert.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

/** Initialise the conformal_state.  */
void IDFileADM_Init (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (CCTK_EQUALS(initial_data, "read from file")) {
    if (CCTK_EQUALS (metric_type, "physical")) {
      *conformal_state = 0;
    } else if (CCTK_EQUALS (metric_type, "static conformal")) {
      if (CCTK_EQUALS (conformal_storage, "factor")) {
        *conformal_state = 1;
      } else if (CCTK_EQUALS (conformal_storage, "factor+derivs")) {
        *conformal_state = 2;
      } else if (CCTK_EQUALS (conformal_storage, "factor+derivs+2nd derivs")) {
        *conformal_state = 3;
      } else {
        CCTK_WARN (0, "Unknown conformal_storage type");
      }
    } else {
      CCTK_WARN (0, "Unknown metric type");
    }
  }
}
