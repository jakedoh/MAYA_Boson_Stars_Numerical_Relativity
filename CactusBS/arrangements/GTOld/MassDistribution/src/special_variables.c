#include <string.h>
#include <assert.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "MassDistribution.h"

/*******************************
 ******** External data ********
 *******************************/
const char *special_variables[NUM_SPECIAL_VARIABLE_TYPES] = SPECIAL_VARIABLE_TYPE_NAMES;
const char *mass_variables[NUM_MASS_VARIABLE_TYPES] = MASS_VARIABLE_TYPE_NAMES;
const char *bin_scalings[NUM_BIN_SCALINGS] = BIN_SCALING_NAMES;

/*******************************
 ****** External Routines ******
 *******************************/
int MassDistribution_resolve_keyword_parameter(const char *name, 
    const char *kw, int num_kws, const char *kws[])
{
  for(int i = 0 ; i < num_kws ; i++)
  {
    if(CCTK_Equals(kw, kws[i]))
      return i;
  }

  CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
             "Invalid keyword '%s' for %s.", kw, name);
  return -1; /* NOTREACHED */
}

void MassDistribution_resolve_special_variable(cGH *cctkGH, const char *varstring, 
    int *vartype, CCTK_REAL **varptr, int num_vartypes, const char *vartypes[])
{
  int ierr;

  // check arguments
  assert(varptr);
  assert(varstring);
  assert(vartype);
  assert(vartypes);
  assert(cctkGH);
  assert(num_vartypes > 0);

  if(strstr(varstring, "::"))
  {
    *vartype = 0;
    *varptr = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH,0,varstring));
    ierr = *varptr == NULL;
  }
  else
  {
    ierr = 1;
    for(int i = 1 ; i < num_vartypes ; i++)
    {
      if(CCTK_Equals(varstring, vartypes[i]))
      {
        *vartype = i;
        *varptr = NULL;
        ierr = 0;
        break;
      }
    }
  }
  if(ierr)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
             "Invalid variable name '%s'.", varstring);
  }

  return;
}

int MassDistribution_validate_special_variable(cGH *cctkGH, const char *varstring, 
    int num_vartypes, const char *vartypes[])
{
  int ierr;

  // check arguments
  assert(varstring);
  assert(vartypes);
  assert(cctkGH);
  assert(num_vartypes > 0);

  if(strstr(varstring, "::"))
  {
    ierr = CCTK_VarIndex(varstring) < 0;
  }
  else
  {
    ierr = 1;
    for(int i = 1 ; i < num_vartypes ; i++)
    {
      if(CCTK_Equals(varstring, vartypes[i]))
      {
        ierr = 0;
        break;
      }
    }
  }

  return ierr;
}
