#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "cctk.h" 
#include "cctk_Arguments.h" 
#include "cctk_Parameters.h" 
#include "cctk_Functions.h"
#include "util_Table.h"
#include "minsearch.h"

//////////////////////////////////////////////////////////////////////
// Internal functions
//////////////////////////////////////////////////////////////////////
static CCTK_INT steer_all_targets(CCTK_POINTER_TO_CONST target, const CCTK_INT target_type, CCTK_STRING targetname, CCTK_POINTER callback_arg);
static size_t strncpyw(char *dst, const char *src, size_t n);

//////////////////////////////////////////////////////////////////////
// Service functions
//////////////////////////////////////////////////////////////////////

extern "C" CCTK_INT MinSearch_TraverseTargetString(cGH *cctkGH, CCTK_STRING targetstring, 
        CCTK_INT (*callback)(CCTK_POINTER_TO_CONST target, const CCTK_INT target_type, CCTK_STRING targetname, CCTK_POINTER callback_arg),
        CCTK_POINTER callback_arg)
{
  char name_buffer[256];
  CCTK_POINTER_TO_CONST target_ptr;
  CCTK_INT ind, ierr;
  char *thorn, *parameter;
  CCTK_INT target_type;
  size_t len_copied;

  // skip whitespace to end or first string
  while(*targetstring != '\0' &&  isspace(*targetstring))
    targetstring++;
  while(targetstring[0] != '\0')
  {
    // extract a single target string
    len_copied = strncpyw(name_buffer, targetstring, sizeof(name_buffer));
    assert(len_copied <= sizeof(name_buffer));

    // try to find out what type the target is
    ind = CCTK_VarIndex(name_buffer);
    if(ind >= 0)
    {
      target_ptr = CCTK_VarDataPtrI(cctkGH,0,ind);
      assert(target_ptr != NULL);

      target_type = grid_scalar_type;
    }
    else
    {
      ierr = CCTK_DecomposeName(name_buffer, &thorn, &parameter);
      assert(ierr == 0);
      target_ptr = CCTK_ParameterData(parameter, thorn);
      free(thorn);
      free(parameter);

      if(target_ptr != NULL)
      {
        target_type = parameter_type;
      }
      else
      {
        CCTK_VWarn(1, __LINE__,__FILE__,CCTK_THORNSTRING, "could not find out target type for '%s'", name_buffer);
        return -1;
      }
    }

    // call callback (any type will do since they are all identical)
    ierr = callback(target_ptr, target_type, name_buffer, callback_arg);
    if(ierr != 0)
      return ierr;

    // advance to possible next string
    while(*targetstring != '\0' && !isspace(*targetstring)) // skip to beyond current target name
      targetstring++;
    while(*targetstring != '\0' &&  isspace(*targetstring)) // skip whitespace to end or next string
      targetstring++;
  }

  return 0;
}

//////////////////////////////////////////////////////////////////////
// Scheduled functions
//////////////////////////////////////////////////////////////////////

extern "C" void MinSearch_SteerPosition(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if (verbose >= 2)
    CCTK_INFO("MinSearch_SteerPosition");

  for(int i = 0 ; i < nminima ; i++)
  {
    CCTK_REAL min_loc[3] = {min_x[i], min_y[i], min_z[i]};
    CCTK_STRING targets[3] = {steer_target_x[i], steer_target_y[i], steer_target_z[i]}; 

    for(int j = 0 ; j < 3 ; j++)
    {
      MinSearch_TraverseTargetString(cctkGH, targets[j], steer_all_targets, &min_loc[j]);
      if(verbose >= 3)
      {
        CCTK_VInfo(CCTK_THORNSTRING, "Steering target '%s' for minimum %d to %g", 
                   targets[j], i, min_loc[j]);
      }
    }
  }
}

// steer target
static CCTK_INT steer_all_targets(CCTK_POINTER_TO_CONST target, const CCTK_INT target_type, CCTK_STRING /*targetname*/, CCTK_POINTER callback_arg)
{
  cParamData *paramdata;
  CCTK_REAL *vardataptr;
  CCTK_REAL min_loc;
  char buffer[32];
  size_t written;
  CCTK_INT ierr;

  min_loc = *(CCTK_REAL*)callback_arg;
  assert(target_type >= 0);

  switch(target_type)
  {
    case parameter_type:
      written = snprintf(buffer, sizeof(buffer), "%.18e", min_loc);
      assert(written < sizeof(buffer));

      paramdata = (cParamData*)target;
      ierr = CCTK_ParameterSet(paramdata->name, paramdata->thorn, buffer);
      break;
    case grid_scalar_type:
      vardataptr = (CCTK_REAL *)target;

      if(vardataptr != NULL)
      {
        *vardataptr = min_loc;
        ierr = 0;
      }
      else
        ierr = -1;
      break;
    case invalid_type: 
      // fall through
    default:
      assert(0 && "internal error");
      break;
  }

  return ierr;
}

// copy at most n characters or until end-of-string or whitespace is encountered
static size_t strncpyw(char *dst, const char *src, size_t n)
{
  char *t;

  // we differ from strncpy in that we always terminate the string (>1 s. >0)
  t = dst;
  while(src[0] != '\0' && !isspace(src[0]) && n > 1)
  {
    *(t++) = *(src++);
    n--;
  }
  *t = '\0'; //terminate string
  
  return (size_t)(t-dst+1);
}
