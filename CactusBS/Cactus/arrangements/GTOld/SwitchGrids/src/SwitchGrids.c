
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#define SQR(x) ((x) * (x))
#define NUMBER_OF_SWITCHES 6
#define DIM(x) (sizeof(x)/sizeof((x)[0]))

/* steer a set of cactus parameters given by name=value pairs in parstring
 * (separated by '\n') */
static void steer_parameters(CCTK_STRING parstring)
{
  DECLARE_CCTK_PARAMETERS;

  char linebuf[1024];   /* holds a copy of the currently processing line, for error output */
  char thornname[256];  /* thorn who owns that parameter */
  char paramname[256];  /* parameter to be set */
  char paramvalue[256]; /* value to be set */
  CCTK_INT ierr;

  do
  {
    /* extract one line of input, trimming whitespace at the front */
    if(sscanf(parstring, " %[^\n]", linebuf) != 1) /* given the pattern, failure means whitespace only line */
      continue;
    /* trim trailing whitespace */
    for(char *p = linebuf + strlen(linebuf) - 1 ; p >= linebuf && isspace(*p) ; --p)
      *p = '\0';

    /* skip commented out and empty lines */
    if(linebuf[0] == '#' || linebuf[0] == '\0' || linebuf[0] == '\n')
      continue;

    /* parse thorn::parameter = value */
    if(sscanf(linebuf, " %[a-zA-Z_]::%[][a-zA-Z_0-9] = %[^\n]", thornname, paramname, paramvalue) != 3)
    {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Malformed parameter assignment in '%s'", linebuf);
      continue;
    }

    /* steer parameter */
    if(verbose >= 3)
    {
      CCTK_VInfo (CCTK_THORNSTRING,
                  "Setting %s::%s to '%s'.", thornname, paramname, paramvalue);
    }
    ierr = CCTK_ParameterSet(paramname, thornname, paramvalue);
    if(ierr)
    {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Could not set %s::%s to '%s': %d",
                 thornname, paramname, paramvalue, ierr);
    }
  } while((parstring = strchr(parstring+1, '\n')) != NULL);
}

/* steer a of cactus grid scalar. Emulates the CCTK_ParameterSet ingerface */
static CCTK_INT variable_set(const cGH * cctkGH, CCTK_STRING variablename, CCTK_STRING variablevalue)
{
  CCTK_INT varindex;
  CCTK_POINTER *vardataptr;

  varindex = CCTK_VarIndex(variablename);
  if(varindex < 0)
    return -2; /* parameter was not found */
  
  vardataptr = CCTK_VarDataPtrI(cctkGH, 0, varindex);
  if(NULL == vardataptr)
    return -3; /* trying to steer a non-steerable parameter */

  switch(CCTK_VarTypeI(varindex))
  {
    case CCTK_VARIABLE_REAL:
      *((CCTK_REAL *)vardataptr) = (CCTK_REAL)atof(variablevalue);
      break;
    case CCTK_VARIABLE_INT:
      *((CCTK_INT *)vardataptr) = (CCTK_INT)atoi(variablevalue);
      break;
    default:
      return -1; /* parameter is out of range */
      break;
  }

  return 0;
}

/* steer a set of cactus grid scalars given by name=value pairs in parstring
 * (separated by '\n') */
static void steer_scalars(const cGH *cctkGH, CCTK_STRING varstring)
{
  DECLARE_CCTK_PARAMETERS;

  char linebuf[1024];   /* holds a copy of the currently processing line, for error output */
  char thornname[256];  /* thorn who owns that parameter */
  char variablename[256];  /* variable to be set */
  char variablevalue[256]; /* value to be set */
  CCTK_INT ierr;
  do
  {
    /* extract one line of input, trimming whitespace at the front */
    if(sscanf(varstring, " %[^\n]", linebuf) != 1) /* given the pattern, failure means whitespace only line */
      continue;
    /* trim trailing whitespace */
    for(char *p = linebuf + strlen(linebuf) - 1 ; p >= linebuf && isspace(*p) ; --p)
      *p = '\0';

    /* skip commented out and empty lines */
    if(linebuf[0] == '#' || linebuf[0] == '\0' || linebuf[0] == '\n')
      continue;

    /* parse thorn::variablename = value */
    if(sscanf(linebuf, " %[a-zA-Z_]::%[][a-zA-Z_0-9] = %[^\n]", thornname, variablename, variablevalue) != 3)
    {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Malformed variable assignment in '%s'", linebuf);
      continue;
    }

    /* steer parameter */
    if(verbose >= 3)
    {
      CCTK_VInfo (CCTK_THORNSTRING,
                  "Setting %s::%s to '%s'.", thornname, variablename, variablevalue);
    }
    ierr = variable_set(cctkGH, variablename, variablevalue);
    if(ierr)
    {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Could not set %s::%s to '%s': %d",
                 thornname, variablename, variablevalue, ierr);
    }
  } while((varstring = strchr(varstring+1, '\n')) != NULL);
}

static CCTK_INT get_sf_centroid(CCTK_ARGUMENTS, int switchnum, int sf, CCTK_REAL xb[3])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT retval = 0;

  if(sf < nsurfaces)
  {
    if(sf >= 0)
    {
      if(sf_valid[sf])
      {
        xb[0] = sf_centroid_x[sf];
        xb[1] = sf_centroid_y[sf];
        xb[2] = sf_centroid_z[sf];
        retval = 1;
      }
      else if(verbose >= 1)
      {
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING, 
                    "Surface index %d specified for switch %d is invalid.", 
                    sf, switchnum);
      }
    }
    else if(verbose >= 2)
      {
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING, 
                    "Ignoring surface index %d specified for switch %d.", 
                    sf, switchnum);
      }

  }
  else if(verbose >= 1)
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING, 
                "Surface index %d specified for switch %d out of bounds.", 
                sf, switchnum);
  }
  if(verbose > 7)
  {
    CCTK_VInfo(CCTK_THORNSTRING, "Using surface %d for switch %d with centre [%g,%g,%g].",
        sf, switchnum, xb[0], xb[1], xb[2]);
  }

  return retval;
}

void SwitchGrids_ParamCheck (CCTK_ARGUMENTS) {

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS

  /* Check that either RegridBoxes or CarpetRegrid2 are active
     and that this matches the regridding_thorn parameter       */

  if ( CCTK_EQUALS( regridding_thorn, "RegridBoxes" ) && !CCTK_IsThornActive("RegridBoxes") ) {
     CCTK_PARAMWARN ( "Parameter chooses RegridBoxes as primary grid control, but RegridBoxes is not active" );
  } 

  if ( CCTK_EQUALS( regridding_thorn, "CarpetRegrid2" ) && !CCTK_IsThornActive("CarpetRegrid2") ) {
     CCTK_PARAMWARN ( "Parameter chooses CarpetRegrid2 as primary grid control, but CarpetRegrid2 is not active" );
  }

}

void SwitchGrids_Switch(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /* only set RegridBoxes parameters once */
  static int already_switched[NUMBER_OF_SWITCHES];
  /* buffers to construct parameter in */
  char parname[] = "refinement_levels[999]";
  char scalarname[] = "CarpetRegrid2::num_levels[99]";
  CCTK_INT *scalarvardataptr;
  char levels[] = "01";
  int ierr, ilen;
  /* coordiantes of tracker and centre for distance calculations */
  CCTK_REAL xb[3], xa[3];
  /* boolean to transfer information out of the loop */
  int do_switch;

  if(debug_fake_motion)
  {
    sf_valid[0] = 1;
    sf_active[0] = 1;
    sf_centroid_x[0] = 20.0-cctk_time;
    sf_centroid_y[0] = 20-1.4*cctk_time;
    sf_centroid_z[0] = 0;
  }

  for(int i = 0 ; i < number_of_switches ; i++)
  {
    do_switch = 0;

    if(centre_index[i] >= 0)
    {
      /* steer centre based on criterion */
      if(CCTK_Equals(switch_criterion[i], "spherical surface"))
      {
        if(get_sf_centroid(CCTK_PASS_CTOC, i, surface_index_b[i], xb)
            && get_sf_centroid(CCTK_PASS_CTOC, i, surface_index_a[i], xa))
        {
          if(SQR(xb[0]-xa[0]) + SQR(xb[1]-xa[1]) + SQR(xb[2]-xa[2]) 
              < SQR(switch_distance[i]))
          {
            do_switch = 1;
          }
        }
      }
      else if(CCTK_Equals(switch_criterion[i], "coordinate"))
      {
        if(get_sf_centroid(CCTK_PASS_CTOC, i, surface_index_a[i], xa))
        {
          if(verbose > 6)
          {
            CCTK_VInfo(CCTK_THORNSTRING, "Testing distance for switch %d: %g.",
                i, sqrt(SQR(coord_x[i]-xa[0]) + SQR(coord_y[i]-xa[1]) 
                              + SQR(coord_z[i]-xa[2])));
          }
          if(SQR(coord_x[i]-xa[0]) + SQR(coord_y[i]-xa[1]) 
              + SQR(coord_z[i]-xa[2]) < SQR(switch_distance[i]))
          {
            do_switch = 1;
          }
        }
      }
      else if(CCTK_Equals(switch_criterion[i], "iteration"))
      {
        if(switch_iteration[i] >= 0 && cctk_iteration >= switch_iteration[i])
        {
          do_switch = 1;
        }
      }
    }
    else
    {
      /* do nothing */
      if(verbose >= 2)
      {
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Ignoring regridboxes centre %d for switch %d.",
                    centre_index[i], i);
      }
    }
       
    if(do_switch && !already_switched[i])
    {
      if(CCTK_Equals(regridding_thorn, "RegridBoxes"))
      {
        if(verbose >= 3)
        {
          CCTK_VInfo (CCTK_THORNSTRING, 
                      "Setting RegridBoxes::refinement_levels[%d] = %d for switch %d.", 
                      centre_index[i], centre_levels[i], i);
        }
        ilen = snprintf(levels, DIM(levels), "%d", centre_levels[i]);
        assert(ilen < DIM(levels));
        ilen = snprintf(parname, DIM(parname), "refinement_levels[%d]", 
                        centre_index[i]);
        assert(ilen < DIM(parname));
        ierr = CCTK_ParameterSet(parname, "RegridBoxes", levels);
        assert(0 == ierr);
      }
      else if(CCTK_Equals(regridding_thorn, "CarpetRegrid2"))
      {
        if(verbose >= 3)
        {
          CCTK_VInfo (CCTK_THORNSTRING, 
                      "Setting CarpetRegrid2::num_levels[%d] = %d for switch %d.", 
                      centre_index[i], centre_levels[i], i);
        }
        ilen = snprintf(scalarname, DIM(scalarname), "num_levels[%d]", 
                        centre_index[i]);
        assert(ilen < DIM(scalarname));
        scalarvardataptr = CCTK_VarDataPtr(cctkGH, 0, scalarname);
        assert(scalarvardataptr);
        *scalarvardataptr = centre_levels[i];
      }
      else
      {
        CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Unknown regridding thorn '%s'", regridding_thorn);
      }

      steer_parameters(steered_parameters[i]);
      steer_scalars(cctkGH, steered_scalars[i]);

      already_switched[i] = 1;
    }
  }
}
