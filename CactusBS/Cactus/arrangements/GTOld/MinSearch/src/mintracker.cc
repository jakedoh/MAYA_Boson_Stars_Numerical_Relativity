
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cctk.h" 
#include "cctk_Arguments.h" 
#include "cctk_Parameters.h" 
#include "cctk_Functions.h"
#include "util_Table.h"
#include "mpi.h"
#include "carpet.hh"
#include "minsearch.h"

//////////////////////////////////////////////////////////////////////
// Global variables
//////////////////////////////////////////////////////////////////////

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define MAX_NUMBER_OF_EXTRA_VARIABLES 30

class SimpleOutputFile;
static const int max_minima = 6;
static SimpleOutputFile *output_files[max_minima];
static CCTK_INT extra_variable_indices[MAX_NUMBER_OF_EXTRA_VARIABLES];
static CCTK_INT number_of_extra_variables;

//////////////////////////////////////////////////////////////////////
// SimpleOutputFile: A class to handle straightforward thorn ascii
// output
//////////////////////////////////////////////////////////////////////

// all inlined so that it can be compiled together with MinTracker (roland)
class SimpleOutputFile
{
public:
SimpleOutputFile(const char *name)
{
  file_written = false;

  CCTK_STRING *out_dir  = (CCTK_STRING *) 
    CCTK_ParameterGet("out_dir", "IOUtil", NULL);

  if (*out_dir == 0)
  {
    CCTK_WARN(1,"Parameter IOUtil::out_dir not found");
    return;
  }

  filename = new char[strlen(*out_dir) + 1 + strlen(name) + 1];
  sprintf(filename, "%s/%s", *out_dir, name);
}

FILE *open()
{
  const char *mode = file_written ? "a" : "w";

  FILE *fp = fopen(filename, mode);

  if (fp == 0)
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Could not write to file \'%s\'", filename);
  }
  else
  {
    file_written = true;
  }
  return fp;
}

private:
  bool file_written;
  char *filename;
};


//////////////////////////////////////////////////////////////////////
// Interpolation functions
//////////////////////////////////////////////////////////////////////

static int interpolate_points(CCTK_ARGUMENTS, 
  const CCTK_REAL *coord_x, const CCTK_REAL *coord_y, const CCTK_REAL *coord_z,
  CCTK_REAL *var_values, CCTK_REAL *extras_values, CCTK_INT *extras_ind, CCTK_INT numextras,
  CCTK_INT npoints)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  int ierr=-1;

  const void* interp_coords[3]
    = { (const void *) coord_x,
        (const void *) coord_y,
        (const void *) coord_z };

  const CCTK_INT num_input_arrays = 1+numextras;
  const CCTK_INT num_output_arrays = 1+numextras;

  CCTK_INT input_array_indices[num_input_arrays];
  input_array_indices[0] = CCTK_VarIndex(var);
  for(int i = 0 ; i < numextras ; i++)
      input_array_indices[i+1] = extras_ind[i];
  for(int i = 0 ; i < num_input_arrays ; i++)
  {
    if (input_array_indices[i] < 0)
    {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
          "MinSearch: couldn't find variable \"%s\"!",
          CCTK_VarName(input_array_indices[i]));
      return 0; /*NOTREACHED*/
    }
  }

  CCTK_INT output_array_types[num_output_arrays];
  for(int i = 0 ; i < num_output_arrays ; i++)
    output_array_types[i] = CCTK_VARIABLE_REAL;

  void * output_arrays[num_output_arrays];
  output_arrays[0] = var_values;
  for(int i = 0 ; i < numextras ; i++)
    output_arrays[i+1] = &extras_values[i*npoints];

  const int operator_handle = CCTK_InterpHandle(interpolator_name);
  if (operator_handle < 0)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
        "MinSearch: couldn't find interpolator \"%s\"!",
        interpolator_name);         /*NOTREACHED*/
    return 0;
  }

  int param_table_handle = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT);
  ierr = Util_TableSetFromString(param_table_handle, interpolator_pars);
  if (ierr < 0)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
        "MinSearch: bad interpolator parameter(s) \"%s\"!",
        interpolator_pars);         /*NOTREACHED*/
    return 0;
  }
  
  // Don't interpolate in time.  
  ierr = Util_TableSetInt(param_table_handle, 1, "InterpNumTimelevels");
  if (ierr < 0)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
        "MinSearch: bad interpolator parameter(s) \"%d\"!",
        1);         /*NOTREACHED*/
    return 0;
  }

  CCTK_INT operand_indices[num_output_arrays];
  CCTK_INT operation_codes[num_output_arrays];
  for(int i = 0 ; i < num_output_arrays ; i++)
  {
     operand_indices[i] = 0;
     operation_codes[i] = 0;
  }

  if (Util_TableSetIntArray(param_table_handle, num_output_arrays, 
                            operand_indices, "operand_indices") < 0)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "MinSearch: cannot set operand indices in parameter table");
  }

  if (Util_TableSetIntArray(param_table_handle, num_output_arrays, 
                            operation_codes, "operation_codes") < 0)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "MinSearch: cannot set operation_codes in parameter table");
  }
  
  const int coord_system_handle = CCTK_CoordSystemHandle(coord_system);
  if (coord_system_handle < 0)
  {
    CCTK_VWarn(-1, __LINE__, __FILE__, CCTK_THORNSTRING,
        "MinSearch: can't get coordinate system handle for\n"
        "                           for coordinate system \"%s\"!",
        coord_system);
    return 0;
  }

  ierr = CCTK_InterpGridArrays(cctkGH,
                               3, // number of dimensions used, hardcoded elsewhere
                               operator_handle,
                               param_table_handle,
                               coord_system_handle,
                               npoints,
                               CCTK_VARIABLE_REAL,
                               interp_coords,
                               num_input_arrays, // Number of input arrays to be interpolated
                               input_array_indices,
                               num_output_arrays, // Number of output arrays
                               output_array_types,
                               output_arrays);
  
  if( verbose >= 2 )
    CCTK_VInfo(CCTK_THORNSTRING,
        "MinSearch_Interp(): status=%d return from interpolator",
        ierr);

  Util_TableDestroy(param_table_handle);
  
  if (ierr < 0)
    return 0;

  return 1;
}

//////////////////////////////////////////////////////////////////////
// Callback routine for paramcheck
//////////////////////////////////////////////////////////////////////
static CCTK_INT check_target(CCTK_POINTER_TO_CONST /*target*/, 
        const CCTK_INT target_type, CCTK_STRING /*targetname*/, 
        CCTK_POINTER /*callback_arg*/)
{
  return target_type == invalid_type ? -1 : 0;
}

//////////////////////////////////////////////////////////////////////
// callback routine to add one extra variable to be output
//////////////////////////////////////////////////////////////////////
static void fill_variable(int idx, const char *optstring, void *callback_arg)
{
  assert(idx >= 0);
  assert(callback_arg);

  struct calldata_t {
      CCTK_INT *extras_ind;
      CCTK_INT num_extras;
  } *calldata = (struct calldata_t *)callback_arg;

  if(calldata->num_extras >= MAX_NUMBER_OF_EXTRA_VARIABLES) { /* no more free slots */
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "too many extra variables, ignoring variable '%s'.",
               CCTK_VarName(idx));
    return;
  }

  calldata->extras_ind[calldata->num_extras++] = idx;
}


//////////////////////////////////////////////////////////////////////
// Scheduled functions
//////////////////////////////////////////////////////////////////////

extern "C" void MinSearch_ParamCheck(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  bool all_valid = true;

  if (verbose >= 2)
    CCTK_INFO("MinSearch_ParamCheck");

  if (verbose >= 1)
  {

    if(newton_gamma != -1.0)
	CCTK_VWarn(1, __LINE__,__FILE__,CCTK_THORNSTRING, "newton_gamma is not used and only available to not break MinTracker compatability");

    if(newton_tolerance != -1.0)
	CCTK_VWarn(1, __LINE__,__FILE__,CCTK_THORNSTRING, "newton_tolerance is not used and only available to not break MinTracker compatability");

  }

  // check whether targets are valid
  for(int i = 0 ; i < nminima ; i++)
  {
    CCTK_STRING targets[3] = {steer_target_x[i], steer_target_y[i], steer_target_z[i]}; 
    
    for(int j = 0 ; j < 3 ; j++)
    {
      if(MinSearch_TraverseTargetString(cctkGH, targets[j], check_target, NULL) != 0)
      {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Invalid target string steer_target_%c[%d]: '%s'", 
                   "xyz"[j],i,targets[i]);
        all_valid = false;
      }
    }
  }

  if(!all_valid)
    CCTK_PARAMWARN("invalid paramters found");
}

extern "C" void MinSearch_SetupOutputFiles(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  char name[256];
  FILE *f;
  struct calldata_t {
      CCTK_INT *extras_ind;
      CCTK_INT num_extras;
  } calldata;

  if (verbose >= 2)
    CCTK_INFO("MinSearch_SetupOutputFiles");

  for (int m = 0; m < nminima; m++)
  {
    /* parse variables string, to be able to write header */
    if(!CCTK_Equals(extra_variables[m], "")) {
      calldata.extras_ind = extra_variable_indices;
      calldata.num_extras = 0;

      const CCTK_INT ierr = 
          CCTK_TraverseString(extra_variables[m], fill_variable, &calldata,
              CCTK_GROUP_OR_VAR);
      assert(ierr > 0);
    
      number_of_extra_variables = calldata.num_extras;
    } else {
      number_of_extra_variables = 0;
    }

    /* open output files and write header */
    if (snprintf(name, sizeof(name), outfile_basename, m) >= (int)sizeof(name))
    {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Filename for MinSearch %d is too long", m); 
    }
    output_files[m] = new SimpleOutputFile(name);

    if (CCTK_MyProc(cctkGH) == 0 && (f = output_files[m]->open()))
    {
      fprintf(f, "# 1:it 2:t 3:x 4:y 5:z 6:min");
      for(int i = 0 ; i < number_of_extra_variables ; i++)
      {
          fprintf(f, " %d:%s", i+7, CCTK_VarName(extra_variable_indices[m*MAX_NUMBER_OF_EXTRA_VARIABLES+i]));
      }
      fputc('\n',f);
      fclose(f);
    }
  }
}

extern "C" void MinSearch_Init(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if (verbose >= 2)
    CCTK_INFO("MinSearch_Init");

  for (int m = 0; m < nminima; m++)
  {
    min_x_new[m] = min_x[m] = x0[m];
    min_y_new[m] = min_y[m] = y0[m];
    min_z_new[m] = min_z[m] = z0[m];
    min_value_new[m] = track_maximum ? -HUGE_VAL : HUGE_VAL;
    
    int sind = surface_index[m];
    if (sind != -1) 
    {
      assert(sind >=0 && sind < nsurfaces);
      sf_valid[sind]=1;
      sf_active[sind]=1;
      sf_centroid_x[sind]=min_x[m];
      sf_centroid_y[sind]=min_y[m];
      sf_centroid_z[sind]=min_z[m];
      if (verbose >= 2)
      {
        CCTK_VInfo(CCTK_THORNSTRING, 
                   "Initializing spherical surface %d\n", sind);
      }
    }
  }
}

extern "C" void MinSearch_Output(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if (cctk_iteration % find_every != 0) // always triggers during INITIALIZE
    return;

  if (verbose >= 2)
    CCTK_INFO("MinSearch_Output");

  for (int m = 0; m < nminima; m++)
  {
    /* output dat to file */
    FILE *fp = 0;
    if (CCTK_MyProc(cctkGH) == 0 && (fp = output_files[m]->open()))
    {
      fprintf(fp, "%d\t%g\t%.19g\t%.19g\t%.19g\t%.19g", cctk_iteration, 
              cctk_time, min_x[m], min_y[m], min_z[m], min_value[m]);
      for(int i = 0 ; i < number_of_extra_variables ; i++)
      {
          fprintf(fp, "\t%.19g", extra_variable_value[m*MAX_NUMBER_OF_EXTRA_VARIABLES+i]);
      }
      fputc('\n',fp);
      fclose(fp);
    }
  }
}

extern "C" void MinSearch_SetSphericalSurface(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if (cctk_iteration % find_every != 0) // always triggers during INITIALIZE
    return;

  if (verbose >= 2)
    CCTK_INFO("MinSearch_SetSphericalSurface");

  for (int m = 0; m < nminima; m++)
  {
    CCTK_REAL dist_moved;

    dist_moved = sqrt(SQR(min_x_new[m]-min_x[m])+SQR(min_y_new[m]-min_y[m])+SQR(min_z_new[m]-min_z[m]));
    if (verbose >= 1 || dist_moved > dist_cutoff )
    {
      CCTK_VInfo(CCTK_THORNSTRING, "min[%d] = [%g,%g,%g] moved distance %g by [%g,%g,%g]", m,
	      min_x_new[m],min_y_new[m],min_z_new[m], dist_moved,
	      min_x_new[m]-min_x[m],min_y_new[m]-min_y[m],min_z_new[m]-min_z[m]);
    }

    /* store final results in grid scalars */
    min_x[m] = min_x_new[m];
    min_y[m] = min_y_new[m];
    min_z[m] = min_z_new[m];
    min_value[m] = min_value_new[m];
    min_value_new[m] = track_maximum ? -HUGE_VAL : HUGE_VAL; // reset for next iteration
    for(int i = 0 ; i < number_of_extra_variables ; i++)
    {
      extra_variable_value[m*MAX_NUMBER_OF_EXTRA_VARIABLES+i] = extra_variable_value_new[m*MAX_NUMBER_OF_EXTRA_VARIABLES+i];
    }

    int sind = surface_index[m];
    if (sind != -1) 
    {
      assert(sind >=0 && sind < nsurfaces);
      sf_valid[sind]=1;
      sf_active[sind]=1;
      sf_centroid_x[sind]=min_x[m];
      sf_centroid_y[sind]=min_y[m];
      sf_centroid_z[sind]=min_z[m];
      if (verbose >= 2)
      {
        CCTK_VInfo(CCTK_THORNSTRING, "Setting spherical surface %d\n", sind);
      }
    }
  }
}

extern "C" void MinSearch_InterpolateMinimum(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  struct calldata_t {
      CCTK_INT *extras_ind;
      CCTK_INT num_extras;
  } calldata;

  if (cctk_iteration % find_every != 0)
    return;

  if (verbose >= 2)
    CCTK_INFO("MinSearch_InterpolateMinimum");

  for (int m = 0; m < nminima; m++)
  {
    CCTK_REAL wmin[3]    = {min_x[m],min_y[m],min_z[m]};
    CCTK_INT ngx = cctk_nghostzones[0], 
	     ngy = cctk_nghostzones[1], ngz = cctk_nghostzones[2];
    CCTK_REAL var_values[2*ngx+1][2*ngy+1][2*ngz+1];
    CCTK_REAL extras_values[MAX_NUMBER_OF_EXTRA_VARIABLES][2*ngx+1][2*ngy+1][2*ngz+1];
    CCTK_REAL coord_x[2*ngx+1][2*ngy+1][2*ngz+1];
    CCTK_REAL coord_y[2*ngx+1][2*ngy+1][2*ngz+1];
    CCTK_REAL coord_z[2*ngx+1][2*ngy+1][2*ngz+1];

    if (verbose >= 1)
    {
      CCTK_VInfo(CCTK_THORNSTRING,
                 "Searching for extremum %d",
                 m);
    }

    /* parse variables string */
    if(!CCTK_Equals(extra_variables[m], "")) {
      calldata.extras_ind = extra_variable_indices;
      calldata.num_extras = 0;

      const CCTK_INT ierr = 
          CCTK_TraverseString(extra_variables[m], fill_variable, &calldata,
              CCTK_GROUP_OR_VAR);
      assert(ierr > 0);
    
      number_of_extra_variables = calldata.num_extras;
    } else {
      number_of_extra_variables = 0;
    }

    // set out grid around the old minimum
    int ii,jj,kk;
    // This will only work for uniform coordinates, but so will Carpet.
    // from Cactus docs:
    // x = CCTK_ORIGIN_SPACE(0) + i*CCTK_DELTA_SPACE(0) [ignoring cctk_lbnd]
    ii = ((int) floor(( wmin[0] - CCTK_ORIGIN_SPACE(0) ) / CCTK_DELTA_SPACE(0)));
    jj = ((int) floor(( wmin[1] - CCTK_ORIGIN_SPACE(1) ) / CCTK_DELTA_SPACE(1)));
    kk = ((int) floor(( wmin[2] - CCTK_ORIGIN_SPACE(2) ) / CCTK_DELTA_SPACE(2)));

    if(ii < ngx || ii >= cctk_gsh[0]-ngx || jj < ngy || jj >= cctk_gsh[1]-ngy || kk < ngz || kk >= cctk_gsh[2]-ngz)
    {
        // might wander outside of current component => skip it
        if (verbose >= 1)
        {
          CCTK_VInfo(CCTK_THORNSTRING,
                     "Minimum (%g,%g,%g) is too close to the edge of the component => skipping",
                     wmin[0],wmin[1],wmin[2]);
        }

        return;
    }
	
    for(int i = -ngx ; i <= ngx ; i++)
    {
      for(int j = -ngy ; j <= ngy ; j++)
      {
        for(int k = -ngz ; k <= ngz ; k++)
	{
	  coord_x[ngx+i][ngy+j][ngz+k] = CCTK_ORIGIN_SPACE(0) + (ii + i)*CCTK_DELTA_SPACE(0); 
	  coord_y[ngx+i][ngy+j][ngz+k] = CCTK_ORIGIN_SPACE(1) + (jj + j)*CCTK_DELTA_SPACE(1); 
	  coord_z[ngx+i][ngy+j][ngz+k] = CCTK_ORIGIN_SPACE(2) + (kk + k)*CCTK_DELTA_SPACE(2); 
	}
      }
    }

    if (verbose >= 2)
    {
      CCTK_VInfo(CCTK_THORNSTRING,
		 "Search volume: %d points in [[%g,%g,%g],[%g,%g,%g]]",
		 (2*ngx+1)*(2*ngy+1)*(2*ngz+1), 
		 coord_x[0][0][0], coord_y[0][0][0], coord_y[0][0][0],
		 coord_x[2*ngx][2*ngy][2*ngz], coord_y[2*ngx][2*ngy][2*ngz], coord_y[2*ngx][2*ngy][2*ngz]);
    }

    // fill array with variable data of the points surrounding the current minimum
    if (!interpolate_points(CCTK_PASS_CTOC, &coord_x[0][0][0],
	  &coord_y[0][0][0], &coord_z[0][0][0], &var_values[0][0][0],
          &extras_values[0][0][0][0], extra_variable_indices, number_of_extra_variables,
          (2*ngx+1)*(2*ngy+1)*(2*ngz+1)))
    {
       exit(1);
    }

    // search through returned values for minimum (HUGE_VAL is max(double))
    CCTK_REAL min = track_maximum ? -HUGE_VAL : HUGE_VAL; 

    int imin=42,jmin=42,kmin=42; // initialize to something distinctive
    for(int i = -ngx ; i <= ngx ; i++)
    {
      for(int j = -ngy ; j <= ngy ; j++)
      {
        for(int k = -ngz ; k <= ngz ; k++)
	{
	  if (verbose >= 4)
	  {
	    CCTK_VInfo(CCTK_THORNSTRING,
	      	 "Searching for extremum at [%d,%d,%d]: %g",
	      	 i,j,k,var_values[ngx+i][ngy+j][ngz+k]);
	  }
	  if(track_maximum ? (var_values[ngx+i][ngy+j][ngz+k] > min) : (var_values[ngx+i][ngy+j][ngz+k] < min))
	  {
	    min = var_values[ngx+i][ngy+j][ngz+k];
	    imin=i;jmin=j;kmin=k;
	  }
	}
      }
    }

    if (verbose >= 1)
    {
	if(imin == -ngx || imin == ngx || jmin == -ngy || jmin == ngy 
	  || kmin == -ngz || kmin == ngz)
	    CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,  
		    "Extremum found on egde of search domain of size %d x %d x %d at face [%d,%d,%d]",
		    ngx,ngy,ngz, imin/ngx,jmin/ngy,kmin/ngz);
    }

    // we cannot set the spherical surface yet since Carpet will drive us from
    // coarse to fine levels, so we wait until we have seen all of them and use
    // the finest value
    
    // calculate coordinate of new minimum
    min_x_new[m] = CCTK_ORIGIN_SPACE(0) + (ii + imin)*CCTK_DELTA_SPACE(0);
    min_y_new[m] = CCTK_ORIGIN_SPACE(1) + (jj + jmin)*CCTK_DELTA_SPACE(1);
    min_z_new[m] = CCTK_ORIGIN_SPACE(2) + (kk + kmin)*CCTK_DELTA_SPACE(2);
    min_value_new[m] = min;
    for(int i = 0 ; i < number_of_extra_variables ; i++)
    {
      extra_variable_value_new[m*MAX_NUMBER_OF_EXTRA_VARIABLES+i] = extras_values[i][ngx+imin][ngy+jmin][ngz+kmin];
    }

    if (verbose >= 2)
    {
      CCTK_VInfo(CCTK_THORNSTRING, "Found extremum at [%g,%g,%g]", 
                 min_x_new[m],min_y_new[m], min_z_new[m]);
    }
    
  }

}

extern "C" void MinSearch_FindGlobalMinimum(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  CCTK_REAL min;      // current best guess for minimum
  CCTK_REAL *varptr;  // points to data of variable 'var'
  CCTK_REAL *weight;
  int ind;            // linear index
  CCTK_REAL *extras_varptr[MAX_NUMBER_OF_EXTRA_VARIABLES];
  struct calldata_t {
      CCTK_INT *extras_ind;
      CCTK_INT num_extras;
  } calldata;

  // very simple minded implementation, will only work for a single minimum or maximum!
  if(nminima == 0)
  {
    return;
  }
  if(nminima > 1)
  {
    CCTK_WARN(0, "Global minimum search currently only supports a single minimum");
    return; // NOTREACHED
  }

  if (cctk_iteration % find_every != 0) // always triggers during INITIALIZE
    return;

  /* parse variables string and allocate temporary memory for them */
  if(!CCTK_Equals(extra_variables[0], "")) {
    calldata.extras_ind = extra_variable_indices;
    calldata.num_extras = 0;

    const CCTK_INT ierr = 
        CCTK_TraverseString(extra_variables[0], fill_variable, &calldata,
            CCTK_GROUP_OR_VAR);
    assert(ierr > 0);
  
    number_of_extra_variables = calldata.num_extras;
  } else {
    number_of_extra_variables = 0;
  }

  // search through returned values for minimum 
  varptr = (CCTK_REAL *)CCTK_VarDataPtr(cctkGH, 0, var); // only look at current (0) timelevel
  if (NULL == varptr)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
        "MinSearch: couldn't find variable \"%s\"!",
        var);
    return; /*NOTREACHED*/
  }
  for(int i = 0 ; i < number_of_extra_variables ; i++)
  {
    extras_varptr[i] = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH, 0, extra_variable_indices[i]); // only look at current (0) timelevel
    if (NULL == extras_varptr[i])
    {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
          "MinSearch: couldn't find variable \"%s\"!",
          CCTK_VarName(i)); 
      return; /*NOTREACHED*/
    }
  }

  if(CCTK_Equals(weight_variable[0], "none"))
  {
    weight = NULL;
  }
  else
  {
    weight = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH,0,weight_variable[0]));
    if (NULL == weight)
    {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
          "MinSearch: couldn't find variable \"%s\"!",
          weight_variable[0]); 
      return; /*NOTREACHED*/
    }
    assert(weight != NULL);
  }
    

  min = min_value_new[0]; // we only support one global minimum so far
  for(int k = cctk_nghostzones[2] ; k < cctk_lsh[2]-cctk_nghostzones[2] ; k++)
  {
    for(int j = cctk_nghostzones[1] ; j < cctk_lsh[1]-cctk_nghostzones[1] ; j++)
    {
      for(int i = cctk_nghostzones[0] ; i < cctk_lsh[0]-cctk_nghostzones[0] ; i++)
      {
        ind = CCTK_GFINDEX3D(cctkGH,i,j,k);
        if ((weight == NULL || weight[ind] > 0.0) && 
            (track_maximum ? varptr[ind] > min : varptr[ind] < min))
        {
          min_value_new[0] = min = varptr[ind];
          min_x_new[0] = x[ind]; //CCTK_ORIGIN_SPACE(0) + (cctk_lbnd[0] + i)*CCTK_DELTA_SPACE(0);
          min_y_new[0] = y[ind]; //CCTK_ORIGIN_SPACE(1) + (cctk_lbnd[1] + j)*CCTK_DELTA_SPACE(1);
          min_z_new[0] = z[ind]; //CCTK_ORIGIN_SPACE(2) + (cctk_lbnd[2] + k)*CCTK_DELTA_SPACE(2);
          for(int e = 0 ; e < number_of_extra_variables ; e++)
          {
            extra_variable_value_new[e] = extras_varptr[e][ind];
          }

          if (verbose >= 4)
          {
            CCTK_VInfo(CCTK_THORNSTRING,
                     "%simum candidate in variable %s [%g,%g,%g]: %g",
                         track_maximum ? "Max" : "Min", var, min_x_new[0], min_y_new[0], min_z_new[0], min_value_new[0]);
          }
        }
      }
    }
  }

  if (verbose >= 3)
  {
    CCTK_VInfo(CCTK_THORNSTRING,
             "Processor local %simum in variable %s [%g,%g,%g]: %g",
		 track_maximum ? "max" : "min", var, min_x_new[0], min_y_new[0], min_z_new[0], min_value_new[0]);
  }

  // more than one minimum is reasonable complicated:
  // * need to make sure that we are not fooled by using the same minimum twice
  //   on coarse and fine levels (use weigth)
  // * need to communicate the found minima across processors
  // * need to keep a stack of found candidates

}

// reduce per-processor results to global result
extern "C" void MinSearch_ReduceGlobalMinimum(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  int rank, nprocs;   // who am I, who else is there?
  MPI_Comm world;
  // four doubles per proc: (minval,x,y,z)
  const size_t sz_myvals = 4+number_of_extra_variables;
  double myvals[sz_myvals], min;
  double *values = NULL; 

  // very simple minded implementation, will only work for a single minimum or maximum!
  if(nminima == 0)
  {
    return;
  }
  if(nminima > 1)
  {
    CCTK_WARN(0, "Global minimum search currently only supports a single minimum");
    return; // NOTREACHED
  }

  if (cctk_iteration % find_every != 0) // always triggers during INITIALIZE
    return;

  // get MPI rank and nprocs
  world = Carpet::CarpetMPIComm();
  MPI_Comm_size(world, &nprocs);
  MPI_Comm_rank(world, &rank);
  // get storage for candidate minima
  values = new double[sz_myvals*nprocs];

  myvals[0] = min_value_new[0];
  myvals[1] = min_x_new[0];
  myvals[2] = min_y_new[0];
  myvals[3] = min_z_new[0];
  for(int i = 0 ; i < number_of_extra_variables ; i++)
  {
    myvals[4+i] = extra_variable_value_new[i];
  }

  // very well, now find the minimum over all processors
  // everything to proc 0
  // TODO: it would likely be better to gather to all since it saves one set of communication
  MPI_Gather(myvals,sz_myvals,MPI_DOUBLE, &values[0],sz_myvals,MPI_DOUBLE, 0, world);
  if(0 == rank)
  {
    // find global minimum
    min = track_maximum ? -HUGE_VAL : HUGE_VAL; 
    for(int n = 0 ; n < nprocs ; n++)
    {
      if (track_maximum ? values[n*sz_myvals+0] > min : values[n*sz_myvals+0] < min)
      {
        myvals[0] = min = values[n*sz_myvals+0];
        myvals[1] = values[n*sz_myvals+1]; 
        myvals[2] = values[n*sz_myvals+2];
        myvals[3] = values[n*sz_myvals+3];
        for(int i = 0 ; i < number_of_extra_variables ; i++)
        {
          myvals[4+i] = values[n*sz_myvals+4+i];
        }
      }
    }
  }
  // let others know of our wisdom
  MPI_Bcast(myvals, sz_myvals, MPI_DOUBLE, 0, world);
  if (verbose >= 2)
  {
    CCTK_VInfo(CCTK_THORNSTRING,
             "Global %simum in variable %s [%g,%g,%g]: %g",
		 track_maximum ? "max" : "min", var, myvals[1], myvals[2], myvals[3], myvals[0]);
  }
  min_value_new[0] = myvals[0];
  min_x_new[0] = myvals[1];
  min_y_new[0] = myvals[2];
  min_z_new[0] = myvals[3];
  for(int i = 0 ; i < number_of_extra_variables ; i++)
  {
    extra_variable_value_new[i] = myvals[4+i];
  }

  delete[] values;
}

extern "C" void MinSearch_FindLocalizedMinimum(CCTK_ARGUMENTS)
{
	DECLARE_CCTK_ARGUMENTS
	DECLARE_CCTK_PARAMETERS

	if (cctk_iteration % find_every != 0) // always triggers during INITIALIZE
		return;
	for(int idx = 0; idx < nminima; idx++)
	{
		CCTK_REAL min;      // current best guess for minimum
		CCTK_REAL *varptr;  // points to data of variable 'var'
		CCTK_REAL *weight;
		int ind;            // linear index
		CCTK_REAL *extras_varptr[MAX_NUMBER_OF_EXTRA_VARIABLES];
		struct calldata_t {
			CCTK_INT *extras_ind;
			CCTK_INT num_extras;
		} calldata;


		if(!CCTK_Equals(extra_variables[idx], "")) {
			calldata.extras_ind = extra_variable_indices;
			calldata.num_extras = 0;

			const CCTK_INT ierr = 
				CCTK_TraverseString(extra_variables[idx], fill_variable, &calldata,
						CCTK_GROUP_OR_VAR);
			assert(ierr > 0);

			number_of_extra_variables = calldata.num_extras;
		} else {
			number_of_extra_variables = 0;
		}

		// search through returned values for minimum 
		varptr = (CCTK_REAL *)CCTK_VarDataPtr(cctkGH, 0, var); // only look at current (0) timelevel
		if (NULL == varptr)
		{
			CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
					"MinSearch: couldn't find variable \"%s\"!",
					var);
			return; /*NOTREACHED*/
		}
		for(int i = 0 ; i < number_of_extra_variables ; i++)
		{
			extras_varptr[i] = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH, 0, extra_variable_indices[i]); // only look at current (0) timelevel
			if (NULL == extras_varptr[i])
			{
				CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
						"MinSearch: couldn't find variable \"%s\"!",
						CCTK_VarName(i)); 
				return; /*NOTREACHED*/
			}
		}
		if(CCTK_Equals(weight_variable[idx], "none"))
		{
			weight = NULL;
		}
		else
		{
			weight = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH,0,weight_variable[idx]));
			if (NULL == weight)
			{
				CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
						"MinSearch: couldn't find variable \"%s\"!",
						weight_variable[idx]); 
				return; /*NOTREACHED*/
			}
			assert(weight != NULL);
		}

		min = min_value_new[idx]; // we only support one global minimum so far
		for(int k = cctk_nghostzones[2] ; k < cctk_lsh[2]-cctk_nghostzones[2] ; k++)
		{
			for(int j = cctk_nghostzones[1] ; j < cctk_lsh[1]-cctk_nghostzones[1] ; j++)
			{
				for(int i = cctk_nghostzones[0] ; i < cctk_lsh[0]-cctk_nghostzones[0] ; i++)
				{
					ind = CCTK_GFINDEX3D(cctkGH,i,j,k);
					CCTK_REAL delx = x[ind] - min_x[idx];
					CCTK_REAL dely = y[ind] - min_y[idx];
					CCTK_REAL delz = z[ind] - min_z[idx];
					CCTK_REAL dist = sqrt(delx*delx + dely*dely + delz*delz);
					if ((weight == NULL || weight[ind] > 0.0) && 
							(track_maximum ? varptr[ind] > min : varptr[ind] < min) && (dist < dist_cutoff))
					{
						//consider putting the loop over minima in here, rather than wrapping this loop over gridpoints
						min_value_new[idx] = min = varptr[ind];
						min_x_new[idx] = x[ind]; //CCTK_ORIGIN_SPACE(0) + (cctk_lbnd[0] + i)*CCTK_DELTA_SPACE(0);
						min_y_new[idx] = y[ind]; //CCTK_ORIGIN_SPACE(1) + (cctk_lbnd[1] + j)*CCTK_DELTA_SPACE(1);
						min_z_new[idx] = z[ind]; //CCTK_ORIGIN_SPACE(2) + (cctk_lbnd[2] + k)*CCTK_DELTA_SPACE(2);
						for(int e = 0 ; e < number_of_extra_variables ; e++)
						{
							extra_variable_value_new[5*idx+e] = extras_varptr[e][ind];
						}

						if (verbose >= 4)
						{
							CCTK_VInfo(CCTK_THORNSTRING,
									"%simum candidate in variable %s [%g,%g,%g]: %g",
									track_maximum ? "Max" : "Min", var, min_x_new[0], min_y_new[0], min_z_new[0], min_value_new[0]);
						}
					}
				}
			}
		}

		if (verbose >= 3)
		{
			CCTK_VInfo(CCTK_THORNSTRING,
					"Processor local %simum #%d in variable %s [%g,%g,%g]: %g",
					track_maximum ? "max" : "min", idx, var, min_x_new[idx], min_y_new[idx], min_z_new[idx], min_value_new[idx]);
		}
	}

}

extern "C" void MinSearch_ReduceLocalizedMinimum(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  if (cctk_iteration % find_every != 0) // always triggers during INITIALIZE
    return;

  int rank, nprocs;   // who am I, who else is there?
  MPI_Comm world;
  world = Carpet::CarpetMPIComm();
  MPI_Comm_size(world, &nprocs);
  MPI_Comm_rank(world, &rank);
  // four doubles per proc: (minval,x,y,z)
  for(int idx = 0; idx < nminima; idx++)
  { const size_t sz_myvals = 4+number_of_extra_variables;
	  double myvals[sz_myvals], min;
	  double *values = NULL; 
	  // get storage for candidate minima
	  values = new double[sz_myvals*nprocs];

	  myvals[0] = min_value_new[idx];
	  myvals[1] = min_x_new[idx];
	  myvals[2] = min_y_new[idx];
	  myvals[3] = min_z_new[idx];
  
	  for(int i = 0 ; i < number_of_extra_variables ; i++)
	  {
		  myvals[4+i] = extra_variable_value_new[5*idx+i];
	  }
	  MPI_Gather(myvals,sz_myvals,MPI_DOUBLE, &values[0],sz_myvals,MPI_DOUBLE, 0, world);
	  if(0 == rank)
	  {
		  // find global minimum
		  min = track_maximum ? -HUGE_VAL : HUGE_VAL; 
		  for(int n = 0 ; n < nprocs ; n++)
		  {
			  if (track_maximum ? values[n*sz_myvals+0] > min : values[n*sz_myvals+0] < min)
			  {
				  myvals[0] = min = values[n*sz_myvals+0];
				  myvals[1] = values[n*sz_myvals+1]; 
				  myvals[2] = values[n*sz_myvals+2];
				  myvals[3] = values[n*sz_myvals+3];
				  for(int i = 0 ; i < number_of_extra_variables ; i++)
				  {
					  myvals[4+i] = values[n*sz_myvals+4+i];
				  }
			  }
		  }
	  }
	  // let others know of our wisdom
	  MPI_Bcast(myvals, sz_myvals, MPI_DOUBLE, 0, world);
	  if (verbose >= 2)
	  {
		  CCTK_VInfo(CCTK_THORNSTRING,
				  "Localized %simum #%d in variable %s [%g,%g,%g]: %g",
				  track_maximum ? "max" : "min", idx, var, myvals[1], myvals[2], myvals[3], myvals[0]);
	  }
	  min_value_new[idx] = myvals[0];
	  min_x_new[idx] = myvals[1];
	  min_y_new[idx] = myvals[2];
	  min_z_new[idx] = myvals[3];
	  for(int i = 0 ; i < number_of_extra_variables ; i++)
	  {
		  extra_variable_value_new[5*idx+i] = myvals[4+i];
	  }
	  delete[] values;
  }
}
