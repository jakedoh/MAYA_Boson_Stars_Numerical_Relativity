
#include <assert.h>
#include <cmath>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "cctk.h" 
#include "cctk_Arguments.h" 
#include "cctk_Parameters.h" 
#include "cctk_Functions.h"

#include "util_Table.h"

//////////////////////////////////////////////////////////////////////
// Global variables
//////////////////////////////////////////////////////////////////////

class SimpleOutputFile;
static const int max_minima = 6;
static SimpleOutputFile *output_files[max_minima];

//////////////////////////////////////////////////////////////////////
// SimpleOutputFile: A class to handle straightforward thorn ascii
// output
//////////////////////////////////////////////////////////////////////

class SimpleOutputFile
{
public:
  SimpleOutputFile(const char *name);
  FILE *open();

private:
  bool file_written;
  char *filename;
};

SimpleOutputFile::SimpleOutputFile(const char *name)
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

FILE *SimpleOutputFile::open()
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

//////////////////////////////////////////////////////////////////////
// 3x3 matrix functions
//////////////////////////////////////////////////////////////////////

static inline CCTK_REAL sqr(CCTK_REAL x)
{
  return x*x;
}

static void invert_matrix_33(const CCTK_REAL m[3][3], CCTK_REAL n[3][3])
{
  CCTK_REAL detm = 
    -(m[0][2]*m[1][1]*m[2][0]) + m[0][1]*m[1][2]*m[2][0] + 
    m[0][2]*m[1][0]*m[2][1] - m[0][0]*m[1][2]*m[2][1] - 
    m[0][1]*m[1][0]*m[2][2] + m[0][0]*m[1][1]*m[2][2];

  CCTK_REAL mInvDetm[3][3] = 
    {-(m[1][2]*m[2][1]) + m[1][1]*m[2][2], m[0][2]*m[2][1] - m[0][1]*m[2][2], 
     -(m[0][2]*m[1][1]) + m[0][1]*m[1][2], m[1][2]*m[2][0] - m[1][0]*m[2][2], 
     -(m[0][2]*m[2][0]) + m[0][0]*m[2][2], m[0][2]*m[1][0] - m[0][0]*m[1][2], 
     -(m[1][1]*m[2][0]) + m[1][0]*m[2][1], m[0][1]*m[2][0] - m[0][0]*m[2][1], 
     -(m[0][1]*m[1][0]) + m[0][0]*m[1][1]};
  
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      n[i][j] = mInvDetm[i][j] / detm;
    }
  }
}

static void multiply_matrix_vector_33(const CCTK_REAL m[3][3], const CCTK_REAL v[3], 
                               CCTK_REAL w[3])
{
  for (int i = 0; i < 3; i++)
  {
    w[i] = 0;
    for (int j = 0; j < 3; j++)
    {
      w[i] += m[i][j] * v[j];
    }
  }
}

static void multiply_matrix_matrix_33(const CCTK_REAL m1[3][3], const CCTK_REAL m2[3][3], CCTK_REAL n[3][3])
{
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      n[i][j] = 0;
      for (int k = 0; k < 3; k++)
      {
        n[i][j] += m1[i][k] * m2[k][j];
      }
    }
  }
}

static void print_matrix(FILE *f, const CCTK_REAL m[3][3])
{
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      fprintf(f, "%g\t", m[i][j]);
    }
    fprintf(f, "\n");
  }
  fprintf(f, "\n");
}

static inline CCTK_REAL vector_norm_33(CCTK_REAL v[3])
{
  return sqrt(sqr(v[0]) + sqr(v[1]) + sqr(v[2]));
}

//////////////////////////////////////////////////////////////////////
// Interpolation functions
//////////////////////////////////////////////////////////////////////

static int interpolate_point(CCTK_ARGUMENTS, const CCTK_REAL p[3], 
                              CCTK_REAL *f,
                              CCTK_REAL df[3], CCTK_REAL ddf[3][3])
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  int ierr=-1;

  const void* interp_coords[3]
    = { (const void *) &p[0],
        (const void *) &p[1],
        (const void *) &p[2] };
  if (verbose >= 2)
    printf("MinTracker: Interpolating point [%g,%g,%g]\n", p[0], p[1], p[2]);

  const CCTK_INT num_input_arrays = 1;
  const CCTK_INT num_output_arrays = 10;
  const CCTK_INT input_array_indices[num_input_arrays]
    = { CCTK_VarIndex(var) };
  if (input_array_indices[0] < 0)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
        "MinTracker: couldn't find variable \"%s\"!",
        var);         /*NOTREACHED*/
    return 0;
  }

  const CCTK_INT output_array_types[num_output_arrays]
    = { CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, 
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL };
  void * output_arrays[num_output_arrays]
    = { (void *) f, 
        (void *) &df[0], (void *) &df[1], (void *) &df[2],
        (void *) &ddf[0][0], (void *) &ddf[0][1], (void *) &ddf[0][2],
        (void *) &ddf[1][1], (void *) &ddf[1][2], (void *) &ddf[2][2] };

  const int operator_handle = CCTK_InterpHandle(interpolator_name);
  if (operator_handle < 0)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
        "MinTracker: couldn't find interpolator \"%s\"!",
        interpolator_name);         /*NOTREACHED*/
    return 0;
  }

  int param_table_handle = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT);
  ierr = Util_TableSetFromString(param_table_handle, interpolator_pars);
  if (ierr < 0)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
        "MinTracker: bad interpolator parameter(s) \"%s\"!",
        interpolator_pars);         /*NOTREACHED*/
    return 0;
  }
  
  // Don't interpolate in time.  
  ierr = Util_TableSetInt(param_table_handle, 1, "InterpNumTimelevels");
  if (ierr < 0)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
        "MinTracker: bad interpolator parameter(s) \"%d\"!",
        1);         /*NOTREACHED*/
    return 0;
  }

  const CCTK_INT operand_indices[num_output_arrays] =
    { 0,0,0,0,0,0,0,0,0,0 };

  const CCTK_INT operation_codes[num_output_arrays] = 
    { 0, 1, 2, 3, 11, 12, 13, 22, 23, 33 }; // What derivative
                                            // operator to use

  if (Util_TableSetIntArray(param_table_handle, num_output_arrays, 
                            operand_indices, "operand_indices") < 0)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "MinTracker: cannot set operand indices in parameter table");
  }

  if (Util_TableSetIntArray(param_table_handle, num_output_arrays, 
                            operation_codes, "operation_codes") < 0)
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "MinTracker: cannot set operation_codes in parameter table");
  }
  
  const int coord_system_handle = CCTK_CoordSystemHandle(coord_system);
  if (coord_system_handle < 0)
  {
    CCTK_VWarn(-1, __LINE__, __FILE__, CCTK_THORNSTRING,
        "MinTracker: can't get coordinate system handle for\n"
        "                           for coordinate system \"%s\"!",
        coord_system);
    return 0;
  }

  ierr = CCTK_InterpGridArrays(cctkGH,
                               3, // number of dimensions used, hardcoded elsewhere
                               operator_handle,
                               param_table_handle,
                               coord_system_handle,
                               1, // One point to interpolate
                               CCTK_VARIABLE_REAL,
                               interp_coords,
                               num_input_arrays, // Number of input arrays to be interpolated
                               input_array_indices,
                               num_output_arrays, // Number of output arrays
                               output_array_types,
                               output_arrays);
  
  if( verbose >= 2 ) {
    {
      CCTK_VInfo(CCTK_THORNSTRING,
                 "MinTracker(): f=%g", *f);
    }
  }

  if( verbose >= 2 )
    CCTK_VInfo(CCTK_THORNSTRING,
        "MinTracker_Interp(): status=%d return from interpolator",
        ierr);

  Util_TableDestroy(param_table_handle);
  
  if (ierr < 0)
    return 0;

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < i; j++)
    {
      ddf[i][j] = ddf[j][i];
    }
  }
  return 1;
}

//////////////////////////////////////////////////////////////////////
// Scheduled functions
//////////////////////////////////////////////////////////////////////

extern "C" void MinTracker_SetupOutputFiles(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if (verbose >= 2)
    CCTK_INFO("MinTracker_SetupOutputFiles");

  for (int m = 0; m < nminima; m++)
  {
    char name[100];
    if (snprintf(name, sizeof(name), "MinTracker%d.asc", m) >= sizeof(name))
    {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Filename for MinTracker %d is too long", m); 
    }
    output_files[m] = new SimpleOutputFile(name);
  }
}

extern "C" void MinTracker_Init(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if (verbose >= 2)
    CCTK_INFO("MinTracker_Init");

  for (int m = 0; m < nminima; m++)
  {
    min_x[m] = x0[m];
    min_y[m] = y0[m];
    min_z[m] = z0[m];

    min_r[m] = sqrt(pow(min_x[m],2) + pow(min_y[m],2) + pow(min_z[m],2));

    min_vx[m] = 0;
    min_vy[m] = 0;
    min_vz[m] = 0;

    min_last_t[m] = 0;
    min_found[m] = 0;
    vel_found[m] = 0;

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

    if (FILE *f = output_files[m]->open())
    {
      fprintf(f, "# 1:it 2:t 3:x 4:y 5:z\n");
      fclose(f);
    }
  }
}

extern "C" void MinTracker_SetSphericalSurface(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if (verbose >= 2)
    CCTK_INFO("MinTracker_SetSphericalSurface");

  for (int m = 0; m < nminima; m++)
  {
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
        CCTK_VInfo(CCTK_THORNSTRING, "Initializing spherical surface %d\n", sind);
      }
    }
  }
}

extern "C" void MinTracker_InterpolateMinimum(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if (verbose >= 2)
    CCTK_INFO("MinTracker_InterpolateMinimum");

  for (int m = 0; m < nminima; m++)
  {
    CCTK_REAL wmin[3]    = {min_x[m],min_y[m],min_z[m]};
    CCTK_REAL f          = 0;
    CCTK_REAL df[3]      = {100,100,100};
    CCTK_REAL ddf[3][3]  = {0,0,0,0,0,0,0,0,0};
    int       iterations = 0;
    const int iteration_limit = 20;
    CCTK_REAL df_norm    = 0;

    if (find_every_distance != -1 && vel_found[m] && (cctk_iteration % find_every != 0))
    {
      // If the estimated distance that the minimum has moved is less
      // than find_every_distance, don't find this minimum on this
      // iteration

      CCTK_REAL dt = cctk_time - min_last_t[m];
      CCTK_REAL dr = dt * sqrt(min_vx[m]*min_vx[m] + min_vy[m]*min_vy[m] + min_vz[m]*min_vz[m]);

      assert(dt >= 0);

      if (verbose >= 2)
      {
       CCTK_VInfo(CCTK_THORNSTRING,
                  "Minimum %d has moved by %g (estimated) in time %g",
                  m, dr, dt);
      }
      if (dr < find_every_distance)
      {
        if (verbose >= 2)
        {
          CCTK_VInfo(CCTK_THORNSTRING, "Skipping minimum find");
        }
        continue;
      }
    }

    if (find_every_distance == -1 && (cctk_iteration % find_every != 0))
      return;

    if (verbose >= 1)
    {
      CCTK_VInfo(CCTK_THORNSTRING,
                 "Finding minimum %d using Newton's method",
                 m);
    }

    if (reflection[m] == -1)
    {
      if (extrapolate_guess && vel_found[m])
      {
        CCTK_REAL dt = cctk_time - min_last_t[m];

        wmin[0] += dt * min_vx[m];
        wmin[1] += dt * min_vy[m];
        wmin[2] += dt * min_vz[m];

        CCTK_VInfo(CCTK_THORNSTRING,
                   "Minimum %d initial guess by extrapolation: \n[%g,%g,%g] ->[%g,%g,%g]",
                   m, min_x[m], min_y[m], min_z[m], wmin[0], wmin[1], wmin[2]);
      }
    
    while((df_norm = vector_norm_33(df))
          > newton_tolerance && iterations < iteration_limit)
    { 
      iterations++;
      if (!interpolate_point(CCTK_PASS_CTOC, wmin, &f, df, ddf))
      {
        exit(1);
      }

      if (verbose >= 2)
      {
        CCTK_VInfo(CCTK_THORNSTRING,
                   "\nMinTracker: At [%f,%f,%f], f = %f, df = [%f,%f,%f]\nddf = \n",  
                   wmin[0], wmin[1],wmin[2], f, df[0], df[1], df[2]);
        print_matrix(stdout, ddf);
      }

      for (int i = 0; i < 3; i++)
      {
        if (df[i] != df[i])
        {
          CCTK_VInfo(CCTK_THORNSTRING, "NaN found when interpolating");
        }
      }

      CCTK_REAL ddfInv[3][3];
      invert_matrix_33(ddf, ddfInv);

      CCTK_REAL id[3][3] = {44,44,44,44,44,44,44,44,44};

      multiply_matrix_matrix_33(ddf, ddfInv, id);

      CCTK_REAL amx[3]; // a minus x
      multiply_matrix_vector_33(ddfInv, df, amx);

      if (verbose >= 2)
      {
        CCTK_VInfo(CCTK_THORNSTRING, "a - x = [%g,%g,%g]", amx[0], amx[1], amx[2]);
      }

      for (int i = 0; i < 3; i++)
      {
        wmin[i] -= newton_gamma * amx[i];
      }

      df_norm = vector_norm_33(df);

      if (verbose >= 2)
      {
        CCTK_VInfo(CCTK_THORNSTRING, "newton_tolerance = %g\n", newton_tolerance);
        CCTK_VInfo(CCTK_THORNSTRING, "x = [%g,%g,%g], f = %g, df = [%g,%g,%g], |df| = %g\n",
               wmin[0], wmin[1],wmin[2], f, df[0], df[1], df[2], df_norm);
      }
      if (verbose >= 1)
      {
        CCTK_VInfo(CCTK_THORNSTRING,
                   "it = %d, |df| = %g, |delta| = %g, min = [%g,%g,%g]",
                   iterations-1, df_norm, vector_norm_33(amx),wmin[0],wmin[1],wmin[2]);
      }
    }

    if (iterations >= iteration_limit)
    {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Minimum not found within tolerance in %d iterations; giving best guess.",
                 iteration_limit); 
    }
    }
    else
    {
      // It is a reflection
      int m2 = reflection[m];
      wmin[0] = -min_x[m2];
      wmin[1] = -min_y[m2];
      wmin[2] = -min_z[m2];
    }

    if (min_found[m])
    {
      CCTK_REAL dt = cctk_time - min_last_t[m];
      CCTK_VInfo(CCTK_THORNSTRING, "dt = %g",
                 dt);

      CCTK_VInfo(CCTK_THORNSTRING, "cctk_time = %g",
                 cctk_time);

      CCTK_VInfo(CCTK_THORNSTRING, "min_last_t[%d] = %g",
                 m, min_last_t[m]);
      assert(dt >= 0);

      if (fabs(dt) > 1e-10)
      {
        min_vx[m] = (wmin[0] - min_x[m])/dt;
        min_vy[m] = (wmin[1] - min_y[m])/dt;
        min_vz[m] = (wmin[2] - min_z[m])/dt;
        vel_found[m] = 1;

        CCTK_REAL dr = dt * sqrt(min_vx[m]*min_vx[m] + min_vy[m]*min_vy[m] + min_vz[m]*min_vz[m]);
        if (verbose >= 1)
        {
          CCTK_VInfo(CCTK_THORNSTRING, "v = [%g,%g,%g]",
                     min_vx[m], min_vy[m], min_vz[m]);
          CCTK_VInfo(CCTK_THORNSTRING, "dr = %g", dr);
        }
      }
      else
      {
        CCTK_VInfo(CCTK_THORNSTRING, "dt = 0, hence cannot compute velocity");
      }
    }

    min_x[m] = wmin[0];
    min_y[m] = wmin[1];
    min_z[m] = wmin[2];
    min_r[m] = sqrt(pow(min_x[m],2) + pow(min_y[m],2) + pow(min_z[m],2));

    if (verbose > 2)
      CCTK_VInfo(CCTK_THORNSTRING, "min_r[%d] == %f", m, min_r[m]);
    
    int min_is_valid = !std::isnan(min_r[m]); // Check for NaN

    min_found[m] = 1;
    min_last_t[m] = cctk_time;

    if (verbose >= 1)
    {
      CCTK_VInfo(CCTK_THORNSTRING, "Converged to (%g,%g,%g)", 
                 min_x[m],min_y[m], min_z[m]);
    }
    
    int sind = surface_index[m];
    if (sind != -1) 
    {
      assert(sind >=0 && sind < nsurfaces);
      if (min_is_valid)
      {
        if (verbose > 2)
          CCTK_VInfo(CCTK_THORNSTRING, "Minimum is valid: setting spherical surface");
        sf_valid[sind]=1;
        sf_active[sind]=1;
        sf_centroid_x[sind]=min_x[m];
        sf_centroid_y[sind]=min_y[m];
        sf_centroid_z[sind]=min_z[m];
      }
      else
      {
        sf_valid[sind]=0;
        CCTK_VInfo(CCTK_THORNSTRING, "Minimum is invalid; setting sf_valid[%d] = 0", sind);
      }
    }
  }
}

extern "C" void MinTracker_Output(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  int output_every_loc = (output_every == -1) ? find_every : output_every;

  if (cctk_iteration % output_every_loc != 0)
    return;

  if (verbose >= 2)
    CCTK_INFO("MinTracker_Output");

  for (int m = 0; m < nminima; m++)
  {
    FILE *fp = 0;

    if (CCTK_MyProc(cctkGH) == 0 && (fp = output_files[m]->open()))
    {
      fprintf(fp, "%d\t%f\t%.19g\t%.19g\t%.19g\n", cctk_iteration, 
              cctk_time, min_x[m], min_y[m], min_z[m]);
      fclose(fp);
    }
  }
}
