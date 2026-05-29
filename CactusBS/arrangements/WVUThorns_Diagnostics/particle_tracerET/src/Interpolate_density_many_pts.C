#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "util_Table.h"
#include "util_String.h"

void Interpolate_density_many_pts(cGH *cctkGH,int interp_num_points,double *particle_x_temp,double *particle_y_temp,double *particle_z_temp, double *particle_density_temp) {
  DECLARE_CCTK_PARAMETERS; 
  int ierr;

  int NUM_INPUT_ARRAYS=1;
  int NUM_OUTPUT_ARRAYS=1;

  /* NEEDED?
  // uni-processor code - only work on CPU 0
  const CCTK_INT myproc= CCTK_MyProc(cctkGH);
  if ( myproc != 0 ) {
    interp_npoints=0;
  }
  */

  CCTK_STRING coord_system = "cart3d";

  // Set up handles
  const int coord_system_handle = CCTK_CoordSystemHandle(coord_system);
  if (coord_system_handle < 0) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
        "can't get coordinate system handle for coordinate system \"%s\"!",
               coord_system);
  }

  const int operator_handle = CCTK_InterpHandle(interpolator_name);
  if (operator_handle < 0)
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "couldn't find interpolator \"%s\"!",
               interpolator_name);

  int param_table_handle = Util_TableCreateFromString(interpolator_pars);
  if (param_table_handle < 0) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "bad interpolator parameter(s) \"%s\"!",
               interpolator_pars);
  }
  
  CCTK_INT operand_indices[1]; //NUM_OUTPUT_ARRAYS + MAX_NUMBER_EXTRAS];
  for(int i = 0 ; i < 1  ; i++) {
    operand_indices[i] = i;
  }
  Util_TableSetIntArray(param_table_handle, NUM_OUTPUT_ARRAYS,
                        operand_indices, "operand_indices");
  

  CCTK_INT operation_codes[1]; //NUM_OUTPUT_ARRAYS + MAX_NUMBER_EXTRAS];
  for(int i = 0 ; i < 1  ; i++) {
    operation_codes[i] = 0;
  }
  Util_TableSetIntArray(param_table_handle, NUM_OUTPUT_ARRAYS, 
                        operation_codes, "operation_codes");

  const void* interp_coords[3] 
    = { (const void *) particle_x_temp,
        (const void *) particle_y_temp,
        (const void *) particle_z_temp };

  // 3d input arrays
  CCTK_STRING input_array_names[1]
    = { "HydroBase::rho" };
  //    = { "IllinoisGRMHD::rho_b" };

  CCTK_INT input_array_indices[1];
  for(int i = 0 ; i < NUM_INPUT_ARRAYS ; i++) {
    input_array_indices[i] = CCTK_VarIndex(input_array_names[i]);
    if(input_array_indices[i] < 0) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
        "COULD NOT FIND VARIABLE '%s'.",
        input_array_names[i]);
      exit(1);
    }
  }

  CCTK_INT output_array_types[1];
  for(int i = 0 ; i < NUM_OUTPUT_ARRAYS ; i++) {
    output_array_types[i] = CCTK_VARIABLE_REAL;
  }

  void * output_arrays[1]
    = { (void *) particle_density_temp };

  // actual interpolation call
  ierr = CCTK_InterpGridArrays(cctkGH,
                               3, // number of dimensions 
                               operator_handle,
                               param_table_handle,
                               coord_system_handle,
                               interp_num_points,
                               CCTK_VARIABLE_REAL,
                               interp_coords,
                               NUM_INPUT_ARRAYS, // Number of input arrays
                               input_array_indices,
                               NUM_OUTPUT_ARRAYS, // Number of output arrays
                               output_array_types,
                               output_arrays);
  if (ierr<0) {
    CCTK_WARN(1,"interpolation screwed up");
    Util_TableDestroy(param_table_handle);
    exit(1);
  }

  ierr = Util_TableDestroy(param_table_handle);
  if (ierr != 0) {
    CCTK_WARN(1,"Could not destroy table");
    exit(1);
  }
  /* - u_t - 1 */
}
