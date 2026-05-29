#include "ShiftTrackerInterpolate.hh"

//ShiftTrackerInterpolate.cc

extern "C" void ShiftTracker_Interp(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // These are hardcoded magic numbers because they have to
  // correspond to the number of input and output arrays.
  // Since that is fixed, these should also be fixed.
  const CCTK_INT num_input_arrays = 3;
  const CCTK_INT num_output_arrays = 3;

  // if you don't know what ierr is you should not be coding.
  int ierr=-1;

  // This is and array of pointers to arrays that contain the 
  // coordinate locations of the points that need to be interpolated
  // to.
  const void* interp_coords[3]
    = { (const void *) st_x,
        (const void *) st_y,
        (const void *) st_z };

  // This is an array of variable indices cooresponding to the
  // physical variables that need interpolation
  const CCTK_INT input_array_indices[num_input_arrays]
    = { CCTK_VarIndex(beta1_var),
        CCTK_VarIndex(beta2_var),
        CCTK_VarIndex(beta3_var) };
  for(int i = 0 ; i < num_input_arrays ; i++) 
  {
    if (input_array_indices[i] < 0)
    {
      const CCTK_STRING varnames[num_input_arrays] 
	  = {beta1_var, beta2_var, beta3_var};
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
          "ShiftTracker: couldn't find variable \"%s\"!",
          varnames[i]);         /*NOTREACHED*/
      return;
    }
  }

  // Another array that keeps the type of the output arrays
  // double, int, ...
  const CCTK_INT output_array_types[num_output_arrays]
    = { CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL };

  // Another one of these arrays of pointers to arrays that will
  // contain the interpolated values of the arrays defined above
  // in input_array_indices at the coordinates defined by interp_coords
  void * output_arrays[num_output_arrays]
    = { (void *) st_vx, (void *) st_vy, (void *) st_vz };

  // Setup the various handles
  // =========================

  // This tells the interpolater what kind of interpolation it should do.
  // This has only been tested with the hermite polynomial interp from
  // the AEI.
  const int operator_handle = CCTK_InterpHandle(interpolator_name);
  if (operator_handle < 0)
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
        "ShiftTracker_Interp(): couldn't find interpolator \"%s\"!",
        interpolator_name);         /*NOTREACHED*/

  // This tells the interpolater what arguments it should use, e.g.
  // interpolation order.
  int param_table_handle = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT);

  // User parameters can be set in the parameter file
  ierr = Util_TableSetFromString(param_table_handle, interpolator_pars);

  if (ierr < 0)
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
        "ShiftTracker_Interp(): bad interpolator parameter(s) \"%s\"!",
        interpolator_pars);         /*NOTREACHED*/
  
  if (interp_num_timelevels != -1)
  {
    // Don't interpolate in time.  I do this because I think that carpet's
    // dertermination of when it should interpolate in time is bad.  It 
    // should never interpolate this in time.
    ierr = Util_TableSetInt(param_table_handle, interp_num_timelevels, 
                            "InterpNumTimelevels");
  }

  if (ierr < 0)
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
        "ShiftTracker_Interp(): bad interpolator parameter(s) \"%d\"!",
        1);         /*NOTREACHED*/
  
  // This tells the interpolator what kind of coordinate system we are on
  // I use cart3d.
  const int coord_system_handle = CCTK_CoordSystemHandle(coord_system);
  if (coord_system_handle < 0)
    CCTK_VWarn(-1, __LINE__, __FILE__, CCTK_THORNSTRING,
        "ShiftTracker_Interp(): can't get coordinate system handle for\n"
        "                           for coordinate system \"%s\"!",
        coord_system);                               /*NOTREACHED*/

  ierr = CCTK_InterpGridArrays(cctkGH,
      3, // number of dimensions used, hardcoded elsewhere
      operator_handle,
      param_table_handle,
      coord_system_handle,
      num_trackers, // One point in each dimension for each tracker
      CCTK_VARIABLE_REAL,
      interp_coords,
      num_input_arrays, // Number of input arrays to be interpolated
      input_array_indices,
      num_output_arrays, // Number of output arrays
      output_array_types,
      output_arrays);

  if( verbose > 1 ) {
    for( int i=0; i < num_trackers; i++ ) {
      CCTK_VInfo(CCTK_THORNSTRING,
          "ShiftTracker_Interp(): Tracker %d RHS: x=%f,y=%f,z=%f",
          i,st_vx[i],st_vy[i],st_vz[i]);
    }
  }
  // because some people really like to see output.
  if( verbose > 1 )
    CCTK_VInfo(CCTK_THORNSTRING,
        "ShifTracker_Interp(): status=%d return from interpolator",
        ierr);

  Util_TableDestroy(param_table_handle);
}


