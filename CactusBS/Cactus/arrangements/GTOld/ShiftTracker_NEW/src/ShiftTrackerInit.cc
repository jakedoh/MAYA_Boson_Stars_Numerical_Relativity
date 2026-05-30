#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

// ShiftTracker_Init:
//   This simply initializes the tracker locations to those
//   specified in the parameter file by the user.
extern "C" void ShiftTracker_Init(CCTK_ARGUMENTS);

extern "C" void ShiftTracker_Init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if( !( CCTK_EQUALS( beta1rhs_var, "DEP") && 
         CCTK_EQUALS( beta2rhs_var, "DEP") && 
         CCTK_EQUALS( beta3rhs_var, "DEP") ) )
  {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
        "ShiftTracker_Init(): You have set depracted parameters beta*rhs, these are not used\n");
  }

  for( int i=0; i < num_trackers; i++ ) {
    st_x[i] = x0[i];
    st_y[i] = y0[i];
    st_z[i] = z0[i];
    vx_last[i] = 0.0;
    vy_last[i] = 0.0;
    vz_last[i] = 0.0;
    if( verbose > 0 ) {
      CCTK_VInfo(CCTK_THORNSTRING,
          "ShiftTracker_Init(): Init tracker %d coords set to x=%f,y=%f,z=%f",
          i,
          (double)st_x[i], 
          (double)st_y[i], 
          (double)st_z[i]);
    }
  }
}
