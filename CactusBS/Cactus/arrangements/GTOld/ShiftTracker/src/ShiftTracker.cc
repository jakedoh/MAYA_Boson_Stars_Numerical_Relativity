#include "ShiftTracker.hh"
#include "ShiftTrackerInterpolate.hh"

//CCTK_FILEVERSION(AnalysisThorns_ShiftTracker_cc)

// ShiftTracker_UpdateLHS:
//   This does a very simple forward euler integration
//   of the ODE dt(x)=-beta(x)
extern "C" void ShiftTracker_UpdateLHS(CCTK_ARGUMENTS);
// ShiftTracker_CalcRHS:
//   This evaluates what the shift should be at the 
//   coordinate location, x.  This requires interpolation.
//   This function simply calls the interpolation and derivatives
//   routines.
extern "C" void ShiftTracker_CalcRHS(CCTK_ARGUMENTS);
  
extern "C" void ShiftTracker_UpdateLHS(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if( verbose > 2)
    CCTK_VInfo(CCTK_THORNSTRING,
        "ShiftTracker_UpdateLHS(): updating the trackers");
  
  // For each tracker evaluate integrate the ODE in time
  // Since v = -beta, and xdot = -beta, xdot = v
  for( int i=0; i < num_trackers; i++ )  {
    // If we have an odd number of intevals, add the last interval as dt*v
    // If we have an even number of intervals add it correctly for simpson's rule
   
    const CCTK_REAL dt = CCTK_DELTA_TIME;

    // This is a bit of a mess so I'll try to explain.  Since I can't do Simpson's 
    // rule util I have 3 points, I do Euler for the first iteration, the trapezoid
    // then I build into simpson's
    // this is equivalent to RK3 for functions depending on t only
    
    // it | rhs terms                    | change from last iter
    // ===|==============================|=======================
    // 0  | h*f0                         | +h*f0
    // 1  | 1/2*h*(f0+f1)                | +.5*h*f1 -.5*h*f0
    // 2  | 1/3*h*(f0+4*f1+f2)           | -1/6*h*f0 +5/6*h*f1 +1/3*h*f2
    // 3  | 1/3*h*(f0+4*f1+f2)+h*f3      | +h*f3
    // 4  | 1/3*h*(f0+4*f1+2*f2+4*f3+f4) | +1/3*h*(f2+f3+f4)

    // As you can see, once I reach step 3, steps three and 4 just repeat.  It is
    // also clear that I need to keep the previous 2 values of the rhs variable.
    if( cctk_iteration == 0 ) {
      st_x[i] += dt*st_vx[i];
      st_y[i] += dt*st_vy[i];
      st_z[i] += dt*st_vz[i];
    }
    else if( cctk_iteration == 1 ) {
      st_x[i] += 0.5*dt*(st_vx[i] - vx_last[i]);
      st_y[i] += 0.5*dt*(st_vy[i] - vy_last[i]);
      st_z[i] += 0.5*dt*(st_vz[i] - vz_last[i]);
    }
    else if( cctk_iteration == 2 ) { // once here, simpson's rule should be in effect
      st_x[i] += dt*(st_vx[i] +2.5*vx_last[i] - 0.5*vx_last2[i])/3.0;
      st_y[i] += dt*(st_vy[i] +2.5*vy_last[i] - 0.5*vy_last2[i])/3.0;
      st_z[i] += dt*(st_vz[i] +2.5*vz_last[i] - 0.5*vz_last2[i])/3.0;
    }
    else if( cctk_iteration%2 != 0 ) {
      st_x[i] += dt*st_vx[i];  // not even #of points so add the last in
      st_y[i] += dt*st_vy[i];  // not even #of points so add the last in
      st_z[i] += dt*st_vz[i];  // not even #of points so add the last in
    }
    else { //cctk_iteration%2 == 0
      st_x[i] += dt*(st_vx[i]+vx_last[i]+vx_last2[i])/3.0;
      st_y[i] += dt*(st_vy[i]+vy_last[i]+vy_last2[i])/3.0;
      st_z[i] += dt*(st_vz[i]+vz_last[i]+vz_last2[i])/3.0;
    }

    vx_last2[i] = vx_last[i];
    vx_last[i] = st_vx[i];
    
    vy_last2[i] = vy_last[i];
    vy_last[i] = st_vy[i];
    
    vz_last2[i] = vz_last[i];
    vz_last[i] = st_vz[i];
    
    if( verbose > 0 ) { 
      CCTK_VInfo(CCTK_THORNSTRING,
        "ShiftTracker_calcRHS(): Tracker %d location: x=%f,y=%f,z=%f",
        i,st_x[i],st_y[i],st_z[i]);
    } 
  }
}

extern "C" void ShiftTracker_CalcRHS(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Have to set the acceleration here to -betarhs
  // and the veolcity to -beta
  // Since they initially come from the interpolator as
  // +betarhs and +beta, I will do this explicity
  CCTK_REAL dt = CCTK_DELTA_TIME;
  for( int i=0; i < num_trackers; i++ ) {
    st_vx[i] *= -1.0; st_vy[i] *= -1.0; st_vz[i] *= -1.0;
    st_ax[i] = (st_vx[i] - vx_last[i])/dt;
    st_ay[i] = (st_vy[i] - vy_last[i])/dt;
    st_az[i] = (st_vz[i] - vz_last[i])/dt;
  }

  
}
