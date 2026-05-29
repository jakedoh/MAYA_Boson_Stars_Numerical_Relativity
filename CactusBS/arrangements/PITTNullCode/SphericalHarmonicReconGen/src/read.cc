#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "util_Table.h"
#include "util_ErrorCodes.h"
#include <stdio.h>
#include <assert.h>


#include "vars.hh"


using namespace std;

namespace
{
  // returns true when the interpolation uses the last valid iteration in the file
  bool load_modes(const CCTK_REAL time, const bool time_derivative_in_file, const bool allow_offcentered_time_stencils)
  {
    // have we just read the last dataset from the files?
    bool read_last_iteration = false;
    // identify iteration which corresponds to the time just before the target time "time"
    int iteration = SHR::C[0]->get_iteration(time + SHR::initial_time_in_file);
    
    if (iteration < 0) {
      stringstream str;
      str << "Error: iteration = -1. time = "<<time << ", SHR::initial_time_in_file="<<SHR::initial_time_in_file;
      CCTK_ERROR(str.str().c_str());
    }
   
    if (iteration < 2) {
      if(allow_offcentered_time_stencils) 
         iteration=2;
      else {
         CCTK_ERROR("Error: there are not enough iterations before interpolation target time to do the interpolation!");
      }
    }
    
    // set all modes from file necessary for 5th-order time interpolation
    for (int k=-2; k <= 3; ++k) {
      for (int i=0; i < NUM_METRIC_COMPONENTS; ++i) {
        SHR::C[i]->set_modes_from_file(5-(k+2), iteration+k);
        SHR::Cr[i]->set_modes_from_file(5-(k+2), iteration+k);
        if (time_derivative_in_file)
          SHR::Ct[i]->set_modes_from_file(5-(k+2), iteration+k);
      }
      for (int i=0; i < SHR::db.size(); ++i)
        read_last_iteration |= (iteration+k == SHR::db[i]->n_timesteps()-1);
    }

    return read_last_iteration;
  }

  void interpolate(const CCTK_REAL time, const int tl, const bool time_derivative_in_file)
  {
    // coefficients
    vector<CCTK_REAL> L(6);
    
    CCTK_REAL x[6] = { SHR::C[0]->get_time_of_coeff_in_file(5),
                       SHR::C[0]->get_time_of_coeff_in_file(4), 
                       SHR::C[0]->get_time_of_coeff_in_file(3),
                       SHR::C[0]->get_time_of_coeff_in_file(2),
                       SHR::C[0]->get_time_of_coeff_in_file(1),
                       SHR::C[0]->get_time_of_coeff_in_file(0) };
    
    // precompute Lagrange coefficients
    for (int i=0; i < 6; ++i)
    {
       L[i] = 1;
       for (int j=0; j < 6; ++j)
       {
          if (i != j)
             L[i] *= (time+SHR::initial_time_in_file - x[j]) / (x[i]-x[j]);
       }
    }
    
    // interpolate and fill internal coeffs[tl]
    for (int i=0; i < NUM_METRIC_COMPONENTS; ++i) {
      SHR::C[i]->interpolate_to_(tl, L);
      SHR::Cr[i]->interpolate_to_(tl, L);
      if (time_derivative_in_file)
        SHR::Ct[i]->interpolate_to_(tl, L);
    }
  }
}


extern "C"
{
  void SphericalHarmonicReconGeneric_SetTimeStep(CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    if (! SHR::initialized)
    {
      CCTK_WARN(CCTK_WARN_ABORT, "Schedule mismatch. \n"
          "      SphericalHarmonicReconGeneric_SetTimeStep must be called after\n"
          "      SphericalHarmonicReconGeneric_Startup");
    }

    CCTK_REAL &dt = cctkGH->cctk_delta_time;
    // no adaptivity for the 1st few time-steps
    if(cctk_iteration<5) {
      if(not time_interpolate)
        CCTK_WARN(1, "cannot adjust time-step when time_interpolate==false");
      return;
    }
    SHR::SPH_database& db=*(SHR::db[0]);
    double current_time_gap_in_file = -1;
    const double time=cctk_time + SHR::initial_time_in_file;
    const double eps=1.e-10*(db.time(1)-db.time(0));
    if(db.n_timesteps()<2) return; // no adjustment
    if(time<db.time(0)-eps) return; // no adjustment
    if(time>eps+db.time(db.n_timesteps()-1)) return; // no adjustment

    // locate by bisection the interval containing time
    for(int l0=0,l1=db.n_timesteps()-1,c=0;;++c) {
      const int lh=static_cast<int>(0.5*(l0+l1));
      if(time<db.time(lh)-eps) l1=lh;
      else l0=lh;
      if(l1==l0+1) {
        current_time_gap_in_file = db.time(l1)-db.time(l0);
        break;
      }
      if(c>2+db.n_timesteps()) CCTK_WARN(0, "run-away bisection algorithm");
    }

    if(current_time_gap_in_file<0)
      CCTK_WARN(0, "cannot identify time gap in data file");

    double target_dt = current_time_gap_in_file/
      adjust_timestep__target_steps_per_time_level;
    // change the time-step gradually
    target_dt = std::min(target_dt,adjust_timestep__maximum_time_step);
    dt += adjust_timestep__rate_of_change*(target_dt-dt);

    if(verbose) CCTK_VInfo(CCTK_THORNSTRING,
        "Time-step has been set to %g",(double)dt);
  }
  
  void SphericalHarmonicReconGeneric_ReadData(CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    if (! SHR::initialized)
    {
      CCTK_WARN(CCTK_WARN_ABORT, "Schedule mismatch. \n"
         "      SphericalHarmonicReconGeneric_ReadData must be called after\n"
         "      SphericalHarmonicReconGeneric_Startup");
    }

    if (!time_interpolate)
    {
	 // read coefficients for current iteration
	 // iteration offset in boundary data file.
	 int it_off = 0;
	 
	 if (!time_derivative_in_file)
	 {
	    if (time_fd_order == 2) it_off = 2;
	    if (time_fd_order == 4) it_off = 4;
	 }
	 
	 it_off += SHR::starting_iteration;
	 
	 const int my_iteration = cctk_iteration - (cctk_time < initial_relaxation_time ? cctk_iteration : initial_relaxation_time / CCTK_DELTA_TIME);
	 
	 {
	    if (my_iteration == 0 || my_iteration*CCTK_DELTA_TIME >= SHR::iteration_in_file*SHR::db[0]->delta_t())
	    {
	       if (verbose)
	          CCTK_INFO("Cycling timesteps...");
	       
	       // cycle timelevels!
	       for (int i=0; i < NUM_METRIC_COMPONENTS; ++i) {
		  SHR::C[i]->cycle_timelevels();
		  SHR::Cr[i]->cycle_timelevels();
		  if (time_derivative_in_file)
		     SHR::Ct[i]->cycle_timelevels();
	       }
	       
	       if (verbose)
	          CCTK_INFO("Loading decomposed variables for current timestep...");
	       
	       for (int i=0; i < 3*NUM_METRIC_COMPONENTS; ++i) {
		  if ((my_iteration+it_off) == SHR::db[i]->n_timesteps()-1) { // will run out in next step
                     static bool reported_done = false;
                     CCTK_TerminateNext(cctkGH);
                     if (!reported_done) {
                        CCTK_INFO("Done reading data from file.");
                        reported_done = true;
                     }
                  }
		  if ((my_iteration+it_off) >= SHR::db[i]->n_timesteps())
		     CCTK_WARN(0, "The input file does not contain more iterations. Stopping.");
	       }
	       
	       // we always read n timesteps ahead (timelevel==0 corresponds to cctk_iteration+it_off) because we need it for time-derivatives!
	       for (int i=0; i < NUM_METRIC_COMPONENTS; ++i) {
		  SHR::C[i]->set_modes_from_file(0, (my_iteration+it_off));
		  SHR::Cr[i]->set_modes_from_file(0, (my_iteration+it_off));
		  if (time_derivative_in_file)
		     SHR::Ct[i]->set_modes_from_file(0, (my_iteration+it_off));
	       }
	       
	       if (initial_relaxation_time < cctk_time)
	          SHR::iteration_in_file += 1;
	    }
	 }
	 
	 // copy internal coefficients array to target coefficient array (since we do not use interpolation)
	 // this is necessary because we keep two internal coeff arrays (interpolated and non-interpolated)
	 for (int i=0; i < NUM_METRIC_COMPONENTS; ++i) {
	    SHR::C[i]->copy_timelevels();
	    SHR::Cr[i]->copy_timelevels();
	    if (time_derivative_in_file)
	       SHR::Ct[i]->copy_timelevels();
	 }
    }
    else
    {
       // get current times. u[2] is target time.
       // Other times are times of stencil points for taking time derivatives
       CCTK_REAL u[5] = {initial_time+cctk_time-2*cctk_delta_time, 
                         initial_time+cctk_time-cctk_delta_time, 
                         initial_time+cctk_time, 
                         initial_time+cctk_time+cctk_delta_time, 
                         initial_time+cctk_time+2*cctk_delta_time };
       bool read_last_iteration = false;
      
       // If the automatically computed time-offset is disabled,
       // we allow for off-centered stencils at the beginning.
       // so that, when the data in the file starts at t=t0,
       // the evolution can also start at t=t0, implying off-centered
       // time-stencils very early in the run. 
       if (!time_derivative_in_file)
       {
          // interpolate coefficients to all stencil points given by u[...]
          
          if (time_fd_order == 2)
             for (int i=0; i < 3; ++i)
	     {
	        read_last_iteration |=
                  load_modes(u[i], false, disable_auto_time_offset);
	        interpolate(u[i], 2-i, false);
	     }
          
	  if (time_fd_order == 4)
	     for (int i=0; i < 5; ++i)
	     {
	        read_last_iteration |=
	          load_modes(u[i], false, disable_auto_time_offset);
	        interpolate(u[i], 4-i, false);
	     }
	  
       }
       else
       {
          // just time interpolate coefficients to u[2]
          
          // load all modes centered around time u[2] to internal "coeff_in_file" array...
          read_last_iteration |=
            load_modes(u[2], true, disable_auto_time_offset);
          
          // ...and use these coeffs to setup internal interpolated coeffs "coeff"
          // and copy obtained value to timelevel 0
          interpolate(u[2], 0, true);
          
       }

       // check if current interpolation used the last iteration in the file, then
       // request termination so as to not run out of data. This may terminate a
       // bit early if the last iteration can be used for more than one
       // interpolation.
       if (read_last_iteration) { // will (likely) run out in next step
          CCTK_TerminateNext(cctkGH);
          CCTK_INFO("Likely done reading data from file.");
       }
       
    }

    SHR::read_data = 1;
  }






  void SphericalHarmonicReconGeneric_PostStep(CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    SHR::read_data = 0;
  }

}


