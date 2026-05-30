 /*@@
   @file      Whisky_Minima.cc
   @date      Tue Aug 29 18:52:10 2006
   @author    
   @desc 
        Sets up the scalars used for the atmosphere, after initial data.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

//#include "Carpet/Carpet/src/carpet.hh"
#include "carpet.hh"

#ifdef HAVE_CARPET
using namespace Carpet;
#endif

#ifdef __cplusplus
  extern "C" {
#endif
    
    /* Scheduled functions */
    void Whisky_Rho_Minima_Setup_Final(CCTK_ARGUMENTS);

    void Whisky_Rho_Minima_Setup_Final_PUGH(CCTK_ARGUMENTS);
    
#ifdef __cplusplus
  } /* extern "C" */
#endif

 /*@@
   @routine    Whisky_Rho_Minima_Setup_Final
   @date       Tue Aug 29 18:54:05 2006
   @author     Luca Baiotti
   @desc 
        After initial data, set up the scalar whisky_rho_min used for checking the atmosphere.
        This is computed taking into account the actual maximum rest-mass density on the grid.
	Version for Carpet.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
@@*/

void Whisky_Rho_Minima_Setup_Final(CCTK_ARGUMENTS)
{

#ifdef HAVE_CARPET

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL  max_rho, max_tau;

  static int flag = true;
  
  if (flag) // Actually run this reduction only the first time the routine is called.
    {
      if (rho_abs_min > 0.0)
	{
	  *whisky_rho_min = rho_abs_min;
	}
      else
	{	  
	  // Go to global mode

          BEGIN_GLOBAL_MODE (cctkGH) {

            const int Reduction_Handle = CCTK_ReductionHandle ("maximum");

            // Find the global maximum of rho
            const CCTK_INT input_array_variable_indices_rho= {CCTK_VarIndex("HydroBase::rho")};
            const int ierrRho = CCTK_Reduce(cctkGH,
					 -1, // target processors; -1 -> all
					 Reduction_Handle,
					 1,  // number of output variables
					 CCTK_VARIABLE_REAL,
					 &max_rho,
					 1,  // number of variables to be reduced
					 input_array_variable_indices_rho);
	    
            if (ierrRho != 0)
	      {
		CCTK_WARN(0, "Failed to compute the global maximum of rho");
	      }
	    
            *whisky_rho_min = max_rho * rho_rel_min;

            // Find the global maximum of tau
            const CCTK_INT input_array_variable_indices_tau= {CCTK_VarIndex("Whisky::tau")};
            const int ierrTau = CCTK_Reduce(cctkGH,
					 -1, // target processors; -1 -> all
					 Reduction_Handle,
					 1,  // number of output variables
					 CCTK_VARIABLE_REAL,
					 &max_tau,
					 1,  // number of variables to be reduced
					 input_array_variable_indices_tau);
	    
            if (ierrTau != 0)
	      {
		CCTK_WARN(0, "Failed to compute the global maximum of tau");
	      }
	    
            *whisky_tau_min = max_tau * tau_rel_min;


	    

	    
            // Go back to local mode
          } END_GLOBAL_MODE;
	  
	}
      // After this has run once, set the flag so that this does not run again
      flag = false;
    } // end if (flag)
  else
    {
      return;
    }

#endif

}


 /*@@
   @routine    Whisky_Rho_Minima_Setup_Final_PUGH
   @date       Tue Aug 29 18:57:42 2006
   @author     Luca Baiotti
   @desc 
         As above, but for PUGH.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
@@*/

void Whisky_Rho_Minima_Setup_Final_PUGH(CCTK_ARGUMENTS)
{ 
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_REAL  max_rho, max_tau;
  
  if (rho_abs_min > 0.0)
    {
      *whisky_rho_min = rho_abs_min;
    }
  else
    {      
      

      // Find the global maximum of rho
      const int Reduction_Handle = CCTK_ReductionHandle("maximum");
      const CCTK_INT input_array_variable_indices={CCTK_VarIndex("HydroBase::rho")};
      const int ierr = CCTK_Reduce(cctkGH,
			 -1, // target processors; -1 -> all
			 Reduction_Handle,
			 1,  // number of output variables   
			 CCTK_VARIABLE_REAL,
			 &max_rho,
			 1,  // number of variables to be reduced
			 input_array_variable_indices);
      if (ierr != 0)
	{
	  CCTK_WARN(0, "Failed to compute the global maximum of rho");
	}
      
      *whisky_rho_min = max_rho * rho_rel_min;    

    }
   
  if (tau_abs_min > 0.0)
    {
     *whisky_tau_min = tau_abs_min;
    }
  else
    {      

      // Find the global maximum of tau
      const int Reduction_Handle = CCTK_ReductionHandle("maximum");
      const CCTK_INT input_array_variable_indices={CCTK_VarIndex("Whisky::tau")};
      const int ierr = CCTK_Reduce(cctkGH,
			 -1, // target processors; -1 -> all
			 Reduction_Handle,
			 1,  // number of output variables   
			 CCTK_VARIABLE_REAL,
			 &max_tau,
			 1,  // number of variables to be reduced
			 input_array_variable_indices);
      if (ierr != 0)
	{
	  CCTK_WARN(0, "Failed to compute the global maximum of tau");
	}
      
      *whisky_tau_min = max_tau * tau_rel_min;    
 
    }
}

