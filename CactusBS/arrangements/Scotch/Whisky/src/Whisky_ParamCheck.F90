 /*@@
   @file      Whisky_ParamCheck.F90
   @date      Sat Feb  9 23:48:01 2002
   @author    
   @desc 
   Parameter checking routine.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

 /*@@
   @routine    Whisky_ParamCheck
   @date       Sat Feb  9 23:48:43 2002
   @author     Ian Hawke
   @desc 
   Checks the parameters.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Whisky_ParamCheck(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer istat

  if (whisky_stencil > minval(cctk_nghostzones)) then
    call CCTK_PARAMWARN("The stencil is larger than the number of ghost zones. Answer will be dependent on processor number...")
  end if

  !! Module selection
  if (CCTK_EQUALS(riemann_solver,"invalid")) then
      call CCTK_PARAMWARN("Invalid value for Riemann solver (riemann_solver) specified.")
  end if
  if (CCTK_EQUALS(recon_method,"invalid")) then
      call CCTK_PARAMWARN("Invalid value for reconstruction method (recon_method) specified.")
  end if

  !! Thorn Activation
  if ( CCTK_EQUALS(evolution_method,"Scotch_QuasiMHD") .or. CCTK_EQUALS(evolution_method,"Scotch_MHD_Bvec") .or. CCTK_EQUALS(evolution_method,"Scotch_MHD_Avec") ) then
      call CCTK_IsFunctionAliased(istat, "Con2PrimMHD")
      if (istat.eq.0) then
          call CCTK_PARAMWARN("Did not find Con2PrimMHD aliased yet activate_mhd is on. Please activate an applicable recovery thorn.")
      end if
  else
      if ( .not.CCTK_EQUALS(evolution_method,"Scotch_HD") ) then
         call CCTK_PARAMWARN("Need to set HydroBase::evolution_method to 'Scotch_HD' to evolve with Maya now");
      end if
      call CCTK_IsFunctionAliased(istat, "Con2PrimPoly")
      if (istat.eq.0) then
          call CCTK_PARAMWARN("Did not find Con2PrimPoly aliased. Please activate an applicable recovery thorn.")
      end if
  end if

  !! EOS
  if (CCTK_EQUALS(whisky_eos_table,"2D_Polytrope").and.&
       (.not.CCTK_EQUALS(whisky_eos_type,"Polytype")) ) then
    call CCTK_PARAMWARN("When using the 2D_Polytrope EOS you need to set eos_type to Polytype")
  end if
  if (CCTK_EQUALS(whisky_eos_table,"Ideal_Fluid").and.&
       (.not.CCTK_EQUALS(whisky_eos_type,"General"))) then
    call CCTK_PARAMWARN("When using the Ideal_Fluid EOS you need to set eos_type to General")
  end if

  !! Metric Type
  if (.not.(CCTK_EQUALS(metric_type, "Physical").or.&
       CCTK_EQUALS(metric_type, "Static Conformal"))) then
    call CCTK_PARAMWARN("Whisky only knows how to deal with physical metric type as other types were unused. ")
  end if

  !! Required parameter
  if (use_mask .eq. 0) then
    call CCTK_PARAMWARN("Whisky requires you to set SpaceMask::use_mask = ""yes""")
  end if

  if (number_of_particles .gt. 0) then
    if (number_of_arrays .ne. 3) then
      call CCTK_PARAMWARN("If tracking particles then you must have number_of_arrays = 3")
    end if
  end if
  
end subroutine Whisky_ParamCheck
