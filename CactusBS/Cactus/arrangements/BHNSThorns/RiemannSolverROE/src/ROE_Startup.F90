 /*@@
   @file      ROE_Startup.F90
   @date      Sun Feb 10 00:02:52 2002
   @author    Ian Hawke
   @desc 
   Startup banner.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

 /*@@
   @routine    ROE_Startup
   @date       Sun Feb 10 00:03:09 2002
   @author     Ian Hawke
   @desc 
   Startup banner.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

integer function ROE_Startup()

  USE ROE_Scalars

  implicit none 

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT :: var_names_param_len, var_gfs_len
  character(len=256) :: eos_name_param, var_names_param, var_gfs
  integer ROE_MakeEOSCall

!!$  Set up EOS calls.

  if (use_eosgeneral .ne. 0) then

    ! names of the extra EOS variables (never changes)
    call CCTK_FortranString(var_names_param_len, eosgeneral_indep_names, &
         var_names_param)

!!$    First, "Plus" boundary extended
    call CCTK_FortranString(var_gfs_len, eosgeneral_indep_plus_gfs, var_gfs)
    if (CCTK_EQUALS(whisky_eos_type, "Polytype")) then
      EOS_RiemannCallPlus = ROE_MakeEOSCall("Rho "//var_names_param(1:var_names_param_len), &
        "whisky::rhoplus "//var_gfs(1:var_gfs_len), 1+eosgeneral_n_indeps, &
        "Pressure DPressureDSpecificInternalEnergy c_s^2", &
        "whisky::pressplus whisky::eos_dpdeps_p whisky::eos_cs2_p", &
        3)
    else
      EOS_RiemannCallPlus = ROE_MakeEOSCall("Rho SpecificInternalEnergy "//var_names_param(1:var_names_param_len), &
        "whisky::rhoplus whisky::epsplus "//var_gfs(1:var_gfs_len), 2+eosgeneral_n_indeps, &
        "Pressure DPressureDSpecificInternalEnergy c_s^2", &
        "whisky::pressplus whisky::eos_dpdeps_p whisky::eos_cs2_p", &
        3)
    end if

!!$    Second, "Minus" boundary extended
    call CCTK_FortranString(var_gfs_len, eosgeneral_indep_minus_gfs, var_gfs)
    if (CCTK_EQUALS(whisky_eos_type, "Polytype")) then
      EOS_RiemannCallMinus = ROE_MakeEOSCall("Rho "//var_names_param(1:var_names_param_len), &
        "whisky::rhominus "//var_gfs(1:var_gfs_len), 1+eosgeneral_n_indeps, &
        "Pressure DPressureDSpecificInternalEnergy c_s^2", &
        "whisky::pressminus whisky::eos_dpdeps_m whisky::eos_cs2_m", &
        3)
    else
      EOS_RiemannCallMinus = ROE_MakeEOSCall("Rho SpecificInternalEnergy "//var_names_param(1:var_names_param_len), &
        "whisky::rhominus whisky::epsminus "//var_gfs(1:var_gfs_len), 2+eosgeneral_n_indeps, &
        "Pressure DPressureDSpecificInternalEnergy c_s^2", &
        "whisky::pressminus whisky::eos_dpdeps_m whisky::eos_cs2_m", &
        3)
    end if

!!$ Third the "Average" state for the Roe solver
    call CCTK_FortranString(var_gfs_len, eosgeneral_indep_gfs, var_gfs)
    if (CCTK_EQUALS(whisky_eos_type, "Polytype")) then
      EOS_RoeAverageCall = ROE_MakeEOSCall("Rho "//var_names_param(1:var_names_param_len), &
        "riemannsolverroe::rho_ave "//var_gfs(1:var_gfs_len), 1+eosgeneral_n_indeps, &
        "Pressure DPressureDSpecificInternalEnergy c_s^2", &
        "riemannsolverroe::press_ave riemannsolverroe::eos_dpdeps_ave riemannsolverroe::eos_cs2_ave", &
        3)
    else
      EOS_RoeAverageCall = ROE_MakeEOSCall("Rho SpecificInternalEnergy "//var_names_param(1:var_names_param_len), &
        "riemannsolverroe::rho_ave riemannsolverroe::eps_ave "//var_gfs(1:var_gfs_len), 2+eosgeneral_n_indeps, &
        "Pressure DPressureDSpecificInternalEnergy c_s^2", &
        "riemannsolverroe::pressminus riemannsolverroe::eos_dpdeps_ave riemannsolverroe::eos_cs2_ave", &
        3)
    end if
    
  end if
  
  ROE_Startup = 0

end function ROE_Startup
