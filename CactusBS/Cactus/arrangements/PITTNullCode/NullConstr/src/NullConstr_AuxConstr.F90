! vim: syntax=fortran

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

#define stereo .true.

subroutine NullConstr_AuxConstr(CCTK_ARGUMENTS)

  use cctk
  use NullGrid_Vars
  use NullInterp

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT, parameter :: izero = 0
  integer, parameter :: ik = kind(izero)

  integer :: l1, l2, l3, i

  call CCTK_ActiveTimeLevels(l1, cctkGH, "NullVars::realcharfuncs")
  call CCTK_ActiveTimeLevels(l2, cctkGH, "NullVars::cmplxcharfuncs_basic")
  call CCTK_ActiveTimeLevels(l3, cctkGH, "NullVars::cmplxcharfuncs_aux")

  if (min(l1,l2,l3).lt.1) then

     call CCTK_WARN(1, "cannot calculate constraints -- not enough allocated time-levels")

     Null_AuxConstr_nucn = 1.e+10
     Null_AuxConstr_ckcn = 1.e+10
     Null_AuxConstr_cbcn = 1.e+10

     if (stereo) then
        Null_AuxConstr_nucs = 1.e+10
        Null_AuxConstr_ckcs = 1.e+10
        Null_AuxConstr_cbcs = 1.e+10
     end if

     return
  end if

  do i = 1, nx

     call NullInterp_d1(Null_AuxConstr_cbcn(:,:,i), dcmplx(bcn(:,:,i)), 0_ik, 1_ik) 
     call NullInterp_d1(Null_AuxConstr_ckcn(:,:,i), sqrt(1. + jcn(:,:,i) * conjg(jcn(:,:,i))), 0_ik, 1_ik) 
     call NullInterp_d1(Null_AuxConstr_nucn(:,:,i), jcn(:,:,i), 2_ik, -1_ik) 

     if (stereo) then

        call NullInterp_d1(Null_AuxConstr_cbcs(:,:,i), dcmplx(bcs(:,:,i)), 0_ik, 1_ik) 
        call NullInterp_d1(Null_AuxConstr_ckcs(:,:,i), sqrt(1. + jcs(:,:,i) * conjg(jcs(:,:,i))), 0_ik, 1_ik) 
        call NullInterp_d1(Null_AuxConstr_nucs(:,:,i), jcs(:,:,i), 2_ik, -1_ik) 

        call NullInterp_3cinterp(cctkGH,&
            tmp_cgfn1, tmp_cgfs1, tmp_cgfn2, tmp_cgfs2, tmp_cgfn3, tmp_cgfs3,&
            Null_AuxConstr_ckcn(:,:,i), Null_AuxConstr_ckcs(:,:,i),&
            Null_AuxConstr_nucn(:,:,i), Null_AuxConstr_nucs(:,:,i),&
            Null_AuxConstr_cbcn(:,:,i), Null_AuxConstr_cbcs(:,:,i), 1_ik, 1_ik, 1_ik)

     end if

  end do

  Null_AuxConstr_nucn = nucn - Null_AuxConstr_nucn
  Null_AuxConstr_ckcn = ckcn - Null_AuxConstr_ckcn
  Null_AuxConstr_cbcn = cbcn - Null_AuxConstr_cbcn

  if (stereo) then
     Null_AuxConstr_nucs = nucs - Null_AuxConstr_nucs
     Null_AuxConstr_ckcs = ckcs - Null_AuxConstr_ckcs
     Null_AuxConstr_cbcs = cbcs - Null_AuxConstr_cbcs
  end if

end subroutine NullConstr_AuxConstr

