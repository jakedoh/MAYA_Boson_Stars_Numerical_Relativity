! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"


 subroutine NullSHRE_AnaSchMetricRecon(CCTK_ARGUMENTS)

!fills the current levels with analytic valuse of the metric, lapse, shift, and their (r,q,p,t) derivatives
  use cctk
  use NullGrid_Vars
  implicit none

    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    DECLARE_CCTK_FUNCTIONS

! zero components

    SHRE_dtgij = 0.d0
    SHRE_dtbeta = 0.d0
    SHRE_dqalpha = 0.d0
    SHRE_dpalpha = 0.d0
    SHRE_dtalpha = 0.d0

! spatial components

    SHRE_gij(:,:,1) = 1+8*mass*qs**2/pp**2/cr
    SHRE_gij(:,:,2) = 1+8*mass*qs**2/pp**2/cr
    SHRE_gij(:,:,3) = 8*mass*qs*ps/cr/pp**2
    SHRE_gij(:,:,4) = -8*mass*qs*ps/cr/pp**2
    SHRE_gij(:,:,5) = -4*mass*(-1+qs**2+ps**2)*qs/cr/pp**2
    SHRE_gij(:,:,6) = 4*mass*(-1+qs**2+ps**2)*qs/cr/pp**2 
    SHRE_gij(:,:,7) = 1+8*mass*ps**2/pp**2/cr
    SHRE_gij(:,:,8) = 1+8*mass*ps**2/pp**2/cr
    SHRE_gij(:,:,9) = -4*mass*(-1+qs**2+ps**2)*ps/cr/pp**2
    SHRE_gij(:,:,10) = -4*mass*(-1+qs**2+ps**2)*ps/cr/pp**2
    SHRE_gij(:,:,11) = 1+2*mass*(-1+qs**2+ps**2)**2/pp**2/cr
    SHRE_gij(:,:,12) = 1+2*mass*(-1+qs**2+ps**2)**2/pp**2/cr

! shift
    SHRE_beta(:,:,1) = -4*mass*qs/(cr+2*mass)/pp
    SHRE_beta(:,:,2) = -4*mass*qs/(cr+2*mass)/pp
    SHRE_beta(:,:,3) = -4*mass*ps/(cr+2*mass)/pp
    SHRE_beta(:,:,4) = 4*mass*ps/(cr+2*mass)/pp
    SHRE_beta(:,:,5) = 2*mass/(cr+2*mass)*(-1+qs**2+ps**2)/pp
    SHRE_beta(:,:,6) = -2*mass/(cr+2*mass)*(-1+qs**2+ps**2)/pp

! lapse
    SHRE_alpha(:,:,1) = 1/sqrt(cr+2*mass)*sqrt(cr)
    SHRE_alpha(:,:,2) = 1/sqrt(cr+2*mass)*sqrt(cr)

! radial derivatives of the metric
    SHRE_drgij(:,:,1) = -8*mass*qs**2/cr**2/pp**2 
    SHRE_drgij(:,:,2) = -8*mass*qs**2/cr**2/pp**2
    SHRE_drgij(:,:,3) = -8*mass*qs*ps/cr**2/pp**2
    SHRE_drgij(:,:,4) = 8*mass*qs*ps/cr**2/pp**2
    SHRE_drgij(:,:,5) = 4*mass*(-1+qs**2+ps**2)*qs/cr**2/pp**2
    SHRE_drgij(:,:,6) = -4*mass*(-1+qs**2+ps**2)*qs/cr**2/pp**2
    SHRE_drgij(:,:,7) = -8*mass*ps**2/cr**2/pp**2
    SHRE_drgij(:,:,8) = -8*mass*ps**2/cr**2/pp**2
    SHRE_drgij(:,:,9) = 4*mass*(-1+qs**2+ps**2)*ps/cr**2/pp**2
    SHRE_drgij(:,:,10) = 4*mass*(-1+qs**2+ps**2)*ps/cr**2/pp**2
    SHRE_drgij(:,:,11) = -2*mass*(-1+qs**2+ps**2)**2/pp**2/cr**2
    SHRE_drgij(:,:,12) = -2*mass*(-1+qs**2+ps**2)**2/pp**2/cr**2

! q derivatives of the metric
    SHRE_dqgij(:,:,1) = -16*mass*qs*(-1+qs**2-ps**2)/pp**3/cr 
    SHRE_dqgij(:,:,2) = -16*mass*qs*(-1+qs**2-ps**2)/pp**3/cr
    SHRE_dqgij(:,:,3) = -8*mass*ps*(-1+3*qs**2-ps**2)/pp**3/cr
    SHRE_dqgij(:,:,4) = 8*mass*ps*(-1+3*qs**2-ps**2)/pp**3/cr
    SHRE_dqgij(:,:,5) = 4*mass*(-6*qs**2+qs**4+1-ps**4)/pp**3/cr
    SHRE_dqgij(:,:,6) = -4*mass*(-6*qs**2+qs**4+1-ps**4)/pp**3/cr
    SHRE_dqgij(:,:,7) = -32*mass*ps**2/pp**3/cr*qs
    SHRE_dqgij(:,:,8) = -32*mass*ps**2/pp**3/cr*qs
    SHRE_dqgij(:,:,9) = 8*mass*qs*ps*(-3+qs**2+ps**2)/pp**3/cr
    SHRE_dqgij(:,:,10) = 8*mass*qs*ps*(-3+qs**2+ps**2)/pp**3/cr
    SHRE_dqgij(:,:,11) = 16*mass*(-1+qs**2+ps**2)*qs/pp**3/cr
    SHRE_dqgij(:,:,12) = 16*mass*(-1+qs**2+ps**2)*qs/pp**3/cr

! p derivatives of the metric
    SHRE_dpgij(:,:,1) = -32*mass*qs**2*ps/cr/pp**3 
    SHRE_dpgij(:,:,2) = -32*mass*qs**2*ps/cr/pp**3
    SHRE_dpgij(:,:,3) = 8*mass*qs*(1+qs**2-3*ps**2)/pp**3/cr
    SHRE_dpgij(:,:,4) = -8*mass*qs*(1+qs**2-3*ps**2)/pp**3/cr
    SHRE_dpgij(:,:,5) = 8*mass*qs*ps*(-3+qs**2+ps**2)/pp**3/cr
    SHRE_dpgij(:,:,6) = -8*mass*qs*ps*(-3+qs**2+ps**2)/pp**3/cr
    SHRE_dpgij(:,:,7) = 16*mass*ps*(1+qs**2-ps**2)/pp**3/cr
    SHRE_dpgij(:,:,8) = 16*mass*ps*(1+qs**2-ps**2)/pp**3/cr
    SHRE_dpgij(:,:,9) = -4*mass*(6*ps**2-ps**4-1+qs**4)/pp**3/cr
    SHRE_dpgij(:,:,10) = -4*mass*(6*ps**2-ps**4-1+qs**4)/pp**3/cr
    SHRE_dpgij(:,:,11) = 16*mass*(-1+qs**2+ps**2)*ps/pp**3/cr
    SHRE_dpgij(:,:,12) = 16*mass*(-1+qs**2+ps**2)*ps/pp**3/cr

! radial derivatives of the shift
    SHRE_drbeta(:,:,1) = 4*mass*qs/(cr+2*mass)**2/pp
    SHRE_drbeta(:,:,2) = 4*mass*qs/(cr+2*mass)**2/pp
    SHRE_drbeta(:,:,3) = 4*mass*ps/(cr+2*mass)**2/pp
    SHRE_drbeta(:,:,4) = -4*mass*ps/(cr+2*mass)**2/pp
    SHRE_drbeta(:,:,5) = -2*mass/(cr+2*mass)**2*(-1+qs**2+ps**2)/pp
    SHRE_drbeta(:,:,6) = 2*mass/(cr+2*mass)**2*(-1+qs**2+ps**2)/pp

! q derivatives of the shift
    SHRE_dqbeta(:,:,1) = 4*mass*(-1+qs**2-ps**2)/(cr+2*mass)/pp**2 
    SHRE_dqbeta(:,:,2) = 4*mass*(-1+qs**2-ps**2)/(cr+2*mass)/pp**2
    SHRE_dqbeta(:,:,3) = 8*mass/(cr+2*mass)*ps/pp**2*qs
    SHRE_dqbeta(:,:,4) = -8*mass/(cr+2*mass)*ps/pp**2*qs
    SHRE_dqbeta(:,:,5) = 8*mass*qs/(cr+2*mass)/pp**2
    SHRE_dqbeta(:,:,6) = -8*mass*qs/(cr+2*mass)/pp**2

! p derivatives of the shift
    SHRE_dpbeta(:,:,1) = 8*mass/(cr+2*mass)*ps/pp**2*qs
    SHRE_dpbeta(:,:,2) = 8*mass/(cr+2*mass)*ps/pp**2*qs
    SHRE_dpbeta(:,:,3) = -4*mass*(1+qs**2-ps**2)/(cr+2*mass)/pp**2
    SHRE_dpbeta(:,:,4) = 4*mass*(1+qs**2-ps**2)/(cr+2*mass)/pp**2
    SHRE_dpbeta(:,:,5) = 8*mass*ps/(cr+2*mass)/pp**2
    SHRE_dpbeta(:,:,6) = -8*mass*ps/(cr+2*mass)/pp**2

! radial derivatives of the lapse
    SHRE_dralpha(:,:,1) = mass/sqrt(cr+2*mass)**3/sqrt(cr)
    SHRE_dralpha(:,:,2) = mass/sqrt(cr+2*mass)**3/sqrt(cr)

!time components of the metric
    SHRE_git(:,:,1) = -4*qs*mass/cr/pp
    SHRE_git(:,:,2) = -4*qs*mass/cr/pp
    SHRE_git(:,:,3) = -4*mass*ps/cr/pp
    SHRE_git(:,:,4) = +4*mass*ps/cr/pp
    SHRE_git(:,:,5) = +2*mass*(-1+qs**2+ps**2)/cr/pp
    SHRE_git(:,:,6) = -2*mass*(-1+qs**2+ps**2)/cr/pp
    SHRE_git(:,:,7) = (2*mass-cr)/cr
    SHRE_git(:,:,8) = (2*mass-cr)/cr

!r derivative for time components of the metric
    SHRE_drgit(:,:,1) = 4*qs*mass/cr**2/pp 
    SHRE_drgit(:,:,2) = 4*qs*mass/cr**2/pp
    SHRE_drgit(:,:,3) = +4*mass*ps/cr**2/pp
    SHRE_drgit(:,:,4) = -4*mass*ps/cr**2/pp
    SHRE_drgit(:,:,5) = -2*mass*(-1+qs**2+ps**2)/cr**2/pp
    SHRE_drgit(:,:,6) = +2*mass*(-1+qs**2+ps**2)/cr**2/pp
    SHRE_drgit(:,:,7) = -2*mass/cr**2
    SHRE_drgit(:,:,8) = -2*mass/cr**2

!q derivative for time components of the metric
    SHRE_dqgit(:,:,1) = 4*mass*(-1+qs**2-ps**2)/pp**2/cr
    SHRE_dqgit(:,:,2) = 4*mass*(-1+qs**2-ps**2)/pp**2/cr
    SHRE_dqgit(:,:,3) = +8*mass/pp**2*qs*ps/cr
    SHRE_dqgit(:,:,4) = -8*mass/pp**2*qs*ps/cr
    SHRE_dqgit(:,:,5) = +8*mass*qs/pp**2/cr
    SHRE_dqgit(:,:,6) = -8*mass*qs/pp**2/cr
    SHRE_dqgit(:,:,7) = 0.d0
    SHRE_dqgit(:,:,8) = 0.d0

!p derivative for time components of the metric
    SHRE_dpgit(:,:,1) = 8*mass/pp**2*ps*qs/cr
    SHRE_dpgit(:,:,2) = 8*mass/pp**2*ps*qs/cr
    SHRE_dpgit(:,:,3) = -4*mass*(1+qs**2-ps**2)/pp**2/cr
    SHRE_dpgit(:,:,4) = +4*mass*(1+qs**2-ps**2)/pp**2/cr 
    SHRE_dpgit(:,:,5) = +8*mass/pp**2*ps/cr
    SHRE_dpgit(:,:,6) = -8*mass/pp**2*ps/cr
    SHRE_dpgit(:,:,7) = 0.d0
    SHRE_dpgit(:,:,8) = 0.d0

!t derivative for time components of the metric
    SHRE_dtgit(:,:,1) = 0.d0
    SHRE_dtgit(:,:,2) = 0.d0 
    SHRE_dtgit(:,:,3) = 0.d0 
    SHRE_dtgit(:,:,4) = 0.d0 
    SHRE_dtgit(:,:,5) = 0.d0 
    SHRE_dtgit(:,:,6) = 0.d0 
    SHRE_dtgit(:,:,7) = 0.d0 
    SHRE_dtgit(:,:,8) = 0.d0 

 end subroutine NullSHRE_AnaSchMetricRecon
