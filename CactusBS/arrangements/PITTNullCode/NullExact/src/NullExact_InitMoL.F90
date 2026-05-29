! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

 subroutine NullExact_InitMoL(CCTK_ARGUMENTS)
 use cctk
 use NullExact_Analytic
 use NullGrid_Vars
 implicit none

 CCTK_COMPLEX, dimension(lsh(1), lsh(2), 2) :: initNews, initdotNews
 CCTK_COMPLEX, dimension(lsh(1), lsh(2), 2) :: pinitNews, pinitdotNews
 CCTK_INT :: patch_ID

 DECLARE_CCTK_ARGUMENTS
 DECLARE_CCTK_FUNCTIONS
 DECLARE_CCTK_PARAMETERS

 !I need to initialize News_MoL, dotNews_MoL 
  patch_ID = ip_n-1

     call NullExact_Analytic_BondiNews_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(N_radial_pts),cctk_time,initNews(:,:,patch_ID+1),initdotNews(:,:,patch_ID+1), Ylm_2)
     call NullExact_Analytic_BondiNews_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(N_radial_pts),cctk_time,initNews(:,:,ip_s),initdotNews(:,:,ip_s), Ylm_2)

     call NullExact_Analytic_BondiNews_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(N_radial_pts),cctk_time-cctk_delta_time,pinitNews(:,:,patch_ID+1),pinitdotNews(:,:,patch_ID+1), Ylm_2)
     call NullExact_Analytic_BondiNews_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(N_radial_pts),cctk_time-cctk_delta_time,pinitNews(:,:,ip_s),pinitdotNews(:,:,ip_s), Ylm_2)

  do patch_ID = 1, 2

     News_MoL(:,:,patch_ID) = EQ_mask * dble(initNews(:,:,patch_ID))
     dotNews_MoL(:,:,patch_ID) = EQ_mask * dble(initdotNews(:,:,patch_ID))

     News_MoL_p(:,:,patch_ID) = EQ_mask * dble(pinitNews(:,:,patch_ID))
     dotNews_MoL_p(:,:,patch_ID) = EQ_mask * dble(pinitdotNews(:,:,patch_ID))

  end do 

 end subroutine NullExact_InitMoL 
