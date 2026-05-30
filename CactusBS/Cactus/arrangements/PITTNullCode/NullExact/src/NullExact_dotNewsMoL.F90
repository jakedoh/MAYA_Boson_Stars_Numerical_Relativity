! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine NullExact_dotNewsMoL(CCTK_ARGUMENTS)
 use cctk
 use NullInterp
 use NullGrid_Vars
 use NullExact_Analytic
 implicit none

  CCTK_INT :: patch_ID
  CCTK_COMPLEX, dimension(lsh(1), lsh(2), 2) :: anaNewsB, anaPsi4

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

!INTEGRATION Int_Psi4 = News
  patch_ID = ip_n-1

     call NullExact_Analytic_BondiNews_2D(patch_ID+1,null_lsh(1),null_lsh(2),&
          stereo_q,stereo_p,null_xb(N_radial_pts),&
          cctk_time,anaNewsB(:,:,patch_ID+1),anaPsi4(:,:,patch_ID+1), Ylm_2)
     call NullExact_Analytic_BondiNews_2D(ip_s,null_lsh(1),null_lsh(2),&
          stereo_q,stereo_p,null_xb(N_radial_pts),cctk_time,&
          anaNewsB(:,:,ip_s),anaPsi4(:,:,ip_s), Ylm_2)

  do patch_ID = 1, 2

     dotNews_MoL(:,:,patch_ID) = EQ_mask * dble(Psi4(:,:,patch_ID))

!the error between the analytic Bondi News evolved and the integrated Psi4
     errNews_MoL(:,:,patch_ID) = EQ_mask * (dble(News_MoL(:,:,patch_ID) - anaNewsB(:,:,patch_ID)))

  end do

end subroutine NullExact_dotNewsMoL 
