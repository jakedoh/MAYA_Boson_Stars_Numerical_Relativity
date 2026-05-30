! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"


subroutine NullExact_AnaFields(CCTK_ARGUMENTS)

  use cctk
  use NullExact_Analytic
  use NullGrid_Vars
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  ! local variables, as needed:
  CCTK_INT :: k, patch_ID
  CCTK_REAL :: fct1, fct2
  CCTK_COMPLEX, dimension(lsh(1), lsh(2), 2) :: anaNews, anaPsi4, anaJ_l, numJ_l
  CCTK_COMPLEX, dimension(lsh(1), lsh(2), N_radial_pts) :: anajcn, anajcs, anaucn, anaucs, anadrucn, anadrucs
  CCTK_REAL, dimension(lsh(1), lsh(2), 2) :: anamu, anabcn, anabcs


  if (verbose.ne.0) call CCTK_INFO("analytic values for main complex evolution variables")

  patch_ID = ip_n-1

  do k = 1, N_radial_pts


     call NullExact_Analytic_J_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,anajcn(:,:,k), Ylm_2)

     call NullExact_Analytic_J_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,anajcs(:,:,k), Ylm_2)

     call NullExact_Analytic_U_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xbh(k),cctk_time,anaucn(:,:,k), Ylm_1)

     call NullExact_Analytic_U_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xbh(k),cctk_time,anaucs(:,:,k), Ylm_1)

     call NullExact_Analytic_drU_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,anadrucn(:,:,k), Ylm_1)

     call NullExact_Analytic_drU_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,anadrucs(:,:,k), Ylm_1)

        call NullExact_Analytic_beta_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,anabcn(:,:,k), Ylm_0)

        call NullExact_Analytic_beta_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,anabcs(:,:,k), Ylm_0)

     ana_qcn(:,:,k) = EQ_mask*(rwt*null_xb(k)/(1-null_xb(k)))**2&
                      *exp(-2.d0*anabcn(:,:,k))*anadrucn(:,:,k)
     ana_qcs(:,:,k) = EQ_mask*(rwt*null_xb(k)/(1-null_xb(k)))**2&
                      *exp(-2.d0*anabcs(:,:,k))*anadrucs(:,:,k)

     ana_jcn(:,:,k) = EQ_mask * (anajcn(:,:,k))
     ana_jcs(:,:,k) = EQ_mask * (anajcs(:,:,k))

     ana_ucn(:,:,k) = EQ_mask * (anaucn(:,:,k))
     ana_ucs(:,:,k) = EQ_mask * (anaucs(:,:,k))

  end do

     call NullExact_Analytic_BondiNews_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(N_radial_pts),cctk_time,anaNews(:,:,patch_ID+1),anaPsi4(:,:,patch_ID+1), Ylm_2)
     call NullExact_Analytic_BondiNews_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(N_radial_pts),cctk_time,anaNews(:,:,ip_s),anaPsi4(:,:,ip_s), Ylm_2)

!compute analytic mu = omega -1, only for l=2, m=0 dynamic minkowsky
     fct1 = 1.D0/375000.D0  
     anamu(:,:,patch_ID+1) = dble(fct1*sin(cctk_time)*Ylm_0(:,:,patch_ID+1))
     anamu(:,:,ip_s) = dble(fct1*sin(cctk_time)*Ylm_0(:,:,ip_s))

!compute analytic J_l, only for l=2, m=0 dynamic minkowsky
     fct2 = 3.D0*dsqrt(6.D0)/2000000.D0  
     anaJ_l(:,:,patch_ID+1) = fct2*Ylm_2(:,:,patch_ID+1)*cos(cctk_time) 
     anaJ_l(:,:,ip_s) = fct2*Ylm_2(:,:,ip_s)*cos(cctk_time)

!compute numeric J_l, by one sided approx, using analytic J 
     numJ_l(:,:,patch_ID+1) = - 0.5d0 * ( 3.0d0 * anajcn(:,:,N_radial_pts) &
- 4.0d0 * anajcn(:,:,N_radial_pts-1) + anajcn(:,:,N_radial_pts-2) ) / (null_dx) * rwt

     numJ_l(:,:,ip_s) = - 0.5d0 * ( 3.0d0 * anajcs(:,:,N_radial_pts) &
- 4.0d0 * anajcs(:,:,N_radial_pts-1) + anajcs(:,:,N_radial_pts-2) ) / (null_dx ) * rwt


  do patch_ID = 1, 2

     ana_Psi4(:,:,patch_ID) = EQ_mask * (anaPsi4(:,:,patch_ID))
     ana_News(:,:,patch_ID) = EQ_mask * (anaNews(:,:,patch_ID))
     ana_mu(:,:,patch_ID) = EQ_mask * dble(anamu(:,:,patch_ID))
     ana_J_l(:,:,patch_ID) = EQ_mask * (anaJ_l(:,:,patch_ID))     
     num_J_l(:,:,patch_ID) = EQ_mask * (numJ_l(:,:,patch_ID))     


  end do 

end subroutine NullExact_AnaFields
