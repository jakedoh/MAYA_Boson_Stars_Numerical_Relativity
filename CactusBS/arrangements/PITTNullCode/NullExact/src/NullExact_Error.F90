! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"


subroutine NullExact_Error_Cmplx(CCTK_ARGUMENTS)

  use cctk
  use NullExact_Analytic
  use NullGrid_Vars
  use NullInterp
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT, parameter :: izero = 0
  integer, parameter :: ik = kind(izero)

  ! local variables, as needed:
  CCTK_INT :: k, patch_ID
  CCTK_COMPLEX, dimension(lsh(1), lsh(2), N_radial_pts) :: anajcn, anajcs, anaucn, anaucs, anadrucn, anadrucs, anaeth2jcn, anaeth2jcs, anaqcn, anaqcs


  if (verbose.ne.0) call CCTK_INFO("checking for error in main complex evolution variables")
  !     write (*,*) "error checking time:", cctk_time

  patch_ID = ip_n-1

  do k = 1, N_radial_pts


!Compute analytic values in complex Characteristic Fields 
     call NullExact_Analytic_J_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,anajcn(:,:,k), Ylm_2)
     call NullExact_Analytic_U_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xbh(k),cctk_time,anaucn(:,:,k), Ylm_1)
     call NullInterp_d2 (anaeth2jcn(:,:,k), anajcn(:,:,k), 2_ik, 1_ik, 1_ik)
     call NullExact_Analytic_drU_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,anadrucn(:,:,k), Ylm_1)

     call NullExact_Analytic_J_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,anajcs(:,:,k), Ylm_2)
     call NullExact_Analytic_U_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xbh(k),cctk_time,anaucs(:,:,k), Ylm_1)
     call NullInterp_d2 (anaeth2jcs(:,:,k), anajcs(:,:,k), 2_ik, 1_ik, 1_ik)
     call NullExact_Analytic_drU_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,anadrucs(:,:,k), Ylm_1)



     jcnr(:,:,k)   = EQ_mask * dble(jcn(:,:,k))
     jcni(:,:,k)   = EQ_mask * dimag(jcn(:,:,k))

!Compute errors in complex Characteristic Fields 
     jcnr_e(:,:,k) = EQ_mask * dble(jcn(:,:,k) - anajcn(:,:,k))
     jcni_e(:,:,k) = EQ_mask * dimag(jcn(:,:,k) - anajcn(:,:,k))

     jcn_e(:,:,k) = EQ_mask * (jcn(:,:,k) - anajcn(:,:,k))
     jcs_e(:,:,k) = EQ_mask * (jcs(:,:,k) - anajcs(:,:,k))

     ucn_e(:,:,k) = EQ_mask * (ucn(:,:,k) - anaucn(:,:,k))
     ucs_e(:,:,k) = EQ_mask * (ucs(:,:,k) - anaucs(:,:,k))

     eth2jcn_e(:,:,k) = EQ_mask * (eth2jcn(:,:,k) - anaeth2jcn(:,:,k))
     eth2jcs_e(:,:,k) = EQ_mask * (eth2jcs(:,:,k) - anaeth2jcs(:,:,k))

     anaqcn(:,:,k) = (rwt*null_xb(k)/(1-null_xb(k)))**2&
                      *exp(-2.d0*bcn(:,:,k))*anadrucn(:,:,k)
     anaqcs(:,:,k) = (rwt*null_xb(k)/(1-null_xb(k)))**2&
                      *exp(-2.d0*bcs(:,:,k))*anadrucs(:,:,k)

     qcn_e(:,:,k) = EQ_mask * (qcn(:,:,k) - anaqcn(:,:,k))
     qcs_e(:,:,k) = EQ_mask * (qcs(:,:,k) - anaqcs(:,:,k))


!Compute masked errors in complex Characteristic Fields
    if (CCTK_EQUALS(error_mask_type, "none")) then

       Mjcn_e(:,:,k) = (jcn(:,:,k) - anajcn(:,:,k))
       Mucn_e(:,:,k) = (ucn(:,:,k) - anaucn(:,:,k))
       Meth2jcn_e(:,:,k) = (eth2jcn(:,:,k) - anaeth2jcn(:,:,k))

    else if (CCTK_EQUALS(error_mask_type, "EG_mask")) then

       Mjcn_e(:,:,k) = EG_mask * (jcn(:,:,k) - anajcn(:,:,k))
       Mucn_e(:,:,k) = EG_mask * (ucn(:,:,k) - anaucn(:,:,k))
       Meth2jcn_e(:,:,k) = EG_mask * (eth2jcn(:,:,k) - anaeth2jcn(:,:,k))

    else if (CCTK_EQUALS(error_mask_type, "EQ_mask")) then

       Mjcn_e(:,:,k) = EQ_mask * (jcn(:,:,k) - anajcn(:,:,k))
       Mucn_e(:,:,k) = EQ_mask * (ucn(:,:,k) - anaucn(:,:,k))
       Meth2jcn_e(:,:,k) = EQ_mask * (eth2jcn(:,:,k) - anaeth2jcn(:,:,k))

    else if (CCTK_EQUALS(error_mask_type, "EV_mask")) then

       Mjcn_e(:,:,k) = EV_mask * (jcn(:,:,k) - anajcn(:,:,k))
       Mucn_e(:,:,k) = EV_mask * (ucn(:,:,k) - anaucn(:,:,k))
       Meth2jcn_e(:,:,k) = EV_mask * (eth2jcn(:,:,k) - anaeth2jcn(:,:,k))

    else if (CCTK_EQUALS(error_mask_type, "guard_mask")) then

       Mjcn_e(:,:,k) = guard_mask * (jcn(:,:,k) - anajcn(:,:,k))
       Mucn_e(:,:,k) = guard_mask * (ucn(:,:,k) - anaucn(:,:,k))
       Meth2jcn_e(:,:,k) = guard_mask * (eth2jcn(:,:,k) - anaeth2jcn(:,:,k))

    else
      call CCTK_WARN(0, "unsupported circular mask type")
    end if

  end do

end subroutine NullExact_Error_Cmplx

subroutine NullExact_Error_Real(CCTK_ARGUMENTS)

  use cctk
  use NullExact_Analytic
  use NullGrid_Vars
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  ! local variables, as needed:

  CCTK_INT :: k, patch_ID

  if (verbose.ne.0) call CCTK_INFO("checking for error in real evolution variables")

  patch_ID = ip_n-1

  do k = 1, N_radial_pts

     call NullExact_Analytic_beta_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,bcn_e(:,:,k), Ylm_0)
     call NullExact_Analytic_W_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,wcn_e(:,:,k), Ylm_0)


     call NullExact_Analytic_beta_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,bcs_e(:,:,k), Ylm_0)
     call NullExact_Analytic_W_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,wcs_e(:,:,k), Ylm_0)

     bcn_e(:,:,k) = EG_mask * (bcn(:,:,k) - bcn_e(:,:,k))
     wcn_e(:,:,k) = EG_mask * (wcn(:,:,k) - wcn_e(:,:,k))

     bcs_e(:,:,k) = EG_mask * (bcs(:,:,k) - bcs_e(:,:,k))
     wcs_e(:,:,k) = EG_mask * (wcs(:,:,k) - wcs_e(:,:,k))

  end do

end subroutine NullExact_Error_Real

subroutine NullExact_Error_CmplxAux(CCTK_ARGUMENTS)

  use cctk
  use NullExact_Analytic
  use NullGrid_Vars
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  ! local variables, as needed:

  CCTK_INT :: k, patch_ID

  if (first_order_scheme.eq.0) then
     if (cctk_iteration.lt.100) call CCTK_WARN(1, "cannot check for error in complex auxiliary variables with NullEvolve::first_order_scheme==0")
     return
  end if

  patch_ID = ip_n-1

  if (verbose.ne.0) call CCTK_INFO("checking for error in auxiliary complex evolution variables")

  do k = 1, N_radial_pts

     call NullExact_Analytic_nu_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,nucn_e(:,:,k), Ylm_1)
     call NullExact_Analytic_ck_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,ckcn_e(:,:,k), Ylm_1)
     call NullExact_Analytic_cb_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,cbcn_e(:,:,k), Ylm_1)

     call NullExact_Analytic_nu_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,nucs_e(:,:,k), Ylm_1)
     call NullExact_Analytic_ck_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,ckcs_e(:,:,k), Ylm_1)
     call NullExact_Analytic_cb_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,cbcs_e(:,:,k), Ylm_1)

     nucn_e(:,:,k) = EG_mask * (nucn(:,:,k) - nucn_e(:,:,k))
     ckcn_e(:,:,k) = EG_mask * (ckcn(:,:,k) - ckcn_e(:,:,k))
     cbcn_e(:,:,k) = EG_mask * (cbcn(:,:,k) - cbcn_e(:,:,k))

     nucs_e(:,:,k) = EG_mask * (nucs(:,:,k) - nucs_e(:,:,k))
     ckcs_e(:,:,k) = EG_mask * (ckcs(:,:,k) - ckcs_e(:,:,k))
     cbcs_e(:,:,k) = EG_mask * (cbcs(:,:,k) - cbcs_e(:,:,k))

  end do

end subroutine NullExact_Error_CmplxAux

subroutine NullExact_Error_NewsB(CCTK_ARGUMENTS)

  use cctk
  use NullExact_Analytic
  use NullGrid_Vars
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL, dimension(lsh(1), lsh(2), 2) :: anamu
  CCTK_COMPLEX, dimension(lsh(1), lsh(2), 2) :: anaJl, numJl
  CCTK_INT :: patch_ID
  CCTK_REAL :: fct1, fct2

  patch_ID = ip_n-1

  if (verbose.ne.0) call CCTK_INFO("checking for error in the Bondi News variables")

  if (linearized_news .eq. 0) then
     call NullExact_Analytic_BondiNews_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(N_radial_pts),cctk_time,NewsB_e(:,:,patch_ID+1),Psi4_e(:,:,patch_ID+1), Ylm_2)
     call NullExact_Analytic_BondiTime_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(N_radial_pts),cctk_time,uBondi_e(:,:,patch_ID+1), Ylm_0)
  else
     call NullExact_Analytic_BondiNews_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(N_radial_pts),time_of_news,NewsB_e(:,:,patch_ID+1),Psi4_e(:,:,patch_ID+1), Ylm_2)
  endif

  if (linearized_news .eq. 0) then
     call NullExact_Analytic_BondiNews_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(N_radial_pts),cctk_time,NewsB_e(:,:,ip_s),Psi4_e(:,:,ip_s), Ylm_2)
     call NullExact_Analytic_BondiTime_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(N_radial_pts),cctk_time,uBondi_e(:,:,ip_s), Ylm_0)
  else
     call NullExact_Analytic_BondiNews_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(N_radial_pts),time_of_news,NewsB_e(:,:,ip_s),Psi4_e(:,:,ip_s), Ylm_2)
  endif

!compute analytic mu = (omega -1), only for dynamic minkowsky, l=2,m=0.
     fct1 = 1.D0/375000.D0  
     anamu(:,:,patch_ID+1) = dble(fct1*sin(cctk_time)*Ylm_0(:,:,patch_ID+1))
     anamu(:,:,ip_s) = dble(fct1*sin(cctk_time)*Ylm_0(:,:,ip_s))

!compute analytic J_l, only for l=2, m=0 dynamic minkowsky
     fct2 = 3.D0*dsqrt(6.D0)/2000000.D0  
     anaJl(:,:,patch_ID+1) = fct2*Ylm_2(:,:,patch_ID+1)*cos(cctk_time) 
     anaJl(:,:,ip_s) = fct2*Ylm_2(:,:,ip_s)*cos(cctk_time)

!compute numeric J_l,using ananlytic J 
!  do k = 1, N_radial_pts
!     call NullExact_Analytic_J_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,anajcn(:,:,k), Ylm_2)
!     call NullExact_Analytic_J_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,anajcs(:,:,k), Ylm_2)
!  end do 
     numJl(:,:,patch_ID+1) = - 0.5d0 * ( 3.0d0 * jcn(:,:,N_radial_pts) &
    - 4.0d0 * jcn(:,:,N_radial_pts-1) + jcn(:,:,N_radial_pts-2) ) / (null_dx) * rwt
     numJl(:,:,patch_ID+1) = - 0.5d0 * ( 3.0d0 * jcs(:,:,N_radial_pts) &
    - 4.0d0 * jcs(:,:,N_radial_pts-1) + jcs(:,:,N_radial_pts-2) ) / (null_dx) * rwt


  do patch_ID = 1, 2

     NewsBr(:,:,patch_ID) = EQ_mask * dble(NewsB(:,:,patch_ID))
     NewsBi(:,:,patch_ID) = EQ_mask * dimag(NewsB(:,:,patch_ID))

     !    In linearized regime, News_exact should be NewsB_exact
     Newsr_e(:,:,patch_ID) = EQ_mask * dble(News(:,:,patch_ID)-NewsB_e(:,:,patch_ID))
     Newsi_e(:,:,patch_ID) = EQ_mask * dimag(News(:,:,patch_ID)-NewsB_e(:,:,patch_ID))

     NewsB_e (:,:,patch_ID) = EQ_mask * (NewsB(:,:,patch_ID) - NewsB_e(:,:,patch_ID))
     NewsBr_e(:,:,patch_ID) = EQ_mask * dble(NewsB_e(:,:,patch_ID))
     NewsBi_e(:,:,patch_ID) = EQ_mask * dimag(NewsB_e(:,:,patch_ID))

     uBondi_e(:,:,patch_ID) = EQ_mask * (uBondi(:,:,patch_ID) - uBondi_e(:,:,patch_ID))
     Psi4r_e(:,:,patch_ID) = EQ_mask * dble(Psi4(:,:,patch_ID)-Psi4_e(:,:,patch_ID))
     Psi4i_e(:,:,patch_ID) = EQ_mask * dimag(Psi4(:,:,patch_ID)-Psi4_e(:,:,patch_ID))

 
!compute mu = omega -1, and error in mu
     muh(:,:,patch_ID) = dble((omegao(:,:,patch_ID)+omegan(:,:,patch_ID))/2.0d0 - 1.0d0)
     EQmuh(:,:,patch_ID) = EQ_mask * dble(muh(:,:,patch_ID))
     mu_e(:,:,patch_ID) = EQ_mask * dble(muh(:,:,patch_ID) - anamu(:,:,patch_ID))
     EGmu_e(:,:,patch_ID) = EG_mask * dble(muh(:,:,patch_ID) - anamu(:,:,patch_ID))

!compute EQJ_l, and error in J_l

     Jh_l(:,:,patch_ID) = (Jo_l(:,:,patch_ID) + Jn_l(:,:,patch_ID))/2.0d0
     EQJh_l(:,:,patch_ID) = EQ_mask * Jh_l(:,:,patch_ID)
     aJ_l_e(:,:,patch_ID) = EQ_mask*(Jh_l(:,:,patch_ID) - anaJl(:,:,patch_ID))
     EGaJ_l_e(:,:,patch_ID) = EG_mask*(Jh_l(:,:,patch_ID) - anaJl(:,:,patch_ID))
     nJ_l_e(:,:,patch_ID) = EQ_mask*(Jh_l(:,:,patch_ID) - numJl(:,:,patch_ID))
     EGnJ_l_e(:,:,patch_ID) = EG_mask*(Jh_l(:,:,patch_ID) - numJl(:,:,patch_ID))

  end do

end subroutine NullExact_Error_NewsB

!=====================================================================
subroutine NullExact_Error_Constr(CCTK_ARGUMENTS)

  use cctk
  use NullExact_Analytic
  use NullGrid_Vars
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  ! local variables, as needed:

  CCTK_INT :: k, patch_ID

  if (verbose.ne.0) call CCTK_INFO("checking for error in constraints")

  patch_ID = ip_n-1


  do k = 1, N_radial_pts

     call NullExact_Analytic_R01_2D(patch_ID+1,null_lsh(1),null_lsh(2), &
          stereo_q,stereo_p,null_xb(k),cctk_time,R01n_e(:,:,k), Ylm_0)

     call NullExact_Analytic_R01_2D(ip_s,null_lsh(1),null_lsh(2), &
          stereo_q,stereo_p,null_xb(k),cctk_time,R01s_e(:,:,k), Ylm_0)

     R01n_e(:,:,k) = EG_mask * (Null_R01(:,:,k) - R01n_e(:,:,k))
     R01s_e(:,:,k) = EG_mask * (Null_R01_south(:,:,k) - R01s_e(:,:,k))

  end do

end subroutine NullExact_Error_Constr
