! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#define stereo .true.
#define full_initial_data .false.

  ! set J on the entire characteristic grid

  subroutine NullExact_Initial(CCTK_ARGUMENTS)

    use cctk
    use NullExact_Analytic
    use NullGrid_Vars
    implicit none

    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    DECLARE_CCTK_FUNCTIONS

    ! local variables, as needed:

    CCTK_INT :: k , patch_ID
    

    if (verbose.ne.0) call CCTK_INFO("providing exact initial data for J")
    !write (*,*) 'initial data time:', cctk_time
    
    patch_ID = ip_n-1

    do k = 1, N_radial_pts

        call NullExact_Analytic_J_2D(patch_ID+1,null_lsh(1),null_lsh(2),&
                      stereo_q,stereo_p,null_xb(k),cctk_time,jcn(:,:,k), Ylm_2)
        call NullExact_AnalytiC_J_2D(patch_ID+1,null_lsh(1),null_lsh(2),&
                      stereo_q,stereo_p,null_xb(k),cctk_time-cctk_delta_time,jcn_p(:,:,k), Ylm_2)

        if (stereo) then
             call NullExact_Analytic_J_2D(ip_s,null_lsh(1),null_lsh(2),&
                        stereo_q,stereo_p,null_xb(k),cctk_time,jcs(:,:,k), Ylm_2)
             call NullExact_Analytic_J_2D(ip_s,null_lsh(1),null_lsh(2),&
                        stereo_q,stereo_p,null_xb(k),cctk_time-cctk_delta_time,jcs_p(:,:,k), Ylm_2)
        endif

        if (full_initial_data) then

            call NullExact_Analytic_U_2D(patch_ID+1,null_lsh(1),null_lsh(2),&
                           stereo_q,stereo_p,null_xbh(k),cctk_time,ucn(:,:,k), Ylm_1)
            call NullExact_Analytic_beta_2D(patch_ID+1,null_lsh(1),null_lsh(2),&
                           stereo_q,stereo_p,null_xb(k),cctk_time,bcn(:,:,k), Ylm_0)
            call NullExact_Analytic_W_2D(patch_ID+1,null_lsh(1),null_lsh(2),&
                           stereo_q,stereo_p,null_xb(k),cctk_time,wcn(:,:,k), Ylm_0)
   
            call NullExact_Analytic_U_2D(patch_ID+1,null_lsh(1),null_lsh(2),&
                        stereo_q,stereo_p,null_xbh(k),cctk_time-cctk_delta_time,ucn_p(:,:,k), Ylm_1)
            call NullExact_Analytic_beta_2D(patch_ID+1,null_lsh(1),null_lsh(2),&
                        stereo_q,stereo_p,null_xb(k),cctk_time-cctk_delta_time,bcn_p(:,:,k), Ylm_0)
            call NullExact_Analytic_W_2D(patch_ID+1,null_lsh(1),null_lsh(2),&
                        stereo_q,stereo_p,null_xb(k),cctk_time-cctk_delta_time,wcn_p(:,:,k), Ylm_0)
   
            if (stereo) then
   
               call NullExact_Analytic_U_2D(ip_s,null_lsh(1),null_lsh(2),&
                        stereo_q,stereo_p,null_xbh(k),cctk_time,ucs(:,:,k), Ylm_1)
               call NullExact_Analytic_beta_2D(ip_s,null_lsh(1),null_lsh(2),&
                        stereo_q,stereo_p,null_xb(k),cctk_time,bcs(:,:,k), Ylm_0)
               call NullExact_Analytic_W_2D(ip_s,null_lsh(1),null_lsh(2),&
                        stereo_q,stereo_p,null_xb(k),cctk_time,wcs(:,:,k), Ylm_0)
   
               call NullExact_Analytic_U_2D(ip_s,null_lsh(1),null_lsh(2),&
                        stereo_q,stereo_p,null_xbh(k),cctk_time-cctk_delta_time,ucs_p(:,:,k), Ylm_1)
               call NullExact_Analytic_beta_2D(ip_s,null_lsh(1),null_lsh(2),&
                        stereo_q,stereo_p,null_xb(k),cctk_time-cctk_delta_time,bcs_p(:,:,k), Ylm_0)
               call NullExact_Analytic_W_2D(ip_s,null_lsh(1),null_lsh(2),&
                        stereo_q,stereo_p,null_xb(k),cctk_time-cctk_delta_time,wcs_p(:,:,k), Ylm_0)
   
            end if
   
            if (first_order_scheme.ne.0) then
   
               call NullExact_Analytic_nu_2D(patch_ID+1,null_lsh(1),null_lsh(2),&
                        stereo_q,stereo_p,null_xb(k),cctk_time,nucn(:,:,k), Ylm_1)
               call NullExact_Analytic_ck_2D(patch_ID+1,null_lsh(1),null_lsh(2),&
                        stereo_q,stereo_p,null_xb(k),cctk_time,ckcn(:,:,k), Ylm_1)
               call NullExact_Analytic_cb_2D(patch_ID+1,null_lsh(1),null_lsh(2),&
                        stereo_q,stereo_p,null_xb(k),cctk_time,cbcn(:,:,k), Ylm_1)
   
               call NullExact_Analytic_nu_2D(patch_ID+1,null_lsh(1),null_lsh(2),&
                        stereo_q,stereo_p,null_xb(k),cctk_time-cctk_delta_time,nucn_p(:,:,k), Ylm_1)
               call NullExact_Analytic_ck_2D(patch_ID+1,null_lsh(1),null_lsh(2),&
                        stereo_q,stereo_p,null_xb(k),cctk_time-cctk_delta_time,ckcn_p(:,:,k), Ylm_1)
               call NullExact_Analytic_cb_2D(patch_ID+1,null_lsh(1),null_lsh(2),&
                        stereo_q,stereo_p,null_xb(k),cctk_time-cctk_delta_time,cbcn_p(:,:,k), Ylm_1)
   
               if (stereo) then
   
                  call NullExact_Analytic_nu_2D(ip_s,null_lsh(1),null_lsh(2),&
                        stereo_q,stereo_p,null_xb(k),cctk_time,nucs(:,:,k), Ylm_1)
                  call NullExact_Analytic_ck_2D(ip_s,null_lsh(1),null_lsh(2),&
                        stereo_q,stereo_p,null_xb(k),cctk_time,ckcs(:,:,k), Ylm_1)
                  call NullExact_Analytic_cb_2D(ip_s,null_lsh(1),null_lsh(2),&
                        stereo_q,stereo_p,null_xb(k),cctk_time,cbcs(:,:,k), Ylm_1)
   
                  call NullExact_Analytic_nu_2D(ip_s,null_lsh(1),null_lsh(2),&
                        stereo_q,stereo_p,null_xb(k),cctk_time-cctk_delta_time,nucs_p(:,:,k), Ylm_1)
                  call NullExact_Analytic_ck_2D(ip_s,null_lsh(1),null_lsh(2),&
                        stereo_q,stereo_p,null_xb(k),cctk_time-cctk_delta_time,ckcs_p(:,:,k), Ylm_1)
                  call NullExact_Analytic_cb_2D(ip_s,null_lsh(1),null_lsh(2),&
                        stereo_q,stereo_p,null_xb(k),cctk_time-cctk_delta_time,cbcs_p(:,:,k), Ylm_1)
   
               end if
            end if
        end if

     end do

  end subroutine NullExact_Initial
