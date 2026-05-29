! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

  ! set everything at the first few radial points

#define stereo .true.

  subroutine NullExact_Boundary(CCTK_ARGUMENTS)

    use cctk
    use NullExact_Analytic
    use NullGrid_Vars
    implicit none

    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    ! local variables, as needed:

    CCTK_INT :: k, patch_ID, i_maskn, i_masks
    CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
    CCTK_REAL :: sine
    CCTK_COMPLEX, dimension(lsh(1), lsh(2), N_radial_pts) :: drucn, drucs   

    if (verbose.ne.0) call CCTK_INFO("providing exact boundary data")

    sine = A_r0*sin(2.d0*pi*F_r0*cctk_time)

!Setting the worldtube mask, north patch    

    WT_r0(:,:,1) = rwt*(1.d0 + sine)
    x_wt(:,:,1) = WT_r0(:,:,1) / (rwt + WT_r0(:,:,1)) 
    WTexact_maskn(:,:) = int ((x_wt(:,:,1)-null_xin)/dx) + 1

    i_maskn = maxval(WTexact_maskn)
    boundary_maskn = i_maskn !WTexact_maskn    
!DEBUG
!    i_maskn = 3
!    boundary_maskn = 3
 
    patch_ID = ip_n-1

    do k = 1, min(i_maskn+3,N_radial_pts)

        call NullExact_Analytic_J_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,jcn(:,:,k), Ylm_2)
        call NullExact_Analytic_beta_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,bcn(:,:,k), Ylm_0)
        call NullExact_Analytic_U_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xbh(k),cctk_time,ucn(:,:,k), Ylm_1)
        call NullExact_Analytic_W_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,wcn(:,:,k), Ylm_0)
        
        if (first_order_scheme.ne.0) then

           call NullExact_Analytic_nu_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,nucn(:,:,k), Ylm_1)
           call NullExact_Analytic_cb_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,cbcn(:,:,k), Ylm_1)
           call NullExact_Analytic_ck_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,ckcn(:,:,k), Ylm_1)

       end if

!The radial derivative for U variable
       call NullExact_Analytic_drU_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,drucn(:,:,k), Ylm_1)
!Filling the Q variable north patch
       qcn(:,:,k) = (rwt*null_xb(k)/(1-null_xb(k)))**2&
                    *exp(-2.d0*bcn(:,:,k))*drucn(:,:,k)

    end do

!Filling the worldtube variables north patch
    j_wt(:,:,1)    = jcn(:,:,i_maskn)
    beta_wt(:,:,1) = bcn(:,:,i_maskn) 
    u_wt(:,:,1)    = ucn(:,:,i_maskn)
    w_wt(:,:,1)    = wcn(:,:,i_maskn)
    q_wt(:,:,1)    = qcn(:,:,i_maskn) 
!    q_wt(:,:,1)    = rwt**2 * drucn(:,:,i_maskn) 

    if (stereo) then

!Setting the worldtube south patch    
       WT_r0(:,:,2) = rwt*(1.d0 + sine)
       x_wt(:,:,2) = WT_r0(:,:,2) / (rwt + WT_r0(:,:,2)) 
       WTexact_masks(:,:) = int ((x_wt(:,:,2)-null_xin)/dx) + 1

       i_masks = maxval(WTexact_masks)
       boundary_masks = i_masks !WTexact_masks    

!DEBUG -- 
!       i_masks = 3
!       boundary_masks = 3

       do k = 1, min(i_masks+3,N_radial_pts)

           call NullExact_Analytic_J_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,jcs(:,:,k), Ylm_2)
           call NullExact_Analytic_beta_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,bcs(:,:,k), Ylm_0)
           call NullExact_Analytic_U_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xbh(k),cctk_time,ucs(:,:,k), Ylm_1)
           call NullExact_Analytic_W_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,wcs(:,:,k), Ylm_0)

           if (first_order_scheme.ne.0) then

              call NullExact_Analytic_nu_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,nucs(:,:,k), Ylm_1)
              call NullExact_Analytic_cb_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,cbcs(:,:,k), Ylm_1)
              call NullExact_Analytic_ck_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,ckcs(:,:,k), Ylm_1)

          end if

!The radial derivative for U variable
          call NullExact_Analytic_drU_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time,drucs(:,:,k), Ylm_1)
!Filling the Q variable south patch
          qcs(:,:,k) = (rwt*null_xb(k)/(1-null_xb(k)))**2&
                       *exp(-2.d0*bcs(:,:,k))*drucs(:,:,k)

       end do

!Filling the worldtube variables, south patch
       j_wt(:,:,2)    = jcs(:,:,i_masks)
       beta_wt(:,:,2) = bcs(:,:,i_masks) 
       u_wt(:,:,2)    = ucs(:,:,i_masks)
       w_wt(:,:,2)    = wcs(:,:,i_masks)
       q_wt(:,:,2)    = qcs(:,:,i_masks) 
!       q_wt(:,:,2)    = rwt**2 * drucs(:,:,i_masks) 

    end if

!Mask security check

!    i_min = min(i_maskn, i_masks)
!    call CCTK_ReduceLocScalar(retval, cctkGH, -1, reduce_handle,&
!                              i_min, gi_min, CCTK_VARIABLE_INT) 
!
!  if (min(minval(x_wt(:,:,1),minval(x_wt(:,:,2)).lt.null_xb(4)) then
!     call CCTK_WARN(0,"Innermost Bondi grid point must be 4 points inside the extraction radius")
!  endif
    

  end subroutine NullExact_Boundary


  subroutine NullExact_BoundaryPast(CCTK_ARGUMENTS)

    use cctk
    use NullExact_Analytic
    use NullGrid_Vars
    implicit none

    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    ! local variables, as needed:

    CCTK_INT :: k, patch_ID, i_maskn, i_masks
    CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
    CCTK_REAL :: sine
    CCTK_COMPLEX, dimension(lsh(1), lsh(2), N_radial_pts) :: drucn_p, drucs_p   

    if (verbose.ne.0) call CCTK_INFO("providing exact initial boundary data (past level)")

    sine = A_r0*sin(2.d0*pi*F_r0*cctk_time)

!Setting the worldtube mask, north patch    
    WT_r0(:,:,1) = rwt*(1.d0 + sine)
    x_wt_p(:,:,1) = WT_r0(:,:,1) / (rwt + WT_r0(:,:,1)) 
    WTexact_maskn(:,:) = int ((x_wt_p(:,:,1)-null_xin)/dx) + 1

    i_maskn = maxval(WTexact_maskn)
    boundary_maskn = i_maskn !WTexact_maskn

!DEBUG -- 
!    i_maskn = 3
!    boundary_maskn = 3

    patch_ID = ip_n-1

    do k = 1, min(i_maskn+3,N_radial_pts)

        call NullExact_Analytic_J_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time-cctk_delta_time,jcn_p(:,:,k), Ylm_2)
        call NullExact_Analytic_beta_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time-cctk_delta_time,bcn_p(:,:,k), Ylm_0)
        call NullExact_Analytic_U_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xbh(k),cctk_time-cctk_delta_time,ucn_p(:,:,k), Ylm_1)
        call NullExact_Analytic_W_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time-cctk_delta_time,wcn_p(:,:,k), Ylm_0)

        if (first_order_scheme.ne.0) then

           call NullExact_Analytic_nu_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time-cctk_delta_time,nucn_p(:,:,k), Ylm_1)
           call NullExact_Analytic_cb_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time-cctk_delta_time,cbcn_p(:,:,k), Ylm_1)
           call NullExact_Analytic_ck_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time-cctk_delta_time,ckcn_p(:,:,k), Ylm_1)

       end if
!The radial derivative for U variable
       call NullExact_Analytic_drU_2D(patch_ID+1,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time-cctk_delta_time,drucn_p(:,:,k), Ylm_1)
!Filling the Q variable north patch
       qcn_p(:,:,k) = (rwt*null_xb(k)/(1-null_xb(k)))**2&
                      *exp(-2.d0*bcn_p(:,:,k))*drucn_p(:,:,k)

    end do

!Filling the worldtube variables, north patch      
    j_wt_p(:,:,1)    = jcn_p(:,:,i_maskn)
    beta_wt_p(:,:,1) = bcn_p(:,:,i_maskn) 
    u_wt_p(:,:,1)    = ucn_p(:,:,i_maskn)
    w_wt_p(:,:,1)    = wcn_p(:,:,i_maskn)
    q_wt_p(:,:,1)    = qcn_p(:,:,i_maskn) 
!    q_wt_p(:,:,1)    = rwt**2 * drucn_p(:,:,i_maskn) 

    if (stereo) then  ! _changed_  ip_n to ip_s in function calls!

!Setting the worldtube mask, south patch    
       WT_r0(:,:,2) = rwt*(1.d0 + sine)
       x_wt_p(:,:,2) = WT_r0(:,:,2) / (rwt + WT_r0(:,:,2)) 
       WTexact_masks(:,:) = int ((x_wt_p(:,:,2)-null_xin)/dx) + 1

       i_masks = maxval(WTexact_masks)
       boundary_masks = i_masks !WTexact_masks

!DEBUG --
!       i_masks = 3
!       boundary_masks = 3

       do k = 1, min(i_masks+3,N_radial_pts)

           call NullExact_Analytic_J_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time-cctk_delta_time,jcs_p(:,:,k), Ylm_2)
           call NullExact_Analytic_beta_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time-cctk_delta_time,bcs_p(:,:,k), Ylm_0)
           call NullExact_Analytic_U_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xbh(k),cctk_time-cctk_delta_time,ucs_p(:,:,k), Ylm_1)
           call NullExact_Analytic_W_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time-cctk_delta_time,wcs_p(:,:,k), Ylm_0)

           if (first_order_scheme.ne.0) then

              call NullExact_Analytic_nu_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time-cctk_delta_time,nucs_p(:,:,k), Ylm_1)
              call NullExact_Analytic_cb_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time-cctk_delta_time,cbcs_p(:,:,k), Ylm_1)
              call NullExact_Analytic_ck_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time-cctk_delta_time,ckcs_p(:,:,k), Ylm_1)

          end if

!The radial derivative for U variable
       call NullExact_Analytic_drU_2D(ip_s,null_lsh(1),null_lsh(2),stereo_q,stereo_p,null_xb(k),cctk_time-cctk_delta_time,drucs_p(:,:,k), Ylm_1)
!Filling the Q variable south patch
       qcs_p(:,:,k) = (rwt*null_xb(k)/(1-null_xb(k)))**2&
                      *exp(-2.d0*bcs_p(:,:,k))*drucs_p(:,:,k)

       end do

!Filling the worldtube variables, south patch
      j_wt_p(:,:,2)    = jcs_p(:,:,i_masks)
      beta_wt_p(:,:,2) = bcs_p(:,:,i_masks) 
      u_wt_p(:,:,2)    = ucs_p(:,:,i_masks)
      w_wt_p(:,:,2)    = wcs_p(:,:,i_masks)
      q_wt_p(:,:,2)    = qcs_p(:,:,i_masks) 
!      q_wt_p(:,:,2)    = rwt**2 * drucs_p(:,:,i_masks) 

    end if


  end subroutine NullExact_BoundaryPast
