! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

subroutine NullEvol_InitialSlice(CCTK_ARGUMENTS)
  use NullEvol_hyper_beta
  use NullEvol_hyper_u
  use NullEvol_hyper_w
  use NullEvol_hyper_q
  use NullEvol_Evol
  use NullEvol_Mask
  use NullInterp
  use NullGrid_Vars
  implicit none

  CCTK_INT :: i, one = 1
  CCTK_INT :: mn, global_mn
  CCTK_INT :: reval
  CCTK_REAL :: SMass
  integer, save :: reduce_handle
  ! character(len=500) :: message

  logical, save :: FirstTime = .true.

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  ! Note: the extraction WT is at rb(char_nx0+1)
  ! start marching at the 2nd point

  if ( FirstTime ) then
     FirstTime = .false.
  
     call CCTK_ReductionArrayHandle(reduce_handle, "minimum");
     if (reduce_handle .lt. 0 ) then
        call CCTK_WARN(0,"Could not get reduction handle")
     endif
     SMass = .5d0*null_rwt
     ! tmp_dt = .25d0*null_start_dt/SMass
     ! dt_fact = (((tmp_dt**3/6.0d0 - tmp_dt**4/24.0d0) - &
     !        .5d0*tmp_dt**2)  + tmp_dt)

     ! it = it_start
     ! null_dt  =  1.0d300
     ! null_dt = null_start_dt
  endif

  call CCTK_INFO("Null Initial Slice");

  mn = min(minval(boundary_masks(1:null_lsh(1),1:null_lsh(2))),&
       minval(boundary_maskn(1:null_lsh(1),1:null_lsh(2))))
  call CCTK_ReduceLocScalar(reval, cctkGH, -1, reduce_handle,&
       mn, global_mn, CCTK_VARIABLE_INT)

  if (reval .ne. 0 ) then
     call CCTK_WARN(0,"Error in obtaining Minimum of mask")
  endif

  ! call CCTK_INFO(trim(message))

#define PRINTVAL(I,NN,FS) \
    if (debug_verbose.ne.0) write (*,'(a,i3,3g10.2)') NN,I,\
        maxval(abs(FS(:,:,I)),(EG_mask.ne.0).and.(guard_mask.eq.0)),\
        maxval(abs(FS(:,:,I)),guard_mask.ne.0),\
        maxval(abs(FS(:,:,I)),EG_mask.eq.0)

  ! Note that dissip_mask is set in the evolution
  ! algorithm before each step.  Therefore we set
  ! it here to ZERO, wich means no dissipation.
  dissip_mask = 0

  if (first_order_scheme.ne.0) then
     do i = 1, N_radial_pts
        call NullEvol_Set_nu (i, jcn, nucn)
        call NullEvol_Set_nu (i, jcs, nucs)
        call NullEvol_Set_ck (i, jcn, ckcn)
        call NullEvol_Set_ck (i, jcs, ckcs)
        call NullEvol_Set_cb (i, bcn, cbcn)
        call NullEvol_Set_cb (i, bcs, cbcs)
        call NullInterp_3cinterp(cctkGH, tmp_cgfn1, tmp_cgfs1, tmp_cgfn2,&
             tmp_cgfs2, tmp_cgfn3, tmp_cgfs3, ckcn(:,:,i), ckcs(:,:,i),&
             nucn(:,:,i), nucs(:,:,i), cbcn(:,:,i), cbcs(:,:,i), one, one, one)
     end do
  end if

  do i = max(3,global_mn-1), N_radial_pts  

     ! If 'B' is the index of the point just inside the WT, 
     ! the hpyersurface equations will only need special algorithm
     ! at the point B+1.


     call NullEvol_remask (i, 1+boundary_masks, evolution_masks)
     call NullEvol_remask (i, 1+boundary_maskn, evolution_maskn)

     if (DEBUG_skip_B_update.eq.0) then

        call NullEvol_beta (i, jcn, bcn, j_wt(:,:,1), &
                            beta_wt(:,:,1), x_wt(:,:,1), evolution_maskn)
        call NullEvol_beta (i, jcs, bcs, j_wt(:,:,2), &
                            beta_wt(:,:,2), x_wt(:,:,2), evolution_masks)

        PRINTVAL(i,"Bn",bcn)
     end if


     if (first_order_scheme.ne.0) then

        call NullEvol_nu (i, jcn, nucn, evolution_maskn)
        call NullEvol_nu (i, jcs, nucs, evolution_masks)

        call NullEvol_ck (i, jcn, ckcn, evolution_maskn)
        call NullEvol_ck (i, jcs, ckcs, evolution_masks)

        call NullEvol_cb (i, bcn, cbcn, evolution_maskn)
        call NullEvol_cb (i, bcs, cbcs, evolution_masks)

        call NullInterp_3cinterp(cctkGH, tmp_cgfn1, tmp_cgfs1, tmp_cgfn2,&
             tmp_cgfs2, tmp_cgfn3, tmp_cgfs3, ckcn(:,:,i), ckcs(:,:,i),&
             nucn(:,:,i), nucs(:,:,i), cbcn(:,:,i), cbcs(:,:,i), one, one, one)

     end if

     if (DEBUG_skip_Q_update.eq.0) then

        call NullEvol_q (i, jcn, nucn, ckcn, bcn, cbcn, qcn, &
                         q_wt(:,:,1), j_wt(:,:,1), beta_wt(:,:,1), x_wt(:,:,1), &
                         evolution_maskn, eth4_maskn, dissip_mask, 0.0d0, first_order_scheme)

        call NullEvol_q (i, jcs, nucs, ckcs, bcs, cbcs, qcs, &
                         q_wt(:,:,2), j_wt(:,:,2), beta_wt(:,:,2), x_wt(:,:,2), &
                         evolution_masks, eth4_masks, dissip_mask, 0.0d0, first_order_scheme)

        call NullInterp_cinterp(cctkGH, tmp_cgfn, tmp_cgfs, qcn(:,:,i),qcs(:,:,i), one)

        PRINTVAL(i,"Qn",qcn)
     end if

     if (DEBUG_skip_U_update.eq.0) then

        ! note that there are no angular derivatives in the update of U
        ! NullEvol_u (i, jns, bns, qns, uns, u_x, u_wt, x_wt, mask)
        call NullEvol_u (i, jcn, bcn, qcn, ucn, u_wt(:,:,1), u_x_wt(:,:,1),&
                         x_wt(:,:,1), evolution_maskn)
        call NullEvol_u (i, jcs, bcs, qcs, ucs, u_wt(:,:,2), u_x_wt(:,:,2),&
                         x_wt(:,:,2), evolution_masks)

        PRINTVAL(i,"Un",ucn)
     end if

     if (DEBUG_skip_W_update.eq.0) then

        call NullEvol_w (i, jcn, nucn, ckcn, bcn, cbcn, ucn, qcn, wcn,&
             w_wt(:,:,1), u_wt(:,:,1), x_wt(:,:,1), &
             evolution_maskn, eth4_maskn, dissip_mask, 0.0d0,first_order_scheme)
        call NullEvol_w (i, jcs, nucs, ckcs, bcs, cbcs, ucs, qcs, wcs,&
             w_wt(:,:,2), u_wt(:,:,2), x_wt(:,:,2), &
             evolution_masks, eth4_masks, dissip_mask, 0.0d0,first_order_scheme)

        call NullInterp_rinterp(cctkGH, tmp_rgfn, tmp_rgfs, wcn(:,:,i), wcs(:,:,i))

        PRINTVAL(i,"Wn",wcn)
     end if
  end do

  call NullEvol_ResetInactive(CCTK_PASS_FTOF,1,N_radial_pts)
  
  write (*,*) 'The Boundary mask is at: ', min(minval(boundary_maskn),minval(boundary_masks))
  write (*,*) 'The Evolution mask is at: ', min(minval(evolution_maskn), minval(evolution_masks))

end subroutine NullEvol_InitialSlice

