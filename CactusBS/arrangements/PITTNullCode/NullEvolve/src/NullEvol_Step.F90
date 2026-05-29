! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

subroutine NullEvol_Step(CCTK_ARGUMENTS)
  use NullEvol_hyper_beta
  use NullEvol_hyper_u
  use NullEvol_hyper_q
  use NullEvol_hyper_w
  use NullEvol_Evol
  use NullEvol_Mask
  use NullEvol_DissipMask
  use NullInterp
  use NullGrid_Vars 
  implicit none

  CCTK_INT, parameter :: izero = 0
  integer, parameter :: ik = kind(izero)

  CCTK_INT :: i, qp, pq, one = 1, two = 2
  CCTK_INT :: mn, global_mn
  CCTK_INT :: reval
  integer, save :: reduce_handle


  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL :: dissip_eps_W, dissip_eps_Q, dissip_eps_J, dissip_eps_Jx

  ! Note: the extraction WT is at rb(char_nx0+1)
  ! start marching at the 2nd point

  call CCTK_ReductionArrayHandle(reduce_handle, "minimum");
  if (reduce_handle .lt. 0 ) then
     call CCTK_WARN(0,"Could not get reduction handle")
  endif

  !   call CCTK_INFO("Null Evolution");

  mn = min(minval(boundary_masks(1:null_lsh(1),1:null_lsh(2))),&
       minval(boundary_maskn(1:null_lsh(1),1:null_lsh(2))))
  call CCTK_ReduceLocScalar(reval, cctkGH, -1, reduce_handle,&
       mn, global_mn, CCTK_VARIABLE_INT)

  if (reval .ne. 0 ) then
     call CCTK_WARN(0,"Error in obtaining Minimum of mask")
  endif

  call NullEvol_set_dissipmask(dissip_mask)
  dissip_eps_Q = dissip_Q * (null_delta(1) * null_delta(2))**2 / null_dx
  dissip_eps_W = dissip_W * (null_delta(1) * null_delta(2))**2 / null_dx
  dissip_eps_Jx = dissip_Jx * (null_delta(1) * null_delta(2))**2 / null_dx
  dissip_eps_J = dissip_J * (null_delta(1) * null_delta(2))**2 / cctk_delta_time

  call NullEvol_ResetInactive(CCTK_PASS_FTOF,1_ik,global_mn)

  do i = global_mn+1, N_radial_pts  

     ! Note that the J evolution algorithm has a modified stencil
     ! for the first TWO points next to the inner boundary.
     ! All other equations have a smaller stencil, needing
     ! special algorithm at the first point next the boundary but
     ! not at the next one.
     ! If 'B' is the index of the point just inside the WT, 
     ! J will have special algorithm at B+1 and B+2 , but
     ! the hpyersurface equations will only need special algorithm
     ! at the point B+1.

     call NullEvol_remask (i, 1+boundary_masks, evolution_masks)
     call NullEvol_remask (i, 1+boundary_maskn, evolution_maskn)

     call NullEvol_remask (i, 2+boundary_masks, auxiliary_masks)
     call NullEvol_remask (i, 2+boundary_maskn, auxiliary_maskn)

     call NullEvol_remask (i, 4+boundary_masks, eth4_masks)
     call NullEvol_remask (i, 4+boundary_maskn, eth4_maskn)

     if (DEBUG_skip_J_update.eq.0) then

        call NullEvol_j (i, jcn, dxjcn, nucn, ckcn, bcn, cbcn, qcn, ucn, wcn,&
             jcn_p, nucn_p, ckcn_p, bcn_p, cbcn_p, qcn_p, ucn_p, wcn_p,&
             x_wt(:,:,1), j_wt(:,:,1), beta_wt(:,:,1), q_wt(:,:,1),&
             u_wt(:,:,1), w_wt(:,:,1), x_wt_p(:,:,1), j_wt_p(:,:,1),&
             beta_wt_p(:,:,1), q_wt_p(:,:,1), u_wt_p(:,:,1), w_wt_p(:,:,1),&
             evolution_maskn, auxiliary_maskn, eth4_maskn, dissip_mask,&
             null_dissip, cctk_delta_time, dissip_fudge, dissip_fudge_maxx,&
             dissip_eps_J, dissip_eps_Jx, first_order_scheme)

        call NullEvol_j (i, jcs, dxjcs, nucs, ckcs, bcs, cbcs, qcs, ucs, wcs,&
             jcs_p, nucs_p, ckcs_p, bcs_p, cbcs_p, qcs_p, ucs_p, wcs_p, &
             x_wt(:,:,2), j_wt(:,:,2), beta_wt(:,:,2), q_wt(:,:,2),&
             u_wt(:,:,2), w_wt(:,:,2), x_wt_p(:,:,2), j_wt_p(:,:,2),&
             beta_wt_p(:,:,2), q_wt_p(:,:,2), u_wt_p(:,:,2), w_wt_p(:,:,2),&
             evolution_masks, auxiliary_masks, eth4_masks, dissip_mask,&
             null_dissip, cctk_delta_time,  dissip_fudge, dissip_fudge_maxx,&
             dissip_eps_J, dissip_eps_Jx, first_order_scheme)

        call NullInterp_cinterp(cctkGH, tmp_cgfn, tmp_cgfs, jcn(:,:,i),jcs(:,:,i), two)

     end if

     if (DEBUG_skip_B_update.eq.0) then

        call NullEvol_beta (i, jcn, bcn, j_wt(:,:,1), &
                            beta_wt(:,:,1), x_wt(:,:,1), evolution_maskn)
        call NullEvol_beta (i, jcs, bcs, j_wt(:,:,2), &
                            beta_wt(:,:,2), x_wt(:,:,2), evolution_masks)

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
                         evolution_maskn, eth4_maskn, dissip_mask, dissip_eps_Q, first_order_scheme)

        call NullEvol_q (i, jcs, nucs, ckcs, bcs, cbcs, qcs, &
                         q_wt(:,:,2), j_wt(:,:,2), beta_wt(:,:,2), x_wt(:,:,2), &
                         evolution_masks, eth4_masks, dissip_mask, dissip_eps_Q, first_order_scheme)

        call NullInterp_cinterp(cctkGH, tmp_cgfn, tmp_cgfs, qcn(:,:,i),qcs(:,:,i), one)

     end if
     if (DEBUG_skip_U_update.eq.0) then
        ! note that there are no angular derivatives in the update of U
        ! NullEvol_u (i, jns, bns, qns, uns, u_x_wt, u_wt, x_wt, mask)
        call NullEvol_u (i, jcn, bcn, qcn, ucn, u_wt(:,:,1), u_x_wt(:,:,1),&
                         x_wt(:,:,1), evolution_maskn)
        call NullEvol_u (i, jcs, bcs, qcs, ucs, u_wt(:,:,2), u_x_wt(:,:,2),&
                         x_wt(:,:,2), evolution_masks)
     end if

     if (DEBUG_skip_W_update.eq.0) then

        call NullEvol_w (i, jcn, nucn, ckcn, bcn, cbcn, ucn, qcn, wcn,&
             w_wt(:,:,1), u_wt(:,:,1), x_wt(:,:,1), &
             evolution_maskn, eth4_maskn, dissip_mask, dissip_eps_W,first_order_scheme)
        call NullEvol_w (i, jcs, nucs, ckcs, bcs, cbcs, ucs, qcs, wcs,&
             w_wt(:,:,2), u_wt(:,:,2), x_wt(:,:,2), &
             evolution_masks, eth4_masks, dissip_mask, dissip_eps_W,first_order_scheme)

        call NullInterp_rinterp(cctkGH, tmp_rgfn, tmp_rgfs, wcn(:,:,i), wcs(:,:,i))

     end if

!Compute eth2_j
   call NullInterp_d2 (eth2jcn(:,:,i), jcn(:,:,i), two, one, one)
   call NullInterp_d2 (eth2jcs(:,:,i), jcs(:,:,i), two, one, one)

!Mask J and eth2_J
      jcn(:,:,i) = EG_mask * jcn(:,:,i)
      jcs(:,:,i) = EG_mask * jcs(:,:,i)

      eth2jcn(:,:,i) = EQ_mask * eth2jcn(:,:,i)
      eth2jcs(:,:,i) = EQ_mask * eth2jcs(:,:,i)

  end do

 if (old_J_xderiv.eq.0) then
!compute dx_J

   dxjcn(:,:,1) = -0.5 * (3.*jcn(:,:,1) - 4.*jcn(:,:,2) + jcn(:,:,3)) / null_dx
   dxjcs(:,:,1) = -0.5 * (3.*jcs(:,:,1) - 4.*jcs(:,:,2) + jcs(:,:,3)) / null_dx

   do i = 2, N_radial_pts-1

      dxjcn(:,:,i) = 0.5 * (jcn(:,:,i+1) - jcn(:,:,i-1)) / null_dx
      dxjcs(:,:,i) = 0.5 * (jcs(:,:,i+1) - jcs(:,:,i-1)) / null_dx

   end do

   dxjcn(:,:,N_radial_pts) = 0.5 * (3.*jcn(:,:,N_radial_pts) &
              - 4.*jcn(:,:,N_radial_pts-1) + jcn(:,:,N_radial_pts-2)) / null_dx
   dxjcs(:,:,N_radial_pts) = 0.5 * (3.*jcs(:,:,N_radial_pts) &
              - 4.*jcs(:,:,N_radial_pts-1) + jcs(:,:,N_radial_pts-2)) / null_dx

 end if

!compute array for the radial profile
  qp = int((null_lsh(1)-1)/2) + 1
  pq = int((null_lsh(2)-1)/2) + 1

  jcn_rad(:) = jcn(qp,pq,:)
  jcs_rad(:) = jcs(qp,pq,:)

  dxjcn_rad(:) = dxjcn(qp,pq,:)
  dxjcs_rad(:) = dxjcs(qp,pq,:)

  call NullEvol_ResetInactive(CCTK_PASS_FTOF,1_ik,N_radial_pts)

!  write (*,*) 'NULL_EVOLSTEP: J_scri N are: ', maxval(abs(jcn(:,:,N_radial_pts)))
!  write (*,*) 'NULL_EVOLSTEP: J_scri S are: ', maxval(abs(jcs(:,:,N_radial_pts)))

end subroutine NullEvol_Step

subroutine NullEvol_CopyBdry(CCTK_ARGUMENTS)
  use NullEvol_Mask
  use NullGrid_Vars 
  implicit none

  CCTK_INT :: i
  CCTK_INT :: mn, global_mn
  CCTK_INT :: reval
  CCTK_INT, save :: reduce_handle

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  ! Note: the extraction WT is at rb(char_nx0+1)
  ! start marching at the 2nd point

  call CCTK_ReductionArrayHandle(reduce_handle, "minimum");
  if (reduce_handle .lt. 0 ) then
     call CCTK_WARN(0,"Could not get reduction handle")
  endif

  call CCTK_INFO("Null Evolution");
  mn = min(minval(boundary_masks(1:null_lsh(1),1:null_lsh(2))),&
       minval(boundary_maskn(1:null_lsh(1),1:null_lsh(2))))
  call CCTK_ReduceLocScalar(reval, cctkGH, -1, reduce_handle,&
       mn, global_mn, CCTK_VARIABLE_INT)

  if (reval .ne. 0 ) then
     call CCTK_WARN(0,"Error in obtaining Minimum of mask")
  endif

  !write(message,*) "Minval of mask is ", global_mn
  !call CCTK_INFO(trim(message))

  do i = global_mn+1, N_radial_pts  

     call NullEvol_remask (i, boundary_masks, evolution_masks)
     call NullEvol_remask (i, boundary_maskn, evolution_maskn)

     jcn(:,:,i) = jcn(:,:,i) * (1 - evolution_maskn) + evolution_maskn * jcn(:,:,i-1) 
     jcs(:,:,i) = jcs(:,:,i) * (1 - evolution_masks) + evolution_masks * jcs(:,:,i-1) 

     bcn(:,:,i) = bcn(:,:,i) * (1 - evolution_maskn) + evolution_maskn * bcn(:,:,i-1) 
     bcs(:,:,i) = bcs(:,:,i) * (1 - evolution_masks) + evolution_masks * bcs(:,:,i-1) 

     qcn(:,:,i) = qcn(:,:,i) * (1 - evolution_maskn) + evolution_maskn * qcn(:,:,i-1) 
     qcs(:,:,i) = qcs(:,:,i) * (1 - evolution_masks) + evolution_masks * qcs(:,:,i-1) 

     ucn(:,:,i) = ucn(:,:,i) * (1 - evolution_maskn) + evolution_maskn * ucn(:,:,i-1) 
     ucs(:,:,i) = ucs(:,:,i) * (1 - evolution_masks) + evolution_masks * ucs(:,:,i-1) 

     wcn(:,:,i) = wcn(:,:,i) * (1 - evolution_maskn) + evolution_maskn * wcn(:,:,i-1) 
     wcs(:,:,i) = wcs(:,:,i) * (1 - evolution_masks) + evolution_masks * wcs(:,:,i-1) 

     if (first_order_scheme.ne.0) then

        nucn(:,:,i) = nucn(:,:,i) * (1 - evolution_maskn) + evolution_maskn * nucn(:,:,i-1) 
        nucs(:,:,i) = nucs(:,:,i) * (1 - evolution_masks) + evolution_masks * nucs(:,:,i-1) 

        ckcn(:,:,i) = ckcn(:,:,i) * (1 - evolution_maskn) + evolution_maskn * ckcn(:,:,i-1) 
        ckcs(:,:,i) = ckcs(:,:,i) * (1 - evolution_masks) + evolution_masks * ckcs(:,:,i-1) 

        cbcn(:,:,i) = cbcn(:,:,i) * (1 - evolution_maskn) + evolution_maskn * cbcn(:,:,i-1) 
        cbcs(:,:,i) = cbcs(:,:,i) * (1 - evolution_masks) + evolution_masks * cbcs(:,:,i-1) 

     end if

  end do


end subroutine NullEvol_CopyBdry

subroutine NullEvol_ResetTop(CCTK_ARGUMENTS)

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  jcn = 1.e+20
  jcs = 1.e+20

  bcn = 1.e+20
  bcs = 1.e+20

  qcn = 1.e+20
  qcs = 1.e+20

  ucn = 1.e+20
  ucs = 1.e+20

  wcn = 1.e+20
  wcs = 1.e+20

  if (first_order_scheme.ne.0) then
     nucn = 1.e+20
     nucs = 1.e+20

     ckcn = 1.e+20
     ckcs = 1.e+20

     cbcn = 1.e+20
     cbcs = 1.e+20

  end if

end subroutine NullEvol_ResetTop



subroutine NullEvol_ResetInactive(CCTK_ARGUMENTS, imin, imax)

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT i, imin, imax

#define reset_inactive(arr) arr = arr * EG_mask

  do i = max(1,imin), min(N_radial_pts, imax)

     if (first_order_scheme.ne.0) then

        reset_inactive(nucn(:,:,i))
        reset_inactive(nucs(:,:,i))

        reset_inactive(ckcn(:,:,i))
        reset_inactive(ckcs(:,:,i))

        reset_inactive(cbcn(:,:,i))
        reset_inactive(cbcs(:,:,i))

     end if

     reset_inactive(jcn(:,:,i))
     reset_inactive(jcs(:,:,i))

     reset_inactive(bcn(:,:,i))
     reset_inactive(bcs(:,:,i))

     reset_inactive(ucn(:,:,i))
     reset_inactive(ucs(:,:,i))

     reset_inactive(wcn(:,:,i))
     reset_inactive(wcs(:,:,i))

  end do

#undef reset_inactive

end subroutine NullEvol_ResetInactive


