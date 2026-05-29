! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#define cr Extract_R

subroutine NullSHRE_ExtractWT (CCTK_ARGUMENTS)

  !---------------------------------------------------------------------------
  ! extraction routine - produces boundary data for the null code
  ! this version takes care of the finite differencing of the time
  ! derivatives within the wt modules.
  !       b_box is the bounding box of the 3+1 spatial grid,
  !       g_    is the metric at level n,
  !       g_nm1 is the metric at level  n-1 etc.
  !---------------------------------------------------------------------------

  use cctk
  use NullSHRE_modVars
  use NullSHRE_modAnaCoord
  use NullSHRE_mod4Metric
  use NullSHRE_modInverse
  use NullSHRE_modNullDir
  use NullSHRE_modEta
  use NullSHRE_modRadius
  use NullSHRE_modBoundary
  use NullInterp
  use NullGrid_Vars

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  integer :: retval
  CCTK_INT :: i_max, gi_max, i_min, gi_min, i
  CCTK_REAL :: rl_min, elld_min
  integer, save :: reduce_handle    
  character(len=500) :: message

#define dt cctk_delta_time
#define dd null_delta(1)

  !---------------------------------------------------------------------------
  ! world-tube, and analytic part of the coordinate transformation
  !---------------------------------------------------------------------------

  if (abs(null_delta(1)-null_delta(2)).gt.1.e-12) then
     call CCTK_WARN(0, "cannot work with gridsteps that are different in the q,p directions")
  end if

  call wt_sigma (lsh(1), lsh(2), qs, ps, pp, ip_n, cr, sigma_n)
  call wt_sigma (lsh(1), lsh(2), qs, ps, pp, ip_s, cr, sigma_s)

  call wt_dsigma (lsh(1), lsh(2), qs, ps, pp, ip_n, cr, dsigma_n)
  call wt_dsigma (lsh(1), lsh(2), qs, ps, pp, ip_s, cr, dsigma_s)

  !---------------------------------------------------------------------------
  ! 4-metric and its derivatives
  !---------------------------------------------------------------------------

  call wt_g(alpha_n, beta_n, g_n)
  call wt_g(alpha_s, beta_s, g_s)

  call wt_dg (alpha_n, beta_n, dalpha_n, dbeta_n, g_n, dg_n)
  call wt_dg (alpha_s, beta_s, dalpha_s, dbeta_s, g_s, dg_s)

  !---------------------------------------------------------------------------
  ! upstairs 3-metric
  !---------------------------------------------------------------------------
  call wt_gup3 (g_n, detg_n, gup_n)
  call wt_gup3 (g_s, detg_s, gup_s)

  !---------------------------------------------------------------------------
  ! normals, null direction, geodesic expansion
  !---------------------------------------------------------------------------

  call wt_na (alpha_n, beta_n, na_n)
  call wt_na (alpha_s, beta_s, na_s)

  call wt_sa (gup_n, sigma_n, sigma2_n, sa_n)
  call wt_sa (gup_s, sigma_s, sigma2_s, sa_s)

 elld_min = elld_min_coef * cctk_delta_time**elld_min_pow

  call wt_ell (g_n, alpha_n, beta_n, sa_n, na_n, elld_n, ell_n, halt_on_negative_elld, elld_min)
  call wt_ell (g_s, alpha_s, beta_s, sa_s, na_s, elld_s, ell_s, halt_on_negative_elld, elld_min)

  call wt_j0inv (lsh(1), lsh(2), qs, ps, pp, ip_n, cr, j0inv_n)
  call wt_j0inv (lsh(1), lsh(2), qs, ps, pp, ip_s, cr, j0inv_s)

  call wt_j0 (lsh(1), lsh(2), qs, ps, pp, ip_n, cr, ell_n, j0_n)
  call wt_j0 (lsh(1), lsh(2), qs, ps, pp, ip_s, cr, ell_s, j0_s)

  call wt_dj0 (lsh(1), lsh(2), qs, ps, pp, ip_n, cr, dj0_n)
  call wt_dj0 (lsh(1), lsh(2), qs, ps, pp, ip_s, cr, dj0_s)

  !---------------------------------------------------------------------------
  ! derivatives of the normals
  !---------------------------------------------------------------------------

  call wt_dna (alpha_n, beta_n, dalpha_n, dbeta_n, dna_n)
  call wt_dna (alpha_s, beta_s, dalpha_s, dbeta_s, dna_s)

  call wt_dsa (sigma2_n, dsigma_n, sa_n, gup_n, dg_n, dsa_n)
  call wt_dsa (sigma2_s, dsigma_s, sa_s, gup_s, dg_s, dsa_s)

  call wt_gup4 (alpha_n, na_n, gup_n)
  call wt_gup4 (alpha_s, na_s, gup_s)

  !-------------------------------------------------------------------------
  ! derivatives with respect to the affine parameter lambda
  !-------------------------------------------------------------------------

  call wt_g1 (dg_n, ell_n, j0inv_n, g1_n)
  call wt_g1 (dg_s, ell_s, j0inv_s, g1_s)

  call wt_j1 (g_n, gup_n, dg_n, dalpha_n, beta_n, dbeta_n, na_n,&
       sa_n, ell_n, dsa_n, dna_n, j0inv_n, elld_n, delld_n, j1_n)
  call wt_j1 (g_s, gup_s, dg_s, dalpha_s, beta_s, dbeta_s, na_s,&
       sa_s, ell_s, dsa_s, dna_s, j0inv_s, elld_s, delld_s, j1_s)

  !---------------------------------------------------------------------------
  ! first coordinate transformation
  !---------------------------------------------------------------------------

  call wt_eta0 (g_n, j0_n, eta0_n)
  call wt_eta0 (g_s, j0_s, eta0_s)

  call wt_eta1 (g_n, g1_n, j0_n, j1_n, eta1_n)
  call wt_eta1 (g_s, g1_s, j0_s, j1_s, eta1_s)

  call wt_etaup0 (eta0_n, temp_n, etaup0_n) !detg is rewritten
  call wt_etaup0 (eta0_s, temp_s, etaup0_s) !detg is rewritten

  call wt_etaup1 (eta1_n, etaup0_n, etaup1_n)
  call wt_etaup1 (eta1_s, etaup0_s, etaup1_s)

  !---------------------------------------------------------------------------
  ! second coordinate transformation
  !---------------------------------------------------------------------------

  call wt_r (lsh(1), lsh(2), pp, eta0_n, r0_n)
  call wt_r (lsh(1), lsh(2), pp, eta0_s, r0_s)

!  write (*,*) 'minimum value of radial coordinate is ', min(minval(r0_n%d),minval(r0_s%d))
!  write (*,*) 'maximum value of radial coordinate is ', max(maxval(r0_n%d),maxval(r0_s%d))

  ! r_{,\lambda} from the determinant condition

  rl_min = rl_min_coef * cctk_delta_time**rl_min_pow


  call wt_rl (etaup0_n, eta1_n, r0_n, temp_n, dr0_n, halt_on_negative_rl, rl_min)
  call wt_rl (etaup0_s, eta1_s, r0_s, temp_s, dr0_s, halt_on_negative_rl, rl_min)

  !---------------------------------------------------------------------------
  ! r_{,A}
  !---------------------------------------------------------------------------

  call wt_ra_long (lsh(1), lsh(2), qs, ps, pp, g_n, dg_n,&
       j0_n, dj0_n, etaup0_n, r0_n, temp_n, deta0_n, dr0_n)
  call wt_ra_long (lsh(1), lsh(2), qs, ps, pp, g_s, dg_s,&
       j0_s, dj0_s, etaup0_s, r0_s, temp_s, deta0_s, dr0_s)

  call wt_ru (dg_n, j0_n, etaup0_n, r0_n, temp_n, dr0_n)
  call wt_ru (dg_s, j0_s, etaup0_s, r0_s, temp_s, dr0_s)

  !---------------------------------------------------------------------------
  ! r_lu, r_lA by finite differences
  !---------------------------------------------------------------------------

  call wt_rlu (dt, dr0_n, rl_old_n, dr1_n)
  call wt_rlu (dt, dr0_s, rl_old_s, dr1_s)

  call NullInterp_wt_da_ethr(cctkGH, tmp_cgfn, tmp_cgfs,&
       dr0_n(1)%d, dr0_s(1)%d, dr1_n(2)%d, dr1_s(2)%d, dr1_n(3)%d, dr1_s(3)%d)

  !----------------------------------------------------------------
  ! compactified x on the world tube
  !----------------------------------------------------------------


  call wt_null_boundary_mask (lsh(1), lsh(2), dx, rwt, null_xin, r0_n, x_wt_n,&
                              boundary_maskn, MinimumDistanceTo2Bp1)
  call wt_null_boundary_mask (lsh(1), lsh(2), dx, rwt, null_xin, r0_s, x_wt_s,&
                              boundary_masks, MinimumDistanceTo2Bp1)

  if (min(minval(x_wt_n%d),minval(x_wt_s%d)).lt.null_xb(4)) then
     call CCTK_WARN(0,"Innermost Bondi grid point must be 4 points inside the extraction radius")
  endif

  call CCTK_ReductionArrayHandle(reduce_handle, "minimum")
  if (reduce_handle .lt. 0 ) then
     call CCTK_WARN(0,"Could not get reduction handle")
  endif

  i_min =min(minval(boundary_maskn),minval(boundary_masks))

  call CCTK_ReduceLocScalar(retval, cctkGH, -1, reduce_handle,&
       i_min, gi_min, CCTK_VARIABLE_INT)

  !---------------------------------------------------------------------------
  ! world-tube data for the null code
  !---------------------------------------------------------------------------

  call wt_jb (lsh(1), lsh(2), pp, eta0_n, eta1_n, r0_n, dr0_n, qa,&
              j_wt_n, j_l_n)
  call wt_jb (lsh(1), lsh(2), pp, eta0_s, eta1_s, r0_s, dr0_s, qa,&
              j_wt_s, j_l_s)

  call wt_bb (r0_n, dr0_n, j_wt_n, j_l_n, beta_wt_n, beta_l_n)
  call wt_bb (r0_s, dr0_s, j_wt_s, j_l_s, beta_wt_s, beta_l_s)

  call wt_rll (beta_l_n, dr0_n, dr1_n)
  call wt_rll (beta_l_s, dr0_s, dr1_s)

  call wt_ub (qs, ps, lsh(1), lsh(2), etaup0_n, etaup1_n, dr0_n, dr1_n,&
              beta_l_n, qa, u_wt_n, u_l_n)
  call wt_ub (qs, ps, lsh(1), lsh(2), etaup0_s, etaup1_s, dr0_s, dr1_s,&
              beta_l_s, qa, u_wt_s, u_l_s)

  call wt_qb(rwt, r0_n, j_wt_n, u_l_n, q_wt_n)
  call wt_qb(rwt, r0_s, j_wt_s, u_l_s, q_wt_s)

  call wt_wb (etaup0_n, etaup1_n, r0_n, dr0_n, dr1_n, beta_l_n, w_wt_n, w_l_n)
  call wt_wb (etaup0_s, etaup1_s, r0_s, dr0_s, dr1_s, beta_l_s, w_wt_s, w_l_s)

  ! apply boundary data to characteristic 3D arrays, fill only first 4 points 

  call CCTK_ReductionArrayHandle(reduce_handle, "maximum")
  if (reduce_handle .lt. 0 ) then
     call CCTK_WARN(0,"Could not get reduction handle")
  endif

  i_max = max(maxval(boundary_maskn), maxval(boundary_masks))
  call CCTK_ReduceLocScalar(retval, cctkGH, -1, reduce_handle,&
       i_max, gi_max, CCTK_VARIABLE_INT)

  if (retval .ne. 0 ) then
     call CCTK_WARN(0,"Error in obtaining the points filled by extraction outside the boundary mask")
  endif


  write(message,*) "CCE world-tube index range is", gi_min, gi_max, null_rb(gi_min), null_rb(gi_max)
  call CCTK_INFO(trim(message))
    

  gi_max = gi_max + 4 ! fill up an extra 4 points outside the world-tube


  ! J

  ! debugging
  gi_max = N_radial_pts

  if (cctk_iteration.le.3) then 

     ! initially give J on all radial points

     call apply_cbdata (N_radial_pts, N_radial_pts, lsh(1), lsh(2), jcn,&
          j_wt_n, j_l_n, dr0_n, rwt, x_wt_n, xb)
     call apply_cbdata (N_radial_pts, N_radial_pts, lsh(1), lsh(2), jcs,&
          j_wt_s, j_l_s, dr0_s, rwt, x_wt_s, xb)

  else

     call apply_cbdata (gi_max, N_radial_pts, lsh(1), lsh(2), jcn,&
          j_wt_n, j_l_n, dr0_n, rwt, x_wt_n, xb)
     call apply_cbdata (gi_max, N_radial_pts, lsh(1), lsh(2), jcs,&
          j_wt_s, j_l_s, dr0_s, rwt, x_wt_s, xb)

  end if

  ! beta

  call apply_rbdata (gi_max, N_radial_pts, lsh(1), lsh(2), bcn,&
       beta_wt_n, beta_l_n, dr0_n, rwt, x_wt_n, xb)
  call apply_rbdata (gi_max, N_radial_pts, lsh(1), lsh(2), bcs,&
       beta_wt_s, beta_l_s, dr0_s, rwt, x_wt_s, xb)

  ! U
  call apply_cbdata (gi_max, N_radial_pts, lsh(1), lsh(2), ucn,&
       u_wt_n, u_l_n, dr0_n, rwt, x_wt_n, xbh)
  call apply_cbdata (gi_max, N_radial_pts, lsh(1), lsh(2), ucs,&
       u_wt_s, u_l_s, dr0_s, rwt, x_wt_s, xbh)

  ! Q
  do i = 1, gi_max 
     qcn(:,:,i) = q_wt_n%d * ((1-null_xb(i))*r0_n%d/rwt/null_xb(i))**2
     qcs(:,:,i) = q_wt_s%d * ((1-null_xb(i))*r0_s%d/rwt/null_xb(i))**2
  end do

  ! U_{,x}

  U_x_wt(:,:,1) = u_l_n%d / dr0_n(1)%d * rwt / (1.d0 - x_wt_n%d)** 2
  U_x_wt(:,:,2) = u_l_s%d / dr0_s(1)%d * rwt / (1.d0 - x_wt_s%d)** 2

  if (first_order_scheme.ne.0) then

     ! compute the following
     ! \nu = \bar \eth J
     ! ck = \eth K
     ! cb = \eth \beta

     ! Take \eth on the Bondi spheres:

     call wt_eth_expand(cctkGH, gi_max, N_radial_pts, lsh(1), lsh(2),&
          tmp_cgfn, tmp_cgfs, bcn, bcs, jcn, jcs, cbcn, cbcs, nucn, nucs, ckcn, ckcs)

  end if ! first_order_scheme

  ! W
  call apply_rbdata (gi_max, N_radial_pts, lsh(1), lsh(2), wcn,&
       w_wt_n, w_l_n, dr0_n, rwt, x_wt_n, xb)
  call apply_rbdata (gi_max, N_radial_pts, lsh(1), lsh(2), wcs,&
       w_wt_s, w_l_s, dr0_s, rwt, x_wt_s, xb)

!  write (*,*) 'minimum detg in ExtractWT', min(minval(detg_n%d),minval(detg_s%d))
!  write (*,*) 'maximum detg in ExtractWT', max(maxval(detg_n%d),maxval(detg_s%d))

end subroutine NullSHRE_ExtractWT
