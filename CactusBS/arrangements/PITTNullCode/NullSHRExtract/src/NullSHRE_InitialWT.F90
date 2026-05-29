! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"


subroutine NullSHRE_InitialWT (CCTK_ARGUMENTS)

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
  CCTK_REAL :: rl_min, elld_min

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

  call wt_etaup0 (eta0_n, temp_n, etaup0_n)
  call wt_etaup0 (eta0_s, temp_s, etaup0_s)

  call wt_etaup1 (eta1_n, etaup0_n, etaup1_n)
  call wt_etaup1 (eta1_s, etaup0_s, etaup1_s)

  !---------------------------------------------------------------------------
  ! second coordinate transformation
  !---------------------------------------------------------------------------

  call wt_r (lsh(1), lsh(2), pp, eta0_n, r0_n)
  call wt_r (lsh(1), lsh(2), pp, eta0_s, r0_s)

  ! r_{,\lambda} from the determinant condition

  rl_min = rl_min_coef * cctk_delta_time**rl_min_pow

  call wt_rl (etaup0_n, eta1_n, r0_n, temp_n, dr0_n, halt_on_negative_rl, rl_min)
  call wt_rl (etaup0_s, eta1_s, r0_s, temp_s, dr0_s, halt_on_negative_rl, rl_min)

  !---------------------------------------------------------------------------
  ! compute J for the sake of having ready on the past time-level
  !---------------------------------------------------------------------------

  call wt_jb_simple (lsh(1), lsh(2), pp, eta0_n, r0_n, qa, j_wt_n)
  call wt_jb_simple (lsh(1), lsh(2), pp, eta0_s, r0_s, qa, j_wt_s)

end subroutine NullSHRE_InitialWT
