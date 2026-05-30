! vim: syntax=fortran
#include "cctk.h"


module NullSHRE_modVars

   use cctk
   use NullSHRE_modGFdef

   implicit none

   !----------------------------------------------------------------------
   ! radius of the worldtube
   ! values at the level "n"
   !----------------------------------------------------------------------

   type (gf2d), save, dimension (4,4)       :: g_n, g_s
   type (gf2d), save                        :: alpha_n, alpha_s
   type (gf2d), save, dimension (3)         :: beta_n, beta_s

   type (gf2d), save, dimension (4,4,4)     :: dg_n, dg_s
   type (gf2d), save, dimension (4)         :: dalpha_n, dalpha_s
   type (gf2d), save, dimension (3,4)       :: dbeta_n, dbeta_s

   type (gf2d), save, dimension (3)         :: sigma_n, sigma_s
   type (gf2d), save, dimension (4,4)       :: g1_n, g1_s
   type (gf2d), save, dimension (4,4)       :: gup_n, gup_s
   type (gf2d), save, dimension (4,4)       :: j0_n, j0_s
   type (gf2d), save, dimension (4,4)       :: j0inv_n, j0inv_s
   type (gf2d), save, dimension (4,4)       :: j1_n, j1_s
   type (gf2d), save                        :: detg_n, detg_s
   type (gf2d), save                        :: temp_n, temp_s
   type (gf2d), save                        :: sigma2_n, sigma2_s
   type (gf2d), save                        :: elld_n, elld_s
   type (gf2d), save                        :: r0_n, r0_s
   type (gf2d), save, dimension (3,2:3)     :: dsigma_n, dsigma_s
   type (gf2d), save, dimension (4)         :: sa_n, sa_s
   type (gf2d), save, dimension (4)         :: na_n, na_s
   type (gf2d), save, dimension (4)         :: ell_n, ell_s
   type (gf2d), save, dimension (4)         :: delld_n, delld_s
   type (gf2d), save, dimension (4)         :: dr0_n, dr0_s
   type (gf2d), save, dimension (2)         :: rl_old_n, rl_old_s
   type (gf2d), save, dimension (4)         :: dr1_n, dr1_s
   type (gf2d), save, dimension (4,2:4)     :: dsa_n, dsa_s
   type (gf2d), save, dimension (4,2:4)     :: dna_n, dna_s
   type (gf2d), save, dimension (3,2:3,2:3) :: dj0_n, dj0_s
   type (gf2d), save, dimension (4,4)       :: eta0_n, eta0_s
   type (gf2d), save, dimension (2:3,2:3,2:3) :: deta0_n, deta0_s
   type (gf2d), save, dimension (4,4)       :: eta1_n, eta1_s
   type (gf2d), save, dimension (4,4)       :: etaup0_n, etaup0_s
   type (gf2d), save, dimension (4,4)       :: etaup1_n, etaup1_s

   type (gf2d), save, dimension (4,2:3)     :: bracket_n, bracket_s
 
   type (gf2dc), save, dimension (2:3) :: qa

   type (gf2d), save  :: x_wt_n, x_wt_s

   type (gf2dc), save :: j_wt_n, j_wt_s
   type (gf2dc), save :: j_l_n, j_l_s

   type (gf2dc), save :: nu_wt_n, nu_wt_s
   type (gf2dc), save :: nu_x_n, nu_x_s

   type (gf2dc), save :: ck_wt_n, ck_wt_s
   type (gf2dc), save :: ck_x_n, ck_x_s

   type (gf2d), save  :: beta_wt_n, beta_wt_s
   type (gf2d), save  :: beta_l_n, beta_l_s

   type (gf2dc), save :: cb_wt_n, cb_wt_s
   type (gf2dc), save :: cb_x_n, cb_x_s

   type (gf2dc), save :: u_wt_n, u_wt_s
   type (gf2dc), save :: u_l_n, u_l_s
   type (gf2dc), save :: q_wt_n, q_wt_s

   type (gf2d), save  :: w_wt_n, w_wt_s
   type (gf2d), save  :: w_l_n, w_l_s

end module NullSHRE_modVars


