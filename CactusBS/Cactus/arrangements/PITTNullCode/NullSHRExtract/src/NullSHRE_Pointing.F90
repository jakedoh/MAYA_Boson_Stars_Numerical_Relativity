! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
 
subroutine NullSHRE_Pointing (CCTK_ARGUMENTS)

  use cctk
  use NullSHRE_modVars
  implicit none

  TARGET :: SHRE_alpha
  TARGET :: SHRE_dralpha, SHRE_dqalpha, SHRE_dpalpha, SHRE_dtalpha 
  TARGET :: WT_detg, WT_temp, WT_sigma2, WT_elld, WT_r0, WT_r1 
  TARGET :: WT_r1_p, WT_r1_p_p 
  TARGET :: qa_2, qa_3, x_wt, j_wt, j_l
  TARGET :: beta_wt, beta_l, u_wt, u_l, q_wt, w_wt, w_l


  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
 
  CCTK_INT, parameter :: izero = 0
  integer, parameter :: ik = kind(izero)

  logical, save :: first_call = .true.
  if (first_call) call nullify_all

  alpha_n%d => SHRE_alpha(:,:,1)
  alpha_s%d => SHRE_alpha(:,:,2)
  dalpha_n(1)%d => SHRE_dralpha(:,:,1)
  dalpha_s(1)%d => SHRE_dralpha(:,:,2)
  dalpha_n(2)%d => SHRE_dqalpha(:,:,1)
  dalpha_s(2)%d => SHRE_dqalpha(:,:,2)
  dalpha_n(3)%d => SHRE_dpalpha(:,:,1)
  dalpha_s(3)%d => SHRE_dpalpha(:,:,2)
  dalpha_n(4)%d => SHRE_dtalpha(:,:,1)
  dalpha_s(4)%d => SHRE_dtalpha(:,:,2)

  call point_v1(null_lsh, 6_ik, 1_ik, 3_ik, beta_n, beta_s, SHRE_beta)
  call point_v1(null_lsh, 6_ik, 1_ik, 3_ik, dbeta_n(:,1), dbeta_s(:,1), SHRE_drbeta)
  call point_v1(null_lsh, 6_ik, 1_ik, 3_ik, dbeta_n(:,2), dbeta_s(:,2), SHRE_dqbeta)
  call point_v1(null_lsh, 6_ik, 1_ik, 3_ik, dbeta_n(:,3), dbeta_s(:,3), SHRE_dpbeta)
  call point_v1(null_lsh, 6_ik, 1_ik, 3_ik, dbeta_n(:,4), dbeta_s(:,4), SHRE_dtbeta)
  
  call point_v2sym(null_lsh, 12_ik, 1_ik, 3_ik,&
                   g_n(1:3_ik,1:3), g_s(1:3_ik,1:3), SHRE_gij)
  call point_v2sym(null_lsh, 12_ik, 1_ik, 3_ik,&
                   dg_n(1:3_ik,1:3_ik,1), dg_s(1:3_ik,1:3_ik,1), SHRE_drgij)
  call point_v2sym(null_lsh, 12_ik, 1_ik, 3_ik,&
                   dg_n(1:3_ik,1:3_ik,2), dg_s(1:3_ik,1:3_ik,2), SHRE_dqgij)
  call point_v2sym(null_lsh, 12_ik, 1_ik, 3_ik,&
                   dg_n(1:3_ik,1:3_ik,3), dg_s(1:3_ik,1:3_ik,3), SHRE_dpgij)
  call point_v2sym(null_lsh, 12_ik, 1_ik, 3_ik,&
                   dg_n(1:3_ik,1:3_ik,4), dg_s(1:3_ik,1:3_ik,4), SHRE_dtgij)

  call point_v1sym(null_lsh, 8_ik, 1_ik, 4_ik, g_n(:,:), g_s(:,:), SHRE_git)
  call point_v1sym(null_lsh, 8_ik, 1_ik, 4_ik, dg_n(:,:,1), dg_s(:,:,1), SHRE_drgit)
  call point_v1sym(null_lsh, 8_ik, 1_ik, 4_ik, dg_n(:,:,2), dg_s(:,:,2), SHRE_dqgit)
  call point_v1sym(null_lsh, 8_ik, 1_ik, 4_ik, dg_n(:,:,3), dg_s(:,:,3), SHRE_dpgit)
  call point_v1sym(null_lsh, 8_ik, 1_ik, 4_ik, dg_n(:,:,4), dg_s(:,:,4), SHRE_dtgit)

  detg_n%d => WT_detg(:,:,1)
  detg_s%d => WT_detg(:,:,2)

  temp_n%d => WT_temp(:,:,1)
  temp_s%d => WT_temp(:,:,2)

  sigma2_n%d => WT_sigma2(:,:,1)
  sigma2_s%d => WT_sigma2(:,:,2)

  elld_n%d =>WT_elld(:,:,1)
  elld_s%d =>WT_elld(:,:,2)

  r0_n%d =>WT_r0(:,:,1)
  r0_s%d =>WT_r0(:,:,2)

  dr0_n(1)%d =>WT_r1(:,:,1)
  dr0_s(1)%d =>WT_r1(:,:,2)
  rl_old_n(1)%d =>WT_r1_p(:,:,1)
  rl_old_s(1)%d =>WT_r1_p(:,:,2)
  rl_old_n(2)%d =>WT_r1_p_p(:,:,1)
  rl_old_s(2)%d =>WT_r1_p_p(:,:,2)

  call point_v2(null_lsh, 12_ik, 1_ik, 3_ik, 2_ik, 3_ik, dsigma_n, dsigma_s, WT_dsigma_vect)

  call point_v2sym(null_lsh, 20_ik, 1_ik, 4_ik, gup_n, gup_s, WT_gup_vect)
  call point_v2sym(null_lsh, 20_ik, 1_ik, 4_ik, g1_n, g1_s, WT_g1_vect)

  call point_v2(null_lsh, 24_ik, 1_ik, 4_ik, 2_ik, 4_ik, dsa_n, dsa_s, WT_DSnorm_vect)
  call point_v2(null_lsh, 24_ik, 1_ik, 4_ik, 2_ik, 4_ik, dna_n, dna_s, WT_DNnorm_vect)

  call point_v1(null_lsh, 6_ik, 2_ik, 4_ik, dr0_n(2:4), dr0_s(2:4), WT_dr0_vect)
  call point_v1(null_lsh, 6_ik, 1_ik, 3_ik, sigma_n, sigma_s, WT_sigma_vect)
  call point_v1(null_lsh, 8_ik, 1_ik, 4_ik, dr1_n, dr1_s, WT_dr1_vect)
  call point_v3_sym23(null_lsh, 18_ik, 1_ik, 3_ik, 2_ik, 3_ik, dj0_n, dj0_s, WT_dj0_vect)

  call point_v2(null_lsh, 32_ik, 1_ik, 4_ik, 1_ik, 4_ik, j0_n, j0_s, WT_j0_vect)
  call point_v2(null_lsh, 32_ik, 1_ik, 4_ik, 1_ik, 4_ik, j0inv_n, j0inv_s, WT_j0inv_vect)
  call point_v2(null_lsh, 32_ik, 1_ik, 4_ik, 1_ik, 4_ik, j1_n, j1_s, WT_j1_vect)

  call point_v1(null_lsh, 8_ik, 1_ik, 4_ik, sa_n, sa_s, WT_sa_vect)
  call point_v1(null_lsh, 8_ik, 1_ik, 4_ik, na_n, na_s, WT_na_vect)
  call point_v1(null_lsh, 8_ik, 1_ik, 4_ik, ell_n, ell_s, WT_ell_vect)
  call point_v1(null_lsh, 8_ik, 1_ik, 4_ik, delld_n, delld_s, WT_delld_vect)

  call point_v2sym(null_lsh, 20_ik, 1_ik, 4_ik, eta0_n, eta0_s, WT_eta0_vect)
  call point_v2sym(null_lsh, 20_ik, 1_ik, 4_ik, eta1_n, eta1_s, WT_eta1_vect)
  call point_v2sym(null_lsh, 20_ik, 1_ik, 4_ik, etaup0_n, etaup0_s, WT_etaup0_vect)
  call point_v2sym(null_lsh, 20_ik, 1_ik, 4_ik, etaup1_n, etaup1_s, WT_etaup1_vect)
  call point_v3_sym12(null_lsh, 12_ik, 2_ik, 3_ik, 2_ik, 3_ik, deta0_n, deta0_s, WT_deta0_vect)

  qa(2)%d => qa_2
  qa(3)%d => qa_3

  x_wt_n%d => x_wt(:,:,1)
  x_wt_s%d => x_wt(:,:,2)

  j_wt_n%d => j_wt(:,:,1)
  j_wt_s%d => j_wt(:,:,2)
  j_l_n%d => j_l(:,:,1)
  j_l_s%d => j_l(:,:,2)

  beta_wt_n%d =>beta_wt(:,:,1)
  beta_wt_s%d =>beta_wt(:,:,2)
  beta_l_n%d => beta_l(:,:,1)
  beta_l_s%d => beta_l(:,:,2)

  u_wt_n%d => u_wt(:,:,1)
  u_wt_s%d => u_wt(:,:,2)
  u_l_n%d => u_l(:,:,1)
  u_l_s%d => u_l(:,:,2)
  q_wt_n%d => q_wt(:,:,1)
  q_wt_s%d => q_wt(:,:,2)

  w_wt_n%d => w_wt(:,:,1)
  w_wt_s%d => w_wt(:,:,2)
  w_l_n%d => w_l(:,:,1)
  w_l_s%d => w_l(:,:,2)

  first_call = .false.

contains

#include "NullSHRE_PointingLib.h"
#include "NullSHRE_PointingNull.h"

end subroutine NullSHRE_Pointing
 

subroutine NullSHRE_Pointing_p (CCTK_ARGUMENTS)

  use cctk
  use NullSHRE_modVars
  implicit none

  CCTK_INT, parameter :: izero = 0
  integer, parameter :: ik = kind(izero)

  TARGET :: SHRE_alpha_p, WT_r1_p
  TARGET :: SHRE_dralpha_p, SHRE_dqalpha_p, SHRE_dpalpha_p, SHRE_dtalpha_p

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS


  interface
     subroutine NullSHRE_Pointing (CCTK_ARGUMENTS)
       implicit none
       TARGET :: SHRE_alpha
       TARGET :: SHRE_dralpha, SHRE_dqalpha, SHRE_dpalpha, SHRE_dtalpha 
       TARGET :: WT_detg, WT_temp, WT_sigma2, WT_elld, WT_r0, WT_r1 
       TARGET :: WT_r1_p, WT_r1_p_p 
       TARGET :: qa_2, qa_3, x_wt, j_wt, j_l
       TARGET :: beta_wt, beta_l, u_wt, u_l, q_wt, w_wt, w_l
       DECLARE_CCTK_ARGUMENTS
     end subroutine NullSHRE_Pointing
  end interface


  call NullSHRE_Pointing(CCTK_PASS_FTOF)
  
  alpha_n%d => SHRE_alpha_p(:,:,1)
  alpha_s%d => SHRE_alpha_p(:,:,2)
  dalpha_n(1)%d => SHRE_dralpha_p(:,:,1)
  dalpha_s(1)%d => SHRE_dralpha_p(:,:,2)
  dalpha_n(2)%d => SHRE_dqalpha_p(:,:,1)
  dalpha_s(2)%d => SHRE_dqalpha_p(:,:,2)
  dalpha_n(3)%d => SHRE_dpalpha_p(:,:,1)
  dalpha_s(3)%d => SHRE_dpalpha_p(:,:,2)
  dalpha_n(4)%d => SHRE_dtalpha_p(:,:,1)
  dalpha_s(4)%d => SHRE_dtalpha_p(:,:,2)

  call point_v1(null_lsh, 6_ik, 1_ik, 3_ik, beta_n, beta_s, SHRE_beta_p)
  call point_v1(null_lsh, 6_ik, 1_ik, 3_ik, dbeta_n(:,1), dbeta_s(:,1), SHRE_drbeta_p)
  call point_v1(null_lsh, 6_ik, 1_ik, 3_ik, dbeta_n(:,2), dbeta_s(:,2), SHRE_dqbeta_p)
  call point_v1(null_lsh, 6_ik, 1_ik, 3_ik, dbeta_n(:,3), dbeta_s(:,3), SHRE_dpbeta_p)
  call point_v1(null_lsh, 6_ik, 1_ik, 3_ik, dbeta_n(:,4), dbeta_s(:,4), SHRE_dtbeta_p)

  call point_v2sym(null_lsh, 12_ik, 1_ik, 3_ik,&
                   g_n(1:3_ik,1:3), g_s(1:3_ik,1:3), SHRE_gij_p)
  call point_v2sym(null_lsh, 12_ik, 1_ik, 3_ik,&
                   dg_n(1:3_ik,1:3_ik,1), dg_s(1:3_ik,1:3_ik,1), SHRE_drgij_p)
  call point_v2sym(null_lsh, 12_ik, 1_ik, 3_ik,&
                   dg_n(1:3_ik,1:3_ik,2), dg_s(1:3_ik,1:3_ik,2), SHRE_dqgij_p)
  call point_v2sym(null_lsh, 12_ik, 1_ik, 3_ik,&
                   dg_n(1:3_ik,1:3_ik,3), dg_s(1:3_ik,1:3_ik,3), SHRE_dpgij_p)
  call point_v2sym(null_lsh, 12_ik, 1_ik, 3_ik,&
                   dg_n(1:3_ik,1:3_ik,4), dg_s(1:3_ik,1:3_ik,4), SHRE_dtgij_p)

  call point_v1(null_lsh, 8_ik, 1_ik, 4_ik, g_n(:,4), g_s(:,4), SHRE_git_p)
  call point_v1(null_lsh, 8_ik, 1_ik, 4_ik, dg_n(:,4_ik,1), dg_s(:,4_ik,1), SHRE_drgit_p)
  call point_v1(null_lsh, 8_ik, 1_ik, 4_ik, dg_n(:,4_ik,2), dg_s(:,4_ik,2), SHRE_dqgit_p)
  call point_v1(null_lsh, 8_ik, 1_ik, 4_ik, dg_n(:,4_ik,3), dg_s(:,4_ik,3), SHRE_dpgit_p)
  call point_v1(null_lsh, 8_ik, 1_ik, 4_ik, dg_n(:,4_ik,4), dg_s(:,4_ik,4), SHRE_dtgit_p)

  dr0_n(1)%d =>WT_r1_p(:,:,1)
  dr0_s(1)%d =>WT_r1_p(:,:,2)

contains

#include "NullSHRE_PointingLib.h"

end subroutine NullSHRE_Pointing_p

 
subroutine NullSHRE_Pointing_p_p (CCTK_ARGUMENTS)

  use cctk
  use NullSHRE_modVars
  implicit none

  TARGET :: SHRE_alpha_p_p, WT_r1_p_p
  TARGET :: SHRE_dralpha_p_p, SHRE_dqalpha_p_p, SHRE_dpalpha_p_p, SHRE_dtalpha_p_p

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
 
  CCTK_INT, parameter :: izero = 0
  integer, parameter :: ik = kind(izero)


  interface
     subroutine NullSHRE_Pointing (CCTK_ARGUMENTS)
       implicit none
       TARGET :: SHRE_alpha
       TARGET :: SHRE_dralpha, SHRE_dqalpha, SHRE_dpalpha, SHRE_dtalpha 
       TARGET :: WT_detg, WT_temp, WT_sigma2, WT_elld, WT_r0, WT_r1 
       TARGET :: WT_r1_p, WT_r1_p_p 
       TARGET :: qa_2, qa_3, x_wt, j_wt, j_l
       TARGET :: beta_wt, beta_l, u_wt, u_l, q_wt, w_wt, w_l
       DECLARE_CCTK_ARGUMENTS
     end subroutine NullSHRE_Pointing
  end interface


  call NullSHRE_Pointing(CCTK_PASS_FTOF)

  alpha_n%d => SHRE_alpha_p_p(:,:,1)
  alpha_s%d => SHRE_alpha_p_p(:,:,2)
  dalpha_n(1)%d => SHRE_dralpha_p_p(:,:,1)
  dalpha_s(1)%d => SHRE_dralpha_p_p(:,:,2)
  dalpha_n(2)%d => SHRE_dqalpha_p_p(:,:,1)
  dalpha_s(2)%d => SHRE_dqalpha_p_p(:,:,2)
  dalpha_n(3)%d => SHRE_dpalpha_p_p(:,:,1)
  dalpha_s(3)%d => SHRE_dpalpha_p_p(:,:,2)
  dalpha_n(4)%d => SHRE_dtalpha_p_p(:,:,1)
  dalpha_s(4)%d => SHRE_dtalpha_p_p(:,:,2)

  call point_v1(null_lsh, 6_ik, 1_ik, 3_ik, beta_n, beta_s, SHRE_beta_p_p)
  call point_v1(null_lsh, 6_ik, 1_ik, 3_ik, dbeta_n(:,1), dbeta_s(:,1), SHRE_drbeta_p_p)
  call point_v1(null_lsh, 6_ik, 1_ik, 3_ik, dbeta_n(:,2), dbeta_s(:,2), SHRE_dqbeta_p_p)
  call point_v1(null_lsh, 6_ik, 1_ik, 3_ik, dbeta_n(:,3), dbeta_s(:,3), SHRE_dpbeta_p_p)
  call point_v1(null_lsh, 6_ik, 1_ik, 3_ik, dbeta_n(:,4), dbeta_s(:,4), SHRE_dtbeta_p_p)

  call point_v2sym(null_lsh, 12_ik, 1_ik, 3_ik,&
                   g_n(1:3_ik,1:3), g_s(1:3_ik,1:3), SHRE_gij_p_p)
  call point_v2sym(null_lsh, 12_ik, 1_ik, 3_ik,&
                   dg_n(1:3_ik,1:3_ik,1), dg_s(1:3_ik,1:3_ik,1), SHRE_drgij_p_p)
  call point_v2sym(null_lsh, 12_ik, 1_ik, 3_ik,&
                   dg_n(1:3_ik,1:3_ik,2), dg_s(1:3_ik,1:3_ik,2), SHRE_dqgij_p_p)
  call point_v2sym(null_lsh, 12_ik, 1_ik, 3_ik,&
                   dg_n(1:3_ik,1:3_ik,3), dg_s(1:3_ik,1:3_ik,3), SHRE_dpgij_p_p)
  call point_v2sym(null_lsh, 12_ik, 1_ik, 3_ik,&
                   dg_n(1:3_ik,1:3_ik,4), dg_s(1:3_ik,1:3_ik,4), SHRE_dtgij_p_p)

  call point_v1(null_lsh, 8_ik, 1_ik, 4_ik, g_n(:,4), g_s(:,4), SHRE_git_p_p)
  call point_v1(null_lsh, 8_ik, 1_ik, 4_ik, dg_n(:,4_ik,1), dg_s(:,4_ik,1), SHRE_drgit_p_p)
  call point_v1(null_lsh, 8_ik, 1_ik, 4_ik, dg_n(:,4_ik,2), dg_s(:,4_ik,2), SHRE_dqgit_p_p)
  call point_v1(null_lsh, 8_ik, 1_ik, 4_ik, dg_n(:,4_ik,3), dg_s(:,4_ik,3), SHRE_dpgit_p_p)
  call point_v1(null_lsh, 8_ik, 1_ik, 4_ik, dg_n(:,4_ik,4), dg_s(:,4_ik,4), SHRE_dtgit_p_p)

  dr0_n(1)%d =>WT_r1_p_p(:,:,1)
  dr0_s(1)%d =>WT_r1_p_p(:,:,2)

contains

#include "NullSHRE_PointingLib.h"

end subroutine NullSHRE_Pointing_p_p
