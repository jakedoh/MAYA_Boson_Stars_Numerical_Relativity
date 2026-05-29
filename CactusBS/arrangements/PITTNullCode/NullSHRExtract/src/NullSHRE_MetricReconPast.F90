! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"


 subroutine NullSHRE_MetricReconPast(CCTK_ARGUMENTS)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!reads the numerical expansion coefficients Clm
!and reconstructs the metric on the worldtube 
!with the spherical harmonics in (p,q) coord.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use cctk
  use NullGrid_Vars
  implicit none

    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    DECLARE_CCTK_FUNCTIONS

    SHRE_gij_p_p     = SHRE_gij
    SHRE_git_p_p     = SHRE_git
    SHRE_drgij_p_p   = SHRE_drgij
    SHRE_dqgij_p_p   = SHRE_dqgij
    SHRE_dpgij_p_p   = SHRE_dpgij
    SHRE_dtgij_p_p   = SHRE_dtgij
    SHRE_drgit_p_p   = SHRE_drgit
    SHRE_dqgit_p_p   = SHRE_dqgit
    SHRE_dpgit_p_p   = SHRE_dpgit
    SHRE_dtgit_p_p   = SHRE_dtgit
    SHRE_beta_p_p    = SHRE_beta
    SHRE_drbeta_p_p  = SHRE_drbeta
    SHRE_dqbeta_p_p  = SHRE_dqbeta
    SHRE_dpbeta_p_p  = SHRE_dpbeta
    SHRE_dtbeta_p_p  = SHRE_dtbeta
    SHRE_alpha_p_p   = SHRE_alpha
    SHRE_dralpha_p_p = SHRE_dralpha
    SHRE_dqalpha_p_p = SHRE_dqalpha
    SHRE_dpalpha_p_p = SHRE_dpalpha
    SHRE_dtalpha_p_p = SHRE_dtalpha

    SHRE_gij_p     = SHRE_gij
    SHRE_git_p     = SHRE_git
    SHRE_drgij_p   = SHRE_drgij
    SHRE_dqgij_p   = SHRE_dqgij
    SHRE_dpgij_p   = SHRE_dpgij
    SHRE_dtgij_p   = SHRE_dtgij
    SHRE_drgit_p   = SHRE_drgit
    SHRE_dqgit_p   = SHRE_dqgit
    SHRE_dpgit_p   = SHRE_dpgit
    SHRE_dtgit_p   = SHRE_dtgit
    SHRE_beta_p    = SHRE_beta
    SHRE_drbeta_p  = SHRE_drbeta
    SHRE_dqbeta_p  = SHRE_dqbeta
    SHRE_dpbeta_p  = SHRE_dpbeta
    SHRE_dtbeta_p  = SHRE_dtbeta
    SHRE_alpha_p   = SHRE_alpha
    SHRE_dralpha_p = SHRE_dralpha
    SHRE_dqalpha_p = SHRE_dqalpha
    SHRE_dpalpha_p = SHRE_dpalpha
    SHRE_dtalpha_p = SHRE_dtalpha

 end subroutine NullSHRE_MetricReconPast

