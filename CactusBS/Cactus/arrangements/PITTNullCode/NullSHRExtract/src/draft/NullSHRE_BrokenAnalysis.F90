!vim: syntax=fortran
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

 subroutine NullSHRE_Analysis(CCTK_ARGUMENTS)

  use cctk
  use NullGrid_Vars
  use NullSHRE_modVars
  use NullSHRE_mod4Metric
  use NullSHRE_modAnalytic
  
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT i,j
 
  call wt_g(alpha_n, beta_n, g_n)
  call wt_g(alpha_s, beta_s, g_s)

!  call modAna_g(null_lsh(1), null_lsh(2), qs, ps, pp, ip_n, cr, mass, errg_n)
!  call modAna_g(null_lsh(1), null_lsh(2), qs, ps, pp, ip_s, cr, mass, errg_s)

!  do i = 1, 4
!     do j = 1, 4
!
!        errg_n(i,j)%d = errg_n(i,j)%d - g_n(i,j)%d
!        errg_s(i,j)%d = errg_s(i,j)%d - g_s(i,j)%d
!
!     end do 
!  end do

!  call modAna_g(null_lsh(1), null_lsh(2), qs, ps, pp, ip_n, cr, mass, errg1_n)
!  call modAna_g(null_lsh(1), null_lsh(2), qs, ps, pp, ip_s, cr, mass, errg1_s)
!
!  do i = 1, 4
!     do j = 1, 4
!
!        errg1_n(i,j)%d = errg1_n(i,j)%d - g1_n(i,j)%d
!        errg1_s(i,j)%d = errg1_s(i,j)%d - g1_s(i,j)%d
!
!     end do 
!  end do

 end subroutine NullSHRE_Analysis 
