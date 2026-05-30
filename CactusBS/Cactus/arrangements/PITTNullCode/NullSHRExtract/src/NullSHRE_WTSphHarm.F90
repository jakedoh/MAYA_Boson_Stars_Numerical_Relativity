! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine NullSHRE_WTSphHarm(CCTK_ARGUMENTS)

  use cctk
  use NullInterp 
  use NullDecomp_Vars, only: lmax
  use NullDecomp_SpinDecomp, only: SpinDecompCoefs
  use NullDecomp_IO

  implicit none
 
  CCTK_INT, parameter :: izero = 0
  integer, parameter :: ik = kind(izero)

!variable of interest on the WorldTube x_wt, j_wt, j_l, beta_wt, beta_l, u_wt, u_l
  CCTK_COMPLEX, dimension (:,:,:), allocatable:: cmplxtmp

  CCTK_COMPLEX, dimension(2:lmax, -lmax:lmax) :: JwtSphCoeff
  CCTK_COMPLEX, dimension(2:lmax, -lmax:lmax) :: Jwt_lSphCoeff

  CCTK_COMPLEX, dimension(1:lmax, -lmax:lmax) :: QwtSphCoeff
  CCTK_COMPLEX, dimension(1:lmax, -lmax:lmax) :: UwtSphCoeff
  CCTK_COMPLEX, dimension(1:lmax, -lmax:lmax) :: Uwt_lSphCoeff

  CCTK_COMPLEX, dimension(0:lmax, -lmax:lmax) :: BwtSphCoeff
  CCTK_COMPLEX, dimension(0:lmax, -lmax:lmax) :: Bwt_lSphCoeff

  CCTK_COMPLEX, dimension(0:lmax, -lmax:lmax) :: WwtSphCoeff
  CCTK_COMPLEX, dimension(0:lmax, -lmax:lmax) :: Wwt_lSphCoeff

  CCTK_COMPLEX, dimension(0:lmax, -lmax:lmax) :: XwtSphCoeff
  CCTK_COMPLEX, dimension(0:lmax, -lmax:lmax) :: detgwtSphCoeff
  CCTK_COMPLEX, dimension(0:lmax, -lmax:lmax) :: r0wtSphCoeff

  logical, save :: first_time = .TRUE.
  logical :: truncate

  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

! do the spherical harmonic decomposition
  
  allocate(cmplxtmp(null_lsh(1), null_lsh(2), 2))     
 
  call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), 2_ik, &
     zeta, j_wt, JwtSphCoeff)
  call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), 2_ik, &
     zeta, j_l, Jwt_lSphCoeff)

  call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), 1_ik, &
     zeta, u_wt, UwtSphCoeff)
  call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), 1_ik, &
     zeta, u_l, Uwt_lSphCoeff)

  call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), 1_ik, &
     zeta, q_wt, QwtSphCoeff)

  cmplxtmp = beta_wt
  call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), 0_ik, &
     zeta, cmplxtmp, BwtSphCoeff)

  cmplxtmp = beta_l
  call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), 0_ik, &
     zeta, cmplxtmp, Bwt_lSphCoeff)

  cmplxtmp = w_wt
  call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), 0_ik, &
     zeta, cmplxtmp, WwtSphCoeff)

  cmplxtmp = w_l
  call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), 0_ik, &
     zeta, cmplxtmp, Wwt_lSphCoeff)

  cmplxtmp = x_wt
  call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), 0_ik, &
     zeta, cmplxtmp, XwtSphCoeff)

  cmplxtmp = WT_r0
  call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), 0_ik, &
     zeta, cmplxtmp, r0wtSphCoeff)

  cmplxtmp = WT_detg
  call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), 0_ik, &
     zeta, cmplxtmp, detgwtSphCoeff)
 
  deallocate(cmplxtmp)


! call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), 1_ik, &
!    zeta, cb_wt, CBwtSphCoeff)
! call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), 1_ik, &
!    zeta, cb_x, CBwt_xSphCoeff)

! call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), 1_ik, &
!    zeta, ck_wt, CKwtSphCoeff)
! call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), 1_ik, &
!    zeta, ck_x, CKwt_xSphCoeff)

! call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), 1_ik, &
!    zeta, nu_wt, NUwtSphCoeff)
! call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), 1_ik, &
!    zeta, nu_x, NUwt_xSphCoeff)

! Only output if this is processor 0
!  if (MyProc /= 0) return 
  if (CCTK_MyProc(cctkGH) .gt. 0) then
     return
  endif

  ! NullDecomp_WriteCoefFile(FieldName, TruncateFile, Lmin, Lmax, YlmCoef, time)

  truncate = (IO_TruncateOutputFiles(cctkGH) .ne. 0) .and. first_time

  call  NullDecomp_WriteCoefFile('J_wt',   truncate, 2_ik, Lmax, JwtSphCoeff, cctk_time)
  call  NullDecomp_WriteCoefFile('Jl_wt',  truncate, 2_ik, Lmax, Jwt_lSphCoeff, cctk_time)

  call  NullDecomp_WriteCoefFile('Q_wt',   truncate, 1_ik, Lmax, QwtSphCoeff, cctk_time)

  call  NullDecomp_WriteCoefFile('U_wt',   truncate, 1_ik, Lmax, UwtSphCoeff, cctk_time)
  call  NullDecomp_WriteCoefFile('Ul_wt',  truncate, 1_ik, Lmax, Uwt_lSphCoeff, cctk_time)

  call  NullDecomp_WriteCoefFile('B_wt',   truncate, 0_ik, Lmax, BwtSphCoeff, cctk_time)
  call  NullDecomp_WriteCoefFile('Bl_wt',  truncate, 0_ik, Lmax, Bwt_lSphCoeff, cctk_time)

  call  NullDecomp_WriteCoefFile('W_wt',   truncate, 0_ik, Lmax, WwtSphCoeff, cctk_time)
  call  NullDecomp_WriteCoefFile('Wl_wt',  truncate, 0_ik, Lmax, Wwt_lSphCoeff, cctk_time)

! if(first_order_scheme.ne.0) then
!   call  NullDecomp_WriteCoefFile('NU_wt',   cctk_iteration.eq.0, 1_ik, Lmax, NUwtSphCoeff, cctk_time)
!   write (*,*) 'NUwtSphCoeff=',maxval(abs(NUwtSphCoeff))
!   call  NullDecomp_WriteCoefFile('NUx_wt',  cctk_iteration.eq.0, 1_ik, Lmax, NUwt_xSphCoeff, cctk_time)

!   call  NullDecomp_WriteCoefFile('CK_wt',   cctk_iteration.eq.0, 1_ik, Lmax, CKwtSphCoeff, cctk_time)
!   call  NullDecomp_WriteCoefFile('CKx_wt',  cctk_iteration.eq.0, 1_ik, Lmax, CKwt_xSphCoeff, cctk_time)

!   call  NullDecomp_WriteCoefFile('CB_wt',   cctk_iteration.eq.0, 1_ik, Lmax, CBwtSphCoeff, cctk_time)
!   call  NullDecomp_WriteCoefFile('CBx_wt',  cctk_iteration.eq.0, 1_ik, Lmax, CBwt_xSphCoeff, cctk_time)
! end if

  call  NullDecomp_WriteCoefFile('X_wt',    truncate, 0_ik, Lmax, XwtSphCoeff, cctk_time)
  call  NullDecomp_WriteCoefFile('r0_wt',   truncate, 0_ik, Lmax, r0wtSphCoeff, cctk_time)
  call  NullDecomp_WriteCoefFile('detg_wt', truncate, 0_ik, Lmax, detgwtSphCoeff, cctk_time)

  first_time = .FALSE.

end subroutine NullSHRE_WTSphHarm



