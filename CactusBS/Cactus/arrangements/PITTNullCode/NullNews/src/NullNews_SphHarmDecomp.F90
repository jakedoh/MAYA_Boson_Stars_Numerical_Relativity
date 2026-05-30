! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"


subroutine NullNews_WriteSphHarm(CCTK_ARGUMENTS)

  use NullNews_Omega
  use NullNews_Bondi
  use NullInterp 
  use NullDecomp_Vars, only: lmax
  use NullDecomp_SpinDecomp, only: SpinDecompCoefs
  use NullDecomp_IO

  implicit none

  CCTK_INT, parameter :: izero = 0
  integer, parameter :: ik = kind(izero)

  CCTK_COMPLEX, dimension (:,:,:), allocatable:: cmplxtmp
  CCTK_COMPLEX, dimension (:,:,:), allocatable:: NewsError
  CCTK_COMPLEX, dimension(2:lmax, -lmax:lmax) :: NewsSphCoeff
  CCTK_COMPLEX, dimension(2:lmax, -lmax:lmax) :: NewsBSphCoeff
  CCTK_COMPLEX, dimension(2:lmax, -lmax:lmax) :: Psi4SphCoeff
  CCTK_COMPLEX, dimension(2:lmax, -lmax:lmax) :: JnSphCoeff
  CCTK_COMPLEX, dimension(2:lmax, -lmax:lmax) :: Jn_lSphCoeff
  CCTK_COMPLEX, dimension(1:lmax, -lmax:lmax) :: QnSphCoeff
  CCTK_COMPLEX, dimension(1:lmax, -lmax:lmax) :: Qn_lSphCoeff
  CCTK_COMPLEX, dimension(1:lmax, -lmax:lmax) :: UnSphCoeff
  CCTK_COMPLEX, dimension(1:lmax, -lmax:lmax) :: Un_lSphCoeff
  CCTK_COMPLEX, dimension(0:lmax, -lmax:lmax) :: omeganSphCoeff
  CCTK_COMPLEX, dimension(0:lmax, -lmax:lmax) :: betanSphCoeff
  CCTK_COMPLEX, dimension(0:lmax, -lmax:lmax) :: WnSphCoeff
  CCTK_COMPLEX, dimension(0:lmax, -lmax:lmax) :: uBondiSphCoeff
  CCTK_COMPLEX, dimension(2:lmax, -lmax:lmax) :: linStrainSphCoeff
  CCTK_REAL :: time

  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  logical, save :: first_time = .TRUE.
  logical :: truncate

  truncate = (IO_TruncateOutputFiles(cctkGH) .ne. 0) .and. first_time

!  write (*,*) "computing ylm coefs: time=", cctk_time, ", J=", maxval(abs(Jn))
! fill the waveforms with Ylm at Scri before decomposition
  if (Ylm_at_Scri == 1) then

     News = YlmScri_2
     NewsB = YlmScri_2
     Psi4 = YlmScri_2
     News_uBondi = YlmScri_2
     NewsB_uBondi = YlmScri_2
     Psi4_uBondi = YlmScri_2
     Jn = YlmScri_2
     Jn_l = YlmScri_2
     Qn = YlmScri_1
     Qn_l = YlmScri_1
     Un = YlmScri_1
     Un_l = YlmScri_1
     omegan = YlmScri_0
     deltan = YlmScri_0
     betan = YlmScri_0
     uBondi = YlmScri_0

  end if


  allocate(NewsError(null_lsh(1), null_lsh(2), 2) )

  ! SpinDecompCoefs(GH, nq, np, s, z, f, sphcoef)

  if (interp_to_constant_uBondi .eq. 0) then
      
      time = cctk_time
      
      call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), NewsSpinWeight, &
            zeta, News, NewsSphCoeff)
      call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), NewsSpinWeight, &
            zeta, NewsB, NewsBSphCoeff)
      call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), NewsSpinWeight, &
            zeta, Psi4, Psi4SphCoeff)
      
      if (compute_lin_strain .ne. 0) then
         call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), NewsSpinWeight, &
            zeta, linStrain, linStrainSphCoeff)
      endif
  
  else
      time = constant_uBondi
  
      call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), NewsSpinWeight, &
            zeta, News_uBondi, NewsSphCoeff)
      call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), NewsSpinWeight, &
            zeta, NewsB_uBondi, NewsBSphCoeff)
      call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), NewsSpinWeight, &
            zeta, Psi4_uBondi, Psi4SphCoeff)
      
      if (compute_lin_strain .ne. 0) then
         call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), NewsSpinWeight, &
            zeta, linStrain_uBondi, linStrainSphCoeff)
      endif
  endif

  if (debug_output == 1) then
      allocate(cmplxtmp(null_lsh(1), null_lsh(2), 2) )
  
      call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), 2_ik, &
         zeta, Jn, JnSphCoeff)
      call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), 2_ik, &
         zeta, Jn_l, Jn_lSphCoeff)
      call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), 1_ik, &
         zeta, Qn, QnSphCoeff)
      call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), 1_ik, &
         zeta, Qn_l, Qn_lSphCoeff)
      call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), 1_ik, &
         zeta, Un, UnSphCoeff)
      call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), 1_ik, &
         zeta, Un_l, Un_lSphCoeff)
      cmplxtmp = omegan
      call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), 0_ik, &
         zeta, cmplxtmp, omeganSphCoeff)
      cmplxtmp = betan
      call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), 0_ik, &
         zeta, cmplxtmp, betanSphCoeff)
      cmplxtmp = Wn
      call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), 0_ik, &
         zeta, cmplxtmp, WnSphCoeff)
      cmplxtmp = uBondi
      call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), 0_ik, &
         zeta, cmplxtmp, uBondiSphCoeff)
      
      deallocate(cmplxtmp)
  endif

  deallocate(NewsError)

  ! Only output if this is processor 0
  if (CCTK_MyProc(cctkGH) .gt. 0) then
     return
  endif
  
  ! subroutine NullDecomp_WriteCoefFile(FieldName, TruncateFile, Lmin, Lmax, YlmCoef, time)
  call NullDecomp_WriteCoefFile('News_scri', truncate, 2_ik, lmax, NewsSphCoeff, time)
  call NullDecomp_WriteCoefFile('NewsB_scri', truncate, 2_ik, lmax, NewsBSphCoeff, time)
  call NullDecomp_WriteCoefFile('Psi4_scri', truncate, 2_ik, lmax, Psi4SphCoeff, time)
  
  if (compute_lin_strain .ne. 0) then
     call NullDecomp_WriteCoefFile('linStrain_scri', truncate, 2_ik, lmax, linStrainSphCoeff, time)
  endif
  
  if (debug_output == 1) then
  
     call NullDecomp_WriteCoefFile('J_scri', truncate, 2_ik, lmax, JnSphCoeff, time)
     call NullDecomp_WriteCoefFile('Jl_scri', truncate, 2_ik, lmax, Jn_lSphCoeff, time)

     call NullDecomp_WriteCoefFile('U_scri', truncate, 1_ik, lmax, UnSphCoeff, time)
     call NullDecomp_WriteCoefFile('Ul_scri', truncate, 1_ik, lmax, Un_lSphCoeff, time)

     call NullDecomp_WriteCoefFile('Q_scri', truncate, 1_ik, lmax, QnSphCoeff, time)
     call NullDecomp_WriteCoefFile('Ql_scri', truncate, 1_ik, lmax, Qn_lSphCoeff, time)

     call NullDecomp_WriteCoefFile('B_scri', truncate, 0_ik, lmax, betanSphCoeff, time)
     call NullDecomp_WriteCoefFile('W_scri', truncate, 0_ik, lmax, WnSphCoeff, time)

     call NullDecomp_WriteCoefFile('omega_scri', truncate, 0_ik, lmax, omeganSphCoeff, time)
     call NullDecomp_WriteCoefFile('uBondi_scri', truncate, 0_ik, lmax, uBondiSphCoeff, time)

  endif


  first_time = .FALSE.

end subroutine NullNews_WriteSphHarm



