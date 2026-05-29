#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"


subroutine wt_WriteSphHarmWT(CCTK_ARGUMENTS)

  use cctk
  use NullInterp 
  use NullDecomp_Vars, only: lmax
  use NullDecomp_SpinDecomp, only: SpinDecompCoefs
  use NullDecomp_IO

  implicit none
 
  CCTK_INT, parameter :: izero = 0
  integer, parameter :: ik = kind(izero)

  CCTK_COMPLEX, dimension (:,:,:), allocatable:: cmplxtmp
  CCTK_COMPLEX, dimension(2:lmax, -lmax:lmax) :: JSphCoeff
  CCTK_COMPLEX, dimension(1:lmax, -lmax:lmax) :: USphCoeff
  CCTK_COMPLEX, dimension(0:lmax, -lmax:lmax) :: betaSphCoeff, WSphCoeff
  CCTK_INT :: Nradial
  CCTK_INT :: i_max, gi_max, i_min, gi_min, retval
  integer, save :: reduce_handle
  logical, save :: first_time = .TRUE.
  logical :: truncate
  
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS


     call CCTK_ReductionArrayHandle(reduce_handle, "minimum")
     if (reduce_handle .lt. 0 ) then
        call CCTK_WARN(0,"Could not get reduction handle")
     endif

     i_min =min(minval(boundary_maskn),minval(boundary_masks))

     call CCTK_ReduceLocScalar(retval, cctkGH, -1, reduce_handle,&
          i_min, gi_min, CCTK_VARIABLE_INT)

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

     Nradial = gi_max


  ! if (decomp == 1) then
      allocate(cmplxtmp(null_lsh(1), null_lsh(2), 2) )
      
      cmplxtmp(:,:,1) = jcn(:,:,Nradial)
      cmplxtmp(:,:,2) = jcs(:,:,Nradial)
      call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), 2_ik, &
         zeta, cmplxtmp, JSphCoeff)

      cmplxtmp(:,:,1) = ucn(:,:,Nradial)
      cmplxtmp(:,:,2) = ucs(:,:,Nradial)
      call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), 1_ik, &
         zeta, cmplxtmp, USphCoeff)
      
      cmplxtmp(:,:,1) = wcn(:,:,Nradial)
      cmplxtmp(:,:,2) = wcs(:,:,Nradial)
      call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), 0_ik, &
         zeta, cmplxtmp, WSphCoeff)

      cmplxtmp(:,:,1) = bcn(:,:,Nradial)
      cmplxtmp(:,:,2) = bcs(:,:,Nradial)
      call SpinDecompCoefs(cctkGH, null_lsh(1), null_lsh(2), 0_ik, &
         zeta, cmplxtmp, betaSphCoeff)
      
      deallocate(cmplxtmp)
  ! endif


  ! Only output if this is processor 0
  !  if (MyProc /= 0) return 
  if (CCTK_MyProc(cctkGH) .gt. 0) then
     return
  endif

  ! NullDecomp_WriteCoefFile(FieldName, TruncateFile, Lmin, Lmax, YlmCoef, time)

  truncate = (IO_TruncateOutputFiles(cctkGH) .ne. 0) .and. first_time

  call  NullDecomp_WriteCoefRadiusFile('nullgridJ_wt',   truncate, 2_ik, lmax, JSphCoeff, cctk_time, null_rb(Nradial))

  call  NullDecomp_WriteCoefRadiusFile('nullgridU_wt',   truncate, 1_ik, lmax, USphCoeff, cctk_time, null_rb(Nradial))

  call  NullDecomp_WriteCoefRadiusFile('nullgridB_wt',   truncate, 0_ik, lmax, betaSphCoeff, cctk_time, null_rb(Nradial))

  call  NullDecomp_WriteCoefRadiusFile('nullgridW_wt',   truncate, 0_ik, lmax, WSphCoeff, cctk_time, null_rb(Nradial))

  first_time = .FALSE.


end subroutine wt_WriteSphHarmWT



