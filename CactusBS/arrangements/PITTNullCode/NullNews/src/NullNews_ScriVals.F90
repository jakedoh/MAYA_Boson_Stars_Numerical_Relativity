! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine NullNews_ScriVals (CCTK_ARGUMENTS)
  use NullNews_ScriUtil
  use NullDecomp_SpinDecomp, only: SpinDecompFilter

  implicit none

  LOGICAL, SAVE :: FirstTime = .true.
  CCTK_COMPLEX, DIMENSION(:,:,:), ALLOCATABLE,  SAVE :: temp

  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, parameter :: izero = 0
  integer, parameter :: ik = kind(izero)

  !   call CCTK_INFO("Null News, obtaining scrivals")

  if (filter_scri_fields.ne.0) then
     if (FirstTime) then
        FirstTime = .false.
        allocate(temp(null_lsh(1), null_lsh(2), 2))
     endif
  endif

  ! get scri quantities from null code variables

  call NullNews_cscrival (null_lsh, N_radial_pts, jcn_p, jcs_p, Jo)
  call NullNews_cscrival (null_lsh, N_radial_pts, jcn, jcs, Jn)

  call NullNews_cscridbydl (Jl_deriv_order, null_lsh, N_radial_pts, null_dx, null_rwt, jcn_p, jcs_p, Jo_l)
  call NullNews_cscridbydl (Jl_deriv_order, null_lsh, N_radial_pts, null_dx, null_rwt, jcn, jcs, Jn_l)

  call NullNews_rscrival (null_lsh, N_radial_pts, bcn_p, bcs_p, betao)
  call NullNews_rscrival (null_lsh, N_radial_pts, bcn, bcs, betan)

  call NullNews_rscrival (null_lsh, N_radial_pts, wcn, wcs, Wn)

  if (first_order_scheme.ne.0) then
    call NullNews_cscrival (null_lsh, N_radial_pts, cbcn_p, cbcs_p, cBo)
    call NullNews_cscrival (null_lsh, N_radial_pts, cbcn, cbcs, cBn)
  end if

  call NullNews_cscrival (null_lsh, N_radial_pts, qcn_p, qcs_p, Qo)
  call NullNews_cscrival (null_lsh, N_radial_pts, qcn, qcs, Qn)

  call NullNews_cscridbydl (Jl_deriv_order, null_lsh, N_radial_pts, null_dx, null_rwt, qcn_p, qcs_p, Qo_l)
  call NullNews_cscridbydl (Jl_deriv_order, null_lsh, N_radial_pts, null_dx, null_rwt, qcn, qcs, Qn_l)

  call NullNews_cscrivalh (null_lsh, N_radial_pts, ucn_p, ucs_p, Uo)
  call NullNews_cscrivalh (null_lsh, N_radial_pts, ucn, ucs, Un)

  ! FIXME -- these can be computed using Q, Q_l
  call NullNews_cscridbydlh (null_lsh, N_radial_pts, null_dx, null_rwt, ucn_p, ucs_p, Uo_l)
  call NullNews_cscridbydlh (null_lsh, N_radial_pts, null_dx, null_rwt, ucn, ucs, Un_l)

  call NullNews_cscridbydl2h (null_lsh, N_radial_pts, null_dx, null_rwt, ucn_p, ucs_p, Uo_l_l)
  call NullNews_cscridbydl2h (null_lsh, N_radial_pts, null_dx, null_rwt, ucn, ucs, Un_l_l)

  if (filter_scri_fields .ne. 0) then
     call SpinDecompFilter(cctkGH, null_lsh(1), null_lsh(2), 2_ik, zeta, Jo)
     call SpinDecompFilter(cctkGH, null_lsh(1), null_lsh(2), 2_ik, zeta, Jn)

     call SpinDecompFilter(cctkGH, null_lsh(1), null_lsh(2), 2_ik, zeta, Jo_l)
     call SpinDecompFilter(cctkGH, null_lsh(1), null_lsh(2), 2_ik, zeta, Jn_l)

     temp = betao
     call SpinDecompFilter(cctkGH, null_lsh(1), null_lsh(2), 0_ik, zeta, temp)
     betao = dble(temp)

     temp = betan
     call SpinDecompFilter(cctkGH, null_lsh(1), null_lsh(2), 0_ik, zeta, temp)
     betan = dble(temp)

     if (first_order_scheme.ne.0) then
       call SpinDecompFilter(cctkGH, null_lsh(1), null_lsh(2), 1_ik, zeta, cBo)
       call SpinDecompFilter(cctkGH, null_lsh(1), null_lsh(2), 1_ik, zeta, cBn)
     end if

     call SpinDecompFilter(cctkGH, null_lsh(1), null_lsh(2), 1_ik, zeta, Qo)
     call SpinDecompFilter(cctkGH, null_lsh(1), null_lsh(2), 1_ik, zeta, Qn)
     call SpinDecompFilter(cctkGH, null_lsh(1), null_lsh(2), 2_ik, zeta, Qo_l)
     call SpinDecompFilter(cctkGH, null_lsh(1), null_lsh(2), 2_ik, zeta, Qn_l)

     call SpinDecompFilter(cctkGH, null_lsh(1), null_lsh(2), 1_ik, zeta, Uo)
     call SpinDecompFilter(cctkGH, null_lsh(1), null_lsh(2), 1_ik, zeta, Un)

     !For psi4
     call SpinDecompFilter(cctkGH, null_lsh(1), null_lsh(2), 2_ik, zeta, Uo_l)
     call SpinDecompFilter(cctkGH, null_lsh(1), null_lsh(2), 2_ik, zeta, Un_l)

     call SpinDecompFilter(cctkGH, null_lsh(1), null_lsh(2), 2_ik, zeta, Uo_l_l)
     call SpinDecompFilter(cctkGH, null_lsh(1), null_lsh(2), 2_ik, zeta, Un_l_l)


  endif

end subroutine NullNews_ScriVals

