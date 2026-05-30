! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine NullNews_Init(CCTK_ARGUMENTS)
  use NullNews_ScriUtil
  use NullDecomp_SpinDecomp, only: SpinDecompFilter

  implicit none

  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, parameter :: izero = 0
  integer, parameter :: ik = kind(izero)

  CCTK_INT :: ip, i, j
  CCTK_COMPLEX :: temp(null_lsh(1),null_lsh(2),2)

  call CCTK_INFO("Null News setup")

  nzeta = 0     ! some dummy value for the non-existant points
  nzetabar = 0  ! some dummy value for the non-existant points

  do ip = 1, 2
     nzeta(1:null_lsh(1),1:null_lsh(2),ip) = zeta
     nzetabar(1:null_lsh(1),1:null_lsh(2),ip) = conjg(zeta)
  end do

  P = 1.0d0 + nzeta * nzetabar

  omegao = 1.0d0
  comegao = (0.0d0, 0.0d0)
  deltao = 0.0d0      ! deltao and mask will be changes 
  uBondio = cctk_time  ! = null_time    ! prior to the first run of news2bondinews
  patch_index = 0          ! see news2bondinews
  !  uBondio =  0.0 ! -- WE DONT UNDERSTAND THIS

  zEvolo =(.0,.0) ! a dummy value for those points not on the
  ! physical grid

  do ip = 1, 2
     zEvolo(1:null_lsh(1),1:null_lsh(2),ip) = zeta
  end do


  circle = 0.0
  do ip = 1, 2
     do j = 1, null_lsh(2)
        do i = 1, null_lsh(1)
           if (stereo_pp(i,j) < 2.0d0) circle(i,j,ip) = 1.0d0
        end do
     end do
  end do

  starting = 1  !tell news_int that this is the first time level 

  NewsB = 0.0d0; uBondi(:,:,1:2) = 0.0d0; time_of_news = 0.0d0
  dMdOmega = 0.0d0; betao = 0.0d0; betan = 0.0d0; omegan = 0.0d0
  redshiftB = 0.0d0; uBondin = 0.0d0; deltan = 0.0d0; comegan = 0.0d0


  News = 0.0d0; zEvoln = 0.0d0; zEvolh = 0.0d0
  sigmaJo = 0.0d0;sigmaJn = 0.0d0;sigmaKo = 0.0d0;sigmaKn = 0.0d0
  sigmauo = 0.0d0;sigmaun = 0.0d0;sigmaro = 0.0d0;sigmarn = 0.0d0
  sigmaruo = 0.0d0;sigmarun = 0.0d0;sigmauuo = 0.0d0
  sigmauun = 0.0d0

  Jo = 1.0d30; Jn = 1.0d30; Jo_l = 1.0d30; Jn_l = 1.0d30; cBo = 1.0d30
  cBn = 1.0d30; Uo = 1.0d30; Un = 1.0d30; Uyo = 1.0d30; Uyn = 1.0d30
  Uo_l=1.0d30; Uo_l_l=1.0d30;Un_l=1.0d30; Un_l_l=1.0d30
  Wn = 1.0d30



  ! get scri quantities from null code variables

  call NullNews_cscrival (null_lsh, N_radial_pts, jcn, jcs, Jo)
  call NullNews_cscrival (null_lsh, N_radial_pts, jcn, jcs, Jn)

  call NullNews_cscridbydl (Jl_deriv_order, null_lsh, N_radial_pts, null_dx, null_rwt, jcn, jcs, Jo_l)
  call NullNews_cscridbydl (Jl_deriv_order, null_lsh, N_radial_pts, null_dx, null_rwt, jcn, jcs, Jn_l)

  call NullNews_rscrival (null_lsh, N_radial_pts, bcn, bcs, betao)
  call NullNews_rscrival (null_lsh, N_radial_pts, bcn, bcs, betan)

  call NullNews_rscrival (null_lsh, N_radial_pts, wcn, wcs, Wn)

  if (first_order_scheme.ne.0) then
    call NullNews_cscrival (null_lsh, N_radial_pts, cbcn, cbcs, cBo)
    call NullNews_cscrival (null_lsh, N_radial_pts, cbcn, cbcs, cBn)
  end if

  call NullNews_cscrivalh (null_lsh, N_radial_pts, qcn, qcs, Qo)
  call NullNews_cscrivalh (null_lsh, N_radial_pts, qcn, qcs, Qn)

  call NullNews_cscridbydlh (null_lsh, N_radial_pts, null_dx, null_rwt, qcn, qcs, Qo_l)
  call NullNews_cscridbydlh (null_lsh, N_radial_pts, null_dx, null_rwt, qcn, qcs, Qn_l)

  call NullNews_cscrivalh (null_lsh, N_radial_pts, ucn, ucs, Uo)
  call NullNews_cscrivalh (null_lsh, N_radial_pts, ucn, ucs, Un)

  call NullNews_cscridbydlh (null_lsh, N_radial_pts, null_dx, null_rwt, ucn, ucs, Uo_l)
  call NullNews_cscridbydlh (null_lsh, N_radial_pts, null_dx, null_rwt, ucn, ucs, Un_l)

  call NullNews_cscridbydl2h (null_lsh, N_radial_pts, null_dx, null_rwt, ucn, ucs, Uo_l_l)
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

     call SpinDecompFilter(cctkGH, null_lsh(1), null_lsh(2), 1_ik, zeta, Uo)
     call SpinDecompFilter(cctkGH, null_lsh(1), null_lsh(2), 1_ik, zeta, Un)

     !For psi4
     call SpinDecompFilter(cctkGH, null_lsh(1), null_lsh(2), 2_ik, zeta, Uo_l)
     call SpinDecompFilter(cctkGH, null_lsh(1), null_lsh(2), 2_ik, zeta, Un_l)

     call SpinDecompFilter(cctkGH, null_lsh(1), null_lsh(2), 2_ik, zeta, Uo_l_l)
     call SpinDecompFilter(cctkGH, null_lsh(1), null_lsh(2), 2_ik, zeta, Un_l_l)


  endif


end subroutine NullNews_Init


