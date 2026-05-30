! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine NullNews_GetNews(CCTK_ARGUMENTS)
  use NullNews_Omega
  use NullNews_Bondi
  use NullNews_CalcPsi4
  use NullNews_ScriUtil, only: NullNews_ResetInactive
  use NullNews_ScriUtil, only: NullNews_ResetInactiveRe
  use NullInterp 
  use NullDecomp_SpinDecomp, only: SpinDecompFilter

  implicit none

  ! local arrays

    CCTK_INT, parameter :: izero = 0
    integer, parameter :: ik = kind(izero)

  CCTK_COMPLEX, dimension (:,:,:), allocatable, save::&
       J, J_u, J_l, J_l_u,  cB, U, comega, eth_J,&
       beth_J, eth_U, beth_U, U_u, U_l, U_l_l, etheth_u0

  CCTK_REAL,    dimension (:,:,:), allocatable, save::&
       beta, omega, beta_u, u0

  logical, save :: FirstTime = .True.
  integer, save :: myproc

  CCTK_INT :: ip

  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  !   call CCTK_INFO("Null News, Calculating News")
  if (FirstTime) then
     FirstTime = .false.
     allocate(   J(null_lsh(1), null_lsh(2), 2),&
          J_u(null_lsh(1), null_lsh(2), 2),&
          J_l(null_lsh(1), null_lsh(2), 2),&
          J_l_u(null_lsh(1), null_lsh(2), 2),&
          cB(null_lsh(1), null_lsh(2), 2),&
          U(null_lsh(1), null_lsh(2), 2),&
          comega(null_lsh(1), null_lsh(2), 2),&
          eth_J(null_lsh(1), null_lsh(2), 2),&
          beth_J(null_lsh(1), null_lsh(2), 2),&
          eth_U(null_lsh(1), null_lsh(2), 2),&
          beth_U(null_lsh(1), null_lsh(2), 2),&
          beta(null_lsh(1), null_lsh(2), 2),&
          beta_u(null_lsh(1), null_lsh(2), 2),&
          U_u(null_lsh(1), null_lsh(2), 2),&
          U_l(null_lsh(1), null_lsh(2), 2),&
          U_l_l(null_lsh(1), null_lsh(2), 2),&
          omega(null_lsh(1), null_lsh(2), 2)    )
     J=0; J_u=0; J_l=0; J_l_u=0; cB=0; U=0; comega=0;
     eth_J=0; beth_J=0; eth_U=0; beth_U=0; beta=0; omega=0
     beta_u = 0; U_u=0; U_l=0; U_l_l=0

     if (compute_lin_strain .ne. 0) then
        allocate(u0(null_lsh(1), null_lsh(2), 2), etheth_u0(null_lsh(1), null_lsh(2), 2))
        u0=0
        etheth_u0=0
     endif
     
     myproc = CCTK_MyProc(cctkGH)     
  endif
  if (use_linearized_omega .eq. 0) then
    call NullNews_integ_omega (cctkGH, null_lsh, cctk_delta_time, tmp_rgfn, tmp_rgfs, Uo, Un, omegao, omegan,betao,betan)
    call NullNews_integ_comega (cctkGH, null_lsh, cctk_delta_time, tmp_cgfn, tmp_cgfs, omegao, omegan, comegao, comegan)
  else 
    call NullNews_lin_omega (cctkGH, null_lsh, zeta, Jn, omegan, comegan)
  endif

  omega = 0.5d0 * (omegan + omegao)
  comega = 0.5d0 * (comegan + comegao)

  call NullNews_ResetInactiveRe(null_lsh, omega, 1.0d0)
  call NullNews_ResetInactive(null_lsh, comega)

  if (filter_omega .ne. 0) then
     J = omega  ! just for temp space
     call SpinDecompFilter(cctkGH, null_lsh(1), null_lsh(2), 0_ik, zeta, J)
     omega = dble(J)
     call SpinDecompFilter(cctkGH, null_lsh(1), null_lsh(2), 1_ik, zeta, comega)
  endif



  ! mid level values, from the characteristic evolution

  J = 0.5d0 * (Jn + Jo)
  J_u = (Jn - Jo) / cctk_delta_time
  J_l = 0.5d0 * (Jn_l + Jo_l)
  J_l_u = (Jn_l - Jo_l) / cctk_delta_time
  beta = 0.5d0 * (betan + betao)
  cB = 0.5d0 * (cBn + cBo)
  U = 0.5d0 * (Un + Uo)

  ! For psi4 calculation
  beta_u = (betan - betao) / cctk_delta_time
  U_u = (Un - Uo) / cctk_delta_time
  U_l = 0.5d0 * (Un_l + Uo_l)
  U_l_l = 0.5d0 * (Un_l_l + Uo_l_l)

  call NullInterp_d1 (eth_U(:,:,1), U(:,:,1), 1_ik, 1_ik)
  call NullInterp_d1 (eth_U(:,:,2), U(:,:,2), 1_ik, 1_ik)

  call NullInterp_d1 (beth_U(:,:,1), U(:,:,1), 1_ik, -1_ik)
  call NullInterp_d1 (beth_U(:,:,2), U(:,:,2), 1_ik, -1_ik)

  call NullInterp_d1 (eth_J(:,:,1), J(:,:,1), 2_ik, 1_ik)
  call NullInterp_d1 (eth_J(:,:,2), J(:,:,2), 2_ik, 1_ik)

  call NullInterp_d1 (beth_J(:,:,1), J(:,:,1), 2_ik, -1_ik)
  call NullInterp_d1 (beth_J(:,:,2), J(:,:,2), 2_ik, -1_ik)


  ! do we really need to do the syncing here .. probably

  call NullInterp_3cinterp(cctkGH,&
             tmp_cgfn1, tmp_cgfs1,&
             tmp_cgfn2, tmp_cgfs2,&
             tmp_cgfn3, tmp_cgfs3,&
             eth_U(:,:,1), eth_U(:,:,2),&
             beth_U(:,:,1), beth_U(:,:,2),&
             eth_J(:,:,1), eth_J(:,:,2),&
             2_ik, 0_ik, 3_ik)

! call NullInterp_cinterp(cctkGH, tmp_cgfn, tmp_cgfs, eth_U(:,:,1), eth_U(:,:,2), 2)
! call NullInterp_cinterp(cctkGH, tmp_cgfn, tmp_cgfs, beth_U(:,:,1), beth_U(:,:,2), 0)
! call NullInterp_cinterp(cctkGH, tmp_cgfn, tmp_cgfs, eth_J(:,:,1), eth_J(:,:,2), 3)
  call NullInterp_cinterp(cctkGH, tmp_cgfn, tmp_cgfs, beth_J(:,:,1), beth_J(:,:,2), 1_ik)

  call NullNews_ResetInactive(null_lsh, eth_U)
  call NullNews_ResetInactive(null_lsh, beth_U)
  call NullNews_ResetInactive(null_lsh, eth_J)
  call NullNews_ResetInactive(null_lsh, beth_J)

  ! assumes 'time' is time on the new level
  time_of_news = cctk_time - .5d0 * cctk_delta_time    

  ! mid level values, from quantities computed above

  sigmaJo = sigmaJn; sigmaKo = sigmaKn; sigmauo = sigmaun
  sigmaro = sigmarn; sigmaruo = sigmarun; sigmauuo = sigmauun

  call NullNews_uframe (cctkGH, null_lsh,&
       tmp_cgfn1, tmp_cgfs1, tmp_cgfn2, tmp_cgfs2, tmp_cgfn3, tmp_cgfs3,&
       J, J_u, J_l, J_l_u,& 
       beta, cB, U, omega, comega, News(:,:,1:2), sigmaJn, sigmaKn, sigmaun,&
       sigmarn, sigmarun, sigmauun,beta_u,U_u,U_l,U_l_l)

  call NullNews_Psi4 (cctkGH, null_lsh, tmp_cgfn, tmp_cgfs, Jo, Jo_l,betao, cBo, Uo, Uo_l,&
       sigmaJn, sigmaKn, sigmaun,sigmarn, sigmarun, sigmauun,&
       sigmaJo, sigmaKo, sigmauo,sigmaro, sigmaruo, sigmauuo,Psi4(:,:,1:2),cctk_delta_time)

  Psi4(:,:,1:2) = Psi4(:,:,1:2)/omega**3

  if (mask_Psi4 .ne. 0) then
     do ip = 1, 2
        Psi4(:,:,ip) = EQ_mask(:,:) * Psi4(:,:,ip) ! EV_mask
     end do
  end if

  if (filter_news .ne. 0) then
     call SpinDecompFilter(cctkGH, null_lsh(1), null_lsh(2), 2_ik, zeta, News(:,:,1:2))
  endif

  if (DEBUG_skip_BondiNews.eq.0) &
       call NullNews_News2BondiNews(cctkGH, null_lsh, cctk_delta_time, qsize,&
       Uo, Un, Uyo, Uyn, zEvolo, zEvoln,&
       zEvolh, News(:,:,1:2), NewsB(:,:,1:2), uBondio, uBondi, uBondin,&
       redshiftB, omega, J, eth_J, beth_J, eth_U, beth_U,&
       J_u, beta, patch_index, deltao, deltan, nzeta,&
       n_tmp1_cgfn, n_tmp1_cgfs, n_tmp2_cgfn, n_tmp2_cgfs,&
       n_tmp3_cgfn, n_tmp3_cgfs, n_tmp1_rgfn, n_tmp1_rgfs,&
       n_tmp2_rgfn, n_tmp2_rgfs, starting, EV_mask)

  dMdOmega = NewsB(:,:,1:2) * conjg(NewsB(:,:,1:2)) * redshiftB * circle

  if (mask_NewsB .ne. 0) then
     do ip = 1, 2
        NewsB(:,:,ip) = EQ_mask(:,:) * NewsB(:,:,ip) ! EV_mask
     end do
  end if

  if (compute_lin_strain .ne. 0) then
     u0 = uBondi
     call NullInterp_d2 (etheth_u0(:,:,1), dcmplx(u0(:,:,1)), 0_ik, 1_ik, 1_ik)
     call NullInterp_d2 (etheth_u0(:,:,2), dcmplx(u0(:,:,2)), 0_ik, 1_ik, 1_ik)
     call NullInterp_cinterp(cctkGH, tmp_cgfn, tmp_cgfs, etheth_u0(:,:,1), etheth_u0(:,:,2), 2_ik)
     call NullNews_ResetInactive(null_lsh, etheth_u0)
     linStrain(:,:,1) = J_l(:,:,1) + tmp_cgfn
     linStrain(:,:,2) = J_l(:,:,2) + tmp_cgfs
  endif
  
  !! at this point its counterproductive to change spin-basis
  !   ! La definicion es z=\tan(\theta/2) e^{i\phi}, asi que aca le 
  !   ! sacamos el z/zb que le pusimos en Spheroid_RJ.
  !
  !   ! Translation: the definition is z=\tan(\theta/2) e^{i\phi}, thus here
  !   ! we take out the factor of zb/z that we stuck J with in Spheroid_RJ
  !
  !   where (nzeta .ne. (0.0, 0.0))
  !      NewsB(:,:,1:2) = NewsB(:,:,1:2) * nzetabar / nzeta
  !   end where

  omegao = omegan
  comegao = comegan


end subroutine NullNews_GetNews



