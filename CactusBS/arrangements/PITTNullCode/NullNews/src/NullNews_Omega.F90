! vim: syntax=fortran
#include "cctk.h"
!#include "cctk_Functions.h"
!#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

module NullNews_Omega
  implicit none
contains
  subroutine NullNews_integ_omega (cctkGH, lsh, dt, tmp_rgfn, tmp_rgfs, Uo, Un, omegao, omegan,&
       betao, betan)
    use NullInterp
    use NullGrid_Vars, only: delta
    implicit none

    CCTK_INT, parameter :: izero = 0
    integer, parameter :: ik = kind(izero)

    CCTK_POINTER,           intent(in) :: cctkGH
    CCTK_INT, dimension(2), intent(in) :: lsh
    CCTK_REAL,              intent(in) :: dt 

    CCTK_REAL,    dimension (lsh(1),lsh(2)),    intent (inout) :: tmp_rgfn, tmp_rgfs 
    CCTK_COMPLEX, dimension (lsh(1),lsh(2), 2), intent (in)    :: Uo, Un
    CCTK_REAL,    dimension (lsh(1),lsh(2), 2), intent (in)    :: omegao,betao, betan
    CCTK_REAL,    dimension (lsh(1),lsh(2), 2), intent (out)   :: omegan
    CCTK_COMPLEX,   dimension (:,:,:), allocatable, save :: eth_ethb_omega

    logical, save :: FirstTime = .true.

    CCTK_COMPLEX, dimension (:,:,:), allocatable, save :: U 
    CCTK_REAL,    dimension (:,:,:), allocatable, save :: omega, omega_u, beta
    DECLARE_CCTK_PARAMETERS

    if (FirstTime) then
       FirstTime = .false.
       allocate( U(lsh(1),lsh(2), 2), omega(lsh(1),lsh(2), 2), omega_u(lsh(1),lsh(2), 2), &
            beta(lsh(1),lsh(2), 2))
       U=0; omega=0; omega_u=0; beta=0
       allocate(eth_ethb_omega(lsh(1),lsh(2), 2))
       eth_ethb_omega=0
    endif

    ! Runge-Kutta for \omega_{,u}

    call NullInterp_d2 (eth_ethb_omega(:,:,1), dcmplx(omegao(:,:,1)), 0_ik, 1_ik,-1_ik)
    call NullInterp_d2 (eth_ethb_omega(:,:,2), dcmplx(omegao(:,:,2)), 0_ik, 1_ik,-1_ik)

    call NullNews_dot_omega (cctkGH, lsh, Uo, omegao, omega_u, betao)
    omegan = omegao + dt * (omega_u - eps_omega * delta(1)*delta(2) * dble(eth_ethb_omega))

    ! we need to sync and interpolate
    call  NullInterp_rinterp(cctkGH, tmp_rgfn, tmp_rgfs, omegan(:,:,1), omegan(:,:,2))
    U = 0.5d0 * (Uo + Un)
    beta = 0.5d0 * (betao+betan)
    omega = 0.5d0 * (omegao + omegan)

    call NullInterp_d2 (eth_ethb_omega(:,:,1), dcmplx(omega(:,:,1)), 0_ik, 1_ik,-1_ik)
    call NullInterp_d2 (eth_ethb_omega(:,:,2), dcmplx(omega(:,:,2)), 0_ik, 1_ik,-1_ik)

    call NullNews_dot_omega (cctkGH, lsh, U, omega, omega_u, beta)

    omegan = omegao + dt * (omega_u - eps_omega * delta(1)*delta(2) * dble(eth_ethb_omega))

    call  NullInterp_rinterp(cctkGH, tmp_rgfn, tmp_rgfs, omegan(:,:,1), omegan(:,:,2))

  end subroutine NullNews_integ_omega

  subroutine NullNews_dot_omega (cctkGH, lsh, U, omega, omega_u, beta)
    use NullInterp
    implicit none

    CCTK_POINTER,               intent(in) :: cctkGH
    CCTK_INT,     dimension(2), intent(in) :: lsh

    CCTK_COMPLEX, dimension (lsh(1),lsh(2),2), intent (in)  :: U
    CCTK_REAL,    dimension (lsh(1),lsh(2),2), intent (in)  :: omega, beta
    CCTK_REAL,    dimension (lsh(1),lsh(2),2), intent (out) :: omega_u

    CCTK_INT, parameter :: izero = 0
    integer, parameter :: ik = kind(izero)

    logical, save :: FirstTime = .true.

    CCTK_COMPLEX,   dimension (:,:,:), allocatable, save :: eth_omega, eth_Ub

    if (FirstTime) then
       FirstTime = .false.
       allocate(eth_omega(lsh(1),lsh(2),2), eth_Ub(lsh(1),lsh(2),2))
       eth_omega=0; eth_Ub=0
    endif

    call NullInterp_d1 (eth_omega(:,:,1), dcmplx(omega(:,:,1)), 0_ik, 1_ik)
    call NullInterp_d1 (eth_omega(:,:,2), dcmplx(omega(:,:,2)), 0_ik, 1_ik)

    call NullInterp_d1 (eth_Ub(:,:,1), conjg(U(:,:,1)), -1_ik, 1_ik)
    call NullInterp_d1 (eth_Ub(:,:,2), conjg(U(:,:,2)), -1_ik, 1_ik)


!!! The term * e^{-2beta} has been added so as to be consistent with HPN. Nigel
    !omega_u = - dble(eth_omega * conjg(U) + 0.5d0 * omega * eth_Ub * exp(-2*beta))

!!! - update .. I thought this may have been wrong so I removed the
!!! e^{-2beta} for testing. Yosef
    omega_u = - dble(eth_omega * conjg(U) + 0.5d0 * omega * eth_Ub)


  end subroutine NullNews_dot_omega

  subroutine NullNews_integ_comega (cctkGH, lsh, dt, tmp_cgfn, tmp_cgfs, omegao, omegan, comegao, comegan)
    use NullInterp
    use NullNews_ScriUtil, only: NullNews_ResetInactive
    implicit none

    CCTK_INT, parameter :: izero = 0
    integer, parameter :: ik = kind(izero)

    CCTK_POINTER,           intent(in) :: cctkGH
    CCTK_INT, dimension(2), intent(in) :: lsh
    CCTK_REAL,              intent(in) :: dt 

    CCTK_COMPLEX, dimension (lsh(1), lsh(2)),    intent (inout)  :: tmp_cgfn, tmp_cgfs 
    CCTK_REAL,    dimension (lsh(1), lsh(2), 2), intent (in)     :: omegao, omegan
    CCTK_COMPLEX, dimension (lsh(1), lsh(2), 2), intent (in)     :: comegao
    CCTK_COMPLEX, dimension (lsh(1), lsh(2), 2), intent (out)    :: comegan

    logical, save :: FirstTime = .true.

    CCTK_COMPLEX, dimension (:,:,:), allocatable, save :: eth_omega_u

    if (FirstTime) then
       FirstTime = .false.
       allocate(eth_omega_u(lsh(1), lsh(2),2))
       eth_omega_u=0
    endif

    ! Mid-point rule for \eth\omega_{,u}

    call NullInterp_d1 (eth_omega_u(:,:,1),&
         dcmplx((omegan(:,:,1) - omegao(:,:,1)) / dt), 0_ik, 1_ik)
    call NullInterp_d1 (eth_omega_u(:,:,2),&
         dcmplx((omegan(:,:,2) - omegao(:,:,2)) / dt), 0_ik, 1_ik)

    call NullInterp_cinterp(cctkGH, tmp_cgfn, tmp_cgfs, eth_omega_u(:,:,1), eth_omega_u(:,:,2), 1_ik)
    call NullNews_ResetInactive(lsh, eth_omega_u)

    comegan = comegao + dt * eth_omega_u

  end subroutine NullNews_integ_comega

  subroutine NullNews_uframe (cctkGH, lsh,&
       tmp_cgfn1, tmp_cgfs1, tmp_cgfn2, tmp_cgfs2, tmp_cgfn3, tmp_cgfs3,&
       J, J_u, J_l, J_l_u, beta, cB, U, omega, comega, News,&
       sigmaJn, sigmaKn, sigmaun, sigmarn, sigmarun, sigmauun,beta_u,U_u,U_l,U_l_l)

    use NullInterp
    use NullNews_ScriUtil, only: NullNews_ResetInactive
    !use cctk
    implicit none

    CCTK_POINTER,               intent(in) :: cctkGH
    CCTK_INT,     dimension(2), intent(in) :: lsh

    CCTK_COMPLEX, dimension (lsh(1), lsh(2)),   intent(inout) :: tmp_cgfn1, tmp_cgfs1
    CCTK_COMPLEX, dimension (lsh(1), lsh(2)),   intent(inout) :: tmp_cgfn2, tmp_cgfs2
    CCTK_COMPLEX, dimension (lsh(1), lsh(2)),   intent(inout) :: tmp_cgfn3, tmp_cgfs3

    CCTK_REAL,    dimension (lsh(1), lsh(2),2), intent(in)    :: beta, omega,beta_u
    CCTK_COMPLEX, dimension (lsh(1), lsh(2),2), intent(in)    :: J, J_u, J_l, J_l_u, &
         cB, U, comega, U_u,U_l,U_l_l

    CCTK_COMPLEX, dimension (lsh(1), lsh(2),2), intent(inout) ::&
         News, sigmaJn, sigmaKn, sigmaun,&
         sigmarn, sigmarun, sigmauun

    logical, save :: FirstTime = .true.

    DECLARE_CCTK_PARAMETERS

    CCTK_INT, parameter :: izero = 0
    integer, parameter :: ik = kind(izero)

    ! temporary arrays

    CCTK_REAL, dimension (:,:,:), allocatable, save ::&
         a, K, K_l, K_u, K_l_u, e2beta
    CCTK_COMPLEX,   dimension (:,:,:), allocatable, save ::&
         Jb, Ub, s1, s2, s3, s4, s5

    CCTK_COMPLEX,   dimension (:,:,:), allocatable, save ::&
         eth_J, ethb_J, eth_J_l, ethb_J_l, &
         eth_K, eth_K_l, eth_beta, eth2_beta, eth_ethb_beta, eth_U, ethb_U, &
         eth_omega, eth2_omega, eth_ethb_omega, &
         eth_a, eth2_a, eth_ethb_a, &
         ethb2_J, eth_ethb_K, ethb2_U, eth_ethb_U, eth_ethb_Ub, ethb_U_u

    if (FirstTime) then
       FirstTime = .false.
       allocate(      a(lsh(1),lsh(2),2),              K(lsh(1),lsh(2),2),&
            K_l(lsh(1),lsh(2),2),            K_u(lsh(1),lsh(2),2),&
            K_l_u(lsh(1),lsh(2),2),             Jb(lsh(1),lsh(2),2),&
            Ub(lsh(1),lsh(2),2),             s1(lsh(1),lsh(2),2),&
            s2(lsh(1),lsh(2),2),             s3(lsh(1),lsh(2),2),&
            s4(lsh(1),lsh(2),2),             s5(lsh(1),lsh(2),2),&
            eth_J(lsh(1),lsh(2),2),         ethb_J(lsh(1),lsh(2),2),&
            eth_J_l(lsh(1),lsh(2),2),       ethb_J_l(lsh(1),lsh(2),2),&
            eth_K(lsh(1),lsh(2),2),        eth_K_l(lsh(1),lsh(2),2),&
            eth_beta(lsh(1),lsh(2),2),      eth2_beta(lsh(1),lsh(2),2),&
            eth_ethb_beta(lsh(1),lsh(2),2),          eth_U(lsh(1),lsh(2),2),&
            ethb_U(lsh(1),lsh(2),2),      eth_omega(lsh(1),lsh(2),2),&
            eth2_omega(lsh(1),lsh(2),2), eth_ethb_omega(lsh(1),lsh(2),2),&
            eth_a(lsh(1),lsh(2),2),         eth2_a(lsh(1),lsh(2),2),&
            ethb2_J(lsh(1),lsh(2),2),         eth_ethb_K(lsh(1),lsh(2),2),&
            ethb2_U(lsh(1),lsh(2),2),         eth_ethb_U(lsh(1),lsh(2),2),&
            eth_ethb_Ub(lsh(1),lsh(2),2),    ethb_U_u(lsh(1),lsh(2),2),&
            eth_ethb_a(lsh(1),lsh(2),2),     e2beta(lsh(1),lsh(2),2)    )
       a=0; K=0; K_l=0; K_u=0; K_l_u=0; Jb=0; Ub=0; s1=0
       s2=0; s3=0; s4=0; s5=0; eth_J=0; ethb_J=0; eth_J_l=0; ethb_J_l=0
       eth_K=0; eth_K_l=0; eth_beta=0; eth2_beta=0
       eth_ethb_beta=0;eth_U=0; ethb_U=0; eth_omega=0
       eth2_omega=0; eth_ethb_omega=0; eth_a=0;eth2_a=0
       eth_ethb_a=0
       ethb2_J=0; eth_ethb_K=0; ethb2_U=0; eth_ethb_U=0
       eth_ethb_Ub =0;ethb_U_u=0;e2beta=0

    endif
    Jb = conjg(J)
    Ub = conjg(U)
    K = sqrt(1.0d0 + J * Jb)

    K_u = dble( J_u * Jb ) / K
    K_l = dble( J_l * Jb ) / K
    K_l_u = dble( J_u * conjg(J_l) + J_l_u * Jb )/ K - K_l * K_u / K 

    call NullInterp_d1 (eth_K(:,:,1),   dcmplx(K(:,:,1)),   0_ik, 1_ik)
    call NullInterp_d1 (eth_K(:,:,2),   dcmplx(K(:,:,2)),   0_ik, 1_ik)

    call NullInterp_d1 (eth_U(:,:,1),  U(:,:,1), 1_ik,  1_ik)
    call NullInterp_d1 (eth_U(:,:,2),  U(:,:,2), 1_ik,  1_ik)

    call NullInterp_d1 (ethb_U(:,:,1), U(:,:,1), 1_ik, -1_ik)
    call NullInterp_d1 (ethb_U(:,:,2), U(:,:,2), 1_ik, -1_ik)

    call NullInterp_d1 (eth_J(:,:,1),  J(:,:,1), 2_ik,  1_ik)
    call NullInterp_d1 (eth_J(:,:,2),  J(:,:,2), 2_ik,  1_ik)

    call NullInterp_d1 (ethb_J(:,:,1), J(:,:,1), 2_ik, -1_ik)
    call NullInterp_d1 (ethb_J(:,:,2), J(:,:,2), 2_ik, -1_ik)

    call NullInterp_d1 (eth_J_l(:,:,1), J_l(:,:,1), 2_ik,  1_ik)
    call NullInterp_d1 (eth_J_l(:,:,2), J_l(:,:,2), 2_ik,  1_ik)

    call NullInterp_d1 (ethb_J_l(:,:,1), J_l(:,:,1),  2_ik, -1_ik)
    call NullInterp_d1 (ethb_J_l(:,:,2), J_l(:,:,2),  2_ik, -1_ik)

    call NullInterp_d1 (eth_K_l(:,:,1), dcmplx(K_l(:,:,1)), 0_ik, 1_ik)
    call NullInterp_d1 (eth_K_l(:,:,2), dcmplx(K_l(:,:,2)), 0_ik, 1_ik)
    if (first_order_scheme.ne.0) then
       eth_omega = comega
       call NullInterp_d1 (eth2_omega(:,:,1), comega(:,:,1), 1_ik,  1_ik)
       call NullInterp_d1 (eth2_omega(:,:,2), comega(:,:,2), 1_ik,  1_ik)
       call NullInterp_d1 (eth_ethb_omega(:,:,1), comega(:,:,1), 1_ik, -1_ik)
       call NullInterp_d1 (eth_ethb_omega(:,:,2), comega(:,:,2), 1_ik, -1_ik)
       eth_beta = cB
       call NullInterp_d1 (eth2_beta(:,:,1), cB(:,:,1), 1_ik,  1_ik)
       call NullInterp_d1 (eth2_beta(:,:,2), cB(:,:,2), 1_ik,  1_ik)
       call NullInterp_d1 (eth_ethb_beta(:,:,1), cB(:,:,1), 1_ik, -1_ik)
       call NullInterp_d1 (eth_ethb_beta(:,:,2), cB(:,:,2), 1_ik, -1_ik)

       call NullInterp_d1 (ethb2_J(:,:,1), ethb_J(:,:,1), 1_ik,  -1_ik) 
       call NullInterp_d1 (ethb2_J(:,:,2), ethb_J(:,:,2), 1_ik,  -1_ik) 
       call NullInterp_d1 (eth_ethb_K(:,:,1), eth_K(:,:,1), 1_ik,  -1_ik) 
       call NullInterp_d1 (eth_ethB_K(:,:,2), eth_K(:,:,2), 1_ik,  -1_ik) 
       call NullInterp_d1 (ethb2_U(:,:,1), ethb_U(:,:,1), 0_ik,  -1_ik) 
       call NullInterp_d1 (ethb2_U(:,:,2), ethb_U(:,:,2), 0_ik,  -1_ik) 
       call NullInterp_d1 (eth_ethb_U(:,:,1), ethb_U(:,:,1), 0_ik,  1_ik) 
       call NullInterp_d1 (eth_ethb_U(:,:,2), ethb_U(:,:,2), 0_ik,  1_ik) 
       call NullInterp_d1 (eth_ethb_Ub(:,:,1), conjg(eth_U(:,:,1)), -2_ik,  1_ik) 
       call NullInterp_d1 (eth_ethb_Ub(:,:,2), conjg(eth_U(:,:,2)), -2_ik,  1_ik) 
    else
       call NullInterp_d1 (eth_omega(:,:,1), dcmplx(omega(:,:,1)), 0_ik,  1_ik)
       call NullInterp_d1 (eth_omega(:,:,2), dcmplx(omega(:,:,2)), 0_ik,  1_ik)
       call NullInterp_d1 (eth_beta(:,:,1), dcmplx(beta(:,:,1)), 0_ik,  1_ik)
       call NullInterp_d1 (eth_beta(:,:,2), dcmplx(beta(:,:,2)), 0_ik,  1_ik)
       call NullInterp_d2 (eth2_omega(:,:,1), dcmplx(omega(:,:,1)), 0_ik, +1_ik, +1_ik)
       call NullInterp_d2 (eth2_omega(:,:,2), dcmplx(omega(:,:,2)), 0_ik, +1_ik, +1_ik)
       call NullInterp_d2 (eth_ethb_omega(:,:,1), dcmplx(omega(:,:,1)), 0_ik, +1_ik, -1_ik)
       call NullInterp_d2 (eth_ethb_omega(:,:,2), dcmplx(omega(:,:,2)), 0_ik, +1_ik, -1_ik)
       call NullInterp_d2 (eth2_beta(:,:,1), dcmplx(beta(:,:,1)), 0_ik, +1_ik,  +1_ik)
       call NullInterp_d2 (eth2_beta(:,:,2), dcmplx(beta(:,:,2)), 0_ik, +1_ik,  +1_ik)
       call NullInterp_d2 (eth_ethb_beta(:,:,1), dcmplx(beta(:,:,1)), 0_ik, +1_ik, -1_ik)
       call NullInterp_d2 (eth_ethb_beta(:,:,2), dcmplx(beta(:,:,2)), 0_ik, +1_ik, -1_ik)
       call NullInterp_d2 (ethb2_J(:,:,1), J(:,:,1), 2_ik, -1_ik,  -1_ik) 
       call NullInterp_d2 (ethb2_J(:,:,2), J(:,:,2), 2_ik, -1_ik,  -1_ik) 
       call NullInterp_d2 (eth_ethb_K(:,:,1), dcmplx(K(:,:,1)), 0_ik, +1_ik,  -1_ik) 
       call NullInterp_d2 (eth_ethB_K(:,:,2), dcmplx(K(:,:,2)), 0_ik, +1_ik,  -1_ik) 
       call NullInterp_d2 (ethb2_U(:,:,1), U(:,:,1), 1_ik, -1_ik,  -1_ik) 
       call NullInterp_d2 (ethb2_U(:,:,2), U(:,:,2), 1_ik, -1_ik,  -1_ik) 
       call NullInterp_d2 (eth_ethb_U(:,:,1), U(:,:,1), 1_ik, +1_ik,  -1_ik) 
       call NullInterp_d2 (eth_ethb_U(:,:,2), U(:,:,2), 1_ik, +1_ik,  -1_ik) 
       call NullInterp_d2 (eth_ethb_Ub(:,:,1), Ub(:,:,1), -1_ik, +1_ik,  -1_ik) 
       call NullInterp_d2 (eth_ethb_Ub(:,:,2), Ub(:,:,2), -1_ik, +1_ik,  -1_ik) 
    end if

    a = omega * exp(2.0d0 * beta)

    ! debug
    !  eth2_omega = (0.0d0, 0.0d0)
    !  eth_ethb_omega = (0.0d0, 0.0d0)
    ! debug

    ! debug
    !  eth2_beta = (0.0d0, 0.0d0)
    !  eth_ethb_beta = (0.0d0, 0.0d0)
    ! debug

    eth_a = exp(2.0d0 * beta) * ( eth_omega + 2.0d0 * omega * eth_beta )

    eth2_a = exp(2.0d0 * beta) * ( 4.0d0 * eth_beta * eth_omega &
         + 4.0d0 * omega * eth_beta**2 &
         + eth2_omega + 2.0d0 * omega * eth2_beta )

    eth_ethb_a = exp(2.0d0 * beta) * ( 4.0d0 * dble(eth_beta * conjg(eth_omega)) &
         + 4.0d0 * omega * eth_beta * conjg(eth_beta) &
         + eth_ethb_omega + 2.0d0 * omega * eth_ethb_beta )

    ! debug
    !  eth2_a = (0.0d0, 0.0d0)
    !  eth_ethb_a = (0.0d0, 0.0d0)
    ! debug

    s1 = ( -2.0d0 * K_l_u * J * (K + 1.0d0) + J_l_u * (K + 1.0d0)**2 &
         + conjg(J_l_u) * J**2 ) / (K + 1.0d0)

    s2 = 0.5d0 / ( K + 1.0d0) * ( &
         (K + 1.0d0)* (eth_J_l *Ub * (K+1.0d0) - 2.0d0* eth_K_l * J *Ub ) &
         + eth_U * (K+1.0d0)* ( -2.0d0 * J * conjg(J_l) + K_l * 2.0d0 * (K+1.0d0) ) &
         + conjg(ethb_U) * (K+1.0d0) * ( -2.0d0* J * K_l + J_l * 2.0d0 * (K+1.0d0) ) &
         + ethb_J_l * U * (K+1.0d0)**2 - conjg(eth_K_l) * 2.0d0 * U * J * (K+1.0d0) &
         + ethb_U * 2.0d0 * J * ( J * conjg(J_l) - (K+1.0d0) * K_l) &
         + J**2 * ( U * conjg(eth_J_l) + conjg(ethb_J_l * U) ) &
         + J * 2.0d0 * conjg(eth_U) * ( J * K_l - J_l * (K+1.0d0) ) )

    s3 = ( J_l * (K + 1.0d0)**2 -2.0d0 * K_l * J * (K + 1.0d0) &
         + conjg(J_l) * J**2) / (K + 1.0d0)

    s4 = 0.5d0 / ( K + 1.0d0) * ( eth_a * eth_omega * (K + 1.0d0)**2 &
         - (K+1.0d0) * J * 2.0d0* dble( eth_a * conjg(eth_omega) ) &
         + J**2 * conjg(eth_a * eth_omega) )

    s5 = 0.25d0 / ( K + 1.0d0) * ( 2.0d0 * eth2_a * (K+1.0d0)**2 &
         + 2.0d0 * J**2 * conjg(eth2_a) &
         - 4.0d0 * eth_ethb_a * J * (K+1.0d0) &
         + Jb * eth_a * eth_J* (K+1.0d0)**2 &
         + J * eth_a * conjg(ethb_J) * (K+1.0d0)**2 &
         - eth_a * eth_K * 2.0d0 * (K+1.0d0) * ( J*Jb + (K+1.0d0) ) & 
         + eth_a * ethb_J * (K+1.0d0) * ( -J*Jb + (K+1.0d0) ) &
         - J**2 * eth_a * conjg(eth_J) * K &
         +  J**2 * Jb * 2.0d0* eth_a * conjg(eth_K) &
         - conjg(eth_a) * eth_J * (K+1.0d0) * ( J*Jb + K+1.0d0 ) &
         - conjg(ethb_J) * conjg(eth_a) * J**2 * ( K + 2.0d0) &
         + J * 2.0d0 * (K+1.0d0)**2 * eth_K * conjg(eth_a)  &
         + J**2 * Jb * ethb_J * conjg(eth_a) &
         + J**3 * conjg(eth_a * eth_J) &
         - 2.0d0* J**2 *K*conjg(eth_K * eth_a) )

    !   News = 0.25d0 * ( s1 + s2 + 0.5d0 * dble(ethb_U) * s3 &
    !   - 4.0d0 * s4 / omega**2 + 2.0d0 * s5 / omega ) / ( omega**2 * exp(2.0d0 * beta) )

    ! change sign of s3 to compensate for a bug in Eqs. 30, 37, and 38 of
    ! HPN
    News = 0.25d0 * ( s1 + s2 - 0.5d0 * dble(ethb_U) * s3 &
         - 4.0d0 * s4 / omega**2 + 2.0d0 * s5 / omega ) / ( omega**2 * exp(2.0d0 * beta) )

    !For testing, the linearized news expression is
    !News =+J_l_u/2 +eth2_beta+eth2_omega/2

    !debug:
    !  s2 = (0.0d0, 0.0d0)
    !  s4 = (0.0d0, 0.0d0)
    !  s5 = (0.0d0, 0.0d0)
    !  News = 0.25d0 * Delta * ( s1 + s2 + 0.5d0 * dble(ethb_U) * s3 &
    !  - 4.0d0 * s4 / omega**2 + 2.0d0 * s5 / omega ) / ( omega**2 * exp(2.0d0 * beta) )
    !debug:

    e2beta = exp(2*beta)
    call NullInterp_d1 (ethb_U_u(:,:,1), U_u(:,:,1), 1_ik,  -1_ik)
    call NullInterp_d1 (ethb_U_u(:,:,2), U_u(:,:,2), 1_ik,  -1_ik)


    s1 = -1/K/16.0
    s4 = 12.0*J*K*K*eth_beta*eth_K*Jb*e2beta-12.0*J*K*&
         eth_beta*conjg(ethb_J)*e2beta+48.0*J*K*K*conjg(eth_beta)*eth_beta*e2beta-8.0&
         *K*J_l_u-4.0*K*U*ethb_J_l-4.0*K*Ub*eth_J_l-8.0*K*K_l*eth_U-8.0*K*J_l*&
         conjg(ethb_U)+8.0*J*K*K*e2beta-16.0*K*eth2_beta*e2beta-32.0*K*eth_beta*&
         eth_beta*e2beta-12.0*conjg(eth2_beta)*J*J*K*e2beta+8.0*K*K*conjg(eth_beta)&
         *eth_J*e2beta+J*conjg(eth_J)*eth_J*e2beta-J*ethb_J*conjg(ethb_J)*e2beta+24.0*&
         J*K*K*eth_ethb_beta*e2beta+2.0*J*conjg(ethb2_J)*K*e2beta+16.0*K*K*&
         eth_beta*eth_K*e2beta-8.0*K*K*eth_beta*ethb_J*e2beta
    s5 = s4-24.0*(conjg(eth_beta)**2)*J*J*K*e2beta&
         -4.0*J*eth_ethb_K*K*e2beta+2.0*J*ethb2_J*K*e2beta-12.0*J*eth2_beta*Jb*K*&
         e2beta-24.0*J*eth_beta*eth_beta*Jb*K*e2beta+12.0*K*K*conjg(eth_beta)*conjg&
         (eth_K)*J*J*e2beta-8.0*K*eth_beta*Jb*eth_J*e2beta-4.0*K*conjg(eth_beta)*J&
         *ethb_J*e2beta+6.0*J*K*K*conjg(eth_beta)*eth_J*Jb*e2beta
    s3 = s5-6.0*eth_beta*Jb*conjg(ethb_J)*J*J*K*&
         e2beta-6.0*conjg(eth_beta)*J*J*ethb_J*Jb*K*e2beta-12.0*conjg(eth_beta)*J*&
         J*eth_K*Jb*K*e2beta-16.0*K*conjg(eth_beta)*J*eth_K*e2beta-6.0*conjg(&
         eth_beta)*J*J*J*conjg(eth_J)*K*e2beta+6.0*K*K*conjg(eth_beta)*J*J*conjg&
         (ethb_J)*e2beta+6.0*K*K*eth_beta*conjg(eth_J)*J*J*e2beta-12.0*eth_beta*Jb*&
         conjg(eth_K)*J*J*K*e2beta-6.0*J*eth_beta*Jb*Jb*eth_J*K*e2beta+6.0*J*K*&
         K*eth_beta*Jb*ethb_J*e2beta
    s4 = 1/e2beta
    s2 = s3*s4
    sigmaJn = s1*s2

    s1 = 1.0/16.0
    s4 = 4.0*K_l*ethb_U+4.0*K_l*conjg(ethb_U)-2.0*ethb2_J*K*e2beta&
         +16.0*eth_ethb_beta*e2beta+4.0*conjg(J_l)*eth_U+8.0*conjg(eth_beta)*J*conjg(&
         ethb_J)*e2beta+4.0*conjg(eth_beta)*K*ethb_J*e2beta+4.0*eth_beta*K*conjg(&
         ethb_J)*e2beta+4.0*J_l*conjg(eth_U)+8.0*K_l_u-conjg(eth_J)*eth_J*e2beta+4.0*U*&
         conjg(eth_K_l)+4.0*Ub*eth_K_l-8.0*K*K*e2beta-2.0*conjg(ethb2_J)*K*e2beta&
         -24.0*K*K*eth_ethb_beta*e2beta+4.0*eth_ethb_K*K*e2beta+12.0*eth2_beta*Jb*K&
         *e2beta
    s3 = s4+24.0*(conjg(eth_beta)**2)*J*K*e2beta+&
         24.0*eth_beta*eth_beta*Jb*K*e2beta+ethb_J*conjg(ethb_J)*e2beta+32.0*conjg(&
         eth_beta)*eth_beta*e2beta+6.0*eth_beta*Jb*conjg(ethb_J)*J*K*e2beta+6.0*conjg&
         (eth_beta)*J*ethb_J*Jb*K*e2beta-6.0*K*K*eth_beta*Jb*ethb_J*e2beta-6.0*K*&
         K*eth_beta*conjg(eth_J)*J*e2beta-12.0*K*K*conjg(eth_beta)*conjg(eth_K)*J*&
         e2beta-6.0*K*K*conjg(eth_beta)*eth_J*Jb*e2beta-12.0*K*K*eth_beta*eth_K*Jb&
         *e2beta-6.0*K*K*conjg(eth_beta)*J*conjg(ethb_J)*e2beta+6.0*eth_beta*Jb*Jb*&
         eth_J*K*e2beta+6.0*conjg(eth_beta)*J*J*conjg(eth_J)*K*e2beta+8.0*eth_beta*&
         Jb*ethb_J*e2beta-48.0*K*K*conjg(eth_beta)*eth_beta*e2beta+12.0*conjg(&
         eth_beta)*J*eth_K*Jb*K*e2beta+12.0*eth_beta*Jb*conjg(eth_K)*J*K*e2beta+&
         12.0*conjg(eth2_beta)*J*K*e2beta
    s4 = 1/e2beta
    s2 = s3*s4
    sigmaKn = s1*s2

    sigmaun = -eth_beta*conjg(ethb_U)/2.0+conjg(ethb2_U)/4.0-eth_beta*ethb_U/2.0+&
         eth_ethb_U/4.0-Ub*sigmaJn-U*sigmaKn

    sigmarn = (K_l*U_l+J_l*conjg(U_l)+K*U_l_l+J*conjg(U_l_l))/e2beta/2.0

    s1 = 1.0/16.0
    s3 = 1/K
    s5 = 2.0*ethb2_J*K*e2beta-4.0*conjg(eth_beta)*K*ethb_J*e2beta&
         -4.0*eth_beta*K*conjg(ethb_J)*e2beta-8.0*Ub*sigmarn*K+conjg(eth_J)*eth_J*&
         e2beta-8.0*U*conjg(sigmarn)*K+8.0*K*K*e2beta+2.0*conjg(ethb2_J)*K*e2beta+8.0*K*K&
         *eth_ethb_beta*e2beta-4.0*eth_ethb_K*K*e2beta-4.0*eth2_beta*Jb*K*e2beta-8.0*&
         (conjg(eth_beta)**2)*J*K*e2beta-8.0*eth_beta*eth_beta*Jb*K*e2beta-&
         ethb_J*conjg(ethb_J)*e2beta
    s4 = s5-2.0*eth_beta*Jb*conjg(ethb_J)*J*K*e2beta&
         -2.0*conjg(eth_beta)*J*ethb_J*Jb*K*e2beta+2.0*K*K*eth_beta*Jb*ethb_J*&
         e2beta+2.0*K*K*eth_beta*conjg(eth_J)*J*e2beta+4.0*K*K*conjg(eth_beta)*&
         conjg(eth_K)*J*e2beta+2.0*K*K*conjg(eth_beta)*eth_J*Jb*e2beta+4.0*K*K*&
         eth_beta*eth_K*Jb*e2beta+2.0*K*K*conjg(eth_beta)*J*conjg(ethb_J)*e2beta-2.0&
         *eth_beta*Jb*Jb*eth_J*K*e2beta-2.0*conjg(eth_beta)*J*J*conjg(eth_J)*K*&
         e2beta+16.0*K*K*conjg(eth_beta)*eth_beta*e2beta-4.0*conjg(eth_beta)*J*eth_K*&
         Jb*K*e2beta-4.0*eth_beta*Jb*conjg(eth_K)*J*K*e2beta-4.0*conjg(eth2_beta)*&
         J*K*e2beta
    s2 = s3*s4
    sigmarun = s1*s2

    sigmauun = U*conjg(eth_beta)*conjg(ethb_U)/4.0+U*Ub*sigmaKn+conjg(ethb_U_u)/&
         4.0+ethb_U_u/4.0+Ub*eth_beta*ethb_U/4.0+U*conjg(eth_beta)*ethb_U/4.0+Ub*&
         eth_beta*conjg(ethb_U)/4.0-Ub*eth_ethb_U/8.0-beta_u*conjg(ethb_U)/2.0-Ub*&
         conjg(ethb2_U)/8.0-U*eth_ethb_Ub/8.0-beta_u*ethb_U/2.0-U*ethb2_U/8.0+U*conjg&
         (ethb_U)/4.0+Ub*Ub*sigmaJn/2.0+U*U*conjg(sigmaJn)/2.0

    ! your mileage may vary...

    call NullInterp_cinterp(cctkGH, tmp_cgfn1, tmp_cgfs1, news(:,:,1), news(:,:,2), 2_ik)

    call NullInterp_3cinterp(cctkGH,&
               tmp_cgfn1, tmp_cgfs1,&
               tmp_cgfn2, tmp_cgfs2,&
               tmp_cgfn3, tmp_cgfs3,&
               sigmaJn(:,:,1), sigmaJn(:,:,2),&
               sigmaKn(:,:,1), sigmaKn(:,:,2),&
               sigmaun(:,:,1), sigmaun(:,:,2),&
               2_ik, 0_ik, 1_ik)

!   call NullInterp_cinterp(cctkGH, tmp_cgfn, tmp_cgfs, sigmaJn(:,:,1), sigmaJn(:,:,2), 2)
!   call NullInterp_cinterp(cctkGH, tmp_cgfn, tmp_cgfs, sigmaKn(:,:,1), sigmaKn(:,:,2), 0)
!   call NullInterp_cinterp(cctkGH, tmp_cgfn, tmp_cgfs, sigmaun(:,:,1), sigmaun(:,:,2), 1)

    call NullInterp_3cinterp(cctkGH,&
               tmp_cgfn1, tmp_cgfs1,&
               tmp_cgfn2, tmp_cgfs2,&
               tmp_cgfn3, tmp_cgfs3,&
               sigmarn(:,:,1), sigmarn(:,:,2),&
               sigmarun(:,:,1), sigmarun(:,:,2),&
               sigmauun(:,:,1), sigmauun(:,:,2),&
               1_ik, 0_ik, 0_ik)

!   call NullInterp_cinterp(cctkGH, tmp_cgfn, tmp_cgfs, sigmarn(:,:,1), sigmarn(:,:,2), 1)
!   call NullInterp_cinterp(cctkGH, tmp_cgfn, tmp_cgfs, sigmarun(:,:,1), sigmarun(:,:,2), 0)
!   call NullInterp_cinterp(cctkGH, tmp_cgfn, tmp_cgfs, sigmauun(:,:,1), sigmauun(:,:,2), 0)

    call NullNews_ResetInactive(lsh, news)
    call NullNews_ResetInactive(lsh, sigmaJn)
    call NullNews_ResetInactive(lsh, sigmaKn)
    call NullNews_ResetInactive(lsh, sigmaun)
    call NullNews_ResetInactive(lsh, sigmarn)
    call NullNews_ResetInactive(lsh, sigmarun)
    call NullNews_ResetInactive(lsh, sigmauun)

  end subroutine NullNews_uframe

  
  subroutine NullNews_lin_omega (cctkGH, lsh, zeta, J, omega, comega)
    use NullDecomp_SpinDecomp, only: SpinDecompCoefs, SpinDecompRecon
    use NullDecomp_Vars, only: lmax
    implicit none

    CCTK_INT, parameter :: izero = 0
    integer, parameter :: ik = kind(izero)

    CCTK_POINTER,           intent(in) :: cctkGH
    CCTK_INT, dimension(2), intent(in) :: lsh
    CCTK_REAL,    dimension (lsh(1),lsh(2), 2), intent (out)   :: omega
    CCTK_COMPLEX, dimension (lsh(1),lsh(2), 2), intent (out)   :: comega
    CCTK_COMPLEX, dimension (lsh(1),lsh(2), 2), intent (in)    :: J
    CCTK_COMPLEX, dimension (lsh(1),lsh(2), 2), intent (in)    :: zeta
    CCTK_COMPLEX, dimension(2:lmax, -lmax:lmax)                :: JnSphCoeff
    CCTK_COMPLEX, dimension(0:lmax, -lmax:lmax)                :: omegaSphCoeff
    CCTK_COMPLEX, dimension(1:lmax, -lmax:lmax)                :: comegaSphCoeff
    CCTK_INT  :: l, m
    CCTK_REAL :: C
    
    ! get harm. coefs of J at scri
    call SpinDecompCoefs(cctkGH, lsh(1), lsh(2), 2_ik, &
                         zeta, J, JnSphCoeff)

    ! get linearzied harm. coefs of omega at scri
    
    omegaSphCoeff = 0.0
    do l = 2,lmax
      do m = -l,l
        omegaSphCoeff(l, m) = -0.25 * sqrt( l*(l+1.0)/((l-1.0)*(l+2.0)) ) * ( JnSphCoeff(l, m) + (-1.0)**m * conjg(JnSphCoeff(l, -m)) )
      end do
    end do

    ! reconstruct omega
    call SpinDecompRecon(cctkGH, lsh(1), lsh(2), 0_ik, zeta, comega, omegaSphCoeff)
    omega = dble(comega) + 1.0

    ! get linearized harm. coefs of comega at scri
    comegaSphCoeff = 0.0
    do l = 2,lmax
      C = l*(l+1.0)   !sqrt(fact(l+1)/fact(l-1))
      do m = -l,l
        comegaSphCoeff(l, m) = C * omegaSphCoeff(l, m)
      end do
    end do

    ! reconstruct comega
    call SpinDecompRecon(cctkGH, lsh(1), lsh(2), 1_ik, zeta, comega, comegaSphCoeff)

  end subroutine NullNews_lin_omega


end module NullNews_Omega
