! vim: syntax=fortran
#include "cctk.h"

#define DMIN  .5
#define DN  2
#define DISMAX  .05

module NullEvol_Evol
    use NullInterp
    use NullGrid_Vars 
    implicit none
    private
    public NullEvol_j

contains

  subroutine NullEvol_j_rhs (i, jns, nuns, ckns, bns, cbns, qns, uns, wns, &
       jos, nuos, ckos, bos, cbos, qos, uos, wos,&
       J, dx_J, dx_Jo, W, dx_beta, h_c, fac, f1, f2, xc, dt,&
       x_wt, j_wt, beta_wt, q_wt, u_wt, w_wt,&
       x_wt_o, j_wt_o, beta_wt_o, q_wt_o, u_wt_o, w_wt_o,&
       first_order_scheme, stencil_type)

    CCTK_REAL,                                    intent (in)    :: dt
    CCTK_INT,                                     intent (in)    :: i, first_order_scheme, stencil_type
    CCTK_COMPLEX, dimension (lsh(1), lsh(2), nx), intent (in)    :: jns
    CCTK_COMPLEX, dimension (lsh(1), lsh(2), nx), intent (in)    :: jos, nuos, ckos, cbos, uos, qos
    CCTK_COMPLEX, dimension (lsh(1), lsh(2), nx), intent (in)    :: nuns, ckns, cbns, uns, qns
    CCTK_REAL,    dimension (lsh(1), lsh(2), nx), intent (in)    :: bns, bos, wns, wos
    CCTK_COMPLEX, dimension (lsh(1), lsh(2)),     intent (out)   :: J, dx_J, dx_Jo, h_c, fac 
    CCTK_REAL,    dimension (lsh(1), lsh(2)),     intent (out)   :: W, dx_beta, f1, f2, xc
    CCTK_COMPLEX, dimension (lsh(1), lsh(2)),     intent (in)    :: j_wt, q_wt, u_wt
    CCTK_COMPLEX, dimension (lsh(1), lsh(2)),     intent (in)    :: j_wt_o, q_wt_o, u_wt_o
    CCTK_REAL,    dimension (lsh(1), lsh(2)),     intent (in)    :: x_wt, beta_wt, w_wt
    CCTK_REAL,    dimension (lsh(1), lsh(2)),     intent (in)    :: x_wt_o, beta_wt_o, w_wt_o

    CCTK_INT, parameter :: izero = 0
    integer, parameter :: ik = kind(izero)

    logical, save :: FirstTime = .true.

    !! probably better to make these into allocated arrays
    CCTK_REAL, dimension (:,:), allocatable, save :: &
         beta, K, dx_K, dx_W, rc, xold, c1, c2, c3

    CCTK_COMPLEX, dimension (:,:), allocatable, save :: &
         U, cB, cK, d2x_jave, dx_U, dx_jave, dx_mu,&
         dx_nu, eth_J, eth_U,  ethb_U, eth_cB, eth_dx_J,&
         eth_dx_U, ethb_dx_U, ethb_cB,jave,&
         mu, nu, J_H, Q, A

    CCTK_COMPLEX :: ii = (0. , 1.)

    CCTK_INT B
    CCTK_INT, parameter :: cran_max = 2

    if (FirstTime) then
       FirstTime = .false.
       allocate(beta(lsh(1),lsh(2)), &
            K(lsh(1),lsh(2)), dx_K(lsh(1),lsh(2)), dx_W(lsh(1),lsh(2)))
       allocate(U(lsh(1),lsh(2)), cB(lsh(1),lsh(2)), cK(lsh(1),lsh(2)),&
            d2x_jave(lsh(1),lsh(2)), dx_U(lsh(1),lsh(2)), &
            dx_jave(lsh(1),lsh(2)), dx_mu(lsh(1),lsh(2)),&
            dx_nu(lsh(1),lsh(2)), eth_J(lsh(1),lsh(2)), eth_U(lsh(1),lsh(2)),&
            ethb_U(lsh(1),lsh(2)), eth_cB(lsh(1),lsh(2)), eth_dx_J(lsh(1),lsh(2)),&
            eth_dx_U(lsh(1),lsh(2)), ethb_dx_U(lsh(1),lsh(2)), ethb_cB(lsh(1),lsh(2)),&
            jave(lsh(1),lsh(2)),&
            mu(lsh(1),lsh(2)), nu(lsh(1),lsh(2)), J_H(lsh(1),lsh(2)),Q(lsh(1),lsh(2)))
       allocate(A(lsh(1),lsh(2)),rc(lsh(1),lsh(2)),xold(lsh(1),lsh(2)),&
                c1(lsh(1),lsh(2)), c2(lsh(1),lsh(2)), c3(lsh(1),lsh(2)))

       beta=0; dx_beta=0; W=0; K=0; dx_K=0; dx_W=0;
       U=0; cB=0; cK=0; d2x_jave=0; dx_U=0; dx_jave=0; dx_mu=0
       dx_nu=0; eth_J=0; eth_U=0; ethb_U=0; eth_cB=0; eth_dx_J=0; eth_dx_U=0;
       ethb_dx_U=0; ethb_cB=0; jave=0; mu=0; nu=0; J_H=0; Q=0;
    end if


    if(stencil_type.eq.0) then ! evolution stencil

       xc = xbh(i-1)

       U    = 0.25 * (uns(:,:,i-1) + uos(:,:,i-1) +&
            uns(:,:,i-2) + uos(:,:,i))

       J    = 0.5 * ( jns(:,:,i-1) +  jos(:,:,i))
       beta = 0.5 * ( bns(:,:,i-1) +  bos(:,:,i))
       W    = 0.5 * ( wns(:,:,i-1) +  wos(:,:,i))
       Q    = 0.5 * ( qns(:,:,i-1) +  qos(:,:,i))

    else

       ! stencil_type=1 => B+1=i
       ! stencil_type=2 => B+2=i
       B=i-stencil_type

       A = 1+rwt*(x_wt)/(1-x_wt)*w_wt ! 1st order accurate
       if(stencil_type.eq.1) then
          rc = 0.5 * ( 0.50*A*dt + rb(B+1) + rwt*(x_wt)/(1-x_wt) )
       else
          rc = 0.5 * ( 0.25*A*dt + rb(B+2) + rb(B+1) )
       end if
       xc = rc / (rwt+rc)
       xold = x_wt + 2*(xc-x_wt)
       ! we interpolate on the old level using wt, B+2, B+3
       ! to compute W at xc
       W  = 0.5 * ( w_wt &
          + w_wt_o*(xold-xb(B+2))*(xold-xb(B+3))/(x_wt_o-xb(B+2))/(x_wt_o-xb(B+3)) &
          + wos(:,:,B+2)*(xold-x_wt_o)*(xold-xb(B+3))/(xb(B+2)-x_wt_o)/(xb(B+2)-xb(B+3)) &
          + wos(:,:,B+3)*(xold-x_wt_o)*(xold-xb(B+2))/(xb(B+3)-x_wt_o)/(xb(B+3)-xb(B+2)) )

       A = 1+rwt*(x_wt)/(1-x_wt)*W ! improved accuracy
       if(stencil_type.eq.1) then
          rc = 0.5 * ( 0.50*A*dt + rb(B+1) + rwt*(x_wt)/(1-x_wt) )
       else
          rc = 0.5*( rb(B+2) + rb(B+1) + 0.25 * A * dt )
       end if
       xc = rc / (rwt+rc)
       xold = x_wt + 2*(xc-x_wt)
       ! we interpolate J,beta,Q,U, and re-interpolate W to this more accurate target point.
       ! Lagrance interpolation coefficients for the half-grid xbh
       c1 = (xold-xbh(B+2))*(xold-xbh(B+3))/(x_wt_o-xbh(B+2))/(x_wt_o-xbh(B+3))
       c2 = (xold-x_wt_o)*(xold-xbh(B+3))/(xbh(B+2)-x_wt_o)/(xbh(B+2)-xbh(B+3))
       c3 = (xold-x_wt_o)*(xold-xbh(B+2))/(xbh(B+3)-x_wt_o)/(xbh(B+3)-xbh(B+2))
      
       U  = 0.5 * ( u_wt + u_wt_o*c1 + uos(:,:,B+2)*c2 + uos(:,:,B+3)*c3 )

       ! Lagrance interpolation coefficients for the xb grid
       c1 = (xold-xb(B+2))*(xold-xb(B+3))/(x_wt_o-xb(B+2))/(x_wt_o-xb(B+3))
       c2 = (xold-x_wt_o)*(xold-xb(B+3))/(xb(B+2)-x_wt_o)/(xb(B+2)-xb(B+3))
       c3 = (xold-x_wt_o)*(xold-xb(B+2))/(xb(B+3)-x_wt_o)/(xb(B+3)-xb(B+2))
      
       J     = 0.5 * ( j_wt    + j_wt_o   *c1 + jos(:,:,B+2)*c2 + jos(:,:,B+3)*c3 )
       beta  = 0.5 * ( beta_wt + beta_wt_o*c1 + bos(:,:,B+2)*c2 + bos(:,:,B+3)*c3 )
       W     = 0.5 * ( w_wt    + w_wt_o   *c1 + wos(:,:,B+2)*c2 + wos(:,:,B+3)*c3 )
       Q     = 0.5 * ( q_wt    + q_wt_o   *c1 + qos(:,:,B+2)*c2 + qos(:,:,B+3)*c3 )

    end if

    K    = sqrt (1.d0 + J * conjg(J))
    dx_U = exp(2. * beta) / rwt / xc**2 &
                 * ( - J * conjg(Q) + K * Q )

    if (first_order_scheme.ne.0 .and. (stencil_type .ne. 1)) then ! all but the B+1 points

       cB   = 0.5 * (cbns(:,:,i-1) + cbos(:,:,i))
       cK   = 0.5 * (ckns(:,:,i-1) + ckos(:,:,i))
       nu   = 0.5 * (nuns(:,:,i-1) + nuos(:,:,i))

    else

       ! nu = ethb_J
       ! ck = eth_K
       ! cb = eth_beta

       call NullInterp_d1(nu, J,            2_ik, -1_ik)
       call NullInterp_d1(ck, dcmplx(K),    0_ik, +1_ik)
       call NullInterp_d1(cb, dcmplx(beta), 0_ik, +1_ik)

    end if

    if (i .ne. nx) then

       if(stencil_type.eq.0) then ! evolution stencil
          dx_J = 0.5 * (jns(:,:,i-1) - jns(:,:,i-2) &
               + jos(:,:,i+1) - jos(:,:,i)) / dx

          dx_beta = 0.5 * (bns(:,:,i-1) - bns(:,:,i-2) &
               + bos(:,:,i+1) - bos(:,:,i)) / dx

          dx_W = 0.5 * (wns(:,:,i-1) - wns(:,:,i-2) &
               + wos(:,:,i+1) - wos(:,:,i)) / dx

!Computed radial derivative of J
          dx_Jo = 0.5 * (jos(:,:,i+1) - jos(:,:,i-1)) / dx
!         dx_Jo = 1/12.d0 * (jos(:,:,i-2) - 8.d0 * jos(:,:,i-1)&
!                   + 8.d0 * jos(:,:,i+1) - jos(:,:,i+2)) / dx
       else ! all boundary stencils

          ! these are 1st order accurate in time
          ! dx_J    = (jos(:,:,i+1) -    j_wt_o) / (xb(i+1)-x_wt_o)
          ! dx_beta = (bos(:,:,i+1) - beta_wt_o) / (xb(i+1)-x_wt_o)
          ! dx_W    = (wos(:,:,i+1) -    w_wt_o) / (xb(i+1)-x_wt_o)
          c1 = ((xc-xb(i+2))/(x_wt_o-xb(i+1))/(x_wt_o-xb(i+2))+(xc-xb(i+1))/(x_wt_o-xb(i+1))/(x_wt_o-xb(i+2)))
          c2 = ((xc-xb(i+2))/(xb(i+1)-x_wt_o)/(xb(i+1)-xb(i+2))+(xc-x_wt_o)/(xb(i+1)-x_wt_o)/(xb(i+1)-xb(i+2)))
          c3 = ((xc-xb(i+1))/(xb(i+2)-x_wt_o)/(xb(i+2)-xb(i+1))+(xc-x_wt_o)/(xb(i+2)-x_wt_o)/(xb(i+2)-xb(i+1)))

          dx_J    = c1 * j_wt_o    + c2 * jos(:,:,i+1) + c3 * jos(:,:,i+2)
          dx_beta = c1 * beta_wt_o + c2 * bos(:,:,i+1) + c3 * bos(:,:,i+2)
          dx_W    = c1 * w_wt_o    + c2 * wos(:,:,i+1) + c3 * wos(:,:,i+2)
          dx_Jo   = dx_J
       end if

       if (first_order_scheme.ne.0 .and. (stencil_type.eq.0)) then

          dx_nu = 0.5 * (nuns(:,:,i-1) - nuns(:,:,i-2) &
               + nuos(:,:,i+1) - nuos(:,:,i)) / dx

       else

          ! nu = ethb_J
          call NullInterp_d1(dx_nu, dx_J, 2_ik, -1_ik)

       end if

    else ! i .eq. nx -- this is the algorithm at scri+

       dx_J = 0.5 * (&
            + 0.5 * (3. * jns(:,:,i-1) - 4. * jns(:,:,i-2) + jns(:,:,i-3))&
            + 0.5 * (3. * jos(:,:,i)   - 4. * jos(:,:,i-1) + jos(:,:,i-2))&
            ) / dx

       dx_beta = 0.5 * (&
            + 0.5 * (3. * bns(:,:,i-1) - 4. * bns(:,:,i-2) + bns(:,:,i-3))&
            + 0.5 * (3. * bos(:,:,i)   - 4. * bos(:,:,i-1) + bos(:,:,i-2))&
            ) / dx

       dx_W = 0.5 * (&
            + 0.5 * (3. * wns(:,:,i-1) - 4. * wns(:,:,i-2) + wns(:,:,i-3))&
            + 0.5 * (3. * wos(:,:,i)   - 4. * wos(:,:,i-1) + wos(:,:,i-2))&
            ) / dx

       dx_Jo = 0.5 * (3. * jos(:,:,i)- 4. *jos(:,:,i-1) + jos(:,:,i-2)) / dx

       if (first_order_scheme.ne.0) then

          dx_nu = 0.5 * (&
               + 0.5 * (3. * nuns(:,:,i-1) - 4. * nuns(:,:,i-2) + nuns(:,:,i-3))&
               + 0.5 * (3. * nuos(:,:,i)   - 4. * nuos(:,:,i-1) + nuos(:,:,i-2))&
               ) / dx

       else

          ! nu = ethb_J
          call NullInterp_d1(dx_nu, dx_J, 2_ik, -1_ik)

       end if

    end if

    call NullInterp_d1(eth_J,  J, 2_ik,  1_ik)
    call NullInterp_d1(eth_U,  U, 1_ik,  1_ik)
    call NullInterp_d1(ethb_U, U, 1_ik, -1_ik)
    call NullInterp_d1(eth_dx_J,  dx_J, 2_ik,  1_ik)
    call NullInterp_d1(eth_dx_U,  dx_U, 1_ik,  1_ik)
    call NullInterp_d1(ethb_dx_U, dx_U, 1_ik, -1_ik)
    call NullInterp_d1(mu, J, 2_ik, 1_ik)
    call NullInterp_d1(dx_mu, dx_J, 2_ik, 1_ik)

    if (first_order_scheme.eq.0) then
       call NullInterp_d2(ethb_cb, dcmplx(beta), 0_ik, -1_ik, +1_ik) ! ethb_eth
       call NullInterp_d2(eth_cb,  dcmplx(beta), 0_ik, +1_ik, +1_ik) ! eth_eth
    else
       call NullInterp_d1(eth_cB,  cB, 1_ik,  1_ik)
       call NullInterp_d1(ethb_cB, cB, 1_ik, -1_ik)
    end if

    dx_K = dble( dx_J * conjg(J) ) / K 
    fac = ( - dx_K * J + K * dx_J )  / K

    ! non-linear terms

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    J_H = (1. - xc) / (rwt * xc) * exp(2. * beta) *&
         ( -K * eth_J * conjg(cB) &
         + (K * nu + conjg(J) * eth_J - 2. * K * ck) * cB&
         + J * ((2. * ck - nu) * conjg(cB) - 2. * K *&
         (ethb_cB + cB * conjg(cB)) &
         + 2. * dble((nu - ck) * conjg(cB) + conjg(J) *&
         (eth_cB + cB * cB)))) &
         + 0.5 * rwt * (1. - xc) * xc**3 * exp(-2. * beta) &
         * ( (K * dx_U + J * conjg(dx_U))**2 &
         - J * dble(conjg(dx_U) * (K * dx_U + J * conjg(dx_U)))) &
         - 0.5 * (  nu * (xc * (1. - xc) * dx_U + 2. * U) &
         + eth_J * conjg(xc * (1. - xc) * dx_U + 2. * U) ) &
         + J * ii * dimag(xc * (1. - xc) * ethb_dx_U + 2. * ethb_U) &
         - xc * (1. - xc) * dx_J * dble(ethb_U) &
         + xc * (1. - xc) * (conjg(U) * eth_J + U * nu) * ii *&
         dimag(J * conjg(dx_J)) &
         - xc * (1. - xc) * (conjg(U) * eth_dx_J + U * dx_nu) &
         - 2. * xc * (1. - xc) * (J * dx_K - K * dx_J) &
         * ( dble(conjg(U) * ck) + ii *&
         dimag(K * ethb_U - conjg(J) * eth_U) )

    ! the rhs ...

    h_c = - K * (xc * (1. - xc) * eth_dx_U + 2. * eth_U )&
         + 2. * (1. - xc) / (rwt * xc) * exp (2. * beta) *&
         (eth_cB + cB * cB)&
         - (xc * (1. - xc) * dx_W + W) * J + J_H

    !these will be used in the main evolution

    f1 = dt * (xc * (1. - xc) * dx_W + W)
    f2 = dt * (1. - xc)**2 * ((1. - xc)/rwt + xc * W)

  end subroutine NullEvol_j_rhs

  subroutine NullEvol_P_u(i, cran, dt, jns, jos, xc, fac, dx_beta, J, W, P_u,&
                          x_wt, x_wt_o, j_wt, j_wt_o, stencil_type)
    CCTK_REAL,                                    intent(in)   :: dt
    CCTK_INT,                                     intent (in)  :: i, cran, stencil_type
    CCTK_COMPLEX, dimension (lsh(1), lsh(2), nx), intent (in)  :: jns, jos
    CCTK_COMPLEX, dimension (lsh(1), lsh(2)),     intent (in)  :: J, fac, j_wt, j_wt_o
    CCTK_REAL,    dimension (lsh(1), lsh(2)),     intent (in)  :: xc, dx_beta, W, x_wt, x_wt_o
    CCTK_COMPLEX, dimension (lsh(1), lsh(2)),     intent (out) :: P_u

    if(stencil_type.ne.1) then ! all but the B+1 points
       if(cran.eq.1) then
          P_u = J * (1. - xc) * (xc * 2. * dble(( jns(:,:,i-1) &
               - jos(:,:,i-1)) &
               / (dt) * conjg(fac)) &
               - 8. * dx_beta * ((1. - xc)/rwt + xc * W))
       else
          P_u = J * (1. - xc) * (xc * 2. *&
               dble((jns(:,:,i) + jns(:,:,i-1)&
               - jos(:,:,i) - jos(:,:,i-1)) &
               / (2. * dt) * conjg(fac)) &
               - 8. * dx_beta * ((1. - xc)/rwt + xc * W))
       end if
    else
       ! we interpolate J on the old level onto x_wt(new)  and then
       ! take a time-derivative along constant x to get a 1st order
       ! accurate expression for P_u

       P_u = J * (1. - xc) * (xc * 2. * dble(( j_wt&
               - j_wt_o*(x_wt-xb(i+1))*(x_wt-xb(i+2))/(x_wt_o-xb(i+1))/(x_wt_o-xb(i+2)) &
               - jos(:,:,i+1)*(x_wt-x_wt_o)*(x_wt-xb(i+2))/(xb(i+1)-x_wt_o)/(xb(i+1)-xb(i+2)) &
               - jos(:,:,i+2)*(x_wt-x_wt_o)*(x_wt-xb(i+1))/(xb(i+2)-x_wt_o)/(xb(i+2)-xb(i+1)) ) &
           / (dt) * conjg(fac)) &
             - 8. * dx_beta * ((1. - xc)/rwt + xc * W))

    end if

  end subroutine NullEvol_P_u

  subroutine NullEvol_j_Bp1(i, j_Bp1, jos, P_u, h_c, xc, W,&
                            j_wt, j_wt_o,  x_wt, x_wt_o, dt)

    CCTK_INT,                                     intent (in)    :: i
    CCTK_REAL,    dimension (lsh(1), lsh(2)),     intent (in)    :: xc, x_wt, x_wt_o, W
    CCTK_COMPLEX, dimension (lsh(1), lsh(2)),     intent (in)    :: j_wt, j_wt_o, P_u, h_c
    CCTK_COMPLEX, dimension (lsh(1), lsh(2)),     intent (out)   :: j_Bp1
    CCTK_COMPLEX, dimension (lsh(1), lsh(2), nx), intent (in)    :: jos
    CCTK_REAL,                                    intent(in)     :: dt

    CCTK_REAL     :: r_0o, r_1o, A, x_0o, x_1o, x_1n, r_1n, r_wt, r_1h, x_1h, eta
    CCTK_COMPLEX  :: j_0o, j_1o, j_1n

    integer i1, i2, B

    B = i-1

    do i2 = 1, lsh(2)
       do i1 = 1, lsh(1)

          A = 1 + W(i1,i2) * rwt * xc(i1,i2) /( 1-xc(i1,i2))

          r_wt = rwt * x_wt(i1,i2) / (1-x_wt(i1,i2))
          r_0o = 0.5 * A * dt + r_wt

          x_1h = 0.5 * (x_wt(i1,i2) + xb(B+1))
          r_1h = rwt * x_1h/(1-x_1h)

          r_1n = -0.25 * A * dt + rb(B+1)
          x_1n = r_1n/(rwt+r_1n)
          if(r_1h .gt. r_1n) then
            r_1n = r_1h
            x_1n = x_1h
            eta = 1
          else
            eta = (xb(B+1)-x_1n)/(x_1n-x_wt(i1,i2))
          end if
          ! r_1o = +0.25 * A * dt + rb(B+1)
          r_1o = r_1n + 0.5 * A * dt

          x_0o = r_0o/(rwt+r_0o)
          x_1o = r_1o/(rwt+r_1o)

          j_0o = j_wt_o(i1,i2)*(x_0o-xb(B+2))*(x_0o-xb(B+3))/(x_wt_o(i1,i2)-xb(B+2))/(x_wt_o(i1,i2)-xb(B+3)) &
              + jos(i1,i2,B+2)*(x_0o-x_wt_o(i1,i2))*(x_0o-xb(B+3))/(xb(B+2)-x_wt_o(i1,i2))/(xb(B+2)-xb(B+3)) &
              + jos(i1,i2,B+3)*(x_0o-x_wt_o(i1,i2))*(x_0o-xb(B+2))/(xb(B+3)-x_wt_o(i1,i2))/(xb(B+3)-xb(B+2))

          j_1o = j_wt_o(i1,i2)*(x_1o-xb(B+2))*(x_1o-xb(B+3))/(x_wt_o(i1,i2)-xb(B+2))/(x_wt_o(i1,i2)-xb(B+3)) &
              + jos(i1,i2,B+2)*(x_1o-x_wt_o(i1,i2))*(x_1o-xb(B+3))/(xb(B+2)-x_wt_o(i1,i2))/(xb(B+2)-xb(B+3)) &
              + jos(i1,i2,B+3)*(x_1o-x_wt_o(i1,i2))*(x_1o-xb(B+2))/(xb(B+3)-x_wt_o(i1,i2))/(xb(B+3)-xb(B+2))

          j_1n = ( r_wt * j_wt(i1,i2) - r_0o * j_0o + r_1o * j_1o&
                  + 0.5 * dt * (r_1n-r_wt) * (P_u(i1,i2)+h_c(i1,i2)) ) / r_1n

          j_Bp1(i1,i2) = j_1n + (j_1n-j_wt(i1,i2)) * eta

       end do
    end do

  end subroutine NullEvol_j_Bp1


  subroutine NullEvol_j_Bp2(i, j_Bp2, jns, jos, P_u, h_c, xc, W,&
                            j_wt, j_wt_o,  x_wt, x_wt_o, dt)

    CCTK_INT,                                     intent (in)    :: i
    CCTK_REAL,    dimension (lsh(1), lsh(2)),     intent (in)    :: xc, x_wt, x_wt_o, W
    CCTK_COMPLEX, dimension (lsh(1), lsh(2)),     intent (in)    :: j_wt, j_wt_o, P_u, h_c
    CCTK_COMPLEX, dimension (lsh(1), lsh(2)),     intent (out)   :: j_Bp2
    CCTK_COMPLEX, dimension (lsh(1), lsh(2), nx), intent (in)    :: jns, jos
    CCTK_REAL,                                    intent(in)     :: dt

    CCTK_REAL,    dimension (lsh(1), lsh(2)) :: r_1o, r_1n, r_2o, r_2n, A,&
              x_1o, x_1n, x_2o, x_2n, r_wt
    CCTK_COMPLEX, dimension (lsh(1), lsh(2)) :: j_1o, j_1n,j_2o, j_2n

    integer B

    B = i-2

    A = 1 + W * rwt * xc /( 1-xc)

    r_wt = rwt * x_wt / (1-x_wt)
    r_1o = rb(B+1) + 0.25*A*dt
    r_1n = rb(B+1) - 0.25*A*dt
    r_2o = rb(B+2) + 0.25*A*dt
    r_2n = rb(B+2) - 0.25*A*dt

    x_1o = r_1o / (rwt+r_1o)
    x_1n = r_1n / (rwt+r_1n)
    x_2o = r_2o / (rwt+r_2o)
    x_2n = r_2n / (rwt+r_2n)

    j_1o = j_wt_o*(x_1o-xb(B+2))*(x_1o-xb(B+3))/(x_wt_o-xb(B+2))/(x_wt_o-xb(B+3)) &
         + jos(:,:,B+2)*(x_1o-x_wt_o)*(x_1o-xb(B+3))/(xb(B+2)-x_wt_o)/(xb(B+2)-xb(B+3)) &
         + jos(:,:,B+3)*(x_1o-x_wt_o)*(x_1o-xb(B+2))/(xb(B+3)-x_wt_o)/(xb(B+3)-xb(B+2))

    j_1n = j_wt*(x_1n-xb(B+1))/(x_wt-xb(B+1))&
         + jns(:,:,B+1)*(x_1n-x_wt)/(xb(B+1)-x_wt);

    j_2o = j_wt_o*(x_2o-xb(B+2))*(x_2o-xb(B+3))/(x_wt_o-xb(B+2))/(x_wt_o-xb(B+3)) &
         + jos(:,:,B+2)*(x_2o-x_wt_o)*(x_2o-xb(B+3))/(xb(B+2)-x_wt_o)/(xb(B+2)-xb(B+3)) &
         + jos(:,:,B+3)*(x_2o-x_wt_o)*(x_2o-xb(B+2))/(xb(B+3)-x_wt_o)/(xb(B+3)-xb(B+2))

    j_2n = ( r_1n * j_1n - r_1o * j_1o + r_2o * j_2o &
           + 0.5 * dt * (r_2n-r_1n) * (P_u+h_c) ) / r_2n

    j_Bp2 = ( j_2n - j_wt*(x_2n-x_1n)*(x_2n-xb(B+2))/(x_wt-x_1n)/(x_wt-xb(B+2)) &
            - j_1n*(x_2n-x_wt)*(x_2n-xb(B+2))/(x_1n-x_wt)/(x_1n-xb(B+2)) ) &
          / ((x_2n-x_wt)*(x_2n-x_1n)/(xb(B+2)-x_wt)/(xb(B+2)-x_1n))

  end subroutine NullEvol_j_Bp2

  subroutine NullEvol_j (i, jns, dxjos, nuns, ckns, bns, cbns, qns, uns, wns,&
       jos, nuos, ckos, bos, cbos, qos, uos, wos, &
       x_wt, j_wt, beta_wt, q_wt, u_wt, w_wt,&
       x_wt_o, j_wt_o, beta_wt_o, q_wt_o, u_wt_o, w_wt_o,&
       mask_Bp1, mask, eth4_mask, dissip_mask, dissip, dt, dissip_fudge,&
       dissip_fudge_maxx, dissip_eps_u, dissip_eps_x, first_order_scheme)

    CCTK_REAL,                                    intent(in)     :: dt
    CCTK_INT,                                     intent (in)    :: i, dissip_fudge, first_order_scheme
    CCTK_COMPLEX, dimension (lsh(1), lsh(2), nx), intent (inout) :: jns, dxjos
    CCTK_COMPLEX, dimension (lsh(1), lsh(2), nx), intent (in)    :: jos, nuos, ckos, cbos, uos, qos
    CCTK_COMPLEX, dimension (lsh(1), lsh(2), nx), intent (in)    :: nuns, ckns, cbns, uns, qns
    CCTK_REAL,    dimension (lsh(1), lsh(2), nx), intent (in)    :: bns, bos, wns, wos
    CCTK_COMPLEX, dimension (lsh(1), lsh(2)),     intent (in)    :: j_wt, q_wt, u_wt
    CCTK_COMPLEX, dimension (lsh(1), lsh(2)),     intent (in)    :: j_wt_o, q_wt_o, u_wt_o
    CCTK_REAL,    dimension (lsh(1), lsh(2)),     intent (in)    :: x_wt, beta_wt, w_wt
    CCTK_REAL,    dimension (lsh(1), lsh(2)),     intent (in)    :: x_wt_o, beta_wt_o, w_wt_o
    CCTK_INT,     dimension (lsh(1), lsh(2)),     intent (in)    :: mask_Bp1, mask, eth4_mask
    CCTK_REAL,    dimension (lsh(1), lsh(2)),     intent (in)    :: dissip_mask
    CCTK_REAL,                                    intent(in)     :: dissip, dissip_fudge_maxx, dissip_eps_u,&
                                                                    dissip_eps_x

    CCTK_INT, parameter :: izero = 0
    integer, parameter :: ik = kind(izero)

    logical, save :: FirstTime = .true.
    CCTK_INT :: ll

    !! probably better to make these into allocated arrays
    CCTK_REAL, dimension (:,:), allocatable, save :: &
         xc, dx_beta, W, f1, f2

    CCTK_COMPLEX, dimension (:,:), allocatable, save :: &
         J, d2x_jave, dx_J, dx_Jo, dx_jave, fac, P_u, h_c, jave, tmp1, tmp2

    CCTK_REAL, dimension(:), allocatable, save :: disarr, tmp_dis

    CCTK_INT cran
    CCTK_INT, parameter :: cran_max = 2

    if (FirstTime) then
       FirstTime = .false.
       allocate(xc(lsh(1),lsh(2)), dx_beta(lsh(1),lsh(2)), W(lsh(1),lsh(2)),&
            f1(lsh(1),lsh(2)), f2(lsh(1),lsh(2)), P_u(lsh(1),lsh(2)),&
            J(lsh(1),lsh(2)), d2x_jave(lsh(1),lsh(2)), dx_J(lsh(1),lsh(2)),&
            dx_Jo(lsh(1),lsh(2)), dx_jave(lsh(1),lsh(2)), tmp1(lsh(1),lsh(2)),&
            tmp2(lsh(1),lsh(2)), fac(lsh(1),lsh(2)), h_c(lsh(1),lsh(2)), jave(lsh(1),lsh(2)))

       xc=0; dx_beta=0; W=0; f1=0; f2=0; P_u=0
       J=0; d2x_jave=0; dx_J=0; dx_Jo=0; dx_jave=0;
       fac=0; h_c=0; jave=0; 

       allocate(disarr(nx), tmp_dis(nx))
       disarr = dissip
       if (dissip_fudge .eq. 1) then
          do ll = 1, nx
             if (xb(ll) .lt. DMIN ) then
                disarr(ll) = DISMAX
             end if
             if (xb(ll) .gt. DMIN .and. xb(ll) .lt. dissip_fudge_maxx) then
                disarr(ll) = dissip + (DISMAX - dissip) *&
                     (1.0d0 - (xb(ll) - DMIN)**3/(dissip_fudge_maxx - DMIN)**3) 
             end if
          end do

          tmp_dis(1) = disarr(1)

          do ll = 2, nx
             tmp_dis(ll) = .5d0 * (disarr(ll) + disarr(ll-1))
          end do
          disarr = tmp_dis
       end if

    end if

    ! Note that mask_Bp1 is set to zero for i <= B+1 and is set to one elsewhere
    ! in a similar way mask is set to zero for i <= B+2 and is set to one elsewhere.

    ! we first apply the B+1 type algorithm

    if(minval(mask_Bp1).eq.0) then

      call NullEvol_j_rhs (i, jns, nuns, ckns, bns, cbns, qns, uns, wns, &
         jos, nuos, ckos, bos, cbos, qos, uos, wos,&
         J, dx_J, dx_Jo, W, dx_beta, h_c, fac, f1, f2, xc, dt,&
         x_wt, j_wt, beta_wt, q_wt, u_wt, w_wt,&
         x_wt_o, j_wt_o, beta_wt_o, q_wt_o, u_wt_o, w_wt_o,&
         first_order_scheme, 1_ik)

      call NullEvol_P_u(i, 0_ik, dt, jns, jos, xc, fac, dx_beta, J, W, P_u,&
                        x_wt, x_wt_o, j_wt, j_wt_o, 1_ik)

      call NullEvol_j_Bp1(i, J, jos, P_u, h_c, xc, W,&
                          j_wt, j_wt_o,  x_wt, x_wt_o, dt)

      ! Here we set jns(:,:,i).   Then later on, with the mask,
      ! the points that are not 'boundary' type will be set by the regular
      ! evolution mask.

      jns(:,:,i) = J

    end if

    if(any(mask_Bp1.eq.1 .and. mask.eq.0)) then

      call NullEvol_j_rhs (i, jns, nuns, ckns, bns, cbns, qns, uns, wns, &
         jos, nuos, ckos, bos, cbos, qos, uos, wos,&
         J, dx_J, dx_Jo, W, dx_beta, h_c, fac, f1, f2, xc, dt,&
         x_wt, j_wt, beta_wt, q_wt, u_wt, w_wt,&
         x_wt_o, j_wt_o, beta_wt_o, q_wt_o, u_wt_o, w_wt_o,&
         first_order_scheme, 2_ik)

      call NullEvol_P_u(i, 0_ik, dt, jns, jos, xc, fac, dx_beta, J, W, P_u,&
                        x_wt, x_wt_o, j_wt, j_wt_o, 2_ik)

      call NullEvol_j_Bp2(i, J, jns, jos, P_u, h_c, xc, W,&
                          j_wt, j_wt_o,  x_wt, x_wt_o, dt)

      ! Here we set jns(:,:,i) for the non-Bp1 points only.  In the next
      ! step the points that are not 'boundary' type will be set by the regular
      ! evolution mask.

      jns(:,:,i) = (1-mask_Bp1)*jns(:,:,i) + mask_Bp1*J

    end if

    ! at this point all points i<=B+1 are filled using the B+1 scheme
    ! next we go through again and see if there are B+2 type points


    call NullEvol_j_rhs (i, jns, nuns, ckns, bns, cbns, qns, uns, wns, &
       jos, nuos, ckos, bos, cbos, qos, uos, wos,&
       J, dx_J, dx_Jo, W, dx_beta, h_c, fac, f1, f2, xc, dt,&
       x_wt, j_wt, beta_wt, q_wt, u_wt, w_wt,&
       x_wt_o, j_wt_o, beta_wt_o, q_wt_o, u_wt_o, w_wt_o,&
       first_order_scheme, 0_ik)

    ! add angular dissipation

    !  h_c -> h_c - eps h^3 * (\eth \ethb)^2 ( J + (1-x) J_{,x} )
    !  note that we use here, temporarily, variables used elsewhere

    if (abs(dissip_eps_u).gt.1.0e-15) then

!Angular dissipation from Psi paper
       tmp1 = J + xc*(1-xc)*dx_J
       call NullInterp_d2 (tmp2,  tmp1, 2_ik,  -1_ik, -1_ik)
       tmp2 = dissip_mask * tmp2
       call NullInterp_d2 (tmp1,  tmp2, 0_ik,  +1_ik, +1_ik)
       h_c = h_c - dissip_eps_u * eth4_mask * tmp1

    end if

    ! add angular dissipation along the radial axis
    ! -- I'll just do this without any factors of r

    !  h_c -> h_c - eps h^3 * (\eth \ethb)^2 J_{,u}
    !  note that we use here, temporarily, variables used elsewhere

    if (abs(dissip_eps_x).gt.1.0e-15) then

!This is the radial dissipation
       tmp1 = xc*(jns(:,:,i-1)-jos(:,:,i-1))/dt
       call NullInterp_d2 (tmp2,  tmp1, 2_ik,  -1_ik, -1_ik)
       tmp2 = dissip_mask * tmp2
       call NullInterp_d2 (tmp1,  tmp2, 0_ik,  +1_ik, +1_ik)
       h_c = h_c - dissip_eps_x * eth4_mask * tmp1
    end if

    ! interpolate the $f = xb j$ values
    !crank initial


    jns(:,:,i) = (1-mask)* jns(:,:,i) + jos(:,:,i) * mask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (i .ne. nx) then

       jave = 0.5 * (xb(i) * jos(:,:,i) + xb(i-1) * jos(:,:,i-1)) &
            + (1. - xc) * ((1. - 0.75 * disarr(i)) &
            * (xb(i) * jos(:,:,i) - xb(i-1) * jos(:,:,i-1)) &
            + 0.25 * disarr(i) &
            * (xb(i+1) * jos(:,:,i+1) - xb(i-2) *&
            jos(:,:,i-2))) / dx 

       d2x_jave = xb(i+1) * jos(:,:,i+1) - 2. * xb(i) * jos(:,:,i) &
            + xb(i-1) * jos(:,:,i-1)

    else

       jave  = 0.5 * (xb(i) * jos(:,:,i) + xb(i-1) * jos(:,:,i-1) ) &
            + (1. - xc) * (xb(i) * jos(:,:,i) - xb(i-1) * jos(:,:,i-1) ) / dx

       d2x_jave = ( 2. * xb(i)   * jos(:,:,i)   -&
            5. * xb(i-1) * jos(:,:,i-1) &
            + 4. * xb(i-2) * jos(:,:,i-2) -&
            xb(i-3) * jos(:,:,i-3) )

    end if

    do cran = 1, cran_max

       call NullEvol_P_u(i, cran, dt, jns, jos, xc, fac, dx_beta, J, W, P_u,&
                         x_wt, x_wt_o, j_wt, j_wt_o, 0_ik)

       jns(:,:,i) = jns(:,:,i) * (1 - mask) &
            + mask * (- xb(i-1) * jns(:,:,i-1) * ( (1. - f1 * 0.25) &
            * (0.5 - (1. - xc) / dx) + 0.5 * f2 / dx**2) &
            + xb(i-2) * jns(:,:,i-2) * f2 * 0.25 / dx**2 &
            + jave * (1. + f1 * 0.25) + f2/dx**2 * d2x_jave*0.25&
            + 0.5 * (h_c + P_u) * dt) &
            / (xb(i) * ((1. - f1 * 0.25)&
            * (0.5 + (1. - xc) / dx) - 0.25 * f2 / dx**2))

    end do

!fill the radial derivative
    dxjos(:,:,i) = dx_Jo

  end subroutine NullEvol_j


end module NullEvol_Evol
