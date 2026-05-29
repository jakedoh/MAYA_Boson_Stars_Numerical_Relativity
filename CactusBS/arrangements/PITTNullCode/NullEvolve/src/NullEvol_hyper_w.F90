! vim: syntax=fortran
#include "cctk.h"

module NullEvol_hyper_w
    use NullInterp
    use NullGrid_Vars
    private
    public NullEvol_w
contains

  subroutine NullEvol_w_rhs(i, jns, nuns, ckns, bns, cbns, uns, qns,&
                            u_wt, x_wt, rhs,&
                            first_order_scheme, compute_at_Bp1)
    implicit none

    CCTK_INT, parameter :: izero = 0
    integer, parameter :: ik = kind(izero)

    CCTK_INT,                                   intent (in)    :: i, first_order_scheme
    CCTK_COMPLEX, dimension (lsh(1),lsh(2),nx), intent (in)    :: jns, nuns, ckns
    CCTK_REAL,    dimension (lsh(1),lsh(2),nx), intent (in)    :: bns 
    CCTK_COMPLEX, dimension (lsh(1),lsh(2),nx), intent (in)    :: cbns, qns, uns
    CCTK_COMPLEX, dimension (lsh(1),lsh(2)),    intent (in)    :: u_wt
    CCTK_REAL,    dimension (lsh(1),lsh(2)),    intent (in)    :: x_wt
    logical,                                    intent (in)    :: compute_at_Bp1
    CCTK_REAL,    dimension (lsh(1),lsh(2)),    intent (out)   :: rhs

    CCTK_REAL xhere
    logical, save :: FirstTime = .true. 
    CCTK_COMPLEX, dimension (:,:), allocatable, save :: &
         ck, cB, nu, J, ethb_U, ethb_dx_U, eth_J, ethb_nu, ethb_ck, &
         eth_cb, ethb_cb, U, dx_U, Q

    CCTK_REAL, dimension (:,:), allocatable, save :: &
         beta, K, Ricci

    if (FirstTime) then
       FirstTime=.false.
       allocate(ck(lsh(1),lsh(2)), cB(lsh(1),lsh(2)),&
            nu(lsh(1),lsh(2)), J(lsh(1),lsh(2)), ethb_U(lsh(1),lsh(2)), ethb_dx_U(lsh(1),lsh(2)),&
            eth_J(lsh(1),lsh(2)), ethb_nu(lsh(1),lsh(2)), ethb_ck(lsh(1),lsh(2)),&
            eth_cb(lsh(1),lsh(2)), ethb_cb(lsh(1),lsh(2)), beta(lsh(1),lsh(2)), K(lsh(1),lsh(2)),&
            Ricci(lsh(1),lsh(2)), U(lsh(1),lsh(2)), dx_U(lsh(1),lsh(2)), Q(lsh(1),lsh(2)))

       ck=0; cB=0; nu=0; J=0; ethb_U=0; ethb_dx_U=0;
       eth_J=0; ethb_nu=0; ethb_ck=0; eth_cb=0; ethb_cb=0; beta=0; K=0;
       Ricci=0; 
    end if

    if(.not. compute_at_Bp1) then
       xhere = xbh(i-1)
       beta = 0.5 * ( bns(:,:,i) +  bns(:,:,i-1))
       J    = 0.5 * ( jns(:,:,i) +  jns(:,:,i-1))
       Q    = 0.5 * ( qns(:,:,i) +  qns(:,:,i-1))
       U    = uns(:,:,i-1)
    else
       xhere = xb(i) ! note that i = B+1
       beta = bns(:,:,i)
       J    = jns(:,:,i)
       Q    = qns(:,:,i)
       ! interpolate onto x(i) using x(i+1/2) and x_wt
       U    = u_wt * (xhere-xbh(i))/(x_wt-xbh(i)) &
            + uns(:,:,i) * (xhere-x_wt)/(xbh(i)-x_wt)
    end if

    K    = sqrt(1. + dble(J * conjg(J)))
    dx_U = exp(2. * beta) / rwt / xhere**2 &
                 * ( - J * conjg(Q) + K * Q )

    if (first_order_scheme.ne.0) then

       if(.not. compute_at_Bp1) then
          ck   = 0.5 * (ckns(:,:,i) + ckns(:,:,i-1))
          cB   = 0.5 * (cbns(:,:,i) + cbns(:,:,i-1))
          nu   = 0.5 * (nuns(:,:,i) + nuns(:,:,i-1))
       else ! compute_at_Bp1
          ck   = ckns(:,:,i)
          cB   = cbns(:,:,i)
          nu   = nuns(:,:,i)
       end if ! compute_at_Bp1

    else ! first_order_scheme

       ! nu = ethb_J
       ! ck = eth_K
       ! cb = eth_beta

       call NullInterp_d1(nu, J,            2_ik, -1_ik)
       call NullInterp_d1(ck, dcmplx(K),    0_ik, +1_ik)
       call NullInterp_d1(cb, dcmplx(beta), 0_ik, +1_ik)

    end if ! first_order_scheme

    call NullInterp_d1(ethb_U,    U ,   1_ik, -1_ik)
    call NullInterp_d1(ethb_dx_U, dx_U, 1_ik, -1_ik)

    call NullInterp_d1(eth_J,   J,  2_ik,  1_ik)

    if (first_order_scheme.ne.0) then

       call NullInterp_d1(ethb_nu, nu, 1_ik, -1_ik)
       call NullInterp_d1(ethb_ck, ck, 1_ik, -1_ik)
       call NullInterp_d1(ethb_cb, cb, 1_ik, -1_ik)
       call NullInterp_d1(eth_cb,  cb, 1_ik,  1_ik)

    else ! first_order_scheme

       ! nu = ethb_J
       ! ck = eth_K
       ! cb = eth_beta

       call NullInterp_d2(ethb_nu, J,            2_ik, -1_ik, -1_ik)
       call NullInterp_d2(ethb_ck, dcmplx(K),    0_ik, -1_ik, +1_ik) ! ethb_eth
       call NullInterp_d2(ethb_cb, dcmplx(beta), 0_ik, -1_ik, +1_ik) ! ethb_eth
       call NullInterp_d2(eth_cb,  dcmplx(beta), 0_ik, +1_ik, +1_ik) ! eth_eth

    end if ! first_order_scheme

    Ricci = dble( 2. * k + ethb_nu - ethb_ck &
         + ( eth_J * conjg(eth_J) - nu * conjg(nu) ) / (4. * k) )

    rhs = dble( &
          (1. - xhere) / (rwt * xhere) * ( exp(2. * beta) * (0.5 * Ricci &
         - K * (ethb_cb + cb * conjg(cb)) + conjg(J) * (eth_cb + cb**2) &
         + conjg(cb) * (nu - ck) ) - 1.) &
         + 2. * ethb_U + 0.5 * xhere * (1. - xhere) * ethb_dx_U &
         - 0.25 * exp(-2. * beta) * rwt * xhere**3 * (1. - xhere) &
         * conjg(dx_U) * (K * dx_U + J * conjg(dx_U)) )

 
  end subroutine NullEvol_w_rhs

  subroutine NullEvol_w (i, jns, nuns, ckns, bns, cbns, uns, qns, wns, &
       w_wt, u_wt, x_wt, mask, eth4_mask, dissip_mask, dissip_eps, first_order_scheme)
    implicit none

    CCTK_INT, parameter :: izero = 0
    integer, parameter :: ik = kind(izero)

    CCTK_INT,                                   intent (in)    :: i, first_order_scheme
    CCTK_COMPLEX, dimension (lsh(1),lsh(2),nx), intent (in)    :: jns, nuns, ckns
    CCTK_REAL,    dimension (lsh(1),lsh(2),nx), intent (in)    :: bns 
    CCTK_COMPLEX, dimension (lsh(1),lsh(2),nx), intent (in)    :: cbns, uns, qns
    CCTK_REAL,    dimension (lsh(1),lsh(2),nx), intent (inout) :: wns
    CCTK_REAL,    dimension (lsh(1),lsh(2)),    intent (in)    :: w_wt, x_wt
    CCTK_COMPLEX, dimension (lsh(1),lsh(2)),    intent (in)    :: u_wt
    CCTK_INT,     dimension (lsh(1),lsh(2)),    intent (in)    :: mask, eth4_mask
    CCTK_REAL,    dimension (lsh(1),lsh(2)),    intent (in)    :: dissip_mask
    CCTK_REAL,                                  intent (in)    :: dissip_eps
    ! angular dissipation

    CCTK_REAL,    dimension (:,:), allocatable, save :: rhs, rhs_Bp1
    CCTK_COMPLEX, dimension (:,:), allocatable, save :: e4_w, cW

    logical, save :: FirstTime = .true. 
    integer :: i1, i2

    if (FirstTime) then
       FirstTime=.false.
       allocate(rhs(lsh(1),lsh(2)), e4_w(lsh(1),lsh(2)), &
                cW(lsh(1),lsh(2)), rhs_Bp1(lsh(1),lsh(2)))
    end if

    ! compute the rhs at x(i-1/2)

    call NullEvol_w_rhs(i, jns, nuns, ckns, bns, cbns, uns, qns, u_wt, x_wt, rhs,&
                        first_order_scheme, .false.)

    if(minval(mask).eq.0) then
       ! compute the rhs at x(i) for sake of interpolation:
       call NullEvol_w_rhs(i, jns, nuns, ckns, bns, cbns, uns, qns, u_wt, x_wt, rhs_Bp1,&
                           first_order_scheme, .true.)
       ! Note that the following expression uses a taylor expansion arount the wt.
       ! The problem with this is that, for data spoied by (constraint violating)
       ! errors, this value can be inconsistent with W_{,x} as provided by the 
       ! Einstein equations on the light-cone.
!      rhs = mask * rhs + (1-mask)*0.5*(rhs_Bp1 + 2*w_wt+x_wt*(1-x_wt)*w_x_wt)
       ! As an alternative, here we compute W_{x} at x(i) and x(i-1/2) and
       ! then we extrapolate onto the wt.  The hope is that this will converge
       ! to the hpyersurface-equation dictated value, even if the WT data is
       ! not consistent.
!      rhs = mask * rhs + (1-mask)*(rhs_Bp1 + (rhs_Bp1-rhs)*(x_wt-xb(i))/(0.5*dx))
       ! As an even safer alternative, we use rhs_Bp1 for the 1st update.  This
       ! will stay away from using values on the other side of the WT except
       ! for the angular derivatives that may reach over in case of a non-spherical
       ! extraction world-tube.
       rhs = mask * rhs + (1-mask)*rhs_Bp1
    end if

    if (abs(dissip_eps).gt.1.0e-15) then

       cW = wns(:,:,i-1)
       call NullInterp_d2(e4_w, cW, 0_ik, 1_ik, -1_ik)
       cW = dissip_mask * e4_w;
       call NullInterp_d2(e4_w, cW, 0_ik, 1_ik, -1_ik)

       rhs = rhs + dissip_eps * eth4_mask *dble(e4_w)

    end if

    wns(:,:,i) = 0.5*rhs + mask * ( wns(:,:,i-1) - 0.5*rhs ) * ((xb(i)-1)*xb(i-1)/xb(i)/(xb(i-1)-1))**2&
                      + (1-mask)* ( w_wt         - 0.5*rhs ) * ( (xb(i)-1)*x_wt/xb(i)/(x_wt-1)  )**2

    ! On points that are very near the boundary copy the boundary value as it will be more accurate.
    if(minval(mask).eq.0) then
      do i2 = 1, lsh(2)
         do i1 = 1, lsh(1)
             if(abs(xb(i)-x_wt(i1,i2)).lt.1.e-5*dx) wns(i1,i2,i) = w_wt(i1,i2)
         end do
      end do
    end if


  end subroutine NullEvol_w

end module NullEvol_hyper_w
