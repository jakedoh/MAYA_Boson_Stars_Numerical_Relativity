! vim: syntax=fortran
#include "cctk.h"

module NullEvol_hyper_q
    use NullInterp
    use NullGrid_Vars
    private
    public NullEvol_q
contains

  subroutine NullEvol_q_rhs(i, jns, nuns, ckns, bns, cbns, j_wt, beta_wt, x_wt, rhs,&
                            first_order_scheme, compute_at_Bp1)
    implicit none

    CCTK_INT, parameter :: izero = 0
    integer, parameter :: ik = kind(izero)

    CCTK_INT,                                   intent (in)    :: i, first_order_scheme
    CCTK_COMPLEX, dimension (lsh(1),lsh(2),nx), intent (in)    :: jns, nuns, ckns, cbns
    CCTK_REAL,    dimension (lsh(1),lsh(2),nx), intent (in)    :: bns 
    logical,                                    intent (in)    :: compute_at_Bp1
    CCTK_COMPLEX, dimension (lsh(1),lsh(2)),    intent (in)    :: j_wt
    CCTK_REAL,    dimension (lsh(1),lsh(2)),    intent (in)    :: beta_wt, x_wt
    CCTK_COMPLEX, dimension (lsh(1),lsh(2)),    intent (inout) :: rhs

    CCTK_REAL :: xhere
    logical, save :: FirstTime = .true. 
    CCTK_COMPLEX, dimension (:,:), allocatable, save :: &
         J, nu, ck, cb, dx_J, dx_nu, dx_ck, dx_cb, eth_J, eth_dx_J, n_Q

    CCTK_REAL, dimension (:,:), allocatable, save :: K, dx_K, beta, dx_beta


    if (FirstTime) then
       FirstTime=.false.
       allocate(J(lsh(1),lsh(2)), nu(lsh(1),lsh(2)), ck(lsh(1),lsh(2)), cb(lsh(1),lsh(2)),&
            dx_J(lsh(1),lsh(2)), dx_nu(lsh(1),lsh(2)), dx_ck(lsh(1),lsh(2)), dx_cb(lsh(1),lsh(2)),&
            eth_J(lsh(1),lsh(2)), eth_dx_J(lsh(1),lsh(2)), &
            n_Q(lsh(1),lsh(2)), K(lsh(1),lsh(2)), dx_K(lsh(1),lsh(2)),&
            beta(lsh(1),lsh(2)), dx_beta(lsh(1),lsh(2)))
       J=0; nu=0; ck=0; cb=0; dx_J=0; dx_nu=0; dx_ck=0; dx_cb=0;
       eth_J=0; eth_dx_J=0; n_Q=0; K=0; dx_K=0;
       beta=0; dx_beta = 0
    end if

    if(.not.compute_at_Bp1) then
       xhere = xbh(i-1)

       J = 0.5 * ( jns(:,:,i-1) + jns(:,:,i) )
       dx_J = (jns(:,:,i) - jns(:,:,i-1)) / dx
    else ! compute_at_Bp1
       xhere = xb(i)

       J = jns(:,:,i)
       dx_J = (jns(:,:,i) - j_wt) / (xhere-x_wt)
    end if ! compute_at_Bp1

    K = sqrt(1. + dble(J * conjg(J)))
    dx_K = dble(dx_J * conjg(J)) / K

    if (first_order_scheme.ne.0) then

       if(.not.compute_at_Bp1) then
          nu = 0.5 * ( nuns(:,:,i-1) + nuns(:,:,i) )
          ck = 0.5 * ( ckns(:,:,i-1) + ckns(:,:,i) )
          cb = 0.5 * ( cbns(:,:,i-1) + cbns(:,:,i) )

          dx_nu = (nuns(:,:,i) - nuns(:,:,i-1)) / dx
          dx_ck = (ckns(:,:,i) - ckns(:,:,i-1)) / dx
          dx_cb = (cbns(:,:,i) - cbns(:,:,i-1)) / dx

       else ! compute_at_Bp1

          nu = nuns(:,:,i)
          ck = ckns(:,:,i)
          cb = cbns(:,:,i)

          call NullInterp_d1(dx_nu, dx_J, 2_ik, -1_ik)
          call NullInterp_d1(dx_ck, dcmplx(dx_K), 0_ik, +1_ik)
          call NullInterp_d1(dx_cb, dcmplx(dx_beta), 0_ik, +1_ik)

       end if ! compute_at_Bp1

    else ! first_order_scheme

       ! nu = ethb_J
       ! ck = eth_K
       ! cb = eth_beta

       if(.not.compute_at_Bp1) then
          beta = 0.5 * ( bns(:,:,i-1) + bns(:,:,i) )
          dx_beta = (bns(:,:,i) - bns(:,:,i-1)) / dx
       else ! compute_at_Bp1
          beta = bns(:,:,i)
          dx_beta = (bns(:,:,i) - beta_wt) / (xhere-x_wt)
       end if ! compute_at_Bp1

       call NullInterp_d1(nu, J, 2_ik, -1_ik) 
       call NullInterp_d1(dx_nu, dx_J, 2_ik, -1_ik) 

       call NullInterp_d1(ck, dcmplx(K), 0_ik, +1_ik) 
       call NullInterp_d1(dx_ck, dcmplx(dx_K), 0_ik, +1_ik) 

       call NullInterp_d1(cb, dcmplx(beta), 0_ik, +1_ik) 
       call NullInterp_d1(dx_cb, dcmplx(dx_beta), 0_ik, +1_ik) 

    end if ! first_order_scheme

    call NullInterp_d1(eth_J, J, 2_ik, 1_ik) 

    call NullInterp_d1(eth_dx_J, dx_J, 2_ik, 1_ik) 

    n_Q = 2. * conjg(J) * conjg(ck) * dx_J * J &
         + 3. * K * conjg(nu) * (dx_J * K - dx_K * J) &
         + 4. * J * ck * dx_K * conjg(J) &
         - conjg(J) * conjg(nu) * dx_J * J &
         - 2. * K * conjg(ck) * dx_K * J &
         - 4. * K * ck * dx_J * conjg(J) &
         - 2. * K * ck * conjg(dx_J) * J &
         + J**2 * conjg(eth_J) * dx_K &
         + J**2 * conjg(nu) * conjg(dx_J) &
         + 2. * conjg(J)**2 * eth_J * dx_J &
         + 2. * K**2 * dx_K * (ck + nu) &
         + K**2 * eth_J * conjg(dx_J) &
         - K * conjg(eth_J) * dx_J * J &
         - conjg(J) * nu * (dx_K * J + dx_J * K) &
         - 3. * conjg(J) * eth_J * dx_K * K &
         + 2. * dx_ck + 2. * dx_nu &
         - 2. * K * (dx_nu + dx_ck) &
         + 2. * J * conjg(dx_ck) &
         + 2. * conjg(J) * eth_dx_J

    rhs = - xhere * (1. - xhere) * (dx_nu + dx_ck - 2. * dx_cb) &
         - 4. * cb + 0.5 * (1. - xhere) * (xhere) * n_Q  

  end subroutine NullEvol_q_rhs

  subroutine NullEvol_q (i, jns, nuns, ckns, bns, cbns, qns, q_wt, j_wt, beta_wt, x_wt, &
                         mask, eth4_mask, dissip_mask, dissip_eps, first_order_scheme)
    implicit none

    CCTK_INT, parameter :: izero = 0
    integer, parameter :: ik = kind(izero)

    CCTK_INT,                                   intent (in)    :: i, first_order_scheme
    CCTK_COMPLEX, dimension (lsh(1),lsh(2),nx), intent (in)    :: jns, nuns, ckns, cbns
    CCTK_REAL,    dimension (lsh(1),lsh(2),nx), intent (in)    :: bns 
    CCTK_COMPLEX, dimension (lsh(1),lsh(2),nx), intent (inout) :: qns
    CCTK_INT,     dimension (lsh(1),lsh(2)),    intent (in)    :: mask, eth4_mask
    CCTK_COMPLEX, dimension (lsh(1),lsh(2)),    intent (in)    :: q_wt, j_wt
    CCTK_REAL,    dimension (lsh(1),lsh(2)),    intent (in)    :: beta_wt, x_wt
    CCTK_REAL,    dimension (lsh(1),lsh(2)),    intent (in)    :: dissip_mask
    CCTK_REAL,                                  intent(in)     :: dissip_eps

    CCTK_COMPLEX, dimension (lsh(1),lsh(2)) :: rhs,tmp
    CCTK_INT i1, i2

    call NullEvol_q_rhs(i, jns, nuns, ckns, bns, cbns, j_wt, beta_wt, x_wt, rhs,&
                        first_order_scheme, .false.)

    if(minval(mask).eq.0) then
       ! this will avoid use of points inside the boundary except for angular derivatives
       ! in case of a non-spherical boundary
       call NullEvol_q_rhs(i, jns, nuns, ckns, bns, cbns, j_wt, beta_wt, x_wt, tmp,&
                           first_order_scheme, .true.)
       rhs = mask*rhs + (1-mask)*tmp

    end if

    ! angular dissipation

    if (abs(dissip_eps).gt.1.0e-15) then

       call NullInterp_d2(tmp, qns(:,:,i-1), 1_ik, -1_ik, -1_ik)
       call NullInterp_d2(tmp, dissip_mask * tmp, -1_ik, +1_ik, +1_ik)
       rhs = rhs + dissip_eps * eth4_mask*tmp

    end if

    qns(:,:,i) = 0.5*rhs + mask * ( qns(:,:,i-1) - 0.5*rhs ) * ((xb(i)-1)*xb(i-1)/xb(i)/(xb(i-1)-1))**2&
                      + (1-mask)* ( q_wt         - 0.5*rhs ) * ( (xb(i)-1)*x_wt/xb(i)/(x_wt-1)  )**2

    ! On points that are very near the boundary, copy the boundary value as it will be more accurate.
    if(minval(mask).eq.0) then
      do i2 = 1, lsh(2)
         do i1 = 1, lsh(1)
             if(abs(xb(i)-x_wt(i1,i2)).lt.1.e-5*dx) qns(i1,i2,i) = q_wt(i1,i2)
         end do
      end do
    end if

  end subroutine NullEvol_q
end module NullEvol_hyper_q   
