! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

module NullExact_Analytic

  use cctk
  implicit none

contains

  !---------------------------------------------------------------------
  ! pointwise defined testbed functions
  ! input parameters:
  !    - patch ID
  !    - two angular coordinates
  !    - compactified radial coordinate
  !    - time
  !---------------------------------------------------------------------

!=====================================================================
!  A particular linearized solution of R01 where b0=0, C2=3/omm^2*C1
!=====================================================================
   function R01_out(s,l_in_Ylm,C1r,C1i,C2r,C2i,b0r,b0i,omm,time) result(R01)

      CCTK_INT,  intent(in)  :: l_in_Ylm
      CCTK_REAL, intent(in)  :: s,C1r,C1i,C2r,C2i,b0r,b0i,omm,time
      CCTK_REAL :: R01, inv_r
      CCTK_COMPLEX  :: C1,C2,b0,einuu,ii, W,W_r,rW_rr, ethbU, ethbU_r
      CCTK_COMPLEX, parameter  :: czero = (0.0d0,0.0d0)
      logical :: first
      data       first /.true./
      save first

      call complexify(C1r,C1i,C2r,C2i,b0r,b0i,omm,time,C1,C2,b0,einuu,ii)
      inv_r = x_to_r_inverse(s)

      if( first .and. (b0.ne.czero .or. abs(C2-3.0d0*C1/omm**2).gt.1d-14) ) then
	 call CCTK_WARN(1, "INVALID PARAMETERS FOR LINEARIZED SOLN OF R01")
         first = .false.
      end if

      select case (l_in_Ylm)
      case (2)
         W    = 3.0d0*C1*(ii*inv_r**3/omm - inv_r**2 + 0.5d0*inv_r**4/omm**2)
	 W_r  = 6.0d0*C1*(inv_r**3 - 1.5d0*ii*inv_r**4/omm - inv_r**5/omm**2)
         rW_rr = 18.0d0*C1*(2*ii*inv_r**4/omm - inv_r**3 + 5.0d0/3.0d0*inv_r**5/omm**2)
         ethbU   = -3.0d0*C1*(inv_r**2 + 2.0d0*ii*inv_r**3/omm + 1.5d0*inv_r**4/omm**2)
         ethbU_r = 6.0d0*C1*(inv_r**3+3.0d0*ii*inv_r**4/omm + 3.0d0*inv_r**5/omm**2)
	 R01 = dble((W*inv_r + 2.0d0*W_r + 0.5d0*rW_rr - inv_r*ethbU - 0.5d0*ethbU_r)*einuu)

      case default
	 call CCTK_WARN(1, "unencoded value for l_in_Ylm")
      end select

   end function R01_out

!.....................................................................
   subroutine NullExact_Analytic_R01_2D(patchID,n1,n2,qs,ps,s,time,R01,Ylm)

      CCTK_INT,  intent(in)  :: n1, n2, patchID
      CCTK_REAL, intent(in)  :: qs(n1,n2), ps(n1,n2), s, time
      CCTK_COMPLEX, intent(in) :: Ylm(n1,n2,2)
      CCTK_REAL, intent(out) :: R01(n1,n2)
      CCTK_REAL :: R01_i
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

      select case (testbed_ID)
      case (2) ! linearized solution
	 R01_i = R01_out(s,l_in_Ylm,Constant_C1r,Constant_C1i,Constant_C2r,Constant_C2i &
			 ,beta0r,beta0i,omm,time)
	 R01 = R01_i*Ylm(:,:,patchID)
      case default
	 R01 = 0
	 call CCTK_WARN(1, "unrecognized value for testbed_ID")
      end select

   end subroutine NullExact_Analytic_R01_2D
!=====================================================================

  function eth2fac(i) result(e2)

    CCTK_INT, intent(in)  :: i
    CCTK_REAL :: e2

      e2 = dsqrt(1.0d0*(i-1)*i*(i+1)*(i+2))

  end function eth2fac

  function eth1fac(i) result(e1)

    CCTK_INT, intent(in)  :: i
    CCTK_REAL :: e1

      e1 = dsqrt(1.0d0*i*(i+1))

  end function eth1fac

  subroutine complexify(C1r,C1i,C2r,C2i,b0r,b0i,omm,time,C1,C2,b0,einuu,ii)

    CCTK_REAL, intent(in)  :: C1r,C1i,C2r,C2i,b0r,b0i,omm,time
    CCTK_COMPLEX, intent(out)  :: C1,C2,b0,einuu,ii
      C1=dcmplx(C1r,C1i)
      C2=dcmplx(C2r,C2i)
      b0=dcmplx(b0r,b0i)
      einuu=dcmplx(cos(omm*time),sin(omm*time))
      ii=dcmplx(0.0d0,1.0d0)

  end subroutine complexify

  function x_to_r_inverse(x) result(i_r)

    use NullGrid_Vars
    CCTK_REAL, intent(in)  :: x
    CCTK_REAL :: i_r

      i_r = (1-x)/(x*rwt)

  end function x_to_r_inverse


  !function x_to_r(x) result(i_r)
  !  use NullGrid_Vars
  !  CCTK_REAL, intent(in)  :: x
  !  CCTK_REAL :: i_r
  !
  !    i_r = (x*rwt)/(1-x)
  !end function


  function jout(s,l_in_Ylm,C1r,C1i,C2r,C2i,b0r,b0i,omm,time) result(J)

    CCTK_INT,  intent(in)  :: l_in_Ylm
    CCTK_REAL, intent(in)  :: s,C1r,C1i,C2r,C2i,b0r,b0i,omm,time
    CCTK_REAL :: J, inv_r
    CCTK_COMPLEX  :: C1,C2,b0,einuu,ii

      call complexify(C1r,C1i,C2r,C2i,b0r,b0i,omm,time,C1,C2,b0,einuu,ii)
      inv_r = x_to_r_inverse(s)

      select case (l_in_Ylm)
      case (2)
        J = dble(einuu*(24*b0+3*ii*omm*C1-omm**3*C2*ii-3*C2*inv_r**3+9*C1*inv_r)/36)
      case (3)
        J = dble(-ii*einuu*(-3*omm*C1+omm**4*C2*ii+60*ii*b0+18*ii*C1*inv_r-45*ii*C2*inv_r**4+30*omm*C2*inv_r**3)/180)
      case default
       call CCTK_WARN(1, "unencoded value for l_in_Ylm")
      end select

  end function jout

  function bout(s,l_in_Ylm,C1r,C1i,C2r,C2i,b0r,b0i,omm,time) result(B)

    CCTK_INT,  intent(in)  :: l_in_Ylm
    CCTK_REAL, intent(in)  :: s,C1r,C1i,C2r,C2i,b0r,b0i,omm,time
    CCTK_REAL :: B
    CCTK_COMPLEX  :: C1,C2,b0,einuu,ii

      call complexify(C1r,C1i,C2r,C2i,b0r,b0i,omm,time,C1,C2,b0,einuu,ii)
      B=dble(einuu*b0)

  end function bout

  function uout(s,l_in_Ylm,C1r,C1i,C2r,C2i,b0r,b0i,omm,time) result(U)

    use NullGrid_Vars
    CCTK_INT,  intent(in)  :: l_in_Ylm
    CCTK_REAL, intent(in)  :: s,C1r,C1i,C2r,C2i,b0r,b0i,omm,time
    CCTK_REAL :: U, inv_r
    CCTK_COMPLEX  :: C1,C2,b0,einuu,ii
    CCTK_INT :: i1

      call complexify(C1r,C1i,C2r,C2i,b0r,b0i,omm,time,C1,C2,b0,einuu,ii)
      inv_r = x_to_r_inverse(s)

!The following is a fix to ensure that if we have to set U at scri in the news routine
!it will be correct. Remember that U uses the half grid and U at scri is evaluated as
!(U(x=1+dx/2) + U(x-1-dx/2))/2
      i1 = 1
      If (inv_r<0) THEN
        inv_r = x_to_r_inverse(s-dx)
        i1 = -1
      end if

      select case (l_in_Ylm)
      case (2)
        U = dble(-einuu*( i1*(-18*C1*inv_r**2-12*ii*omm*C2*inv_r**3-9*C2*inv_r**4-72*b0*inv_r) &
                 +24*ii*omm*b0-3*omm**2*C1+omm**4*C2)/36)
        !U = dble(-einuu*( i1*(-18*C1*inv_r**2-12*ii*omm*C2*inv_r**3-9*C2*inv_r**4-72*b0*inv_r) &
        !         +24*ii*omm*b0-3*omm**2*C1+omm**4*C2)/36)
      case (3)
        U = dble(-einuu*( i1*(-180*C2*inv_r**5-225*ii*omm*C2*inv_r**4+120*omm**2*C2*inv_r**3 &
                -90*C1*inv_r**2-360*b0*inv_r) &
                 -3*omm**2*C1+omm**5*C2*ii+60*ii*omm*b0)/180)
      case default
       call CCTK_WARN(1, "unencoded value for l_in_Ylm")
      end select

  end function uout

  function druout(s,l_in_Ylm,C1r,C1i,C2r,C2i,b0r,b0i,omm,time) result(drU)

    use NullGrid_Vars
    CCTK_INT,  intent(in)  :: l_in_Ylm
    CCTK_REAL, intent(in)  :: s,C1r,C1i,C2r,C2i,b0r,b0i,omm,time
    CCTK_REAL :: drU, inv_r
    CCTK_COMPLEX  :: C1,C2,b0,einuu,ii
    CCTK_INT :: i1

      call complexify(C1r,C1i,C2r,C2i,b0r,b0i,omm,time,C1,C2,b0,einuu,ii)
      inv_r = x_to_r_inverse(s)

      i1 = 1
      If (inv_r<0) THEN
        inv_r = x_to_r_inverse(s-dx)
        i1 = -1
      end if

      select case (l_in_Ylm)
      case (2)
        drU = dble(-einuu*(i1*(2*b0*inv_r**2+C1*inv_r**3+ii*omm*C2*inv_r**4+C2*inv_r**5)))
      case (3)
        drU = dble(-einuu*(i1*(2*b0*inv_r**2+C1*inv_r**3-2*omm**2*C2*inv_r**4&
                              +5*ii*omm*C2*inv_r**5+5*C2*inv_r**6)))
      case default
       call CCTK_WARN(1, "unencoded value for l_in_Ylm")
      end select

  end function druout

  function wout(s,l_in_Ylm,C1r,C1i,C2r,C2i,b0r,b0i,omm,time) result(W)

    CCTK_INT,  intent(in)  :: l_in_Ylm
    CCTK_REAL, intent(in)  :: s,C1r,C1i,C2r,C2i,b0r,b0i,omm,time
    CCTK_REAL :: W, inv_r!, my_r
    CCTK_COMPLEX  :: C1,C2,b0,einuu,ii

      call complexify(C1r,C1i,C2r,C2i,b0r,b0i,omm,time,C1,C2,b0,einuu,ii)
      inv_r = x_to_r_inverse(s)
      !my_r = x_to_r(s)

      select case (l_in_Ylm)
      case (2)
        W = dble(einuu*(-12*b0*inv_r+3*C2*inv_r**4+6*ii*omm*C2*inv_r**3+24*ii*omm*b0 &
                  -3*omm**2*C1+omm**4*C2+6*ii*omm*C1*inv_r &
                  -2*ii*omm**3*C2*inv_r-6*omm**2*C2*inv_r**2)/6)
        !W = dble(einuu*(-12*b0*inv_r+3*C2*inv_r**4+6*ii*omm*C2*inv_r**3+24*ii*omm*b0 &
        !          -3*omm**2*C1+omm**4*C2+6*ii*omm*C1*inv_r &
        !          -2*ii*omm**3*C2*inv_r-6*omm**2*C2*inv_r**2)/6)
      case (3)
        W = dble(einuu*(-30*b0*inv_r+45*C2*inv_r**5+75*ii*omm*C2*inv_r**4-60*omm**2*C2*inv_r**3 &
                 +60*ii*omm*b0-3*omm**2*C1+ii*omm**5*C2 &
                 +15*ii*omm*C1*inv_r+5*omm**4*C2*inv_r-30*ii*omm**3*C2*inv_r**2)/15)
      case default
       call CCTK_WARN(1, "unencoded value for l_in_Ylm")
      end select

  end function wout


  function NullExact_Analytic_J_pt(patchID,q,p,s,time) result(J)

    CCTK_INT,  intent(in)  :: patchID
    CCTK_REAL, intent(in)  :: q, p, s, time
    CCTK_COMPLEX :: J
    DECLARE_CCTK_PARAMETERS

    ! local variables, as needed:

    CCTK_REAL :: rand(2)

       call random_number(rand)
       J = amplitude * dcmplx(rand(1)-0.5, rand(2)-0.5)

  end function NullExact_Analytic_J_pt

  function NullExact_Analytic_beta_pt(patchID,q,p,s,time) result(beta)

    CCTK_INT,  intent(in)  :: patchID
    CCTK_REAL, intent(in)  :: q, p, s, time
    CCTK_REAL :: beta
    DECLARE_CCTK_PARAMETERS

    ! local variables, as needed:

       call random_number(beta)
       beta = amplitude * (beta-0.5)

  end function NullExact_Analytic_beta_pt

  function NullExact_Analytic_U_pt(patchID,q,p,s,time) result(U)

    CCTK_INT,  intent(in)  :: patchID
    CCTK_REAL, intent(in)  :: q, p, s, time
    CCTK_COMPLEX :: U

    DECLARE_CCTK_PARAMETERS

    ! local variables, as needed:

    CCTK_REAL :: rand(2)

    select case (testbed_ID)
    case (0) ! minkowski
       U = 0
    case (1) ! random
       call random_number(rand)
       U = amplitude * dcmplx(rand(1)-0.5, rand(2)-0.5)
    case default
       U = 0
       call CCTK_WARN(1, "unrecognized value for testbed_ID")
    end select

  end function NullExact_Analytic_U_pt

  function NullExact_Analytic_W_pt(patchID,q,p,s,time) result(W)

    CCTK_INT,  intent(in)  :: patchID
    CCTK_REAL, intent(in)  :: q, p, s, time
    CCTK_REAL :: W

    DECLARE_CCTK_PARAMETERS

    ! local variables, as needed:

       call random_number(W)
       W = amplitude * (W-0.5)

  end function NullExact_Analytic_W_pt

  function NullExact_Analytic_nu_pt(patchID,q,p,s,time) result(nu)

    CCTK_INT,  intent(in)  :: patchID
    CCTK_REAL, intent(in)  :: q, p, s, time
    CCTK_COMPLEX :: nu

    DECLARE_CCTK_PARAMETERS

    ! local variables, as needed:

    CCTK_REAL :: rand(2)

       call random_number(rand)
       nu = amplitude * dcmplx(rand(1)-0.5, rand(2)-0.5)

  end function NullExact_Analytic_nu_pt

  function NullExact_Analytic_cb_pt(patchID,q,p,s,time) result(cb)

    CCTK_INT,  intent(in)  :: patchID
    CCTK_REAL, intent(in)  :: q, p, s, time
    CCTK_COMPLEX :: cb

    DECLARE_CCTK_PARAMETERS

    ! local variables, as needed:

    CCTK_REAL :: rand(2)

       call random_number(rand)
       cb = amplitude * dcmplx(rand(1)-0.5, rand(2)-0.5)

  end function NullExact_Analytic_cb_pt

  function NullExact_Analytic_ck_pt(patchID,q,p,s,time) result(ck)

    CCTK_INT,  intent(in)  :: patchID
    CCTK_REAL, intent(in)  :: q, p, s, time
    CCTK_COMPLEX :: ck

    DECLARE_CCTK_PARAMETERS

    ! local variables, as needed:

    CCTK_REAL :: rand(2)

       call random_number(rand)
       ck = amplitude * dcmplx(rand(1)-0.5, rand(2)-0.5)

  end function NullExact_Analytic_ck_pt


  !---------------------------------------------------------------------
  ! 2D interface routines to the pointwise defined testbed
  !---------------------------------------------------------------------

  subroutine NullExact_Analytic_J_2D(patchID,n1,n2,qs,ps,s,time,J,Ylm_2)

    CCTK_INT,  intent(in)  :: n1, n2,patchID
    CCTK_REAL, intent(in)  :: qs(n1,n2), ps(n1,n2), s, time
    CCTK_COMPLEX, intent(in) :: Ylm_2(n1,n2,2)
    CCTK_COMPLEX, intent(out) :: J(n1,n2)
    CCTK_INT :: i1, i2
    CCTK_REAL :: j1
    DECLARE_CCTK_PARAMETERS
    select case (testbed_ID)
    case (0) ! minkowski
       J = 0
    case (1) ! random

    do i2 = 1, n2
       do i1 = 1, n1
          J(i1,i2) = NullExact_Analytic_J_pt(patchID,qs(i1,i2),ps(i1,i2),s,time)
       end do
    end do

    case (2) ! linearized solution
       j1 = jout(s,l_in_Ylm,Constant_C1r,Constant_C1i,Constant_C2r,Constant_C2i &
                ,beta0r,beta0i,omm,time)
       J = j1*Ylm_2(:,:,patchID)*eth2fac(l_in_Ylm)
    case (3) ! eth, eth_eth, spin_trafo testbed
       J = Ylm_2(:,:,patchID)
    case default
       J = 0
       call CCTK_WARN(1, "unrecognized value for testbed_ID")
    end select

  end subroutine NullExact_Analytic_J_2D


  subroutine NullExact_Analytic_beta_2D(patchID,n1,n2,qs,ps,s,time,beta,Ylm)

    CCTK_INT,  intent(in)  :: n1, n2, patchID
    CCTK_REAL, intent(in)  :: qs(n1,n2), ps(n1,n2), s, time
    CCTK_COMPLEX, intent(in) :: Ylm(n1,n2,2)
    CCTK_REAL, intent(out) :: beta(n1,n2)
    CCTK_REAL :: b1
    CCTK_INT :: i1, i2
    DECLARE_CCTK_PARAMETERS

    select case (testbed_ID)
    case (0) ! minkowski
       beta = 0
    case (1) ! random
    do i2 = 1, n2
       do i1 = 1, n1
          beta(i1,i2) = NullExact_Analytic_beta_pt(patchID,qs(i1,i2),ps(i1,i2),s,time)
       end do
    end do
    case (2) ! linearized solution
       b1 = bout(s,l_in_Ylm,Constant_C1r,Constant_C1i,Constant_C2r,Constant_C2i &
                ,beta0r,beta0i,omm,time)
       beta = b1*Ylm(:,:,patchID)
    case (3) ! eth, eth_eth, spin_trafo testbed
       beta = Ylm(:,:,patchID)
    case default
       beta = 0
       call CCTK_WARN(1, "unrecognized value for testbed_ID")
    end select

  end subroutine NullExact_Analytic_beta_2D

  subroutine NullExact_Analytic_U_2D(patchID,n1,n2,qs,ps,s,time,U,Ylm_1)

    CCTK_INT,  intent(in)  :: n1, n2, patchID
    CCTK_REAL, intent(in)  :: qs(n1,n2), ps(n1,n2), s, time
    CCTK_COMPLEX, intent(in) :: Ylm_1(n1,n2,2)
    CCTK_COMPLEX, intent(out) :: U(n1,n2)
    CCTK_REAL :: u1
    CCTK_INT :: i1, i2
    DECLARE_CCTK_PARAMETERS

    select case (testbed_ID)
    case (0) ! minkowski
       U = 0
    case (1) ! random
    do i2 = 1, n2
       do i1 = 1, n1
          U(i1,i2) = NullExact_Analytic_U_pt(patchID,qs(i1,i2),ps(i1,i2),s,time)
       end do
    end do
    case (2) ! linearized solution
       u1 = uout(s,l_in_Ylm,Constant_C1r,Constant_C1i,Constant_C2r,Constant_C2i &
                ,beta0r,beta0i,omm,time)
       U = u1*Ylm_1(:,:,patchID)*eth1fac(l_in_Ylm)
    case (3) ! eth, eth_eth, spin_trafo testbed
       U = Ylm_1(:,:,patchID)
    case default
       U = 0
       call CCTK_WARN(1, "unrecognized value for testbed_ID")
    end select

  end subroutine NullExact_Analytic_U_2D

  subroutine NullExact_Analytic_drU_2D(patchID,n1,n2,qs,ps,s,time,drU,Ylm_1)

    CCTK_INT,  intent(in)  :: n1, n2, patchID
    CCTK_REAL, intent(in)  :: qs(n1,n2), ps(n1,n2), s, time
    CCTK_COMPLEX, intent(in) :: Ylm_1(n1,n2,2)
    CCTK_COMPLEX, intent(out) :: drU(n1,n2)
    CCTK_REAL :: dru1
    CCTK_INT :: i1, i2
    DECLARE_CCTK_PARAMETERS

    select case (testbed_ID)
    case (0) ! minkowski
       drU = 0
    case (1) ! random
    do i2 = 1, n2
       do i1 = 1, n1
          drU(i1,i2) = NullExact_Analytic_U_pt(patchID,qs(i1,i2),ps(i1,i2),s,time)
       end do
    end do
    case (2) ! linearized solution
       dru1 = druout(s,l_in_Ylm,Constant_C1r,Constant_C1i,Constant_C2r,Constant_C2i &
                ,beta0r,beta0i,omm,time)
       drU = dru1*Ylm_1(:,:,patchID)*eth1fac(l_in_Ylm)
    case (3) ! eth, eth_eth, spin_trafo testbed
       drU = Ylm_1(:,:,patchID)
    case default
       drU = 0
       call CCTK_WARN(1, "unrecognized value for testbed_ID")
    end select

  end subroutine NullExact_Analytic_drU_2D

  subroutine NullExact_Analytic_W_2D(patchID,n1,n2,qs,ps,s,time,W,Ylm)

    CCTK_INT,  intent(in)  :: n1, n2, patchID
    CCTK_REAL, intent(in)  :: qs(n1,n2), ps(n1,n2), s, time
    CCTK_COMPLEX, intent(in) :: Ylm(n1,n2,2)
    CCTK_REAL, intent(out) :: W(n1,n2)
    CCTK_REAL :: w1
    CCTK_INT :: i1, i2
    DECLARE_CCTK_PARAMETERS

    select case (testbed_ID)
    case (0) ! minkowski
       W = 0
    case (1) ! random
    do i2 = 1, n2
       do i1 = 1, n1
          W(i1,i2) = NullExact_Analytic_W_pt(patchID,qs(i1,i2),ps(i1,i2),s,time)
       end do
    end do
    case (2) ! linearized solution
       w1 = wout(s,l_in_Ylm,Constant_C1r,Constant_C1i,Constant_C2r,Constant_C2i &
                ,beta0r,beta0i,omm,time)
       W = w1*Ylm(:,:,patchID)
    case (3) ! eth, eth_eth, spin_trafo testbed
       W = Ylm(:,:,patchID)
    case default
       W = 0
       call CCTK_WARN(1, "unrecognized value for testbed_ID")
    end select

  end subroutine NullExact_Analytic_W_2D

  subroutine NullExact_Analytic_nu_2D(patchID,n1,n2,qs,ps,s,time,nu,Ylm_1)

    CCTK_INT,  intent(in)  :: n1, n2, patchID
    CCTK_REAL, intent(in)  :: qs(n1,n2), ps(n1,n2), s, time
    CCTK_COMPLEX, intent(in) :: Ylm_1(n1,n2,2)
    CCTK_COMPLEX, intent(out) :: nu(n1,n2)
    CCTK_REAL :: j1
    CCTK_INT :: i1, i2
    DECLARE_CCTK_PARAMETERS

    select case (testbed_ID)
    case (0) ! minkowski
       nu = 0
    case (1) ! random
    do i2 = 1, n2
       do i1 = 1, n1
          nu(i1,i2) = NullExact_Analytic_nu_pt(patchID,qs(i1,i2),ps(i1,i2),s,time)
       end do
    end do
    case (2) ! linearized solution
       j1 = jout(s,l_in_Ylm,Constant_C1r,Constant_C1i,Constant_C2r,Constant_C2i &
                ,beta0r,beta0i,omm,time)
       nu = j1*Ylm_1(:,:,patchID)*eth1fac(l_in_Ylm)*(2-l_in_Ylm*(l_in_Ylm+1))
    case (3) ! eth, eth_eth, spin_trafo testbed
       nu = Ylm_1(:,:,patchID)
    case default
       nu = 0
       call CCTK_WARN(1, "unrecognized value for testbed_ID")
    end select

  end subroutine NullExact_Analytic_nu_2D

  subroutine NullExact_Analytic_cb_2D(patchID,n1,n2,qs,ps,s,time,cb,Ylm_1)

    CCTK_INT,  intent(in)  :: n1, n2, patchID
    CCTK_REAL, intent(in)  :: qs(n1,n2), ps(n1,n2), s, time
    CCTK_COMPLEX, intent(in) :: Ylm_1(n1,n2,2)
    CCTK_COMPLEX, intent(out) :: cb(n1,n2)
    CCTK_REAL :: b1
    CCTK_INT :: i1, i2
    DECLARE_CCTK_PARAMETERS

    select case (testbed_ID)
    case (0) ! minkowski
       cb = 0
    case (1) ! random
    do i2 = 1, n2
       do i1 = 1, n1
          cb(i1,i2) = NullExact_Analytic_cb_pt(patchID,qs(i1,i2),ps(i1,i2),s,time)
       end do
    end do
    case (2) ! linearized solution
       b1 = bout(s,l_in_Ylm,Constant_C1r,Constant_C1i,Constant_C2r,Constant_C2i &
                ,beta0r,beta0i,omm,time)
       cb = b1*Ylm_1(:,:,patchID)*eth1fac(l_in_Ylm)
    case (3) ! eth, eth_eth, spin_trafo testbed
       cb = Ylm_1(:,:,patchID)
    case default
       cb = 0
       call CCTK_WARN(1, "unrecognized value for testbed_ID")
    end select


  end subroutine NullExact_Analytic_cb_2D

  subroutine NullExact_Analytic_ck_2D(patchID,n1,n2,qs,ps,s,time,ck,Ylm_1)

    CCTK_INT,  intent(in)  :: n1, n2, patchID
    CCTK_REAL, intent(in)  :: qs(n1,n2), ps(n1,n2), s, time
    CCTK_COMPLEX, intent(in) :: Ylm_1(n1,n2,2)
    CCTK_COMPLEX, intent(out) :: ck(n1,n2)

    CCTK_INT :: i1, i2
    DECLARE_CCTK_PARAMETERS

    select case (testbed_ID)
    case (0) ! minkowski
       ck = 0
    case (1) ! random
    do i2 = 1, n2
       do i1 = 1, n1
          ck(i1,i2) = NullExact_Analytic_ck_pt(patchID,qs(i1,i2),ps(i1,i2),s,time)
       end do
    end do
    case (2) ! linearized solution
       ck = 0
    case (3) ! eth, eth_eth, spin_trafo testbed
       ck = Ylm_1(:,:,patchID)
    case default
       ck = 0
       call CCTK_WARN(1, "unrecognized value for testbed_ID")
    end select


  end subroutine NullExact_Analytic_ck_2D


!  subroutine NullExact_Analytic_mu_2D(patchID,n1,n2,qs,ps,s,time,mu,Ylm_0)
!
!    CCTK_INT,  intent(in)  :: n1, n2,patchID
!    CCTK_REAL, intent(in)  :: qs(n1,n2), ps(n1,n2), s, time
!    CCTK_COMPLEX, intent(in) :: Ylm_0(n1,n2,2)
!    CCTK_COMPLEX, intent(out) :: mu(n1,n2)
!    CCTK_REAL :: fact
!
!    DECLARE_CCTK_PARAMETERS
!    select case (testbed_ID)
!    case (0) ! minkowski
!       mu = 0
!    case (1) ! random
!       mu = 0
!    case (2) ! linearized solution
!       fact = 1/375000*sin(t)
!       mu = Ylm_0(:,:,patchID)*fact
!    case (3) ! eth, eth_eth, spin_trafo testbed
!       mu = Ylm_0(:,:,patchID)
!    case default
!       mu = 0
!       call CCTK_WARN(1, "unrecognized value for testbed_ID")
!    end select
!
!  end subroutine NullExact_Analytic_mu_2D


  subroutine NullExact_Analytic_BondiNews_2D(patchID,n1,n2,qs,ps,s,time,News,Psi4e,Ylm_2)

    CCTK_INT,  intent(in)  :: n1, n2, patchID
    CCTK_REAL, intent(in)  :: qs(n1,n2), ps(n1,n2), s, time
    CCTK_COMPLEX, intent(in) :: Ylm_2(n1,n2,2)
    CCTK_COMPLEX, intent(out) :: News(n1,n2),Psi4e(n1,n2)
    CCTK_REAL :: nn,npsi
    CCTK_COMPLEX  :: C1,C2,b0,einuu,ii

    DECLARE_CCTK_PARAMETERS

    select case (testbed_ID)
    case (0) ! minkowski
       News = 0
    case (1) ! random
       News = 1.e+10
    case (2) ! linearized solution
      call complexify(Constant_C1r,Constant_C1i,Constant_C2r,Constant_C2i,beta0r,beta0i,omm,time,C1,C2,b0,einuu,ii)
      select case (l_in_Ylm)
        case (2)
          nn = dble(ii*einuu*omm**3*C2/24)
          npsi = dble(-einuu*omm**4*C2/24)
        case (3)
          nn = dble(-einuu*omm**4*C2/60)
          npsi = dble(-ii*einuu*omm**5*C2/60)
        case default
          call CCTK_WARN(1, "unencoded vale for l_in_Ylm")
          nn = 0
      end select
      News = nn*Ylm_2(:,:,patchID)*eth2fac(l_in_Ylm)
      Psi4e = npsi*Ylm_2(:,:,patchID)*eth2fac(l_in_Ylm)
    case (3) ! eth, eth_eth, spin_trafo testbed
       News = Ylm_2(:,:,patchID)
    case default
       News = 0
       Psi4e = 0
       call CCTK_WARN(1, "unrecognized value for testbed_ID")
    end select

  end subroutine NullExact_Analytic_BondiNews_2D

  subroutine NullExact_Analytic_BondiTime_2D(patchID,n1,n2,qs,ps,s,time,BondiTime,Ylm_0)

    CCTK_INT,  intent(in)  :: n1, n2, patchID
    CCTK_REAL, intent(in)  :: qs(n1,n2), ps(n1,n2), s, time
    CCTK_COMPLEX, intent(in) :: Ylm_0(n1,n2,2)
    CCTK_REAL,    intent(out) :: BondiTime(n1,n2)

    DECLARE_CCTK_PARAMETERS

!    call CCTK_WARN(1, "this is to be implemented")
    BondiTime = 1.e+10

  end subroutine NullExact_Analytic_BondiTime_2D


end module NullExact_Analytic
