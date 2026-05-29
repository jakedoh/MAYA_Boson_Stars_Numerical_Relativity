/*
  vim: syntax=fortran
  module dummy
  end module dummy
*/
#include "cctk.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Module SpinWeightedYlm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Purpose: Provide the spin-weighted spherical harmonics
!!        on the north/south stereographic patchs
!!
!! Conventions:  sYlm = Sqrt[(l-s)!/(l+s)!] \eth^s Ylm (for s>0)
!!                    = (-1)^s  Sqrt[(l+s)!/(l-s)!] \bar\eth^(-s) Ylm (for s<0)
!!
!!               \eth f = PP/2 (\partial_q f + I \partial_p f) for spin zero f
!!                    where PP = 1 + q^2 + PP^2, and 
!!                Zn = q+ip = Tan(theta/2) e^{I phi} in the north patch, and
!!                Zs =q+ip = Cot(theta/2) e^{-I phi} in the south patch.
!!
!!                    Note: our Zs is the standard \bar\zeta used in 
!!                       e.g. Goldberg et al, J. Math. Phys. 8, 2155, 1967
!!                           Stewart 1990 "Advanced General Relativity"
!!                    Also Note: Our \eth corresponds to Stewart's \beth
!!                           Our s to Stewart's -s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! 12/05/02  -- created
! status :  Confirmed that 2Y20 is spin2, Confirmed that eth_eth 0Ylm 
!   yeilds Sqrt[(l+s)!/(l-s)!)] 2Ylm, for s=2, l=2, 3, and 4 (with all m)
!   also verfied, via the nullcodes interpolators, that the above 2Ylm are
!   spin-weight 2. 
! - still need to show that the 0Ylm are the spherical harmonics.
! - make sure there are not gotchas with negative spins etc.
! 
! 12/06/02 -- tests
!   confirmed that I get the correct Y98, Y9-8, Y97, Y9-7, Y90,
!                                    Y87, Y8-7, Y86, Y8-6, Y80 
!     so it seems to work for both even and odd l and m and neg/pos m
!
!   confirmed that I get the correct 3Ylm, -3Ylm, 2Ylm, -2Ylm
!     for the above l and m combos. So it seems to work for even
!     and odd s with both positve and negative spin
!   Note the sYl0 seem to have an error ~10x of the others

module NullDecomp_sYlm
  implicit none
  private
  CCTK_REAL :: dummy_var
  integer, parameter :: wp = kind(dummy_var)
  CCTK_REAL, parameter :: pi = 3.1415926535897932384626433832795_wp
  logical, save :: FirstTimeCalled  = .true.
  CCTK_REAL, dimension(0:20), save :: fact_
  CCTK_COMPLEX, dimension(:,:,:), Allocatable, save :: zp
  CCTK_COMPLEX, dimension(:,:), allocatable,  save :: tm
  CCTK_REAL, dimension(:,:,:),  allocatable, save :: PP, P_d

  public rsYlm
  public sYlm
  public sYlm_north_opt
  public sYlm_north
  public sYlm_south_opt
  public sYlm_south

contains

  subroutine rsYlm(s, l, m, nq, np, zeta, out)
    implicit none
    CCTK_INT, intent(in) :: s, l, m, nq, np
    CCTK_COMPLEX, dimension(nq,np), intent(in) :: zeta
    CCTK_COMPLEX, dimension(nq,np,2), intent(out) :: out

    CCTK_COMPLEX, dimension(nq,np,2) :: tmp
    CCTK_COMPLEX :: ii = (0, 1.0)

    select case (m)
    case(:-1)
      call sYlm(s,l,-m,nq,np,zeta,out)
      call sYlm(s,l,m,nq,np,zeta,tmp)
      out = -ii/ dsqrt(2.0d0)*(out -(-1)**m * tmp)
    case(0)
      call sYlm(s,l,m,nq,np,zeta,out)
    case(1:)
      call sYlm(s,l,m,nq,np,zeta,out)
      call sYlm(s,l,-m,nq,np,zeta,tmp)
      out = 1.0d0/ dsqrt(2.0d0)*(out +(-1)**m * tmp)
    case DEFAULT
      call CCTK_WARN(0,"bad m")
    end select
    call sYlm(s,l,m,nq,np,zeta,tmp)
  end subroutine rsYlm

  subroutine sYlm(s, l, m, nq,np, zeta, out)
    implicit none
    CCTK_INT, intent(in) :: s, l, m, nq,np
    CCTK_COMPLEX, dimension(nq,np), intent(in) :: zeta
    CCTK_COMPLEX, dimension(nq,np,2), intent(out) :: out

    CCTK_INT :: ss, mm, n

    if (l < 0 ) then
       stop "sYlm ERROR: l < 0 "
    endif
    if (abs(m) > l)  then
      stop "sYlm ERROR: abs(m) > l"
    endif
    if (abs(s) > l)  then
      stop "sYlm ERROR: abs(s) > l"
    endif

    if (FirstTimeCalled) then
      FirstTimeCalled = .false.
      fact_(0) = 1.0_wp
      do n = 1, 20
        fact_(n) = n*fact_(n-1)
      end do
      allocate(zp(nq,np,0:20), PP(nq,np,0:20), P_d(nq,np,0:20), tm(nq,np))

      zp(:,:,0) = 1.0_wp   
      PP(:,:,0)  = 1.0_wp
      P_d(:,:,0)  = 1.0_wp
      do n = 1, 20
        zp(:,:,n) = zeta(:,:)**n
        PP(:,:,n) = (1.0_wp + abs(zeta(:,:)*conjg(zeta(:,:))))**n
        P_d(:,:,n) = 1.0_wp / (1.0_wp + abs(zeta(:,:)*conjg(zeta(:,:))))**n
      end do
    endif 

    ss = s
    mm = m
    if (ss< 0) ss= -ss
    if (mm< 0) mm= -mm

    if ( (l + ss) .le. 20 .and. (l+mm) .le. 20 ) then

      call sYlm_north_opt(s, l, m, nq,np,out(:,:,1))
      call sYlm_south_opt(s, l, m, nq,np,out(:,:,2))
    else
      !write(*,*) "using unoptimized routines"
      call sYlm_north(s, l, m, nq,np, zeta(:,:), out(:,:,1))
      call sYlm_south(s, l, m, nq,np, zeta(:,:), out(:,:,2))
    endif
  end subroutine sYlm

  subroutine sYlm_north_opt(s, l, m, nq,np, north)
    implicit none  
    CCTK_INT, intent(in) :: s, l, m, nq,np
    CCTK_COMPLEX, dimension(nq,np), intent(out) :: north
    CCTK_INT :: n

    north = 0.0

    tm = 0.0 
    do n  = max(0, s+m), min(l+s, l+m)
      tm = tm + (-1)**n * zp(:,:,n) * conjg(zp(:,:,n-m-s)) / &
              (fact_(n)*fact_(l+m-n)*fact_(l+s-n)*fact_(n-m-s))
    end do
    north = dsqrt(1.0_wp + 2.0_wp * l)*P_d(:,:,l)  * &
            dsqrt(fact_(l-m)*fact_(l+m)*fact_(l-s)*fact_(l+s)) * tm /&
            (2.0_wp * dsqrt(pi))
  end subroutine sYlm_north_opt


  subroutine sYlm_north(s, l, m, nq,np, zeta, north)
    implicit none  
    CCTK_INT, intent(in) :: s, l, m, nq,np
    CCTK_COMPLEX, dimension(nq,np), intent(in) :: zeta
    CCTK_COMPLEX, dimension(nq,np), intent(out) :: north

    CCTK_COMPLEX, dimension(nq,np) :: zetabar, tmp
    CCTK_INT :: n

    zetabar = conjg(zeta)

    if (abs(m) > l)  then
      stop "sYlm_north ERROR: abs(m) > l"
    endif
    if (abs(s) > l)  then
      stop "sYlm_north ERROR: abs(s) > l"
    endif

    north = 0.0

    tmp = 0.0 
    do n  = max(0, s+m), min(l+s, l+m)
      tmp = tmp + (-1)**n * zeta**n * zetabar**(n-m-s) / &
              (fact(n)*fact(l+m-n)*fact(l+s-n)*fact(n-m-s))
    end do
    north = dsqrt(1.0_wp + 2.0_wp * l) / (  &
                      (1.0_wp + zeta * zetabar )**l ) * &
            dsqrt(fact(l-m)*fact(l+m)*fact(l-s)*fact(l+s)) * tmp /&
            (2.0_wp * dsqrt(pi))
  end subroutine sYlm_north 

  subroutine sYlm_south_opt(s, l, m, nq,np, south)
    implicit none  
    CCTK_INT, intent(in) :: s, l, m, nq,np
    CCTK_COMPLEX, dimension(nq,np), intent(out) :: south
    CCTK_INT :: n

    south = 0.0

    tm = 0.0 
    do n  = max(0, s+m), min(l+s, l+m)
      tm = tm  + (-1)**n * zp(:,:,l+s-n) * conjg(zp(:,:,l+m-n)) / &
              (fact_(n)*fact_(l+m-n)*fact_(l+s-n)*fact_(n-m-s))
    end do

    south = (-1)**s * dsqrt(1.0_wp + 2.0_wp * l) * P_d(:,:,l) * &
            dsqrt(fact_(l-m)*fact_(l+m)*fact_(l-s)*fact_(l+s)) * tm /&
            (2.0_wp * dsqrt(pi))


  end subroutine sYlm_south_opt


  subroutine sYlm_south(s, l, m, nq,np, zeta, south)
    implicit none  
    CCTK_INT, intent(in) :: s, l, m, nq,np
    CCTK_COMPLEX, dimension(nq,np), intent(in) :: zeta
    CCTK_COMPLEX, dimension(nq,np), intent(out) :: south

    CCTK_COMPLEX, dimension(nq,np) :: zetabar, tmp
    CCTK_INT :: n

    zetabar = conjg(zeta)

    if (abs(m) > l)  then
      stop "sYlm_south ERROR: abs(m) > l"
    endif
    if (abs(s) > l)  then
      stop "sYlm_south ERROR: abs(s) > l"
    endif

    south = 0.0

    tmp = 0.0 
    do n  = max(0, s+m), min(l+s, l+m)
      tmp = tmp  + (-1)**n * zeta**(l+s-n) * zetabar**(l+m-n) / &
              (fact(n)*fact(l+m-n)*fact(l+s-n)*fact(n-m-s))
    end do

    south = (-1)**s * dsqrt(1.0_wp + 2.0_wp * l) / &
                      (1.0_wp + zeta * zetabar )**l * &
            dsqrt(fact(l-m)*fact(l+m)*fact(l-s)*fact(l+s)) * tmp /&
            (2.0_wp * dsqrt(pi))

  end subroutine sYlm_south 

  CCTK_REAL function fact(j) 
    implicit none
    CCTK_INT, intent(in) :: j

    CCTK_INT :: i
    CCTK_REAL :: res

    if (j < 0 ) then 
      call CCTK_WARN(0,"fact called with a negative argument")
    endif
    if (j > 100 ) then
      call CCTK_WARN(0,"you gotta be kidding")
    endif
    if (j == 0 .OR. j == 1) then
      fact = 1.0_wp
    else
      res = 1.0_wp
      do i = 2, j
        res = res * i
      end do
      fact = res
    end if
    return
  end function fact 
 

end module NullDecomp_sYlm
