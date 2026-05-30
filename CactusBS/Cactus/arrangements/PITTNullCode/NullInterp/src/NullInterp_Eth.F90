! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

module NullInterp_Eth
  use cctk
implicit none

contains

subroutine NullInterp_eth1(cctkGH, tmp_cgfn, tmp_cgfs, F_out, F_in, spin, e1)
  use NullGrid_Vars
  use NullInterp_Interp

  implicit none
  CCTK_POINTER, intent(in) :: cctkGH
  CCTK_INT, intent(in) :: spin, e1
  CCTK_COMPLEX, dimension(lsh(1),lsh(2)), intent(inout) :: tmp_cgfn, tmp_cgfs
  CCTK_COMPLEX, dimension(:,:,:), intent(inout) :: F_out
  CCTK_COMPLEX, dimension(:,:,:), intent(in) :: F_in

  integer i, j

  DECLARE_CCTK_PARAMETERS

  if (poison_test.ne.0) then

     tmp_cgfn = F_in(:,:,1)
     tmp_cgfs = F_in(:,:,2)
     do j = 1, lsh(2)
        do i = 1, lsh(1)
           if (EG(i,j).eq.0) then
              tmp_cgfn(i,j) = poison_value * (1.0, 1.0)
              tmp_cgfs(i,j) = poison_value * (1.0, 1.0)
           end if
        end do
     end do

     call NullInterp_d1(F_out(:,:,1), tmp_cgfn, spin, e1) 
     call NullInterp_d1(F_out(:,:,2), tmp_cgfs, spin, e1) 
    
  else

     call NullInterp_d1(F_out(:,:,1), F_in(:,:,1), spin, e1) 
     call NullInterp_d1(F_out(:,:,2), F_in(:,:,2), spin, e1) 

  end if

  call NullInterp_cinterp(cctkGH, tmp_cgfn, tmp_cgfs, F_out(:,:,1), F_out(:,:,2), spin+e1)

end subroutine NullInterp_eth1

subroutine NullInterp_eth2(cctkGH, tmp_cgfn, tmp_cgfs, F_out, F_in, spin, e1, e2)

  use NullGrid_Vars
  use NullInterp_Interp

  implicit none
  CCTK_POINTER, intent(in) :: cctkGH
  CCTK_INT, intent(in) :: spin, e1, e2
  CCTK_COMPLEX, dimension(lsh(1),lsh(2)), intent(inout) :: tmp_cgfn, tmp_cgfs
  CCTK_COMPLEX, dimension(:,:,:), intent(inout) :: F_out
  CCTK_COMPLEX, dimension(:,:,:), intent(in) :: F_in

  integer i, j

  DECLARE_CCTK_PARAMETERS

  if (poison_test.ne.0) then

     tmp_cgfn = F_in(:,:,1)
     tmp_cgfs = F_in(:,:,2)
     do j = 1, lsh(2)
        do i = 1, lsh(1)
           if (EG(i,j).eq.0) then
              tmp_cgfn(i,j) = poison_value * (1.0, 1.0)
              tmp_cgfs(i,j) = poison_value * (1.0, 1.0)
           end if
        end do
     end do

     call NullInterp_d2(F_out(:,:,1), tmp_cgfn, spin, e1, e2) 
     call NullInterp_d2(F_out(:,:,2), tmp_cgfs, spin, e1, e2) 
    
  else

     call NullInterp_d2(F_out(:,:,1), F_in(:,:,1), spin, e1, e2) 
     call NullInterp_d2(F_out(:,:,2), F_in(:,:,2), spin, e1, e2) 
  
  end if

  call NullInterp_cinterp(cctkGH, tmp_cgfn, tmp_cgfs, F_out(:,:,1), F_out(:,:,2), spin+e1+e2)

end subroutine NullInterp_eth2

subroutine NullInterp_eth3(cctkGH, tmp_cgfn, tmp_cgfs, F_out, F_in, spin, e1, e2, e3)

  use NullGrid_Vars
  use NullInterp_Interp

  implicit none
  CCTK_POINTER, intent(in) :: cctkGH
  CCTK_INT, intent(in) :: spin, e1, e2, e3
  CCTK_COMPLEX, dimension(lsh(1),lsh(2)), intent(inout) :: tmp_cgfn, tmp_cgfs
  CCTK_COMPLEX, dimension(:,:,:), intent(inout) :: F_out
  CCTK_COMPLEX, dimension(:,:,:), intent(in) :: F_in

  integer i, j

  DECLARE_CCTK_PARAMETERS

  if (poison_test.ne.0) then

     tmp_cgfn = F_in(:,:,1)
     tmp_cgfs = F_in(:,:,2)
     do j = 1, lsh(2)
        do i = 1, lsh(1)
           if (EG(i,j).eq.0) then
              tmp_cgfn(i,j) = poison_value * (1.0, 1.0)
              tmp_cgfs(i,j) = poison_value * (1.0, 1.0)
           end if
        end do
     end do

     call NullInterp_d3(F_out(:,:,1), tmp_cgfn, spin, e1, e2, e3)
     call NullInterp_d3(F_out(:,:,2), tmp_cgfs, spin, e1, e2, e3)
    
  else

     call NullInterp_d3(F_out(:,:,1), F_in(:,:,1), spin, e1, e2, e3) 
     call NullInterp_d3(F_out(:,:,2), F_in(:,:,2), spin, e1, e2, e3) 
  
  end if

  call NullInterp_cinterp(cctkGH, tmp_cgfn, tmp_cgfs, F_out(:,:,1), F_out(:,:,2), spin+e1+e2+e3)

end subroutine NullInterp_eth3

! the following implements the 3rd eth operator, defined in Maple as follows:
!Eth3 := proc(f, s, e1, e2, e3)
!         return Eth1(Eth1(Eth1(f, s, e1), s+e1, e2), s+e1+e2, e3);
!end proc:

subroutine NullInterp_d3 (F_out, F_in, spin, e1, e2, e3)
use NullGrid_Vars
   implicit none

   CCTK_COMPLEX,   dimension (:,:), intent (out) :: F_out
   CCTK_COMPLEX,   dimension (:,:), intent (in)  :: F_in
   CCTK_INT,                        intent (in)  :: spin, e1,e2, e3
   CCTK_INT :: nq, np
   CCTK_REAL :: dq, dp
   CCTK_REAL :: h1q, h1p, h2q, h2p, h3q, h3p
   CCTK_COMPLEX :: ii

   logical, save :: FirstTime = .true.
   CCTK_COMPLEX, dimension (:,:), allocatable, save ::&
       p_x, p_y, p_xx, p_xy, p_yy, p_xxx, p_xxy, p_xyy, p_yyy

  DECLARE_CCTK_PARAMETERS

   nq = lsh(1)
   np = lsh(2)
   dq = dble(delta(1))
   dp = dble(delta(2))
   ii = dcmplx(0.,1.)

  if (FirstTime) then
     FirstTime = .false.
     allocate(p_x(nq,np), p_y(nq,np),&
              p_xx(nq,np), p_xy(nq,np),&
              p_yy(nq,np), p_xxx(nq,np),&
              p_xxy(nq,np), p_xyy(nq,np),&
              p_yyy(nq,np))

     p_x = 0; p_y = 0
     p_xx = 0; p_xy = 0; p_yy = 0
     p_xxx = 0; p_xxy = 0; p_xyy = 0; p_yyy = 0

  endif

if (deriv_accuracy.eq.2) then

   h1q = dble(1. / (2. * dq))
   h2q = dble(1. / (dq * dq))
   h3q = dble(1. / (dq * dq * dq))
   h1p = dble(1. / (2. * dp))
   h2p = dble(1. / (dp * dp))
   h3p = dble(1. / (dp * dp * dp))

   p_x(1,1:np) = (-3. * F_in(1,1:np) + 4. * F_in(2,1:np) - F_in(3,1:np)) * h1q
   p_x(2:nq-1,1:np) = (F_in(3:nq,1:np) - F_in(1:nq-2,1:np)) * h1q
   p_x(nq,1:np) = (3. * F_in(nq,1:np) - 4. * F_in(nq-1,1:np) + F_in(nq-2,1:np)) * h1q

   p_y(1:nq,1) = (-3. * F_in(1:nq,1) + 4. * F_in(1:nq,2) - F_in(1:nq,3)) * h1p
   p_y(1:nq,2:np-1) = (F_in(1:nq,3:np) - F_in(1:nq,1:np-2)) * h1p
   p_y(1:nq,np) = (3. * F_in(1:nq,np) - 4. * F_in(1:nq,np-1) + F_in(1:nq,np-2)) * h1p

   p_xx(2:nq-1,1:np) = (F_in(3:nq,1:np)- 2.*F_in(2:nq-1,1:np)+F_in(1:nq-2,1:np)) * h2q
   p_yy(1:nq,2:np-1) = (F_in(1:nq,3:np)- 2.*F_in(1:nq,2:np-1)+F_in(1:nq,1:np-2)) * h2p

   p_xy(2:nq-1,2:np-1) = ( F_in(3:nq,3:np)   - F_in(3:nq,1:np-2)&
                         - F_in(1:nq-2,3:np) + F_in(1:nq-2,1:np-2))&
                       / (4. * dq * dp)

   p_xy(1,2:np-1) = (-3. * F_in(1,3:np)   + 4. * F_in(2,3:np)   - F_in(3,3:np)&
                    + 3. * F_in(1,1:np-2) - 4. * F_in(2,1:np-2) +&
                    & F_in(3,1:np-2)) * h1q * h1p


   !! p_xxx = 0; p_xxy = 0; p_xyy = 0; p_yyy = 0

   p_xxy(2:nq-1,1:np) = (p_y(3:nq,1:np)- 2.*p_y(2:nq-1,1:np)+p_y(1:nq-2,1:np)) * h2q
   p_xyy(1:nq,2:np-1) = (p_x(1:nq,3:np)- 2.*p_x(1:nq,2:np-1)+p_x(1:nq,1:np-2)) * h2p

   !> simplify(subs(x=x1,dx^3*diff(pol5pt,x,x,x)));
   !                              5 f1                         3 f5
   !                            - ---- + 9 f2 - 12 f3 + 7 f4 - ----
   !                               2                            2

   p_xxx(1,:) = ( - 2.5 * F_in(1,:) + 9 * F_in(2,:) - 12 * F_in(3,:) + 7 * F_in(4,:) - 1.5 * F_in(5,:) ) * h3q
   p_yyy(:,1) = ( - 2.5 * F_in(:,1) + 9 * F_in(:,2) - 12 * F_in(:,3) + 7 * F_in(:,4) - 1.5 * F_in(:,5) ) * h3p

   !> simplify(subs(x=x2,dx^3*diff(pol5pt,x,x,x)));
   !                               3 f1                         f5
   !                             - ---- + 5 f2 - 6 f3 + 3 f4 - ----
   !                                2                           2
 
   p_xxx(2,:) = ( - 1.5 * F_in(1,:) + 5 * F_in(2,:) - 6 * F_in(3,:) + 3 * F_in(4,:) - 0.5 * F_in(5,:) ) * h3q
   p_yyy(:,2) = ( - 1.5 * F_in(:,1) + 5 * F_in(:,2) - 6 * F_in(:,3) + 3 * F_in(:,4) - 0.5 * F_in(:,5) ) * h3p

   !> simplify(subs(x=x3,dx^3*diff(pol5pt,x,x,x)));
   !                                     f1               f5
   !                                  - ---- + f2 - f4 + ----
   !                                     2                2
 
   p_xxx(3:nq-2,:) = ( - 0.5 * F_in(1:nq-4,:) + F_in(2:nq-3,:) - F_in(4:nq-1,:) + 0.5 * F_in(5:nq,:) ) * h3q
   p_yyy(:,3:nq-2) = ( - 0.5 * F_in(:,1:nq-4) + F_in(:,2:nq-3) - F_in(:,4:nq-1) + 0.5 * F_in(:,5:nq) ) * h3p

   !> simplify(subs(x=x4,dx^3*diff(pol5pt,x,x,x)));
   !                               f1                         3 f5
   !                              ---- - 3 f2 + 6 f3 - 5 f4 + ----
   !                               2                           2
 
   p_xxx(nq-1,:) = ( + 0.5 * F_in(nq-4,:) - 3 * F_in(nq-3,:) + 6 * F_in(nq-2,:) - 5 * F_in(nq-1,:) + 1.5 * F_in(nq,:) ) * h3q
   p_yyy(:,nq-1) = ( + 0.5 * F_in(:,nq-4) - 3 * F_in(:,nq-3) + 6 * F_in(:,nq-2) - 5 * F_in(:,nq-1) + 1.5 * F_in(:,nq) ) * h3p

   !> simplify(subs(x=x5,dx^3*diff(pol5pt,x,x,x)));
   !                             3 f1                         5 f5
   !                             ---- - 7 f2 + 12 f3 - 9 f4 + ----
   !                              2                            2

   p_xxx(nq,:) = ( + 1.5 * F_in(nq-4,:) - 7 * F_in(nq-3,:) + 12 * F_in(nq-2,:) - 9 * F_in(nq-1,:) + 2.5 * F_in(nq,:) ) * h3q
   p_yyy(:,nq) = ( + 1.5 * F_in(:,nq-4) - 7 * F_in(:,nq-3) + 12 * F_in(:,nq-2) - 9 * F_in(:,nq-1) + 2.5 * F_in(:,nq) ) * h3p

end if

if (deriv_accuracy.eq.4) then

   h1q = dble(1./(dq))
   h2q = dble(1./(dq*dq))
   h3q = dble(1./(dq*dq*dq))
   h1p = dble(1./(dp))
   h2p = dble(1./(dp*dp))
   h3p = dble(1./(dp*dp*dp))
!--------------------------D1---------------------------------------
!2nd order
 !  p_x(1,:) = -1./2.*(3.*F_in(1,:) - 4.*F_in(2,:) + F_in(3,:))*h1q
 !  p_x(2,:) = -1./2.*(F_in(1,:) - F_in(3,:))*h1q
 !  p_x(3:nq-2,:) = -1./2.*(F_in(2:nq-3,:) - F_in(4:nq-1,:))*h1q !2nd order centered
 !  p_x(nq-1,:) = -1./2.*(F_in(nq-2,:) - F_in(nq,:))*h1q
 !  p_x(nq,:) = 1./2.*(3.*F_in(nq,:) - 4.*F_in(nq-1,:) + F_in(nq-2,:))*h1q

 !  p_y(:,1) = -1./2.*(3.*F_in(:,1) - 4.*F_in(:,2) + F_in(:,3))*h1p
 !  p_y(:,2) = -1./2.*(F_in(:,1) - F_in(:,3))*h1p
 !  p_y(:,3:np-2) = -1./2.*(F_in(:,2:np-3) - F_in(:,4:np-1))*h1p !2nd order centered
 !  p_y(:,np-1) = -1./2.*(F_in(:,np-2) - F_in(:,np))*h1p
 !  p_y(:,np) = 1./2.*(3.*F_in(:,np) - 4.*F_in(:,np-1) + F_in(:,np-2))*h1p

!4th order
   p_x(1,:) = -1./2.*(3.*F_in(1,:) - 4.*F_in(2,:) + F_in(3,:))*h1q
   p_x(2,:) = -1./12.*(3.*F_in(1,:) + 10.*F_in(2,:) - 18.*F_in(3,:)&
                       + 6.*F_in(4,:) - F_in(5,:))*h1q
   p_x(3:nq-2,:) = 1./12.*(F_in(1:nq-4,:) - 8.*F_in(2:nq-3,:)&
                           + 8.*F_in(4:nq-1,:) - F_in(5:nq,:))*h1q !4th order centered
   p_x(nq-1,:) = -1./12.*(F_in(nq-4,:) - 6.*F_in(nq-3,:) + 18.*F_in(nq-2,:)&
                          - 10.*F_in(nq-1,:) - 3.*F_in(nq,:))*h1q
   p_x(nq,:) = 1./2.*(3.*F_in(nq,:) - 4.*F_in(nq-1,:) + F_in(nq-2,:))*h1q

   p_y(:,1) = -1./2.*(3.*F_in(:,1) - 4.*F_in(:,2) + F_in(:,3))*h1p
   p_y(:,2) = -1./12.*(3.*F_in(:,1) + 10.*F_in(:,2) - 18.*F_in(:,3)&
                       + 6.*F_in(:,4) - F_in(:,5))*h1p
   p_y(:,3:np-2) = 1./12.*(F_in(:,1:np-4) - 8.*F_in(:,2:np-3)&
                           + 8.*F_in(:,4:np-1) - F_in(:,5:np))*h1p !4th order centered
   p_y(:,np-1) = -1./12.*(F_in(:,np-4) - 6.*F_in(:,np-3) + 18.*F_in(:,np-2)&
                          - 10.*F_in(:,np-1) - 3.*F_in(:,np))*h1p
   p_y(:,np) = 1./2.*(3.*F_in(:,np) - 4.*F_in(:,np-1) + F_in(:,np-2))*h1p


!--------------------------Mixed D2---------------------------------------   
!2nd order
!   p_xy(1,:) = -1./2.*(3.*p_y(1,:) - 4.*p_y(2,:) + p_y(3,:))*h1q
!   p_xy(2,:) = -1./2.*(p_y(1,:) - p_y(3,:))*h1q
!   p_xy(3:nq-2,:) = -1./2.*(p_y(2:nq-3,:) - p_y(4:nq-1,:))*h1q !2nd order centered
!   p_xy(nq-1,:) = -1./2.*(p_y(nq-2,:) - p_y(nq,:))*h1q
!   p_xy(nq,:) = 1./2.*(3.*p_y(nq,:) - 4.*p_y(nq-1,:) + p_y(nq-2,:))*h1q

!4th order   
   p_xy(1,:) = -1./2.*(3.*p_y(1,:) - 4.*p_y(2,:) + p_y(3,:))*h1q
   p_xy(2,:) = -1./12.*(3.*p_y(1,:) + 10.*p_y(2,:) - 18.*p_y(3,:)&
                       + 6.*p_y(4,:) - p_y(5,:))*h1q
   p_xy(3:nq-2,:) = 1./12.*(p_y(1:nq-4,:) - 8.*p_y(2:nq-3,:)&
                           + 8.*p_y(4:nq-1,:) - p_y(5:nq,:))*h1q !4th order centered
   p_xy(nq-1,:) = -1./12.*(p_y(nq-4,:) - 6.*p_y(nq-3,:) + 18.*p_y(nq-2,:)&
                          - 10.*p_y(nq-1,:) - 3.*p_y(nq,:))*h1q
   p_xy(nq,:) = 1./2.*(3.*p_y(nq,:) - 4.*p_y(nq-1,:) + p_y(nq-2,:))*h1q
!------------------------------------------------------------------

!--------------------------D2---------------------------------------
!2nd order
!   p_xx(1,:) = (F_in(1,:) - 2.*F_in(2,:) + F_in(3,:))*h2q
!   p_xx(2,:) = (F_in(1,:) - 2.*F_in(2,:) + F_in(3,:))*h2q
!   p_xx(3:nq-2,:) = (F_in(2:nq-3,:) - 2.*F_in(3:nq-2,:) + F_in(4:nq-1,:))*h2q !2nd order centered
!   p_xx(nq-1,:) = (F_in(nq-2,:) - 2.*F_in(nq-1,:) + F_in(nq,:))*h2q
!   p_xx(nq,:) = (F_in(nq-2,:) - 2.*F_in(nq-1,:) + F_in(nq,:))*h2q

!   p_yy(:,1) = (F_in(:,1) - 2.*F_in(:,2) + F_in(:,3))*h2p
!   p_yy(:,2) = (F_in(:,1) - 2.*F_in(:,2) + F_in(:,3))*h2p
!   p_yy(:,3:np-2) = (F_in(:,2:np-3) - 2.*F_in(:,3:np-2) + F_in(:,4:np-1))*h2p !2nd order centered
!   p_yy(:,np-1) = (F_in(:,np-2) - 2.*F_in(:,np-1) + F_in(:,np))*h2p
!   p_yy(:,np) = (F_in(:,np-2) - 2.*F_in(:,np-1) + F_in(:,np))*h2p

!4th order
   p_xx(1,:) = (F_in(1,:) - 2.*F_in(2,:) + F_in(3,:))*h2q
   p_xx(2,:) = 1./12.*(11.*F_in(1,:) - 20.*F_in(2,:) + 6.*F_in(3,:)&
                       + 4.*F_in(4,:) - F_in(5,:))*h2q !4th order upwind
   p_xx(3:nq-2,:) = -1./12.*(F_in(1:nq-4,:) - 16.*F_in(2:nq-3,:)&
                             + 30.*F_in(3:nq-2,:) - 16.*F_in(4:nq-1,:)&
                             + F_in(5:nq,:))*h2q !4th order centered
   p_xx(nq-1,:) = -1./12.*(F_in(nq-4,:) - 4.*F_in(nq-3,:) - 6.*F_in(nq-2,:)&
                           + 20.*F_in(nq-1,:) - 11.*F_in(nq,:))*h2q !4th order downwind
   p_xx(nq,:) = (F_in(nq-2,:) - 2.*F_in(nq-1,:) + F_in(nq,:))*h2q

   p_yy(:,1) = (F_in(:,1) - 2.*F_in(:,2) + F_in(:,3))*h2p
   p_yy(:,2) = 1./12.*(11.*F_in(:,1) - 20.*F_in(:,2) + 6.*F_in(:,3)&
                       + 4.*F_in(:,4) - F_in(:,5))*h2p !4th order upwind
   p_yy(:,3:np-2) = -1./12.*(F_in(:,1:np-4) - 16.*F_in(:,2:np-3)&
                             + 30.*F_in(:,3:np-2) - 16.*F_in(:,4:np-1)&
                             + F_in(:,5:np))*h2p !4th order centered
   p_yy(:,np-1) = -1./12.*(F_in(:,np-4) - 4.*F_in(:,np-3) - 6.*F_in(:,np-2)&
                           + 20.*F_in(:,np-1) - 11.*F_in(:,np))*h2p !4th order downwind
   p_yy(:,np) = (F_in(:,np-2) - 2.*F_in(:,np-1) + F_in(:,np))*h2p

!--------------------------Mixed D3---------------------------------------
!2nd order
!   p_xyy(1,:) = -1./2.*(3.*p_yy(1,:) - 4.*p_yy(2,:) + p_yy(3,:))*h1q
!   p_xyy(2,:) = -1./2.*(p_yy(1,:) - p_yy(3,:))*h1q
!   p_xyy(3:nq-2,:) = -1./2.*(p_yy(2:nq-3,:) - p_yy(4:nq-1,:))*h1q !2nd order centered
!   p_xyy(nq-1,:) = -1./2.*(p_yy(nq-2,:) - p_yy(nq,:))*h1q
!   p_xyy(nq,:) = 1./2.*(3.*p_yy(nq,:) - 4.*p_yy(nq-1,:) + p_yy(nq-2,:))*h1q

!   p_xxy(:,1) = -1./2.*(3.*p_xx(:,1) - 4.*p_xx(:,2) + p_xx(:,3))*h1p
!   p_xxy(:,2) = -1./2.*(p_xx(:,1) - p_xx(:,3))*h1p
!   p_xxy(:,3:np-2) = -1./2.*(p_xx(:,2:np-3) - p_xx(:,4:np-1))*h1p !2nd order centered
!   p_xxy(:,np-1) = -1./2.*(p_xx(:,np-2) - p_xx(:,np))*h1p
!   p_xxy(:,np) = 1./2.*(3.*p_xx(:,np) - 4.*p_xx(:,np-1) + p_xx(:,np-2))*h1p

!4th order
   p_xyy(1,:) = -1./2.*(3.*p_yy(1,:) - 4.*p_yy(2,:) + p_yy(3,:))*h1q
   p_xyy(2,:) = -1./12.*(3.*p_yy(1,:) + 10.*p_yy(2,:) - 18.*p_yy(3,:)&
                         + 6.*p_yy(4,:) - p_yy(5,:))*h1q
   p_xyy(3:nq-2,:) = 1./12.*(p_yy(1:nq-4,:) - 8.*p_yy(2:nq-3,:)&
                           + 8.*p_yy(4:nq-1,:) - p_yy(5:nq,:))*h1q !4th order centered
   p_xyy(nq-1,:) = -1./12.*(p_yy(nq-4,:) - 6.*p_yy(nq-3,:) + 18.*p_yy(nq-2,:)&
                          - 10.*p_yy(nq-1,:) - 3.*p_yy(nq,:))*h1q
   p_xyy(nq,:) = 1./2.*(3.*p_yy(nq,:) - 4.*p_yy(nq-1,:) + p_yy(nq-2,:))*h1q

   p_xxy(:,1) = -1./2.*(3.*p_xx(:,1) - 4.*p_xx(:,2) + p_xx(:,3))*h1p
   p_xxy(:,2) = -1./12.*(3.*p_xx(:,1) + 10.*p_xx(:,2) - 18.*p_xx(:,3)&
                         + 6.*p_xx(:,4) - p_xx(:,5))*h1p
   p_xxy(:,3:np-2) = 1./12.*(p_xx(:,1:np-4) - 8.*p_xx(:,2:np-3)&
                           + 8.*p_xx(:,4:np-1) - p_xx(:,5:np))*h1p !4th order centered
   p_xxy(:,np-1) = -1./12.*(p_xx(:,np-4) - 6.*p_xx(:,np-3) + 18.*p_xx(:,np-2)&
                          - 10.*p_xx(:,np-1) - 3.*p_xx(:,np))*h1p
   p_xxy(:,np) = 1./2.*(3.*p_xx(:,np) - 4.*p_xx(:,np-1) + p_xx(:,np-2))*h1p

!--------------------------D3---------------------------------------
! second order
!   p_xxx(1,:) = -1./2.*(5.*F_in(1,:) - 18.*F_in(2,:) + 24.*F_in(3,:)&
!                        - 14.*F_in(4,:) + 3.*F_in(5,:))*h3q
!   p_xxx(2,:) = -1./2.*(3.*F_in(1,:) - 10.*F_in(2,:) + 12.*F_in(3,:)&
!                        - 6.*F_in(4,:) + F_in(5,:))*h3q	
!   p_xxx(3,:) = -1./2.*(F_in(1,:) - 2.*F_in(2,:)&
!                         + 2.*F_in(4,:) - F_in(5,:))*h3q !2nd order centered
!   p_xxx(4:nq-3,:) = -1./2.*(F_in(2:nq-5,:) - 2.*F_in(3:nq-4,:)&
!                             + 2.*F_in(5:nq-2,:) - F_in(6:nq-1,:))*h3q !2nd order centered
!   p_xxx(nq-2,:) = -1./2.*(F_in(nq-4,:) - 2.*F_in(nq-3,:)&
!                           + 2.*F_in(nq-1,:) - F_in(nq,:))*h3q !2nd order centered
!   p_xxx(nq-1,:) = 1./2.*(3.*F_in(nq,:) - 10.*F_in(nq-1,:) + 12.*F_in(nq-2,:)&
!                        - 6.*F_in(nq-3,:) + F_in(nq-4,:))*h3q 
!   p_xxx(nq,:) = 1./2.*(5.*F_in(nq,:) - 18.*F_in(nq-1,:) + 24.*F_in(nq-2,:)&
!                        - 14.*F_in(nq-3,:) + 3.*F_in(nq-4,:))*h3q					

!   p_yyy(:,1) = -1./2.*(5.*F_in(:,1) - 18.*F_in(:,2) + 24.*F_in(:,3)&
!                        - 14.*F_in(:,4) + 3.*F_in(:,5))*h3p
!   p_yyy(:,2) = -1./2.*(3.*F_in(:,1) - 10.*F_in(:,2) + 12.*F_in(:,3)&
!                        - 6.*F_in(:,4) + F_in(:,5))*h3p
!   p_yyy(:,3) = -1./2.*(F_in(:,1) - 2.*F_in(:,2)&
!					    + 2.*F_in(:,4) - F_in(:,5))*h3p !2nd order centered
!   p_yyy(:,4:np-3) = -1./2.*(F_in(:,2:np-5) - 2.*F_in(:,3:np-4)&
!                             + 2.*F_in(:,5:np-2) - F_in(:,6:np-1))*h3p !2nd order centered
!   p_yyy(:,np-2) = -1./2.*(F_in(:,np-4) - 2.*F_in(:,np-3)&
!                             + 2.*F_in(:,np-1) - F_in(:,np))*h3p !2nd order centered
!   p_yyy(:,np-1) = 1./2.*(3.*F_in(:,np) - 10.*F_in(:,np-1) + 12.*F_in(:,np-2)&
!                        - 6.*F_in(:,np-3) + F_in(:,np-4))*h3p
!   p_yyy(:,np) = 1./2.*(5.*F_in(:,np) - 18.*F_in(:,np-1) + 24.*F_in(:,np-2)&
!                        - 14.*F_in(:,np-3) + 3.*F_in(:,np-4))*h3p	

!4th order
   p_xxx(1,:) = -1./2.*(5.*F_in(1,:) - 18.*F_in(2,:) + 24.*F_in(3,:)&
                        - 14.*F_in(4,:) + 3.*F_in(5,:))*h3q
   p_xxx(2,:) = -1./8.*(15.*F_in(1,:) - 56.*F_in(2,:) + 83.*F_in(3,:)&
                        - 64.*F_in(4,:) + 29.*F_in(5,:)&
                         - 8.*F_in(6,:) + F_in(7,:))*h3q !4th order upwind
   p_xxx(3,:) = -1./8.*(F_in(1,:) + 8.*F_in(2,:) - 35.*F_in(3,:)&
                                  + 48.*F_in(4,:) - 29.*F_in(5,:)&
                                  + 8.*F_in(6,:) - F_in(7,:))*h3q !4th order upwind
   p_xxx(4:nq-3,:) = 1./8.*(F_in(1:nq-6,:) - 8.*F_in(2:nq-5,:) + 13.*F_in(3:nq-4,:)&
                   - 13.*F_in(5:nq-2,:) + 8.*F_in(6:nq-1,:) - F_in(7:nq,:))*h3q !4th order centered
   p_xxx(nq-2,:) = -1./8.*(F_in(nq-6,:) - 8.*F_in(nq-5,:) + 29.*F_in(nq-4,:)&
                  - 48.*F_in(nq-3,:) + 35.*F_in(nq-2,:)&
                  - 8.*F_in(nq-1,:) - F_in(nq,:))*h3q !4th order downwind
   p_xxx(nq-1,:) = 1./2.*(3.*F_in(nq,:) - 10.*F_in(nq-1,:) + 12.*F_in(nq-2,:)&
                        - 6.*F_in(nq-3,:) + F_in(nq-4,:))*h3q 
   p_xxx(nq,:) = 1./2.*(5.*F_in(nq,:) - 18.*F_in(nq-1,:) + 24.*F_in(nq-2,:)&
                        - 14.*F_in(nq-3,:) + 3.*F_in(nq-4,:))*h3q

   p_yyy(:,1) = -1./2.*(5.*F_in(:,1) - 18.*F_in(:,2) + 24.*F_in(:,3)&
                        - 14.*F_in(:,4) + 3.*F_in(:,5))*h3p
   p_yyy(:,2) = -1./8.*(15.*F_in(:,1) - 56.*F_in(:,2) + 83.*F_in(:,3)&
                - 64.*F_in(:,4) + 29.*F_in(:,5)&
                - 8.*F_in(:,6) + F_in(:,7))*h3p !4th order upwind
   p_yyy(:,3) = -1./8.*(F_in(:,1) + 8.*F_in(:,2) - 35.*F_in(:,3)&
                                  + 48.*F_in(:,4) - 29.*F_in(:,5)&
                                  + 8.*F_in(:,6) - F_in(:,7))*h3p !4th order upwind
   p_yyy(:,4:np-3) = 1./8.*(F_in(:,1:np-6) - 8.*F_in(:,2:np-5) + 13.*F_in(:,3:np-4)&
                   - 13.*F_in(:,5:np-2) + 8.*F_in(:,6:np-1) - F_in(:,7:np))*h3p !4th order centered
   p_yyy(:,np-2) = -1./8.*(F_in(:,np-6) - 8.*F_in(:,np-5) + 29.*F_in(:,np-4)&
                  - 48.*F_in(:,np-3) + 35.*F_in(:,np-2)&
                  - 8.*F_in(:,np-1) - F_in(:,np))*h3p !4th order downwind
   p_yyy(:,np-1) = 1./2.*(3.*F_in(:,np) - 10.*F_in(:,np-1) + 12.*F_in(:,np-2)&
                        - 6.*F_in(:,np-3) + F_in(:,np-4))*h3p
   p_yyy(:,np) = 1./2.*(5.*F_in(:,np) - 18.*F_in(:,np-1) + 24.*F_in(:,np-2)&
                        - 14.*F_in(:,np-3) + 3.*F_in(:,np-4))*h3p
!------------------------------------------------------------------

end if

   F_out = pp * (qs * (qs * (p_x + ii * e1 * p_y) + pp * (p_xx + ii * e1 *&
        & p_xy) / 0.2D1 + e1 * spin * F_in + (e1 * qs + ii * ps) * spin * p_x + ii&
        & * e2 * (ps * (p_x + ii * e1 * p_y) + pp * (p_xy + ii * e1 *&
        & p_yy) / 0.2D1 + ii * spin * F_in + (e1 * qs + ii * ps) * spin * p_y)) +&
        & pp * (p_x + ii * e1 * p_y + 0.2D1 * qs * (p_xx + ii * e1 * p_xy)&
        & + pp * (p_xxx + ii * e1 * p_xxy) / 0.2D1 + 0.2D1 * e1 * spin *&
        & p_x + (e1 * qs + ii * ps) * spin * p_xx + ii * e2 * (ps * (p_xx + ii&
        & * e1 * p_xy) + qs * (p_xy + ii * e1 * p_yy) + pp * (p_xxy + ii *&
        & e1 * p_xyy) / 0.2D1 + ii * spin * p_x + e1 * spin * p_y + (e1 * qs +&
        & ii * ps) * spin * p_xy)) / 0.2D1 + e2 * (spin + e1) * (pp * (p_x +&
        & ii * e1 * p_y) / 0.2D1 + (e1 * qs + ii * ps) * spin * F_in) + (e2 *&
        & qs + ii * ps) * (spin + e1) * (qs * (p_x + ii * e1 * p_y) + pp *&
        & (p_xx + ii * e1 * p_xy) / 0.2D1 + e1 * spin * F_in + (e1 * qs + ii *&
        & ps) * spin * p_x) + ii * e3 * (ps * (qs * (p_x + ii * e1 * p_y) + pp&
        & * (p_xx + ii * e1 * p_xy) / 0.2D1 + e1 * spin * F_in + (e1 * qs + ii&
        & * ps) * spin * p_x + ii * e2 * (ps * (p_x + ii * e1 * p_y) + pp *&
        & (p_xy + ii * e1 * p_yy) / 0.2D1 + ii * spin * F_in + (e1 * qs + ii *&
        & ps) * spin * p_y)) + pp * (qs * (p_xy + ii * e1 * p_yy) + ps * (p_xx &
        &+ ii * e1 * p_xy) + pp * (p_xxy + ii * e1 * p_xyy) / 0.2D1 + e1&
        & * spin * p_y + ii * spin * p_x + (e1 * qs + ii * ps) * spin * p_xy + ii *&
        & e2 * (p_x + ii * e1 * p_y + 0.2D1 * ps * (p_xy + ii * e1 * p_yy)&
        & + pp * (p_xyy + ii * e1 * p_yyy) / 0.2D1 + 0.2D1 * ii * spin *&
        & p_y + (e1 * qs + ii * ps) * spin * p_yy)) / 0.2D1 + ii * (spin + e1) &
        &* (pp * (p_x + ii * e1 * p_y) / 0.2D1 + (e1 * qs + ii * ps) * spin &
        &* F_in) + (e2 * qs + ii * ps) * (spin + e1) * (ps * (p_x + ii * e1 *&
        & p_y) + pp * (p_xy + ii * e1 * p_yy) / 0.2D1 + ii * spin * F_in + (e1&
        & * qs + ii * ps) * spin * p_y))) / 0.2D1 + (e3 * qs + ii * ps) * (spin &
        &+ e1 + e2) * (pp * (qs * (p_x + ii * e1 * p_y) + pp * (p_xx + ii&
        & * e1 * p_xy) / 0.2D1 + e1 * spin * F_in + (e1 * qs + ii * ps) * spin *&
        & p_x + ii * e2 * (ps * (p_x + ii * e1 * p_y) + pp * (p_xy + ii *&
        & e1 * p_yy) / 0.2D1 + ii * spin * F_in + (e1 * qs + ii * ps) * spin *&
        & p_y)) / 0.2D1 + (e2 * qs + ii * ps) * (spin + e1) * (pp * (p_x +&
        & ii * e1 * p_y) / 0.2D1 + (e1 * qs + ii * ps) * spin * F_in))


   if (use_edge_check.ne.0) then
      ! deliberate error
      F_out(1,:) = 1.0d50 * (1.,1.)
      F_out(nq,:) = 1.0d55 * (1.,1.)
      F_out(:,1) = 1.0d60 * (1.,1.)
      F_out(:,np) = 1.0d65 * (1.,1.)
   end if

end subroutine NullInterp_d3



subroutine NullInterp_d2 (F_out, F_in, spin, e1, e2)
use NullGrid_Vars
   implicit none

   CCTK_COMPLEX,   dimension (:,:), intent (out) :: F_out
   CCTK_COMPLEX,   dimension (:,:), intent (in)  :: F_in
   CCTK_INT,                        intent (in)  :: spin, e1,e2
   CCTK_INT :: nq, np
   CCTK_REAL :: dq, dp
   CCTK_COMPLEX :: ii
   CCTK_REAL :: h1q, h1p, h2q, h2p, c1, c2, c3, c4, c5, c6

   logical, save :: FirstTime = .true.
   CCTK_COMPLEX, dimension (:,:), allocatable, save ::&
       p_x, p_y, p_xx, p_xy, p_yy

   DECLARE_CCTK_PARAMETERS

   nq = lsh(1)
   np = lsh(2)
   dq = dble(delta(1))
   dp = dble(delta(2))
   ii = dcmplx(0.,1.)

  if (FirstTime) then
     FirstTime = .false.
     allocate(p_x(nq,np), p_y(nq,np),&
              p_xx(nq,np), p_xy(nq,np),&
              p_yy(nq,np))

     p_x = 0; p_y = 0
     p_xx = 0; p_xy = 0
     p_yy = 0
  endif
  
   c1 = dble(1 + e1 * e2 + spin * (e1 + e2))
   c2 = dble(2 * spin + e1 + e2)
   c3 = dble(2 * spin * e1 * e2 + e1 + e2)
   c4 = dble(e1 * e2)
   c5 = dble(e1 + e2)
   c6 = dble(e2 - e1)

if (deriv_accuracy.eq.2) then

   h1q = dble(1. / (2. * dq))
   h2q = dble(1. / (dq * dq))
   h1p = dble(1. / (2. * dp))
   h2p = dble(1. / (dp * dp))

!  p_xx(1,1:np)  = (0., 0.)
!  p_xx(nq,1:np) = (0., 0.)
!  p_xx(1:nq,1)  = (0., 0.)
!  p_xx(1:nq,np) = (0., 0.)

!  p_xy(1,1:np)  = (0., 0.)
!  p_xy(nq,1:np) = (0., 0.)
!  p_xy(1:nq,1)  = (0., 0.)
!  p_xy(1:nq,np) = (0., 0.)

!  p_yy(1,1:np)  = (0., 0.)
!  p_yy(nq,1:np) = (0., 0.)
!  p_yy(1:nq,1)  = (0., 0.)
!  p_yy(1:nq,np) = (0., 0.)

   p_x(1,1:np) = (-3. * F_in(1,1:np) + 4. * F_in(2,1:np) - F_in(3,1:np)) * h1q
   p_x(2:nq-1,1:np) = (F_in(3:nq,1:np) - F_in(1:nq-2,1:np)) * h1q
   p_x(nq,1:np) = (3. * F_in(nq,1:np) - 4. * F_in(nq-1,1:np) + F_in(nq-2,1:np)) * h1q

   p_y(1:nq,1) = (-3. * F_in(1:nq,1) + 4. * F_in(1:nq,2) - F_in(1:nq,3)) * h1p
   p_y(1:nq,2:np-1) = (F_in(1:nq,3:np) - F_in(1:nq,1:np-2)) * h1p
   p_y(1:nq,np) = (3. * F_in(1:nq,np) - 4. * F_in(1:nq,np-1) + F_in(1:nq,np-2)) * h1p

   p_xx(2:nq-1,1:np) = (F_in(3:nq,1:np)- 2.*F_in(2:nq-1,1:np)+F_in(1:nq-2,1:np)) * h2q
   p_yy(1:nq,2:np-1) = (F_in(1:nq,3:np)- 2.*F_in(1:nq,2:np-1)+F_in(1:nq,1:np-2)) * h2p

   p_xy(2:nq-1,2:np-1) = ( F_in(3:nq,3:np)   - F_in(3:nq,1:np-2)&
                         - F_in(1:nq-2,3:np) + F_in(1:nq-2,1:np-2))&
                       / (4. * dq * dp)

   p_xy(1,2:np-1) = (-3. * F_in(1,3:np)   + 4. * F_in(2,3:np)   - F_in(3,3:np)&
                    + 3. * F_in(1,1:np-2) - 4. * F_in(2,1:np-2) + F_in(3,1:np-2)) * h1q * h1p
					
end if

if (deriv_accuracy.eq.4) then

   h1q = dble(1./(dq))
   h2q = dble(1./(dq*dq))
   h1p = dble(1./(dp))
   h2p = dble(1./(dp*dp))

!--------------------------D1---------------------------------------
!2nd order
 !  p_x(1,:) = -1./2.*(3.*F_in(1,:) - 4.*F_in(2,:) + F_in(3,:))*h1q
 !  p_x(2,:) = -1./2.*(F_in(1,:) - F_in(3,:))*h1q
 !  p_x(3:nq-2,:) = -1./2.*(F_in(2:nq-3,:) - F_in(4:nq-1,:))*h1q !2nd order centered
 !  p_x(nq-1,:) = -1./2.*(F_in(nq-2,:) - F_in(nq,:))*h1q
 !  p_x(nq,:) = 1./2.*(3.*F_in(nq,:) - 4.*F_in(nq-1,:) + F_in(nq-2,:))*h1q

 !  p_y(:,1) = -1./2.*(3.*F_in(:,1) - 4.*F_in(:,2) + F_in(:,3))*h1p
 !  p_y(:,2) = -1./2.*(F_in(:,1) - F_in(:,3))*h1p
 !  p_y(:,3:np-2) = -1./2.*(F_in(:,2:np-3) - F_in(:,4:np-1))*h1p !2nd order centered
 !  p_y(:,np-1) = -1./2.*(F_in(:,np-2) - F_in(:,np))*h1p
 !  p_y(:,np) = 1./2.*(3.*F_in(:,np) - 4.*F_in(:,np-1) + F_in(:,np-2))*h1p

!4th order
   p_x(1,:) = -1./2.*(3.*F_in(1,:) - 4.*F_in(2,:) + F_in(3,:))*h1q
   p_x(2,:) = -1./12.*(3.*F_in(1,:) + 10.*F_in(2,:) - 18.*F_in(3,:)&
                       + 6.*F_in(4,:) - F_in(5,:))*h1q
   p_x(3:nq-2,:) = 1./12.*(F_in(1:nq-4,:) - 8.*F_in(2:nq-3,:)&
                           + 8.*F_in(4:nq-1,:) - F_in(5:nq,:))*h1q !4th order centered
   p_x(nq-1,:) = -1./12.*(F_in(nq-4,:) - 6.*F_in(nq-3,:) + 18.*F_in(nq-2,:)&
                          - 10.*F_in(nq-1,:) - 3.*F_in(nq,:))*h1q
   p_x(nq,:) = 1./2.*(3.*F_in(nq,:) - 4.*F_in(nq-1,:) + F_in(nq-2,:))*h1q

   p_y(:,1) = -1./2.*(3.*F_in(:,1) - 4.*F_in(:,2) + F_in(:,3))*h1p
   p_y(:,2) = -1./12.*(3.*F_in(:,1) + 10.*F_in(:,2) - 18.*F_in(:,3)&
                       + 6.*F_in(:,4) - F_in(:,5))*h1p
   p_y(:,3:np-2) = 1./12.*(F_in(:,1:np-4) - 8.*F_in(:,2:np-3)&
                           + 8.*F_in(:,4:np-1) - F_in(:,5:np))*h1p !4th order centered
   p_y(:,np-1) = -1./12.*(F_in(:,np-4) - 6.*F_in(:,np-3) + 18.*F_in(:,np-2)&
                          - 10.*F_in(:,np-1) - 3.*F_in(:,np))*h1p
   p_y(:,np) = 1./2.*(3.*F_in(:,np) - 4.*F_in(:,np-1) + F_in(:,np-2))*h1p


!--------------------------Mixed D2---------------------------------------   
!2nd order
!   p_xy(1,:) = -1./2.*(3.*p_y(1,:) - 4.*p_y(2,:) + p_y(3,:))*h1q
!   p_xy(2,:) = -1./2.*(p_y(1,:) - p_y(3,:))*h1q
!   p_xy(3:nq-2,:) = -1./2.*(p_y(2:nq-3,:) - p_y(4:nq-1,:))*h1q !2nd order centered
!   p_xy(nq-1,:) = -1./2.*(p_y(nq-2,:) - p_y(nq,:))*h1q
!   p_xy(nq,:) = 1./2.*(3.*p_y(nq,:) - 4.*p_y(nq-1,:) + p_y(nq-2,:))*h1q

!4th order   
   p_xy(1,:) = -1./2.*(3.*p_y(1,:) - 4.*p_y(2,:) + p_y(3,:))*h1q
   p_xy(2,:) = -1./12.*(3.*p_y(1,:) + 10.*p_y(2,:) - 18.*p_y(3,:)&
                       + 6.*p_y(4,:) - p_y(5,:))*h1q
   p_xy(3:nq-2,:) = 1./12.*(p_y(1:nq-4,:) - 8.*p_y(2:nq-3,:)&
                           + 8.*p_y(4:nq-1,:) - p_y(5:nq,:))*h1q !4th order centered
   p_xy(nq-1,:) = -1./12.*(p_y(nq-4,:) - 6.*p_y(nq-3,:) + 18.*p_y(nq-2,:)&
                          - 10.*p_y(nq-1,:) - 3.*p_y(nq,:))*h1q
   p_xy(nq,:) = 1./2.*(3.*p_y(nq,:) - 4.*p_y(nq-1,:) + p_y(nq-2,:))*h1q
!------------------------------------------------------------------

!--------------------------D2---------------------------------------
!2nd order
!   p_xx(1,:) = (F_in(1,:) - 2.*F_in(2,:) + F_in(3,:))*h2q
!   p_xx(2,:) = (F_in(1,:) - 2.*F_in(2,:) + F_in(3,:))*h2q
!   p_xx(3:nq-2,:) = (F_in(2:nq-3,:) - 2.*F_in(3:nq-2,:) + F_in(4:nq-1,:))*h2q !2nd order centered
!   p_xx(nq-1,:) = (F_in(nq-2,:) - 2.*F_in(nq-1,:) + F_in(nq,:))*h2q
!   p_xx(nq,:) = (F_in(nq-2,:) - 2.*F_in(nq-1,:) + F_in(nq,:))*h2q

!   p_yy(:,1) = (F_in(:,1) - 2.*F_in(:,2) + F_in(:,3))*h2p
!   p_yy(:,2) = (F_in(:,1) - 2.*F_in(:,2) + F_in(:,3))*h2p
!   p_yy(:,3:np-2) = (F_in(:,2:np-3) - 2.*F_in(:,3:np-2) + F_in(:,4:np-1))*h2p !2nd order centered
!   p_yy(:,np-1) = (F_in(:,np-2) - 2.*F_in(:,np-1) + F_in(:,np))*h2p
!   p_yy(:,np) = (F_in(:,np-2) - 2.*F_in(:,np-1) + F_in(:,np))*h2p

!4th order
   p_xx(1,:) = (F_in(1,:) - 2.*F_in(2,:) + F_in(3,:))*h2q
   p_xx(2,:) = 1./12.*(11.*F_in(1,:) - 20.*F_in(2,:) + 6.*F_in(3,:)&
                       + 4.*F_in(4,:) - F_in(5,:))*h2q !4th order upwind
   p_xx(3:nq-2,:) = -1./12.*(F_in(1:nq-4,:) - 16.*F_in(2:nq-3,:)&
                             + 30.*F_in(3:nq-2,:) - 16.*F_in(4:nq-1,:)&
                             + F_in(5:nq,:))*h2q !4th order centered
   p_xx(nq-1,:) = -1./12.*(F_in(nq-4,:) - 4.*F_in(nq-3,:) - 6.*F_in(nq-2,:)&
                           + 20.*F_in(nq-1,:) - 11.*F_in(nq,:))*h2q !4th order downwind
   p_xx(nq,:) = (F_in(nq-2,:) - 2.*F_in(nq-1,:) + F_in(nq,:))*h2q

   p_yy(:,1) = (F_in(:,1) - 2.*F_in(:,2) + F_in(:,3))*h2p
   p_yy(:,2) = 1./12.*(11.*F_in(:,1) - 20.*F_in(:,2) + 6.*F_in(:,3)&
                       + 4.*F_in(:,4) - F_in(:,5))*h2p !4th order upwind
   p_yy(:,3:np-2) = -1./12.*(F_in(:,1:np-4) - 16.*F_in(:,2:np-3)&
                             + 30.*F_in(:,3:np-2) - 16.*F_in(:,4:np-1)&
                             + F_in(:,5:np))*h2p !4th order centered
   p_yy(:,np-1) = -1./12.*(F_in(:,np-4) - 4.*F_in(:,np-3) - 6.*F_in(:,np-2)&
                           + 20.*F_in(:,np-1) - 11.*F_in(:,np))*h2p !4th order downwind
   p_yy(:,np) = (F_in(:,np-2) - 2.*F_in(:,np-1) + F_in(:,np))*h2p

end if


   F_out(1:nq,1:np) = 0.25 * pp(1:nq,1:np) * pp(1:nq,1:np) * (p_xx(1:nq,1:np)&
                    - c4 * p_yy(1:nq,1:np) + ii * c5 * p_xy(1:nq,1:np))&
        + 0.5 * pp(1:nq,1:np) * ((c1 * qs(1:nq,1:np) + ii * c2 * ps(1:nq,1:np)) * p_x(1:nq,1:np)&
        + (c2 * qs(1:nq,1:np) + ii * c1 * ps(1:nq,1:np)) * ii * c4 * p_y(1:nq,1:np))&
        + spin * F_in(1:nq,1:np) * 0.5 *&
        (c6 + c3 * (qs(1:nq,1:np) * qs(1:nq,1:np) - c4 * ps(1:nq,1:np) * ps(1:nq,1:np))&
        + 2. * c1 * ii * qs(1:nq,1:np) * ps(1:nq,1:np))

   if (use_edge_check.ne.0) then
      ! deliberate error
      F_out(1,:) = 1.0d50 * (1.,1.)
      F_out(nq,:) = 1.0d55 * (1.,1.)
      F_out(:,1) = 1.0d60 * (1.,1.)
      F_out(:,np) = 1.0d65 * (1.,1.)
   end if

end subroutine NullInterp_d2

subroutine NullInterp_d1(F_out, F_in, spin, e1)
  use NullGrid_Vars 
   implicit none
   CCTK_COMPLEX, dimension (:,:), intent (in) :: F_in
   CCTK_COMPLEX, dimension (:,:), intent (out) :: F_out
   CCTK_INT, intent (in)  :: spin, e1

   logical, save :: FirstTime = .true.
   CCTK_COMPLEX, dimension (:,:), allocatable, save :: p_x, p_y
   CCTK_INT :: nq, np
   CCTK_REAL :: dq, dp
   CCTK_COMPLEX :: ii 
   CCTK_REAL :: h1q, h1p

   DECLARE_CCTK_PARAMETERS

   nq = lsh(1)
   np = lsh(2)
   dq = dble(delta(1))
   dp = dble(delta(2))
   ii = dcmplx(0.,1.)

   if (FirstTime) then
     FirstTime = .false.
     allocate(p_x(nq,np), p_y(nq, np))
     p_x = 0; p_y = 0
   endif

 if (deriv_accuracy.eq.2) then

   h1q = dble(.5 / dq)
   h1p = dble(.5 / dp)

   p_x(1,1:np) = (-3. * F_in(1,1:np) + 4. * F_in(2,1:np) - F_in(3,1:np)) * h1q
   p_x(2:nq-1,1:np) = (F_in(3:nq,1:np) - F_in(1:nq-2,1:np)) * h1q
   p_x(nq,1:np) = (3. * F_in(nq,1:np) - 4. * F_in(nq-1,1:np) + F_in(nq-2,1:np)) * h1q

   p_y(1:nq,1) = (-3. * F_in(1:nq,1) + 4. * F_in(1:nq,2) - F_in(1:nq,3)) * h1p
   p_y(1:nq,2:np-1) = (F_in(1:nq,3:np) - F_in(1:nq,1:np-2)) *h1p
   p_y(1:nq,np) = (3. * F_in(1:nq,np) - 4. * F_in(1:nq,np-1) + F_in(1:nq,np-2)) * h1p

end if 

if (deriv_accuracy.eq.4) then

   h1q = dble(1./(dq))
   h1p = dble(1./(dp))

!--------------------------D1---------------------------------------
!2nd order
 !  p_x(1,:) = -1./2.*(3.*F_in(1,:) - 4.*F_in(2,:) + F_in(3,:))*h1q
 !  p_x(2,:) = -1./2.*(F_in(1,:) - F_in(3,:))*h1q
 !  p_x(3:nq-2,:) = -1./2.*(F_in(2:nq-3,:) - F_in(4:nq-1,:))*h1q !2nd order centered
 !  p_x(nq-1,:) = -1./2.*(F_in(nq-2,:) - F_in(nq,:))*h1q
 !  p_x(nq,:) = 1./2.*(3.*F_in(nq,:) - 4.*F_in(nq-1,:) + F_in(nq-2,:))*h1q

 !  p_y(:,1) = -1./2.*(3.*F_in(:,1) - 4.*F_in(:,2) + F_in(:,3))*h1p
 !  p_y(:,2) = -1./2.*(F_in(:,1) - F_in(:,3))*h1p
 !  p_y(:,3:np-2) = -1./2.*(F_in(:,2:np-3) - F_in(:,4:np-1))*h1p !2nd order centered
 !  p_y(:,np-1) = -1./2.*(F_in(:,np-2) - F_in(:,np))*h1p
 !  p_y(:,np) = 1./2.*(3.*F_in(:,np) - 4.*F_in(:,np-1) + F_in(:,np-2))*h1p

!4th order
   p_x(1,:) = -1./2.*(3.*F_in(1,:) - 4.*F_in(2,:) + F_in(3,:))*h1q
   p_x(2,:) = -1./12.*(3.*F_in(1,:) + 10.*F_in(2,:) - 18.*F_in(3,:)&
                       + 6.*F_in(4,:) - F_in(5,:))*h1q
   p_x(3:nq-2,:) = 1./12.*(F_in(1:nq-4,:) - 8.*F_in(2:nq-3,:)&
                           + 8.*F_in(4:nq-1,:) - F_in(5:nq,:))*h1q !4th order centered
   p_x(nq-1,:) = -1./12.*(F_in(nq-4,:) - 6.*F_in(nq-3,:) + 18.*F_in(nq-2,:)&
                          - 10.*F_in(nq-1,:) - 3.*F_in(nq,:))*h1q
   p_x(nq,:) = 1./2.*(3.*F_in(nq,:) - 4.*F_in(nq-1,:) + F_in(nq-2,:))*h1q

   p_y(:,1) = -1./2.*(3.*F_in(:,1) - 4.*F_in(:,2) + F_in(:,3))*h1p
   p_y(:,2) = -1./12.*(3.*F_in(:,1) + 10.*F_in(:,2) - 18.*F_in(:,3)&
                       + 6.*F_in(:,4) - F_in(:,5))*h1p
   p_y(:,3:np-2) = 1./12.*(F_in(:,1:np-4) - 8.*F_in(:,2:np-3)&
                           + 8.*F_in(:,4:np-1) - F_in(:,5:np))*h1p !4th order centered
   p_y(:,np-1) = -1./12.*(F_in(:,np-4) - 6.*F_in(:,np-3) + 18.*F_in(:,np-2)&
                          - 10.*F_in(:,np-1) - 3.*F_in(:,np))*h1p
   p_y(:,np) = 1./2.*(3.*F_in(:,np) - 4.*F_in(:,np-1) + F_in(:,np-2))*h1p


end if


   F_out(1:nq,1:np) = 0.5 * pp(1:nq,1:np) * (p_x(1:nq,1:np) + ii * e1 * p_y(1:nq,1:np)) +&
          spin * (e1 * qs(1:nq,1:np) + ii * ps(1:nq,1:np)) * F_in(1:nq,1:np)

   if (use_edge_check.ne.0) then
      ! deliberate error
      F_out(1,:) = 1.0d50 * (1.,1.)
      F_out(nq:,:) = 1.0d55 * (1.,1.)
      F_out(:,1) = 1.0d60 * (1.,1.)
      F_out(:,np:) = 1.0d65 * (1.,1.)
   end if
end subroutine NullInterp_d1

end module NullInterp_Eth
