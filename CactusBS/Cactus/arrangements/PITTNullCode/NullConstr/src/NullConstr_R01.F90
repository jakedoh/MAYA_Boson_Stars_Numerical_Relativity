! vim: syntax=fortran

#include "cctk.h"

module NullConstr_R01

contains

subroutine NullConstr_R01_calc (R01,nq,np,r, &
  w_00,w_01,w_02,w_03,w_04,w_11,w_12,w_13,w_14, &
  w_22,w_23,w_24,w_33,w_34, &
  b_00,b_01,b_02,b_03,b_04,b_11,b_12,b_13,b_14, &
  b_22,b_23,b_24,b_33,b_34, &
  k_00,k_01,k_02,k_03,k_04,k_11,k_12,k_13,k_14, &
  k_22,k_23,k_24,k_33,k_34, &
  j_00,j_01,j_02,j_03,j_04,j_11,j_12,j_13,j_14, &
  j_22,j_23,j_24,j_33,j_34, &
  jb_00,jb_01,jb_02,jb_03,jb_04,jb_11,jb_12,jb_13,jb_14, &
  jb_22,jb_23,jb_24,jb_33,jb_34, &
  u_00,u_01,u_02,u_03,u_04,u_11,u_12,u_13,u_14, &
  u_22,u_23,u_24,u_33,u_34, &
  ub_00,ub_01,ub_02,ub_03,ub_04,ub_11,ub_12,ub_13,ub_14, &
  ub_22,ub_23,ub_24,ub_33,ub_34)
! ub_22,ub_23,ub_24,ub_33,ub_34,omm,cctk_time,Ylm_0,patchID)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compilation on Intel 18 is very slow.
!! This forces a lower level of optimization (-O2)
!DIR$ OPTIMIZE: 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  implicit none
  
  CCTK_INT, intent(in) :: nq,np
  CCTK_REAL, intent(in) :: r
  CCTK_REAL, dimension(nq,np), intent(out) :: R01
  CCTK_REAL, dimension(nq,np), intent(in) :: w_00,w_01,w_04,w_11,w_14
  CCTK_COMPLEX, dimension(nq,np), intent(in) :: w_02,w_03,w_12,w_13,w_22, &
  w_23,w_24,w_33,w_34
  CCTK_REAL, dimension(nq,np), intent(in) :: b_00,b_01,b_04,b_11,b_14
  CCTK_COMPLEX, dimension(nq,np), intent(in) :: b_02,b_03,b_12,b_13,b_22, &
  b_23,b_24,b_33,b_34
  CCTK_REAL, dimension(nq,np), intent(in) :: k_00,k_01,k_04,k_11,k_14
  CCTK_COMPLEX, dimension(nq,np), intent(in) :: k_02,k_03,k_12,k_13,k_22, &
  k_23,k_24,k_33,k_34
  CCTK_COMPLEX, dimension(nq,np), intent(in) :: j_00,j_01,j_04,j_11,j_14
  CCTK_COMPLEX, dimension(nq,np), intent(in) :: j_02,j_03,j_12,j_13,j_22, &
  j_23,j_24,j_33,j_34
  CCTK_COMPLEX, dimension(nq,np), intent(in) :: jb_00,jb_01,jb_04,jb_11,jb_14
  CCTK_COMPLEX, dimension(nq,np), intent(in) :: jb_02,jb_03,jb_12,jb_13,jb_22, &
  jb_23,jb_24,jb_33,jb_34
  CCTK_COMPLEX, dimension(nq,np), intent(in) :: u_00,u_01,u_04,u_11,u_14
  CCTK_COMPLEX, dimension(nq,np), intent(in) :: u_02,u_03,u_12,u_13,u_22, &
  u_23,u_24,u_33,u_34
  CCTK_COMPLEX, dimension(nq,np), intent(in) :: ub_00,ub_01,ub_04,ub_11,ub_14
  CCTK_COMPLEX, dimension(nq,np), intent(in) :: ub_02,ub_03,ub_12,ub_13,ub_22, &
  ub_23,ub_24,ub_33,ub_34
  CCTK_COMPLEX :: s2,s4,s5,s6,s7
  CCTK_REAL :: s1,s3,e2bi,e2b

  integer :: p,q

! logical, parameter :: debug = .true.
! CCTK_COMPLEX  :: einuu,ii
! CCTK_REAL, intent(in) :: omm,cctk_time
! CCTK_COMPLEX, intent(in) :: Ylm_0(nq,np,2)
! CCTK_INT, intent(in)  :: patchID

  do p=1,np
    do q=1,nq
  
  e2b = exp(2*b_00(q,p))
  e2bi = 1.0d0/e2b

      s1 = -1/(r*r)/16.0
      s3 = 1/(k_00(q,p)*k_00(q,p)*k_00(q,p))
      s6 = 4.0*j_00(q,p)*r*r*r*r*ub_00(q,p)*ub_11(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-2.0*r*r*u_00(q,p)*jb_02(q,p)* &
k_00(q,p)*k_01(q,p)+2.0*r*r*u_00(q,p)*jb_12(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+4.0*jb_00(q,p)*r*r*u_00(q,p)*k_00(q,p)*k_00(q,p)* &
k_01(q,p)*k_02(q,p)-2.0*r*r*u_00(q,p)*jb_01(q,p)*k_00(q,p)*k_02(q,p)+2.0*r*r*u_00(q,p)*jb_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_02(q,p) &
-6.0*r*r*u_00(q,p)*jb_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)+3.0*j_00(q,p)*r*r*u_00(q,p)*jb_01(q,p)*jb_02(q,p)*k_00(q,p)*k_00(q,p) &
-2.0*j_00(q,p)*jb_00(q,p)*r*r*u_00(q,p)*jb_01(q,p)*k_00(q,p)*k_02(q,p)-2.0*r*r*u_00(q,p)*k_00(q,p)*k_01(q,p)*j_00(q,p)*jb_00(q,p)* &
jb_02(q,p)+2.0*r*r*u_00(q,p)*jb_12(q,p)*k_00(q,p)*k_00(q,p)+r*r*u_00(q,p)*jb_01(q,p)*j_00(q,p)*jb_02(q,p)-2.0*jb_00(q,p)*jb_00(q,p)*r* &
r*u_00(q,p)*j_12(q,p)*k_00(q,p)*k_00(q,p)+2.0*j_00(q,p)*jb_01(q,p)*jb_00(q,p)*r*r*u_02(q,p)*k_00(q,p)*k_00(q,p)-4.0*jb_00(q,p)*r*r* &
u_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)+jb_00(q,p)*r*r*u_00(q,p)*j_00(q,p)*j_00(q,p)*jb_01(q,p)*jb_02(q,p)+2.0*r*r*jb_01(q,p)* &
u_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+2.0*r*r*jb_01(q,p)*u_02(q,p)*k_00(q,p)*k_00(q,p)+4.0*r*r*u_13(q,p)*k_00(q,p)*k_00(q,p)* &
k_00(q,p)+8.0*r*u_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+4.0*r*r*ub_12(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+8.0*r*ub_02(q,p)*k_00(q,p)* &
k_00(q,p)*k_00(q,p)-32.0*r*b_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-32.0*r*r*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_01(q,p)
      s5 = 32.0*b_14(q,p)*r*r*k_00(q,p)*k_00(q,p)*k_00(q,p)-2.0*r*r*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*j_03(q,p)*jb_01(q,p) &
-16.0*r*r*r*b_11(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)-48.0*r*r*b_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)-16.0*r*r &
*r*b_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_01(q,p)+4.0*r*r*k_00(q,p)*k_00(q,p)*k_01(q,p)*j_00(q,p)*jb_04(q,p)+2.0*r*r*k_00(q,p)*k_00(q,p) &
*k_00(q,p)*j_01(q,p)*jb_04(q,p)+4.0*j_00(q,p)*jb_01(q,p)*r*r*k_00(q,p)*k_00(q,p)*k_04(q,p)+2.0*r*r*k_00(q,p)*j_04(q,p)*jb_01(q,p)-4.0* &
j_00(q,p)*j_00(q,p)*r*r*k_00(q,p)*jb_01(q,p)*jb_04(q,p)+8.0*r*r*u_00(q,p)*b_13(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+8.0*r*r*u_01(q,p)* &
b_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+s6+16.0*b_03(q,p)*r*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+2.0*r*r*u_00(q,p)*k_00(q,p)*k_00(q,p)* &
k_00(q,p)*j_01(q,p)*jb_03(q,p)-2.0*j_00(q,p)*r*r*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*jb_13(q,p)+2.0*jb_00(q,p)*r*r*u_00(q,p)*k_00(q,p) &
*k_00(q,p)*k_00(q,p)*j_13(q,p)-4.0*j_00(q,p)*r*r*jb_01(q,p)*u_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-j_00(q,p)*r*r*ub_00(q,p)*j_01(q,p)* &
jb_03(q,p)+2.0*r*r*j_01(q,p)*ub_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+2.0*r*r*ub_00(q,p)*j_13(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)* &
k_00(q,p)-4.0*j_00(q,p)*r*r*ub_00(q,p)*j_03(q,p)*jb_01(q,p)*k_00(q,p)*k_00(q,p)+2.0*r*r*ub_00(q,p)*j_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p) &
*k_01(q,p)-2.0*r*r*ub_00(q,p)*j_03(q,p)*k_00(q,p)*k_01(q,p)+r*r*j_00(q,p)*j_00(q,p)*j_00(q,p)*ub_00(q,p)*jb_01(q,p)*jb_03(q,p)
      s6 = s5-2.0*j_00(q,p)*j_00(q,p)*r*r*ub_03(q,p)*jb_01(q,p)*k_00(q,p)*k_00(q,p)-2.0*j_00(q,p)*j_00(q,p)*r*r*ub_00(q,p)* &
jb_13(q,p)*k_00(q,p)*k_00(q,p)+j_00(q,p)*r*r*ub_00(q,p)*j_01(q,p)*jb_03(q,p)*k_00(q,p)*k_00(q,p)+2.0*r*r*ub_00(q,p)*j_13(q,p)*k_00(q,p)* &
k_00(q,p)+2.0*r*r*j_01(q,p)*ub_03(q,p)*k_00(q,p)*k_00(q,p)+4.0*r*r*u_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)+8.0*r* &
r*ub_00(q,p)*b_12(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+2.0*j_00(q,p)*r*r*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*jb_12(q,p)+4.0*j_00(q,p)*r* &
r*ub_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*jb_01(q,p)-2.0*r*r*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*j_01(q,p)*jb_02(q,p)-2.0*jb_00(q,p)* &
r*r*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*j_12(q,p)+2.0*r*r*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*j_02(q,p)*jb_01(q,p)-8.0*r*r* &
k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*k_04(q,p)+2.0*r*r*k_00(q,p)*k_00(q,p)*k_00(q,p)*j_04(q,p)*jb_01(q,p)+2.0*r*r*k_00(q,p)*j_01(q,p)* &
jb_04(q,p)-16.0*r*r*b_11(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-8.0*r*r*r*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_11(q,p)-16.0*r*k_00(q,p)* &
k_00(q,p)*k_00(q,p)*w_00(q,p)-4.0*r*r*ub_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)+16.0*b_02(q,p)*r*ub_00(q,p)*k_00(q,p)* &
k_00(q,p)*k_00(q,p)+8.0*r*r*ub_01(q,p)*b_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+8.0*j_00(q,p)*b_33(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2b+ &
8.0*b_03(q,p)*j_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2b-8.0*b_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_03(q,p)*e2b
      s7 = s6-16.0*b_23(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2b-32.0*b_03(q,p)*b_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)* &
k_00(q,p)*e2b+8.0*jb_00(q,p)*b_22(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2b+16.0*j_00(q,p)*b_03(q,p)*b_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)* &
e2b+16.0*jb_00(q,p)*b_02(q,p)*b_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2b-8.0*b_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_02(q,p)*e2b+ &
8.0*b_02(q,p)*jb_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2b+4.0*jb_00(q,p)*r*r*r*r*u_00(q,p)*u_11(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)* &
e2bi+4.0*r*r*r*r*jb_01(q,p)*u_00(q,p)*u_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+16.0*jb_00(q,p)*r*r*r*u_00(q,p)*u_01(q,p) &
*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+16.0*j_00(q,p)*r*r*r*ub_00(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+4.0*r*r*r* &
r*ub_00(q,p)*u_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi
      s4 = s7-8.0*jb_00(q,p)*b_01(q,p)*r*r*r*r*u_00(q,p)*u_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+4.0*j_00(q,p)*r*r &
*r*r*ub_01(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+8.0*r*r*r*r*ub_01(q,p)*u_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p) &
*e2bi+16.0*r*r*r*ub_00(q,p)*u_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+16.0*r*r*r*u_00(q,p)*ub_01(q,p)*k_00(q,p) &
*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+4.0*jb_00(q,p)*r*r*r*r*u_01(q,p)*u_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+4.0*r*r*r* &
r*ub_00(q,p)*u_11(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-8.0*b_01(q,p)*r*r*r*r*ub_00(q,p)*u_01(q,p)*k_00(q,p)*k_00(q,p)* &
k_00(q,p)*k_00(q,p)*e2bi-8.0*b_01(q,p)*r*r*r*r*u_00(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+4.0*r*r*r* &
r*u_00(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi+4.0*r*r*r*r*u_00(q,p)*ub_11(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)* &
k_00(q,p)*e2bi+4.0*r*r*r*r*j_01(q,p)*ub_00(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-8.0*j_00(q,p)*b_01(q,p)*r*r*r &
*r*ub_00(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi
      s2 = s3*s4
      R01(q,p) = s1*s2

    end do
  end do

!s4=(4*b_23-4*r*u_03-2*r**2*u_13+4*r*w_00+8*r**2*w_01+2*r**3*w_11)/(4*r**2)
!R01=real(s4)
!R01=real(s4-s1*s2)

!  if (debug) then
!!    write(0,*) 'In routine: NullConstr_R01_calc, cctk_time,omm = ', cctk_time, omm
!
!     ii=dcmplx(0.0d0,1.0d0)
!     einuu=dcmplx(cos(omm*cctk_time),sin(omm*cctk_time))
!     R01 = dble((w_00/r + 2.0d0*w_01 + 0.5d0*r*w_11 - 0.5d0/r*u_03 - &
!           0.5d0/r*ub_02 - 0.25d0*u_13 - 0.25d0*ub_12)*einuu)*Ylm_0(:,:,patchID)
!  end if

end subroutine NullConstr_R01_calc

end module NullConstr_R01
