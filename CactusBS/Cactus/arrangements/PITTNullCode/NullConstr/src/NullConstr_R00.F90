! vim: syntax=fortran

#include "cctk.h"

module NullConstr_R00

contains

subroutine NullConstr_R00_calc (R00,nq,np,r, &
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compilation on Intel 18 is very slow.
!! This forces a lower level of optimization (-O2)
!DIR$ OPTIMIZE: 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  implicit none
  
  CCTK_INT, intent(in) :: nq,np
  CCTK_REAL, intent(in) :: r
  CCTK_REAL, dimension(nq,np), intent(out) :: R00
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
  CCTK_COMPLEX :: s2,s4,s5,s6,s7,s8,s9,s10
  CCTK_REAL :: s1,s3,e2bi,e2b

  integer :: p,q

  do p=1,np
    do q=1,nq
 
  e2b = exp(2*b_00(q,p))
  e2bi = 1.0d0/e2b
      s1 = 1/(r*r)/16.0
      s3 = 1/(k_00(q,p)*k_00(q,p)*k_00(q,p))
      s9 = 8.0*r*r*u_00(q,p)*u_00(q,p)*b_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_03(q,p)+16.0*r*r*b_04(q,p)*u_03(q,p)*&
k_00(q,p)*k_00(q,p)*k_00(q,p)+16.0*r*r*u_00(q,p)*b_34(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-8.0*r*r*u_00(q,p)*u_00(q,p)*b_03(q,p)*b_03(q,p)*&
k_00(q,p)*k_00(q,p)*k_00(q,p)+32.0*b_03(q,p)*b_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*r*e2b+16.0*b_23(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*r*e2b-8.0*b_02(q,p)*jb_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*r*e2b-16.0*j_00(q,p)*&
b_03(q,p)*b_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*r*e2b+12.0*r*r*ub_02(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*&
k_02(q,p)-8.0*b_03(q,p)*j_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*r*e2b-8.0*j_00(q,p)*b_33(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*&
r*e2b-8.0*r*r*ub_00(q,p)*ub_00(q,p)*b_02(q,p)*b_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-8.0*j_00(q,p)*r*r*r*r*ub_00(q,p)*ub_11(q,p)&
*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+8.0*b_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*r*k_03(q,p)*e2b+4.0*r*w_02(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)*k_03(q,p)*e2b+8.0*b_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*r*k_02(q,p)*e2b+4.0*r*w_03(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)*k_02(q,p)*e2b+8.0*r*w_02(q,p)*b_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2b
      s8 = s9-8.0*jb_00(q,p)*b_02(q,p)*r*w_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2b-4.0*jb_00(q,p)*r*w_22(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)*e2b-4.0*r*w_02(q,p)*jb_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2b-8.0*j_00(q,p)*r*w_03(q,p)*b_03(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)*e2b-4.0*r*w_03(q,p)*j_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2b-4.0*j_00(q,p)*r*w_33(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)*e2b-4.0*j_00(q,p)*r*r*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*jb_24(q,p)+4.0*r*r*r*r*r*r*u_00(q,p)*ub_00(q,p)*&
u_01(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi*e2bi-2.0*j_00(q,p)*j_00(q,p)*r*r*r*r*r*r*ub_00(q,p)*ub_00(q,p)*ub_01(q,p)&
*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi*e2bi-2.0*jb_00(q,p)*jb_00(q,p)*r*r*r*r*r*r*u_00(q,p)*u_00(q,p)*u_01(q,p)*u_01(q,p)&
*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi*e2bi-2.0*r*r*r*r*r*r*ub_00(q,p)*ub_00(q,p)*u_01(q,p)*u_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)*e2bi*e2bi-4.0*jb_00(q,p)*r*r*r*r*r*r*ub_00(q,p)*u_00(q,p)*u_01(q,p)*u_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)*e2bi*e2bi-4.0*j_00(q,p)*r*r*r*r*r*r*u_00(q,p)*ub_00(q,p)*ub_01(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)&
*e2bi*e2bi-4.0*j_00(q,p)*r*r*r*r*r*r*ub_00(q,p)*ub_00(q,p)*u_01(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi&
*e2bi+32.0*r*b_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+32.0*r*r*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_01(q,p)-32.0*b_14(q,p)*r*r*&
k_00(q,p)*k_00(q,p)*k_00(q,p)+32.0*r*r*r*b_11(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)
      s9 = 80.0*r*r*b_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)+16.0*r*r*r*b_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*&
w_01(q,p)-8.0*r*r*r*r*r*r*u_00(q,p)*ub_00(q,p)*u_01(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi*e2bi&
-4.0*jb_00(q,p)*r*r*r*r*r*r*u_00(q,p)*u_00(q,p)*ub_01(q,p)*u_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi*e2bi-2.0*r&
*r*r*r*r*r*u_00(q,p)*u_00(q,p)*ub_01(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi*e2bi-16.0*r*r*&
u_00(q,p)*b_13(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-8.0*r*r*u_01(q,p)*b_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-16.0*b_03(q,p)*r*u_00(q,p)*k_00(q,p)&
*k_00(q,p)*k_00(q,p)-16.0*r*r*ub_00(q,p)*b_12(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+16.0*r*r*b_11(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+8.0*&
r*r*r*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_11(q,p)+16.0*r*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)-16.0*b_02(q,p)*r*ub_00(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)-8.0*r*r*ub_01(q,p)*b_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-8.0*j_00(q,p)*b_33(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2b&
-16.0*jb_00(q,p)*b_02(q,p)*b_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*r*e2b+8.0*r*w_03(q,p)*b_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)*e2b-8.0*b_03(q,p)*j_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2b+8.0*b_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_03(q,p)*e2b
      s7 = s9+16.0*b_23(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2b+32.0*b_03(q,p)*b_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)*e2b-8.0*jb_00(q,p)*b_22(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2b-16.0*j_00(q,p)*b_03(q,p)*b_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*&
e2b-16.0*jb_00(q,p)*b_02(q,p)*b_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2b+8.0*b_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_02(q,p)*e2b&
-8.0*b_02(q,p)*jb_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2b-8.0*jb_00(q,p)*r*r*r*r*u_00(q,p)*u_11(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*&
e2bi-8.0*r*r*r*r*jb_01(q,p)*u_00(q,p)*u_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-32.0*jb_00(q,p)*r*r*r*u_00(q,p)*u_01(q,p)&
*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-32.0*j_00(q,p)*r*r*r*ub_00(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-8.0*r*r*r*&
r*ub_00(q,p)*u_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi+16.0*jb_00(q,p)*b_01(q,p)*r*r*r*r*u_00(q,p)*u_01(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)*e2bi-4.0*j_00(q,p)*r*r*r*r*ub_01(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-8.0*r*r*r*r*&
ub_01(q,p)*u_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-32.0*r*r*r*ub_00(q,p)*u_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*&
e2bi-32.0*r*r*r*u_00(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-4.0*jb_00(q,p)*r*r*r*r*u_01(q,p)*u_01(q,p)&
*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+s8
      s9 = s7-8.0*r*r*r*r*ub_00(q,p)*u_11(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+16.0*b_01(q,p)*r*r*r*r&
*ub_00(q,p)*u_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+16.0*b_01(q,p)*r*r*r*r*u_00(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)*e2bi-8.0*r*r*r*r*u_00(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi-8.0*r*r*r*r*&
u_00(q,p)*ub_11(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-8.0*r*r*r*r*j_01(q,p)*ub_00(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)&
*e2bi+16.0*j_00(q,p)*b_01(q,p)*r*r*r*r*ub_00(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-32.0*r*b_04(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)+16.0*r*r*w_04(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+16.0*r*r*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*w_00(q,p)+48.0*&
r*r*r*b_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*w_00(q,p)+16.0*r*r*r*r*b_11(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*w_00(q,p)+&
8.0*r*r*r*r*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*w_11(q,p)+32.0*r*r*r*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*w_01(q,p)-32.0*r&
*r*b_04(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)-32.0*r*r*r*b_14(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)+16.0*r*r*r*r*&
b_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*w_01(q,p)
      s8 = s9+4.0*r*r*u_00(q,p)*u_33(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-2.0*r*r*u_00(q,p)*u_00(q,p)*k_00(q,p)&
*k_00(q,p)*k_00(q,p)*k_03(q,p)*k_03(q,p)+4.0*r*r*u_00(q,p)*u_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_03(q,p)-16.0*b_01(q,p)*u_00(q,p)*&
r*r*r*w_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-8.0*r*r*u_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_04(q,p)-4.0*r*r*jb_04(q,p)*&
j_04(q,p)*k_00(q,p)-4.0*r*r*jb_04(q,p)*j_04(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+4.0*j_00(q,p)*j_00(q,p)*r*r*k_00(q,p)*jb_04(q,p)*jb_04(q,p)+&
4.0*r*r*u_00(q,p)*jb_04(q,p)*j_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-8.0*j_00(q,p)*jb_04(q,p)*r*r*k_00(q,p)*k_00(q,p)*k_04(q,p)+r*r*&
jb_03(q,p)*u_00(q,p)*u_00(q,p)*j_03(q,p)*k_00(q,p)+r*r*jb_03(q,p)*u_00(q,p)*u_00(q,p)*j_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-4.0*jb_00(q,p)*r*r&
*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*j_34(q,p)-4.0*j_00(q,p)*jb_03(q,p)*r*r*u_00(q,p)*u_00(q,p)*b_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-8.0*&
j_00(q,p)*r*r*r*jb_01(q,p)*u_00(q,p)*w_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-j_00(q,p)*j_00(q,p)*r*r*jb_03(q,p)*jb_03(q,p)*u_00(q,p)*u_00(q,p)*&
k_00(q,p)-4.0*r*r*u_00(q,p)*u_33(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-4.0*r*r*u_02(q,p)*ub_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)+4.0*r*r*ub_00(q,p)*ub_22(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)
      s9 = 4.0*j_00(q,p)*r*r*u_00(q,p)*jb_34(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+2.0*j_00(q,p)*jb_03(q,p)*r*r*u_00(q,p)*u_00(q,p)*&
k_00(q,p)*k_00(q,p)*k_03(q,p)+4.0*j_00(q,p)*r*r*u_03(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*jb_03(q,p)+8.0*j_00(q,p)*r*r*jb_04(q,p)*&
u_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+2.0*j_00(q,p)*r*r*u_00(q,p)*u_00(q,p)*jb_33(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-4.0*r*r*u_00(q,p)*&
k_00(q,p)*k_00(q,p)*k_00(q,p)*j_04(q,p)*jb_03(q,p)+2.0*jb_00(q,p)*jb_00(q,p)*r*r*u_02(q,p)*u_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+8.0*&
jb_00(q,p)*r*r*jb_02(q,p)*u_00(q,p)*u_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+2.0*jb_00(q,p)*r*r*u_00(q,p)*u_00(q,p)*jb_22(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)+8.0*r*r*ub_00(q,p)*u_02(q,p)*b_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-4.0*r*r*u_02(q,p)*ub_03(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)-16.0*r*r*u_00(q,p)*ub_00(q,p)*b_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_03(q,p)+4.0*jb_00(q,p)*jb_00(q,p)*r*r*&
u_00(q,p)*u_22(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-16.0*r*r*r*b_01(q,p)*ub_00(q,p)*w_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+4.0*j_00(q,p)*&
j_00(q,p)*r*r*ub_00(q,p)*ub_33(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-4.0*r*r*jb_04(q,p)*u_02(q,p)*k_00(q,p)*k_00(q,p)+2.0*r*r*jb_02(q,p)*&
u_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_03(q,p)+8.0*r*r*r*w_03(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)+4.0*&
j_00(q,p)*jb_00(q,p)*jb_04(q,p)*r*r*u_00(q,p)*k_00(q,p)*k_02(q,p)
      s6 = s9+4.0*jb_00(q,p)*jb_00(q,p)*r*r*u_00(q,p)*j_24(q,p)*k_00(q,p)*k_00(q,p)-2.0*jb_00(q,p)*r*r*u_00(q,p)*j_00(q,p)*&
j_00(q,p)*jb_02(q,p)*jb_04(q,p)+2.0*j_00(q,p)*jb_00(q,p)*jb_02(q,p)*r*r*u_00(q,p)*u_03(q,p)*k_00(q,p)*k_00(q,p)+4.0*j_00(q,p)*jb_00(q,p)*&
jb_03(q,p)*r*r*u_00(q,p)*u_00(q,p)*b_02(q,p)*k_00(q,p)*k_00(q,p)-2.0*j_00(q,p)*jb_04(q,p)*r*r*u_00(q,p)*jb_02(q,p)+12.0*r*r*u_00(q,p)*&
jb_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_04(q,p)+4.0*r*r*u_00(q,p)*jb_02(q,p)*k_00(q,p)*k_04(q,p)+2.0*r*r*u_03(q,p)*jb_02(q,p)*u_00(q,p)*&
k_00(q,p)*k_00(q,p)-4.0*r*r*u_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)+4.0*j_00(q,p)*jb_00(q,p)*jb_01(q,p)*r*r*r*w_02(q,p)*&
u_00(q,p)*k_00(q,p)*k_00(q,p)-4.0*j_00(q,p)*jb_00(q,p)*jb_04(q,p)*r*r*u_02(q,p)*k_00(q,p)*k_00(q,p)-4.0*r*r*jb_04(q,p)*u_02(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)*k_00(q,p)+8.0*r*w_23(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2b-4.0*u_00(q,p)*u_00(q,p)*jb_02(q,p)*r*r*b_03(q,p)&
*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-2.0*r*r*jb_03(q,p)*u_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_02(q,p)-6.0*j_00(q,p)*jb_04(q,p)&
*r*r*u_00(q,p)*jb_02(q,p)*k_00(q,p)*k_00(q,p)+4.0*j_00(q,p)*jb_00(q,p)*r*r*u_00(q,p)*jb_02(q,p)*k_00(q,p)*k_04(q,p)-6.0*r*r*u_03(q,p)*&
jb_02(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-4.0*j_00(q,p)*j_00(q,p)*r*r*ub_00(q,p)*ub_00(q,p)*b_03(q,p)*jb_02(q,p)*k_00(q,p)*&
k_00(q,p)+s8
      s9 = s6+j_00(q,p)*r*r*u_00(q,p)*ub_00(q,p)*j_03(q,p)*jb_03(q,p)+4.0*jb_02(q,p)*r*r*u_00(q,p)*ub_00(q,p)*k_00(q,p)*k_02(q,p)&
+6.0*r*r*ub_00(q,p)*u_02(q,p)*jb_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-8.0*r*r*r*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*&
w_13(q,p)+4.0*r*r*r*w_02(q,p)*jb_01(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)-8.0*jb_00(q,p)*r*r*r*w_02(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)*k_01(q,p)-8.0*jb_00(q,p)*r*r*u_00(q,p)*u_00(q,p)*b_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_03(q,p)-4.0*jb_00(q,p)*r*r*u_03(q,p)*&
u_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-8.0*jb_00(q,p)*r*r*u_00(q,p)*u_23(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-4.0*jb_00(q,p)*&
r*r*u_00(q,p)*u_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_02(q,p)-8.0*jb_00(q,p)*r*r*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_02(q,p)*k_04(q,p)+8.0*&
jb_00(q,p)*r*r*u_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_04(q,p)+16.0*r*r*b_04(q,p)*ub_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+8.0*r*r*&
k_00(q,p)*k_00(q,p)*k_00(q,p)*k_04(q,p)*k_04(q,p)-8.0*r*r*u_34(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+2.0*r*r*u_03(q,p)*u_03(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-4.0*r*r*u_03(q,p)*u_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)
      s8 = s9-4.0*r*r*ub_02(q,p)*ub_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+2.0*r*r*ub_02(q,p)*ub_02(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)*k_00(q,p)-8.0*r*r*ub_24(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+4.0*r*r*ub_02(q,p)*jb_02(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)&
*k_00(q,p)*k_00(q,p)-8.0*r*r*r*b_01(q,p)*u_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)-16.0*r*r*u_00(q,p)*b_03(q,p)*k_00(q,p)*k_00(q,p)&
*k_00(q,p)*w_00(q,p)+2.0*j_03(q,p)*r*r*u_00(q,p)*ub_00(q,p)*k_00(q,p)*k_03(q,p)-2.0*r*r*u_00(q,p)*ub_00(q,p)*j_33(q,p)*k_00(q,p)*k_00(q,p)+&
6.0*r*r*u_00(q,p)*ub_03(q,p)*j_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+4.0*j_00(q,p)*r*r*u_03(q,p)*ub_03(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)-16.0*r*r*r*u_00(q,p)*b_13(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)+4.0*r*r*r*u_01(q,p)*w_03(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)+4.0*jb_00(q,p)*r*r*ub_02(q,p)*u_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-2.0*j_03(q,p)*r*r*u_00(q,p)*ub_00(q,p)&
*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_03(q,p)+8.0*r*r*ub_00(q,p)*b_02(q,p)*jb_02(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+4.0*j_00(q,p)&
*jb_00(q,p)*jb_02(q,p)*r*r*u_00(q,p)*ub_00(q,p)*k_00(q,p)*k_02(q,p)-4.0*j_00(q,p)*jb_00(q,p)*jb_02(q,p)*r*r*u_00(q,p)*ub_02(q,p)*k_00(q,p)*&
k_00(q,p)-2.0*j_00(q,p)*jb_00(q,p)*jb_02(q,p)*r*r*ub_00(q,p)*u_02(q,p)*k_00(q,p)*k_00(q,p)-4.0*r*r*ub_02(q,p)*jb_02(q,p)*u_00(q,p)*&
k_00(q,p)*k_00(q,p)
      s9 = s8+4.0*jb_02(q,p)*r*r*u_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_02(q,p)-3.0*j_00(q,p)*r*r*u_00(q,p)*&
ub_00(q,p)*jb_02(q,p)*jb_02(q,p)*k_00(q,p)*k_00(q,p)+4.0*jb_00(q,p)*r*r*u_00(q,p)*ub_22(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+4.0*&
jb_00(q,p)*r*r*ub_00(q,p)*u_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_02(q,p)-4.0*jb_00(q,p)*r*r*u_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*&
k_02(q,p)*k_02(q,p)+8.0*jb_00(q,p)*r*r*u_00(q,p)*ub_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_02(q,p)-2.0*r*r*u_00(q,p)*ub_03(q,p)*j_03(q,p)*&
k_00(q,p)*k_00(q,p)+8.0*r*r*u_00(q,p)*u_03(q,p)*b_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+3.0*j_00(q,p)*r*r*u_00(q,p)*ub_00(q,p)*j_03(q,p)*&
jb_03(q,p)*k_00(q,p)*k_00(q,p)+2.0*j_00(q,p)*j_00(q,p)*r*r*u_00(q,p)*ub_00(q,p)*jb_33(q,p)*k_00(q,p)*k_00(q,p)-2.0*r*r*u_00(q,p)*ub_00(q,p)&
*jb_22(q,p)*k_00(q,p)*k_00(q,p)-j_00(q,p)*r*r*u_00(q,p)*ub_00(q,p)*jb_02(q,p)*jb_02(q,p)-2.0*r*r*ub_00(q,p)*u_02(q,p)*jb_02(q,p)*k_00(q,p)*&
k_00(q,p)-j_00(q,p)*j_00(q,p)*j_00(q,p)*r*r*u_00(q,p)*ub_00(q,p)*jb_03(q,p)*jb_03(q,p)+2.0*jb_00(q,p)*jb_00(q,p)*r*r*u_00(q,p)*ub_00(q,p)*&
j_22(q,p)*k_00(q,p)*k_00(q,p)-8.0*r*r*b_01(q,p)*u_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+4.0*jb_00(q,p)*r*r*ub_00(q,p)*u_22(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)*k_00(q,p)+2.0*r*r*u_00(q,p)*ub_00(q,p)*jb_22(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)
      s7 = s9-jb_00(q,p)*r*r*u_00(q,p)*ub_00(q,p)*j_00(q,p)*j_00(q,p)*jb_02(q,p)*jb_02(q,p)+6.0*j_00(q,p)*ub_00(q,p)*r*r*&
j_03(q,p)*jb_04(q,p)*k_00(q,p)*k_00(q,p)-2.0*j_00(q,p)*j_00(q,p)*r*r*ub_00(q,p)*ub_02(q,p)*jb_03(q,p)*k_00(q,p)*k_00(q,p)-j_00(q,p)*r*r*&
j_02(q,p)*ub_00(q,p)*ub_00(q,p)*jb_03(q,p)*k_00(q,p)*k_00(q,p)+2.0*ub_00(q,p)*r*r*j_00(q,p)*j_03(q,p)*jb_04(q,p)+j_00(q,p)*jb_02(q,p)*r*r*&
j_03(q,p)*ub_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)+4.0*j_00(q,p)*j_00(q,p)*r*r*ub_03(q,p)*jb_04(q,p)*k_00(q,p)*k_00(q,p)+4.0*j_00(q,p)*&
j_00(q,p)*ub_00(q,p)*r*r*jb_34(q,p)*k_00(q,p)*k_00(q,p)+4.0*j_00(q,p)*j_00(q,p)*r*r*ub_00(q,p)*u_03(q,p)*jb_03(q,p)*k_00(q,p)*k_00(q,p)+2.0&
*j_00(q,p)*j_00(q,p)*r*r*ub_03(q,p)*u_00(q,p)*jb_03(q,p)*k_00(q,p)*k_00(q,p)+4.0*j_00(q,p)*r*r*ub_00(q,p)*u_33(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)+8.0*r*r*u_00(q,p)*j_03(q,p)*ub_00(q,p)*b_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-4.0*r*r*ub_00(q,p)*u_03(q,p)*&
j_03(q,p)*k_00(q,p)*k_00(q,p)+4.0*j_00(q,p)*r*r*u_00(q,p)*ub_33(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+4.0*r*r*ub_00(q,p)*u_03(q,p)*&
j_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+2.0*r*r*u_00(q,p)*ub_00(q,p)*j_33(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-4.0*j_00(q,p)*&
j_00(q,p)*r*r*r*w_03(q,p)*ub_00(q,p)*jb_01(q,p)*k_00(q,p)*k_00(q,p)-2.0*ub_00(q,p)*r*r*j_00(q,p)*j_00(q,p)*j_00(q,p)*jb_03(q,p)*jb_04(q,p)+&
4.0*r*r*r*w_03(q,p)*j_01(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)
      s9 = s7+4.0*r*r*r*w_03(q,p)*j_01(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)-4.0*j_00(q,p)*r*r*ub_02(q,p)*ub_03(q,p)*k_00(q,p)&
*k_00(q,p)*k_00(q,p)*k_00(q,p)-8.0*r*r*ub_00(q,p)*ub_03(q,p)*j_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-4.0*r*r*r*u_03(q,p)*&
k_00(q,p)*k_00(q,p)*k_00(q,p)*w_01(q,p)-8.0*r*r*r*u_01(q,p)*b_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)+4.0*r*r*j_02(q,p)*ub_00(q,p)&
*ub_00(q,p)*b_03(q,p)*k_00(q,p)*k_00(q,p)+4.0*ub_00(q,p)*r*r*j_04(q,p)*k_00(q,p)*k_03(q,p)-4.0*ub_00(q,p)*r*r*j_04(q,p)*k_00(q,p)*k_00(q,p)&
*k_00(q,p)*k_03(q,p)-8.0*j_00(q,p)*r*r*ub_00(q,p)*ub_23(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-4.0*r*r*j_03(q,p)*ub_00(q,p)*&
ub_00(q,p)*b_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-6.0*r*r*ub_00(q,p)*ub_02(q,p)*j_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+2.0&
*r*r*ub_00(q,p)*ub_02(q,p)*j_03(q,p)*k_00(q,p)*k_00(q,p)+4.0*r*r*u_00(q,p)*jb_04(q,p)*k_00(q,p)*k_02(q,p)+4.0*r*r*jb_03(q,p)*u_00(q,p)&
*u_00(q,p)*b_02(q,p)*k_00(q,p)*k_00(q,p)-4.0*r*r*u_00(q,p)*jb_04(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_02(q,p)+4.0*j_00(q,p)*r*r*ub_00(q,p)*&
ub_00(q,p)*b_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*jb_02(q,p)-8.0*r*r*ub_00(q,p)*u_23(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)
      s8 = s9+2.0*j_00(q,p)*j_00(q,p)*r*r*ub_03(q,p)*ub_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+8.0*j_00(q,p)*r*r*r*jb_01(q,p)&
*ub_00(q,p)*w_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+2.0*j_00(q,p)*jb_02(q,p)*r*r*ub_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_02(q,p)-4.0*r*&
r*u_03(q,p)*ub_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-12.0*r*r*u_00(q,p)*ub_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*&
k_02(q,p)-8.0*r*r*ub_00(q,p)*u_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_02(q,p)-4.0*r*r*ub_00(q,p)*u_02(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)*k_03(q,p)+4.0*r*r*u_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_02(q,p)*k_03(q,p)-8.0*r*r*u_00(q,p)*ub_02(q,p)*&
k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_03(q,p)+r*r*j_02(q,p)*ub_00(q,p)*ub_00(q,p)*jb_02(q,p)*k_00(q,p)+8.0*r*r*u_00(q,p)*ub_03(q,p)*&
b_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+r*r*j_02(q,p)*ub_00(q,p)*ub_00(q,p)*jb_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-4.0*j_00(q,p)*jb_02(q,p)*r*&
r*ub_00(q,p)*ub_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-4.0*r*r*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*jb_04(q,p)*j_02(q,p)+4.0*jb_00(q,p)*&
r*r*ub_00(q,p)*j_24(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+2.0*jb_00(q,p)*r*r*ub_00(q,p)*ub_00(q,p)*j_22(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+4.0&
*r*r*ub_00(q,p)*j_04(q,p)*jb_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-j_00(q,p)*j_00(q,p)*r*r*ub_00(q,p)*ub_00(q,p)*k_00(q,p)*jb_02(q,p)*jb_02(q,p)&
-8.0*r*r*u_00(q,p)*ub_23(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)
      s9 = s8+4.0*r*r*r*w_02(q,p)*jb_01(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-8.0*r*r*u_00(q,p)*u_02(q,p)*&
jb_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-8.0*r*r*r*ub_01(q,p)*b_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)-16.0*r*r*&
ub_00(q,p)*b_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)-4.0*r*r*r*ub_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_01(q,p)-8.0*r*r*r*&
w_02(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)-8.0*r*r*r*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_12(q,p)-4.0*r*r&
*ub_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)-8.0*r*r*b_01(q,p)*ub_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+8.0*r*r*ub_00(q,p)*&
ub_02(q,p)*b_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+4.0*r*r*r*ub_01(q,p)*w_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-16.0*r*r*r*ub_00(q,p)*&
b_12(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)+2.0*j_00(q,p)*r*r*ub_00(q,p)*ub_00(q,p)*j_33(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-16.0*r*r*&
ub_00(q,p)*b_02(q,p)*u_00(q,p)*b_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-2.0*jb_00(q,p)*r*r*r*r*u_00(q,p)*u_00(q,p)*u_02(q,p)*jb_01(q,p)*k_00(q,p)&
*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-4.0*r*r*r*r*jb_04(q,p)*u_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi-16.0*r*r*r&
*r*u_00(q,p)*b_04(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+4.0*jb_01(q,p)*r*r*r*r*u_00(q,p)*u_00(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi
      s5 = s9+16.0*r*r*r*ub_00(q,p)*u_03(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-4.0*jb_00(q,p)*r*r&
*r*r*r*u_01(q,p)*u_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*e2bi-16.0*r*r*r*r*ub_00(q,p)*b_04(q,p)*u_01(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+16.0*j_00(q,p)*r*r*r*r*r*ub_00(q,p)*b_01(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*&
e2bi+4.0*jb_00(q,p)*jb_00(q,p)*r*r*r*r*u_00(q,p)*u_00(q,p)*u_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi+8.0*r*r*r*&
u_00(q,p)*u_00(q,p)*ub_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-8.0*r*r*r*r*ub_00(q,p)*u_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*&
k_01(q,p)*k_03(q,p)*e2bi-2.0*jb_01(q,p)*r*r*r*r*u_00(q,p)*u_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_02(q,p)*e2bi+8.0*r*r*r*&
r*ub_00(q,p)*u_14(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-8.0*r*r*r*b_01(q,p)*ub_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)+&
8.0*r*r*ub_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_04(q,p)+8.0*j_00(q,p)*r*r*j_03(q,p)*ub_00(q,p)*ub_03(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)-4.0*r*r*ub_00(q,p)*ub_22(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+16.0*r*r*ub_00(q,p)*b_24(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)+4.0&
*jb_00(q,p)*r*r*jb_02(q,p)*u_00(q,p)*u_00(q,p)*b_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-j_00(q,p)*jb_01(q,p)*jb_01(q,p)*r*r*r*r*u_00(q,p)*&
u_00(q,p)*k_00(q,p)*e2bi+2.0*j_00(q,p)*jb_04(q,p)*r*r*r*r*u_00(q,p)*u_00(q,p)*jb_01(q,p)*k_00(q,p)*e2bi-16.0*jb_00(q,p)*r*r*&
r*r*u_00(q,p)*b_04(q,p)*u_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-u_00(q,p)*u_00(q,p)*r*r*r*r*k_00(q,p)*j_00(q,p)*j_00(q,p)*jb_00(q,p)*&
jb_01(q,p)*jb_01(q,p)*e2bi+2.0*j_00(q,p)*j_00(q,p)*j_00(q,p)*r*r*r*r*ub_00(q,p)*ub_00(q,p)*k_00(q,p)*jb_01(q,p)*jb_04(q,p)*e2bi
      s9 = s5-4.0*r*r*r*r*u_00(q,p)*ub_00(q,p)*j_01(q,p)*ub_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+&
4.0*r*r*r*r*u_00(q,p)*u_01(q,p)*ub_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+4.0*r*r*r*r*ub_00(q,p)*u_00(q,p)*&
u_00(q,p)*jb_12(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-4.0*r*r*r*r*jb_04(q,p)*u_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*&
k_01(q,p)*e2bi+4.0*jb_01(q,p)*u_00(q,p)*u_00(q,p)*r*r*r*r*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi+8.0*jb_00(q,p)*r*r*r*r*&
u_00(q,p)*u_00(q,p)*u_13(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-16.0*jb_00(q,p)*r*r*r*u_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*&
w_00(q,p)*e2bi+8.0*r*r*r*r*ub_00(q,p)*ub_02(q,p)*u_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-8.0*jb_00(q,p)*r*r*r&
*r*u_00(q,p)*u_00(q,p)*u_01(q,p)*b_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-8.0*j_00(q,p)*r*r*r*r*ub_00(q,p)*ub_00(q,p)*ub_01(q,p)*&
b_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-12.0*jb_01(q,p)*r*r*r*r*u_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*e2bi+&
4.0*j_00(q,p)*jb_00(q,p)*jb_01(q,p)*u_00(q,p)*u_00(q,p)*r*r*r*r*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi-2.0*r*r*r*r*u_00(q,p)*&
u_00(q,p)*u_00(q,p)*jb_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_03(q,p)*e2bi+4.0*r*r*r*u_00(q,p)*u_00(q,p)*u_00(q,p)*jb_03(q,p)*&
k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-j_00(q,p)*jb_01(q,p)*jb_01(q,p)*r*r*r*r*r*u_00(q,p)*u_00(q,p)*k_00(q,p)*w_00(q,p)*e2bi+4.0*&
jb_01(q,p)*u_00(q,p)*u_00(q,p)*r*r*r*r*r*k_00(q,p)*k_00(q,p)*k_01(q,p)*w_00(q,p)*e2bi-2.0*j_00(q,p)*jb_00(q,p)*jb_01(q,p)*r*r*r*r&
*u_00(q,p)*u_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_03(q,p)*e2bi
      s8 = s9-4.0*jb_01(q,p)*r*r*r*r*r*u_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_01(q,p)*e2bi+8.0*r*r*r*&
r*jb_14(q,p)*u_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+4.0*r*r*r*r*ub_00(q,p)*u_00(q,p)*jb_01(q,p)*u_02(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)*e2bi+4.0*jb_00(q,p)*r*r*r*u_00(q,p)*u_00(q,p)*ub_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-3.0*j_00(q,p)*&
jb_01(q,p)*jb_01(q,p)*r*r*r*r*r*u_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*e2bi+8.0*jb_00(q,p)*r*r*r*r*u_00(q,p)&
*u_14(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+6.0*r*r*r*r*jb_03(q,p)*u_00(q,p)*u_00(q,p)*u_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-&
u_00(q,p)*u_00(q,p)*r*r*r*r*r*k_00(q,p)*j_00(q,p)*j_00(q,p)*jb_00(q,p)*jb_01(q,p)*jb_01(q,p)*w_00(q,p)*e2bi-2.0*j_00(q,p)*jb_00(q,p)*&
jb_03(q,p)*r*r*r*r*u_00(q,p)*u_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi+4.0*j_00(q,p)*jb_00(q,p)*jb_01(q,p)*u_00(q,p)*u_00(q,p)&
*r*r*r*r*r*k_00(q,p)*k_00(q,p)*k_01(q,p)*w_00(q,p)*e2bi-4.0*j_00(q,p)*r*r*r*r*jb_01(q,p)*j_04(q,p)*ub_00(q,p)*ub_00(q,p)*&
k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-2.0*j_00(q,p)*r*r*r*r*jb_01(q,p)*j_02(q,p)*ub_00(q,p)*ub_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)*e2bi+4.0*j_00(q,p)*ub_00(q,p)*ub_00(q,p)*j_03(q,p)*r*r*b_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-2.0*r*r*ub_00(q,p)*&
ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_02(q,p)*k_02(q,p)-4.0*jb_01(q,p)*r*r*r*r*u_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_04(q,p)*e2bi+&
8.0*jb_00(q,p)*r*r*r*r*u_00(q,p)*u_01(q,p)*u_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+2.0*j_00(q,p)*jb_00(q,p)*jb_01(q,p)*r*r*r&
*r*u_00(q,p)*u_00(q,p)*u_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-4.0*r*r*r*r*u_00(q,p)*u_00(q,p)*jb_01(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)*k_04(q,p)*e2bi+4.0*r*r*r*r*u_00(q,p)*ub_00(q,p)*ub_00(q,p)*j_13(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi
      s9 = s8+4.0*jb_01(q,p)*r*r*r*r*r*u_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*w_00(q,p)*e2bi&
-8.0*r*r*r*r*u_00(q,p)*ub_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*k_02(q,p)*e2bi-4.0*r*r*r*r*u_00(q,p)*u_00(q,p)*&
jb_11(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-2.0*r*r*r*r*jb_03(q,p)*u_00(q,p)*u_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi&
-2.0*jb_01(q,p)*r*r*r*r*u_00(q,p)*u_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_03(q,p)*e2bi-8.0*jb_01(q,p)*r*r*r*r*r*u_00(q,p)*&
u_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*e2bi+16.0*jb_00(q,p)*r*r*r*r*r*u_00(q,p)*b_01(q,p)*u_01(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)*w_00(q,p)*e2bi-32.0*jb_00(q,p)*r*r*r*r*u_00(q,p)*u_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*e2bi+2.0*r*r*r&
*r*ub_01(q,p)*jb_02(q,p)*u_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+4.0*r*r*r*r*ub_00(q,p)*u_01(q,p)*u_03(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-2.0*r*r*r*r*jb_03(q,p)*u_00(q,p)*u_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*&
e2bi-2.0*j_00(q,p)*jb_03(q,p)*r*r*u_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_02(q,p)+12.0*jb_00(q,p)*r*r*r*u_00(q,p)*u_00(q,p)*&
u_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-8.0*j_00(q,p)*jb_01(q,p)*r*r*r*r*r*u_00(q,p)*ub_00(q,p)*k_00(q,p)*k_01(q,p)*w_00(q,p)*&
e2bi+16.0*r*r*r*r*u_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*k_04(q,p)*e2bi-2.0*r*r*r*r*&
u_00(q,p)*ub_00(q,p)*ub_00(q,p)*j_01(q,p)*jb_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-4.0*j_00(q,p)*jb_00(q,p)*jb_01(q,p)*r*r*r&
*r*u_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_04(q,p)*e2bi-2.0*r*r*r*r*u_00(q,p)*ub_00(q,p)*ub_00(q,p)*j_02(q,p)*jb_01(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi
      s7 = s9-8.0*r*r*r*r*r*u_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*w_01(q,p)*e2bi-24.0*r*r*r&
*r*u_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*w_00(q,p)*e2bi-4.0*r*r*r*r*u_00(q,p)*ub_00(q,p)*j_04(q,p)*jb_01(q,p)*&
k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+4.0*jb_00(q,p)*r*r*r*r*u_00(q,p)*u_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*&
k_03(q,p)*e2bi-2.0*j_00(q,p)*jb_00(q,p)*jb_00(q,p)*jb_01(q,p)*r*r*r*r*u_00(q,p)*u_00(q,p)*u_02(q,p)*k_00(q,p)*k_00(q,p)*e2bi-8.0*&
j_00(q,p)*jb_01(q,p)*r*r*r*r*u_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_04(q,p)*e2bi-8.0*jb_00(q,p)*r*r*r*r*r*u_00(q,p)&
*u_11(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*e2bi+8.0*j_00(q,p)*r*r*r*r*u_00(q,p)*u_03(q,p)*ub_00(q,p)*jb_01(q,p)*k_00(q,p)*k_00(q,p)&
*k_00(q,p)*k_00(q,p)*e2bi+16.0*r*r*r*u_00(q,p)*ub_00(q,p)*ub_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-8.0*r*r*r*&
r*ub_00(q,p)*ub_00(q,p)*u_01(q,p)*b_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-8.0*jb_00(q,p)*r*r*r*r*ub_00(q,p)*u_00(q,p)*&
u_01(q,p)*b_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-4.0*r*r*u_00(q,p)*u_00(q,p)*jb_23(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-4.0*r&
*r*u_00(q,p)*jb_24(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-4.0*r*r*u_00(q,p)*jb_24(q,p)*k_00(q,p)*k_00(q,p)-4.0*r*r*ub_00(q,p)*&
ub_00(q,p)*j_23(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-4.0*r*r*ub_00(q,p)*j_34(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-4.0*r*r*&
j_04(q,p)*ub_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-4.0*r*r*j_04(q,p)*ub_03(q,p)*k_00(q,p)*k_00(q,p)-4.0*r*r*ub_00(q,p)*j_34(q,p)*&
k_00(q,p)*k_00(q,p)
      s9 = s7-2.0*r*r*r*r*ub_00(q,p)*u_00(q,p)*u_00(q,p)*jb_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_02(q,p)*e2bi+&
2.0*j_00(q,p)*j_00(q,p)*j_00(q,p)*r*r*r*r*ub_00(q,p)*ub_00(q,p)*ub_03(q,p)*jb_01(q,p)*k_00(q,p)*k_00(q,p)*e2bi+8.0*r*r*r*r*&
u_00(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_04(q,p)*e2bi+6.0*j_00(q,p)*r*r*r*r*u_00(q,p)*u_00(q,p)*jb_01(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)*jb_04(q,p)*e2bi-4.0*r*r*r*r*r*u_00(q,p)*u_00(q,p)*jb_11(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*e2bi-4.0*&
jb_00(q,p)*u_00(q,p)*u_00(q,p)*r*r*r*r*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*k_01(q,p)*e2bi-8.0*r*r*r*r*r*ub_00(q,p)*u_01(q,p)*&
k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*w_00(q,p)*e2bi-4.0*j_00(q,p)*jb_00(q,p)*jb_04(q,p)*r*r*r*r*u_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)&
*k_01(q,p)*e2bi+4.0*r*r*r*r*u_00(q,p)*u_00(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_03(q,p)*e2bi-2.0*jb_00(q,p)*r*r*r&
*r*u_00(q,p)*u_00(q,p)*u_02(q,p)*jb_01(q,p)*k_00(q,p)*k_00(q,p)*e2bi-16.0*r*r*r*r*u_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*&
k_04(q,p)*e2bi+4.0*r*r*r*r*u_00(q,p)*ub_00(q,p)*j_04(q,p)*jb_01(q,p)*e2bi-4.0*jb_00(q,p)*r*r*r*r*u_00(q,p)*ub_00(q,p)*&
j_11(q,p)*k_00(q,p)*k_00(q,p)*e2bi+8.0*r*r*r*r*u_00(q,p)*ub_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*k_02(q,p)*&
e2bi-2.0*r*r*r*r*ub_00(q,p)*u_00(q,p)*u_00(q,p)*j_03(q,p)*jb_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-4.0*j_00(q,p)*&
j_00(q,p)*r*r*r*r*ub_00(q,p)*u_00(q,p)*u_00(q,p)*jb_01(q,p)*jb_03(q,p)*e2bi+4.0*r*r*r*r*ub_00(q,p)*ub_00(q,p)*u_02(q,p)*k_00(q,p)&
*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi
      s8 = s9+2.0*r*r*r*r*ub_00(q,p)*ub_00(q,p)*j_01(q,p)*ub_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi&
-4.0*j_00(q,p)*j_00(q,p)*r*r*r*r*u_00(q,p)*ub_00(q,p)*ub_00(q,p)*jb_01(q,p)*jb_02(q,p)*e2bi-2.0*r*r*r*r*ub_00(q,p)*&
ub_00(q,p)*u_03(q,p)*j_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+8.0*r*r*r*r*ub_00(q,p)*u_00(q,p)*ub_12(q,p)*k_00(q,p)&
*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+4.0*r*r*r*r*u_00(q,p)*ub_00(q,p)*ub_01(q,p)*j_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-32.0&
*r*r*r*r*u_01(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*e2bi-4.0*j_00(q,p)*jb_02(q,p)*r*r*r*r*u_00(q,p)*&
ub_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi+4.0*j_00(q,p)*jb_02(q,p)*r*r*r*r*u_00(q,p)*ub_00(q,p)*ub_00(q,p)*&
k_00(q,p)*k_01(q,p)*e2bi-4.0*j_00(q,p)*r*r*r*r*u_00(q,p)*ub_00(q,p)*jb_11(q,p)*k_00(q,p)*k_00(q,p)*e2bi+4.0*r*r*r*r*&
u_00(q,p)*ub_00(q,p)*j_01(q,p)*ub_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-4.0*j_01(q,p)*jb_01(q,p)*r*r*r*r*r*u_00(q,p)*ub_00(q,p)*&
w_00(q,p)*e2bi-16.0*r*r*u_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+8.0*j_00(q,p)*r*r*r*u_00(q,p)*&
ub_00(q,p)*ub_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+4.0*j_00(q,p)*r*r*r*r*ub_00(q,p)*ub_01(q,p)*u_03(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)*e2bi-4.0*r*r*r*r*u_00(q,p)*ub_00(q,p)*j_01(q,p)*jb_04(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+4.0*r*r*r&
*u_00(q,p)*ub_00(q,p)*ub_00(q,p)*j_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-4.0*j_00(q,p)*r*r*ub_00(q,p)*u_02(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)*jb_03(q,p)+8.0*jb_00(q,p)*r*r*r*r*u_00(q,p)*ub_00(q,p)*j_14(q,p)*k_00(q,p)*k_00(q,p)*e2bi+4.0*r*r*r*r*u_00(q,p)*&
ub_00(q,p)*u_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi
      s9 = s8+4.0*j_00(q,p)*r*r*r*r*u_00(q,p)*ub_01(q,p)*ub_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+4.0*j_00(q,p)*&
j_00(q,p)*jb_01(q,p)*jb_01(q,p)*r*r*r*r*r*u_00(q,p)*ub_00(q,p)*w_00(q,p)*e2bi+8.0*r*r*r*r*ub_00(q,p)*ub_00(q,p)*u_12(q,p)*&
k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+4.0*j_00(q,p)*j_00(q,p)*r*r*r*r*ub_00(q,p)*u_00(q,p)*u_00(q,p)*jb_01(q,p)*jb_03(q,p)*&
k_00(q,p)*k_00(q,p)*e2bi+8.0*j_00(q,p)*jb_01(q,p)*r*r*r*r*u_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi-4.0*&
jb_00(q,p)*r*r*r*r*r*u_00(q,p)*ub_00(q,p)*j_11(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*e2bi+8.0*r*r*r*r*u_01(q,p)*ub_00(q,p)*k_00(q,p)&
*k_00(q,p)*k_00(q,p)*k_04(q,p)*e2bi+4.0*jb_00(q,p)*r*r*r*r*ub_00(q,p)*u_01(q,p)*u_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-8.0*&
j_00(q,p)*r*r*r*r*u_00(q,p)*ub_00(q,p)*ub_01(q,p)*b_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+3.0*j_00(q,p)*j_01(q,p)*jb_01(q,p)*r*r*&
r*r*ub_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-8.0*r*r*r*r*u_00(q,p)*u_00(q,p)*ub_01(q,p)*b_03(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)*e2bi+8.0*r*r*r*r*u_00(q,p)*ub_00(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_02(q,p)*e2bi-8.0*r*r*r&
*r*r*u_00(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*w_00(q,p)*e2bi-8.0*r*r*r*r*u_00(q,p)*ub_00(q,p)*u_03(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi+16.0*r*r*r*u_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_04(q,p)*e2bi+&
8.0*r*r*r*u_00(q,p)*ub_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_02(q,p)*e2bi+4.0*r*r*r*r*ub_00(q,p)*ub_00(q,p)*&
u_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_02(q,p)*e2bi+4.0*j_01(q,p)*jb_01(q,p)*r*r*r*r*u_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)*e2bi
      s10 = s9+4.0*j_00(q,p)*r*r*r*r*u_00(q,p)*u_00(q,p)*ub_03(q,p)*jb_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+&
4.0*j_00(q,p)*jb_03(q,p)*r*r*r*r*ub_00(q,p)*u_00(q,p)*u_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi+8.0*r*r*r*u_00(q,p)*u_00(q,p)*ub_00(q,p)&
*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_03(q,p)*e2bi+8.0*r*r*r*r*u_00(q,p)*u_03(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+&
4.0*r*r*r*r*u_00(q,p)*ub_00(q,p)*ub_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi+4.0*j_00(q,p)*r*r*r*r*u_00(q,p)*&
ub_00(q,p)*ub_00(q,p)*jb_12(q,p)*k_00(q,p)*k_00(q,p)*e2bi-16.0*j_00(q,p)*r*r*r*ub_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*&
w_00(q,p)*e2bi-4.0*r*r*r*r*ub_00(q,p)*u_00(q,p)*jb_01(q,p)*u_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-2.0*&
j_00(q,p)*r*r*r*r*ub_00(q,p)*ub_00(q,p)*ub_03(q,p)*j_01(q,p)*k_00(q,p)*k_00(q,p)*e2bi
      s6 = s10-2.0*j_00(q,p)*r*r*r*r*ub_00(q,p)*ub_00(q,p)*ub_03(q,p)*j_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi&
+8.0*k_00(q,p)*k_00(q,p)*k_01(q,p)*k_01(q,p)*r*r*r*r*r*u_00(q,p)*ub_00(q,p)*w_00(q,p)*e2bi-4.0*j_00(q,p)*r*r*r*r*r*&
ub_01(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*e2bi+j_00(q,p)*j_01(q,p)*jb_01(q,p)*r*r*r*r*r*ub_00(q,p)*ub_00(q,p)*k_00(q,p)&
*w_00(q,p)*e2bi+4.0*r*r*r*r*ub_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*j_04(q,p)*k_01(q,p)*e2bi+2.0*r*r*r&
*r*j_02(q,p)*ub_00(q,p)*ub_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi+8.0*jb_00(q,p)*r*r*r*r*ub_00(q,p)&
*u_00(q,p)*u_12(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-32.0*j_00(q,p)*r*r*r*r*ub_00(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*&
w_00(q,p)*e2bi-8.0*jb_00(q,p)*r*r*r*r*u_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_01(q,p)*e2bi+8.0*r*r*r*r*u_00(q,p)&
*u_00(q,p)*ub_13(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-8.0*j_00(q,p)*r*r*r*r*r*ub_00(q,p)*ub_11(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)*w_00(q,p)*e2bi
      s9 = s6+j_00(q,p)*j_01(q,p)*jb_01(q,p)*r*r*r*r*ub_00(q,p)*ub_00(q,p)*k_00(q,p)*e2bi+8.0*r*r*r*r*ub_00(q,p)*&
u_00(q,p)*u_13(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-8.0*r*r*r*r*u_00(q,p)*ub_00(q,p)*u_01(q,p)*b_03(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)*e2bi+r*r*r*r*u_00(q,p)*u_00(q,p)*u_00(q,p)*jb_01(q,p)*k_00(q,p)*j_00(q,p)*j_00(q,p)*jb_00(q,p)*jb_03(q,p)*e2bi&
-8.0*j_01(q,p)*r*r*r*ub_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-8.0*j_00(q,p)*r*r*r*r*ub_00(q,p)*ub_00(q,p)*&
k_00(q,p)*k_00(q,p)*k_00(q,p)*w_01(q,p)*e2bi+3.0*j_00(q,p)*j_01(q,p)*jb_01(q,p)*r*r*r*r*r*ub_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)*w_00(q,p)*e2bi+j_00(q,p)*jb_03(q,p)*r*r*r*r*u_00(q,p)*u_00(q,p)*u_00(q,p)*jb_01(q,p)*k_00(q,p)*e2bi-16.0*r*r*r*&
u_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi-8.0*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*k_01(q,p)*r*r*r*r*&
u_00(q,p)*ub_00(q,p)*e2bi+8.0*k_00(q,p)*k_00(q,p)*k_01(q,p)*k_01(q,p)*r*r*r*r*u_00(q,p)*ub_00(q,p)*e2bi+2.0*j_00(q,p)*j_00(q,p)*&
jb_03(q,p)*r*r*u_00(q,p)*ub_00(q,p)*k_00(q,p)*jb_02(q,p)-4.0*r*r*r*r*u_00(q,p)*u_00(q,p)*ub_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)&
*k_00(q,p)*k_01(q,p)*e2bi+16.0*r*r*r*r*r*ub_00(q,p)*b_01(q,p)*u_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*e2bi&
-32.0*r*r*r*u_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*e2bi-16.0*r*r*r*r*u_00(q,p)*ub_00(q,p)*&
k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_01(q,p)*e2bi-8.0*r*r*r*r*r*u_01(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*&
w_00(q,p)*e2bi
      s10 = s9-8.0*r*r*r*r*r*ub_00(q,p)*u_11(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*e2bi-8.0*k_00(q,p)*&
k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*k_01(q,p)*r*r*r*r*r*u_00(q,p)*ub_00(q,p)*w_00(q,p)*e2bi-4.0*j_00(q,p)*j_00(q,p)*jb_01(q,p)*&
jb_01(q,p)*r*r*r*r*r*u_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*e2bi+4.0*j_01(q,p)*jb_01(q,p)*r*r*r*r*r*u_00(q,p)*&
ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*e2bi+8.0*r*r*r*r*u_00(q,p)*ub_00(q,p)*ub_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)&
*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi+16.0*r*r*r*r*r*u_00(q,p)*b_01(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*&
e2bi+r*r*u_00(q,p)*j_03(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*jb_02(q,p)+4.0*jb_00(q,p)*r*r*r*r*ub_00(q,p)*u_00(q,p)*u_00(q,p)&
*j_13(q,p)*k_00(q,p)*k_00(q,p)*e2bi+8.0*j_00(q,p)*j_00(q,p)*r*r*r*r*ub_00(q,p)*u_00(q,p)*jb_01(q,p)*jb_04(q,p)*k_00(q,p)*k_00(q,p)*&
e2bi
      s8 = s10+4.0*r*r*r*r*ub_00(q,p)*ub_00(q,p)*u_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi+4.0*j_00(q,p)*&
j_00(q,p)*r*r*r*r*u_00(q,p)*ub_00(q,p)*ub_00(q,p)*jb_01(q,p)*jb_02(q,p)*k_00(q,p)*k_00(q,p)*e2bi-8.0*j_00(q,p)*j_00(q,p)*r*r*r*r*&
ub_00(q,p)*u_00(q,p)*jb_01(q,p)*jb_04(q,p)*e2bi+8.0*j_00(q,p)*r*r*r*r*u_00(q,p)*ub_00(q,p)*jb_14(q,p)*k_00(q,p)*k_00(q,p)*e2bi&
-4.0*j_00(q,p)*jb_01(q,p)*r*r*r*r*ub_00(q,p)*u_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_03(q,p)*e2bi+2.0*r*r*r*r*&
u_00(q,p)*ub_00(q,p)*ub_00(q,p)*j_01(q,p)*jb_02(q,p)*e2bi+8.0*j_00(q,p)*r*r*r*r*u_00(q,p)*ub_00(q,p)*jb_01(q,p)*k_00(q,p)*k_04(q,p)*&
e2bi+4.0*j_00(q,p)*r*r*r*r*ub_00(q,p)*u_00(q,p)*u_00(q,p)*jb_13(q,p)*k_00(q,p)*k_00(q,p)*e2bi-2.0*r*r*r*r*ub_00(q,p)*&
u_00(q,p)*u_00(q,p)*j_01(q,p)*jb_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-4.0*j_00(q,p)*r*r*r*r*r*u_00(q,p)*ub_00(q,p)*&
jb_11(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*e2bi
      s9 = s8+4.0*r*r*r*r*u_00(q,p)*u_00(q,p)*ub_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi-4.0*j_00(q,p)*&
jb_03(q,p)*r*r*r*r*ub_00(q,p)*u_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi+2.0*r*r*r*r*ub_00(q,p)*u_00(q,p)*&
u_00(q,p)*j_01(q,p)*jb_03(q,p)*e2bi+4.0*j_00(q,p)*r*r*r*r*ub_00(q,p)*u_00(q,p)*u_00(q,p)*jb_01(q,p)*k_00(q,p)*k_03(q,p)*e2bi+8.0*&
j_00(q,p)*jb_01(q,p)*r*r*r*r*r*u_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*w_00(q,p)*e2bi-8.0*j_00(q,p)*r*r*r*r*&
jb_01(q,p)*ub_00(q,p)*u_00(q,p)*ub_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-4.0*j_00(q,p)*r*r*r*r*jb_01(q,p)*ub_00(q,p)*&
ub_00(q,p)*u_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+8.0*j_00(q,p)*jb_04(q,p)*r*r*r*r*u_00(q,p)*ub_00(q,p)*k_00(q,p)*k_01(q,p)&
*e2bi-8.0*r*r*r*r*r*u_00(q,p)*ub_11(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*e2bi+4.0*j_00(q,p)*r*r*r*r*&
u_00(q,p)*ub_00(q,p)*ub_00(q,p)*jb_01(q,p)*k_00(q,p)*k_02(q,p)*e2bi+4.0*j_00(q,p)*j_00(q,p)*jb_01(q,p)*jb_01(q,p)*r*r*r*r*u_00(q,p)*&
ub_00(q,p)*e2bi-4.0*j_01(q,p)*jb_01(q,p)*r*r*r*r*u_00(q,p)*ub_00(q,p)*e2bi-4.0*j_00(q,p)*j_00(q,p)*jb_01(q,p)*jb_01(q,p)*r*r&
*r*r*u_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-4.0*j_00(q,p)*jb_01(q,p)*r*r*r*r*u_00(q,p)*ub_00(q,p)*ub_00(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)*k_02(q,p)*e2bi-8.0*j_00(q,p)*jb_01(q,p)*r*r*r*r*u_00(q,p)*ub_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi+4.0*jb_00(q,p)&
*r*r*r*r*u_00(q,p)*ub_00(q,p)*ub_00(q,p)*j_12(q,p)*k_00(q,p)*k_00(q,p)*e2bi+2.0*r*r*r*r*u_00(q,p)*ub_00(q,p)*ub_00(q,p)*j_02(q,p)&
*jb_01(q,p)*e2bi-32.0*r*r*r*r*u_00(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*e2bi
      s10 = s9-8.0*j_00(q,p)*jb_04(q,p)*r*r*r*r*u_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi+8.0*r&
*r*r*r*u_00(q,p)*ub_14(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+4.0*r*r*r*r*u_00(q,p)*ub_00(q,p)*j_01(q,p)*jb_04(q,p)*&
e2bi+2.0*r*r*r*r*ub_00(q,p)*u_00(q,p)*u_00(q,p)*j_03(q,p)*jb_01(q,p)*e2bi-8.0*jb_00(q,p)*b_22(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*&
w_00(q,p)*r*e2b-j_00(q,p)*jb_02(q,p)*r*r*r*r*ub_00(q,p)*ub_00(q,p)*ub_00(q,p)*j_01(q,p)*k_00(q,p)*e2bi+4.0*jb_00(q,p)*r*r*r*&
r*u_00(q,p)*u_01(q,p)*ub_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-3.0*j_00(q,p)*jb_01(q,p)*jb_01(q,p)*r*r*r*r*u_00(q,p)*u_00(q,p)*&
k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+8.0*r*r*r*r*jb_04(q,p)*u_00(q,p)*u_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi
      s7 = s10-2.0*j_00(q,p)*r*r*r*r*jb_01(q,p)*u_00(q,p)*j_03(q,p)*ub_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi&
-8.0*r*r*r*r*ub_00(q,p)*u_00(q,p)*ub_01(q,p)*b_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-j_00(q,p)*j_00(q,p)*j_00(q,p)*&
ub_00(q,p)*ub_00(q,p)*r*r*r*r*k_00(q,p)*jb_01(q,p)*jb_01(q,p)*e2bi+4.0*r*r*r*r*ub_00(q,p)*ub_01(q,p)*u_02(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-2.0*j_00(q,p)*jb_00(q,p)*jb_01(q,p)*r*r*r*r*u_00(q,p)*u_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*&
k_02(q,p)*e2bi-4.0*j_00(q,p)*jb_00(q,p)*jb_01(q,p)*r*r*r*r*u_00(q,p)*u_02(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+8.0*&
r*r*r*r*ub_00(q,p)*u_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*k_03(q,p)*e2bi+j_00(q,p)*j_00(q,p)*j_00(q,p)*r*r*r&
*r*ub_00(q,p)*ub_00(q,p)*ub_00(q,p)*k_00(q,p)*jb_01(q,p)*jb_02(q,p)*e2bi+4.0*jb_00(q,p)*r*r*r*r*u_00(q,p)*u_00(q,p)*ub_00(q,p)*&
k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*k_02(q,p)*e2bi+3.0*j_00(q,p)*jb_02(q,p)*r*r*r*r*u_00(q,p)*u_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)&
*k_00(q,p)*jb_01(q,p)*e2bi
      s9 = s7+8.0*r*r*r*ub_00(q,p)*ub_00(q,p)*u_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+4.0*j_00(q,p)*j_00(q,p)&
*r*r*r*r*ub_00(q,p)*ub_03(q,p)*u_00(q,p)*jb_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-2.0*j_01(q,p)*r*r*r*r*ub_00(q,p)*&
ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi+2.0*j_00(q,p)*j_00(q,p)*r*r*r*r*ub_00(q,p)*ub_00(q,p)*u_03(q,p)*k_00(q,p)&
*k_00(q,p)*k_00(q,p)*jb_01(q,p)*e2bi-2.0*r*r*r*r*ub_00(q,p)*jb_02(q,p)*u_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi-2.0&
*r*r*r*r*ub_00(q,p)*jb_02(q,p)*u_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi+8.0*jb_00(q,p)*r*r*r*r*&
u_00(q,p)*u_02(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi+8.0*r*r*r*ub_00(q,p)*ub_00(q,p)*j_04(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)*e2bi-4.0*r*r*r*r*ub_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*j_04(q,p)*k_01(q,p)*e2bi+8.0*r*r*r*r*&
j_04(q,p)*ub_00(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-4.0*r*r*r*r*ub_00(q,p)*ub_00(q,p)*j_11(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)*e2bi+8.0*jb_00(q,p)*r*r*r*ub_00(q,p)*u_00(q,p)*u_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-j_00(q,p)*jb_03(q,p)*r*r*r&
*r*u_00(q,p)*ub_00(q,p)*ub_00(q,p)*j_01(q,p)*k_00(q,p)*e2bi+3.0*j_00(q,p)*r*r*r*r*jb_03(q,p)*u_00(q,p)*u_00(q,p)*u_00(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)*jb_01(q,p)*e2bi+8.0*r*r*r*r*u_00(q,p)*ub_00(q,p)*u_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_03(q,p)*e2bi-2.0*&
r*r*r*r*j_02(q,p)*ub_00(q,p)*ub_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi+8.0*jb_00(q,p)*r*r*r*r*u_00(q,p)*u_00(q,p)*&
k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*k_04(q,p)*e2bi+8.0*j_00(q,p)*r*r*r*r*ub_00(q,p)*ub_00(q,p)*ub_12(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)*e2bi
      s8 = s9+12.0*j_00(q,p)*r*r*r*ub_00(q,p)*ub_00(q,p)*ub_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-2.0*j_00(q,p)*&
jb_00(q,p)*jb_01(q,p)*r*r*r*r*u_00(q,p)*u_00(q,p)*ub_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+2.0*j_01(q,p)*r*r*r*r*ub_00(q,p)*&
ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi-8.0*j_01(q,p)*r*r*r*r*r*ub_00(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*&
e2bi-4.0*j_01(q,p)*r*r*r*r*r*ub_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_01(q,p)*e2bi-2.0*j_01(q,p)*r*r*r*r*r&
*ub_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*w_00(q,p)*e2bi+8.0*j_00(q,p)*r*r*r*r*ub_00(q,p)*ub_01(q,p)*&
ub_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-8.0*j_00(q,p)*r*r*ub_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+8.0*j_00(q,p)&
*jb_03(q,p)*r*r*u_00(q,p)*ub_00(q,p)*b_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)-j_00(q,p)*j_00(q,p)*j_00(q,p)*ub_00(q,p)*ub_00(q,p)*r*r*r*r*r*&
k_00(q,p)*jb_01(q,p)*jb_01(q,p)*w_00(q,p)*e2bi+2.0*j_00(q,p)*j_00(q,p)*jb_00(q,p)*jb_04(q,p)*r*r*r*r*u_00(q,p)*u_00(q,p)*jb_01(q,p)*&
k_00(q,p)*e2bi-4.0*r*r*r*r*r*ub_00(q,p)*ub_00(q,p)*j_11(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*e2bi-12.0*j_01(q,p)*r*&
r*r*r*ub_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*w_00(q,p)*e2bi-j_00(q,p)*jb_02(q,p)*r*r*r*r*ub_00(q,p)*ub_00(q,p)*ub_00(q,p)&
*j_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-8.0*jb_00(q,p)*r*r*u_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-2.0*j_00(q,p)*&
jb_02(q,p)*r*r*u_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_03(q,p)+2.0*j_01(q,p)*r*r*r*r*r*ub_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*&
k_01(q,p)*w_00(q,p)*e2bi-4.0*j_00(q,p)*r*r*ub_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*jb_23(q,p)+4.0*j_00(q,p)*r*r*u_00(q,p)*&
ub_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*jb_02(q,p)
      s9 = s8-j_00(q,p)*jb_03(q,p)*r*r*r*r*u_00(q,p)*ub_00(q,p)*ub_00(q,p)*j_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+j_00(q,p)&
*jb_02(q,p)*r*r*r*r*ub_00(q,p)*u_00(q,p)*u_00(q,p)*jb_01(q,p)*k_00(q,p)*e2bi+2.0*r*r*r*r*u_00(q,p)*u_00(q,p)*jb_01(q,p)*u_03(q,p)&
*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-r*r*ub_00(q,p)*u_00(q,p)*k_00(q,p)*j_02(q,p)*jb_03(q,p)+r*r*r*r*ub_00(q,p)*&
u_00(q,p)*u_00(q,p)*jb_01(q,p)*k_00(q,p)*j_00(q,p)*j_00(q,p)*jb_00(q,p)*jb_02(q,p)*e2bi+j_00(q,p)*j_00(q,p)*j_00(q,p)*r*r*r*r*u_00(q,p)*&
ub_00(q,p)*ub_00(q,p)*k_00(q,p)*jb_01(q,p)*jb_03(q,p)*e2bi+4.0*r*r*r*r*u_00(q,p)*ub_01(q,p)*ub_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)*e2bi-2.0*r*r*r*r*u_00(q,p)*u_00(q,p)*jb_01(q,p)*ub_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-2.0*&
j_00(q,p)*j_00(q,p)*r*r*r*r*jb_01(q,p)*ub_00(q,p)*ub_00(q,p)*ub_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-2.0*j_00(q,p)*jb_00(q,p)*&
jb_02(q,p)*r*r*r*r*u_00(q,p)*u_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi+4.0*r*r*r*ub_00(q,p)*u_00(q,p)*u_00(q,p)*&
jb_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+4.0*r*r*r*r*ub_00(q,p)*ub_00(q,p)*j_01(q,p)*ub_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*&
e2bi+4.0*r*r*r*ub_00(q,p)*ub_00(q,p)*ub_00(q,p)*j_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+8.0*j_00(q,p)*r*r*r*r*&
ub_00(q,p)*ub_14(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+6.0*r*r*r*r*ub_00(q,p)*ub_00(q,p)*ub_01(q,p)*j_02(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)*e2bi+4.0*r*r*r*r*ub_00(q,p)*ub_00(q,p)*ub_00(q,p)*j_12(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-8.0*jb_01(q,p)*r*r&
*r*u_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+2.0*r*r*r*r*j_03(q,p)*ub_00(q,p)*ub_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi
      s4 = s9+2.0*r*r*r*r*u_01(q,p)*j_03(q,p)*ub_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+8.0*r*r*r*&
u_00(q,p)*u_00(q,p)*jb_04(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+4.0*r*r*r*r*u_00(q,p)*u_00(q,p)*jb_01(q,p)*u_03(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)*e2bi+8.0*r*r*r*r*j_14(q,p)*ub_00(q,p)*ub_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+8.0*j_00(q,p)*r*r*r*r*&
u_00(q,p)*ub_00(q,p)*ub_13(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-2.0*r*r*r*r*j_03(q,p)*ub_00(q,p)*ub_00(q,p)*u_00(q,p)*k_00(q,p)*&
k_00(q,p)*k_01(q,p)*e2bi-4.0*jb_00(q,p)*u_00(q,p)*u_00(q,p)*r*r*r*r*r*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*k_01(q,p)*w_00(q,p)*&
e2bi+4.0*j_00(q,p)*r*r*r*ub_00(q,p)*ub_00(q,p)*u_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi-2.0*j_00(q,p)*r*r*r*r*jb_04(q,p)&
*ub_00(q,p)*ub_00(q,p)*j_01(q,p)*k_00(q,p)*e2bi-2.0*j_00(q,p)*r*r*r*r*jb_04(q,p)*ub_00(q,p)*ub_00(q,p)*j_01(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)*e2bi-16.0*j_00(q,p)*r*r*r*r*ub_00(q,p)*b_04(q,p)*ub_01(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+4.0*r*r*r*r*&
u_00(q,p)*ub_00(q,p)*u_01(q,p)*jb_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*e2bi+4.0*jb_00(q,p)*r*r*r*r*u_00(q,p)*u_00(q,p)*ub_02(q,p)*&
k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi+4.0*r*r*r*r*u_00(q,p)*u_00(q,p)*u_00(q,p)*jb_13(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*&
e2bi-4.0*jb_00(q,p)*r*r*r*r*u_00(q,p)*u_00(q,p)*u_03(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*k_01(q,p)*e2bi-3.0*r*r*&
ub_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*j_02(q,p)*jb_03(q,p)-r*r*u_00(q,p)*j_03(q,p)*ub_00(q,p)*k_00(q,p)*jb_02(q,p)-4.0*jb_00(q,p)*r&
*r*ub_00(q,p)*u_00(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*j_23(q,p)-8.0*j_00(q,p)*r*r*ub_00(q,p)*jb_02(q,p)*u_00(q,p)*b_03(q,p)*k_00(q,p)*k_00(q,p)*&
k_00(q,p)-8.0*j_00(q,p)*r*r*ub_02(q,p)*k_00(q,p)*k_00(q,p)*k_00(q,p)*jb_04(q,p)
      s2 = s3*s4
      R00(q,p) = real(s1*s2)

    end do
  end do

#ifdef DEBUG_LARGE_VALUE
  if (maxval(abs(R00)).gt.DEBUG_LARGE_VALUE) then
     write (*,*) 'R00 = ', maxval(abs(R00))
     write (*,*) 'big: input:', maxval(abs(w_11))

     write (*,*) 'input1:', maxval(abs(w_00)),maxval(abs(w_01)),maxval(abs(w_02)),maxval(abs(w_03)),maxval(abs(w_04)),&
               maxval(abs(w_11)),maxval(abs(w_12)),maxval(abs(w_13)),maxval(abs(w_14))
     write (*,*) 'input2:',maxval(abs(w_22)),maxval(abs(w_23)),maxval(abs(w_24)),maxval(abs(w_33)),maxval(abs(w_34))
     write (*,*) 'input3:',maxval(abs(b_00)),maxval(abs(b_01)),maxval(abs(b_02)),maxval(abs(b_03)),maxval(abs(b_04)),&
               maxval(abs(b_11)),maxval(abs(b_12)),maxval(abs(b_13)),maxval(abs(b_14))
     write (*,*) 'input4:',maxval(abs(b_22)),maxval(abs(b_23)),maxval(abs(b_24)),maxval(abs(b_33)),maxval(abs(b_34))
     write (*,*) 'input5:',maxval(abs(k_00)),maxval(abs(k_01)),maxval(abs(k_02)),maxval(abs(k_03)),maxval(abs(k_04)),&
               maxval(abs(k_11)),maxval(abs(k_12)),maxval(abs(k_13)),maxval(abs(k_14))
     write (*,*) 'input6:',maxval(abs(k_22)),maxval(abs(k_23)),maxval(abs(k_24)),maxval(abs(k_33)),maxval(abs(k_34))
     write (*,*) 'input7:',maxval(abs(j_00)),maxval(abs(j_01)),maxval(abs(j_02)),maxval(abs(j_03)),maxval(abs(j_04)),&
               maxval(abs(j_11)),maxval(abs(j_12)),maxval(abs(j_13)),maxval(abs(j_14))
     write (*,*) 'input8:',maxval(abs(j_22)),maxval(abs(j_23)),maxval(abs(j_24)),maxval(abs(j_33)),maxval(abs(j_34))
     write (*,*) 'input9:',maxval(abs(jb_00)),maxval(abs(jb_01)),maxval(abs(jb_02)),maxval(abs(jb_03)),maxval(abs(jb_04)),&
               maxval(abs(jb_11)),maxval(abs(jb_12)),maxval(abs(jb_13)),maxval(abs(jb_14))
     write (*,*) 'input10:',maxval(abs(jb_22)),maxval(abs(jb_23)),maxval(abs(jb_24)),maxval(abs(jb_33)),maxval(abs(jb_34))
     write (*,*) 'input11:',maxval(abs(u_00)),maxval(abs(u_01)),maxval(abs(u_02)),maxval(abs(u_03)),maxval(abs(u_04)),&
               maxval(abs(u_11)),maxval(abs(u_12)),maxval(abs(u_13)),maxval(abs(u_14))
     write (*,*) 'input12:',maxval(abs(u_22)),maxval(abs(u_23)),maxval(abs(u_24)),maxval(abs(u_33)),maxval(abs(u_34))
     write (*,*) 'input13:',maxval(abs(ub_00)),maxval(abs(ub_01)),maxval(abs(ub_02)),maxval(abs(ub_03)),maxval(abs(ub_04)),&
               maxval(abs(ub_11)),maxval(abs(ub_12)),maxval(abs(ub_13)),maxval(abs(ub_14))
     write (*,*) 'input14:',maxval(abs(ub_22)),maxval(abs(ub_23)),maxval(abs(ub_24)),maxval(abs(ub_33)),maxval(abs(ub_34))
     stop
  end if
#endif
  
end subroutine NullConstr_R00_calc

end module NullConstr_R00
