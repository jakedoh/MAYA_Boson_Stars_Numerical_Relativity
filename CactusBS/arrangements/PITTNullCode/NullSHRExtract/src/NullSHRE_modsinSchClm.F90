! vim: syntax=fortran
#include "cctk.h"

module NullSHRE_modsinSchClm

  ! gives the numerical spherical harmonic coefficients Clm for the IEF sinSchwarzschild metric 

  use cctk
  implicit none

contains

  subroutine sinSchAlm_00(R,M,T,F,A, A00)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: A00 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: sine 

      sine = A*sin(2.d0*pi*F*T)

      A00 = 2.d0/dsqrt(R+2.d0*M*sine)*dsqrt(R)*dsqrt(pi)

  end subroutine sinSchAlm_00  


  subroutine sinSchClm_00(R,M,T,F,A,&
             C00_11, C00_12, C00_13, C00_22, C00_23, C00_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: C00_11, C00_12, C00_13, C00_22, C00_23, C00_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: sine 

      sine = A*sin(2.d0*pi*F*T)

      C00_11 = 2.d0/3.d0*dsqrt(pi)*(2.d0*M*sine+3.d0*R)/R
      C00_12 = 0.d0
      C00_13 = 0.d0
      C00_22 = 2.d0/3.d0*dsqrt(pi)*(2.d0*M*sine+3.d0*R)/R
      C00_23 = 0.d0
      C00_33 = 2.d0/3.d0*dsqrt(pi)*(2.d0*M*sine+3.d0*R)/R

  end subroutine sinSchClm_00


  subroutine sinSchBlm_1m1(R,M,T,F,A, B1m1_1, B1m1_2, B1m1_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: B1m1_1, B1m1_2, B1m1_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: sine 

      sine = A*sin(2.d0*pi*F*T)

      B1m1_1 = -2.d0/3.d0*dsqrt(pi)*M*sine/(R+2.d0*M*sine)*dsqrt(6.d0)
      B1m1_2 = dcmplx(0.d0,-2.d0/3.d0)*dsqrt(pi)*M*sine/(R+2.d0*M*sine)*dsqrt(6.d0)
      B1m1_3 = 0.d0

  end subroutine sinSchBlm_1m1


  subroutine sinSchBlm_10(R,M,T,F,A, B10_1, B10_2, B10_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: B10_1, B10_2, B10_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: sine 

      sine = A*sin(2.d0*pi*F*T)

      B10_1 = 0.d0
      B10_2 = 0.d0
      B10_3 = -4.d0/3.d0*M*sine/(R+2.d0*M*sine)*dsqrt(3.d0)*dsqrt(pi)

  end subroutine sinSchBlm_10


  subroutine sinSchBlm_11(R,M,T,F,A, B11_1, B11_2, B11_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: B11_1, B11_2, B11_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: sine 

      sine = A*sin(2.d0*pi*F*T)

      B11_1 = 2.d0/3.d0*dsqrt(pi)*M*sine/(R+2.d0*M*sine)*dsqrt(6.d0)
      B11_2 = dcmplx(0.d0,-2.d0/3.d0)*dsqrt(pi)*M*sine/(R+2.d0*M*sine)*dsqrt(6.d0)
      B11_3 = 0.d0

  end subroutine sinSchBlm_11


  subroutine sinSchClm_2m2(R,M,T,F,A,&
             C2m2_11, C2m2_12, C2m2_13, C2m2_22, C2m2_23, C2m2_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: C2m2_11, C2m2_12, C2m2_13, C2m2_22, C2m2_23, C2m2_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: sine 

      sine = A*sin(2.d0*pi*F*T)

      C2m2_11 = 2.d0/15.d0*dsqrt(30.d0)/R*dsqrt(pi)*M*sine
      C2m2_12 = dcmplx(0.d0,2.d0/15.d0)*dsqrt(pi)*M*sine/R*dsqrt(30.d0)
      C2m2_13 = 0.d0
      C2m2_22 = -2.d0/15.d0*dsqrt(30.d0)/R*dsqrt(pi)*M*sine
      C2m2_23 = 0.d0
      C2m2_33 = 0.d0

  end subroutine sinSchClm_2m2


  subroutine sinSchClm_2m1(R,M,T,F,A,&
             C2m1_11, C2m1_12, C2m1_13, C2m1_22, C2m1_23, C2m1_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: C2m1_11, C2m1_12, C2m1_13, C2m1_22, C2m1_23, C2m1_33 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: sine 

      sine = A*(2.d0*pi*F*T)

      C2m1_11 = 0.d0
      C2m1_12 = 0.d0
      C2m1_13 = 2.d0/15.d0*dsqrt(30.d0)/R*dsqrt(pi)*M*sine
      C2m1_22 = 0.d0
      C2m1_23 = dcmplx(0.d0,2.d0/15.d0)*dsqrt(pi)*M*sine/R*dsqrt(30.d0)
      C2m1_33 = 0.d0

  end subroutine sinSchClm_2m1


  subroutine sinSchClm_20(R,M,T,F,A,&
             C20_11, C20_12, C20_13, C20_22, C20_23, C20_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: C20_11, C20_12, C20_13, C20_22, C20_23, C20_33 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: sine 

      sine = A*sin(2.d0*pi*F*T)

      C20_11 = -4.d0/15.d0*dsqrt(pi)*M*sine*dsqrt(5.d0)/R
      C20_12 = 0.d0
      C20_13 = 0.d0
      C20_22 = -4.d0/15.d0*dsqrt(pi)*M*sine*dsqrt(5.d0)/R
      C20_23 = 0.d0
      C20_33 = 8.d0/15.d0*dsqrt(pi)*M*sine*dsqrt(5.d0)/R

  end subroutine sinSchClm_20


  subroutine sinSchClm_21(R,M,T,F,A,&
             C21_11, C21_12, C21_13, C21_22, C21_23, C21_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: C21_11, C21_12, C21_13, C21_22, C21_23, C21_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: sine 

      sine = A*sin(2.d0*pi*F*T)

      C21_11 = 0.d0
      C21_12 = 0.d0
      C21_13 = -2.d0/15.d0*dsqrt(30.d0)/R*dsqrt(pi)*M*sine
      C21_22 = 0.d0
      C21_23 = dcmplx(0.d0,2.d0/15.d0)*dsqrt(pi)*M*sine/R*dsqrt(30.d0)
      C21_33 = 0.d0

  end subroutine sinSchClm_21


  subroutine sinSchClm_22(R,M,T,F,A,&
             C22_11, C22_12, C22_13, C22_22, C22_23, C22_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: C22_11, C22_12, C22_13, C22_22, C22_23, C22_33 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: sine 

      sine = A*sin(2.d0*pi*F*T)

      C22_11 = 2.d0/15.d0*dsqrt(30.d0)/R*dsqrt(pi)*M*sine
      C22_12 = dcmplx(0.d0,-2.d0/15.d0)*dsqrt(pi)*M*sine/R*dsqrt(30.d0)
      C22_13 = 0.d0
      C22_22 = -2.d0/15.d0*dsqrt(30.d0)/R*dsqrt(pi)*M*sine
      C22_23 = 0.d0
      C22_33 = 0.d0

  end subroutine sinSchClm_22

!time derivatives of spherical harmonic coefficients

  subroutine sinSchdtAlm_00(R,M,T,F,A, dtA00)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: dtA00 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: sine, ocos 

      sine = A*sin(2.d0*pi*F*T)
      ocos = 2.d0*A*pi*F*cos(2.d0*pi*F*T)

      dtA00 = -2.d0/dsqrt(R+2.d0*M*sine)**3*dsqrt(R)*dsqrt(pi)*M*ocos

  end subroutine sinSchdtAlm_00  


  subroutine sinSchdtClm_00(R,M,T,F,A,&
             dtC00_11, dtC00_12, dtC00_13, dtC00_22, dtC00_23, dtC00_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: dtC00_11, dtC00_12, dtC00_13, dtC00_22, dtC00_23, dtC00_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: sine, ocos 

      sine = A*sin(2.d0*pi*F*T)
      ocos = 2.d0*A*pi*F*cos(2.d0*pi*F*T)

      dtC00_11 = 4.d0/3.d0*dsqrt(pi)*M*ocos/R
      dtC00_12 = 0.d0
      dtC00_13 = 0.d0
      dtC00_22 = 4.d0/3.d0*dsqrt(pi)*M*ocos/R
      dtC00_23 = 0.d0
      dtC00_33 = 4.d0/3.d0*dsqrt(pi)*M*ocos/R

  end subroutine sinSchdtClm_00


  subroutine sinSchdtBlm_1m1(R,M,T,F,A, dtB1m1_1, dtB1m1_2, dtB1m1_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: dtB1m1_1, dtB1m1_2, dtB1m1_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: sine, ocos 

      sine = A*sin(2.d0*pi*F*T)
      ocos = 2.d0*A*pi*F*cos(2.d0*pi*F*T)

      dtB1m1_1 = -2.d0/3.d0*dsqrt(pi)*M*ocos*R/(R+2.d0*M*sine)**2*dsqrt(6.d0)
      dtB1m1_2 = dcmplx(0.d0,-2.d0/3.d0)*dsqrt(pi)*M*ocos*R/(R+2.d0*M*sine)**2*dsqrt(6.d0)
      dtB1m1_3 = 0.d0

  end subroutine sinSchdtBlm_1m1


  subroutine sinSchdtBlm_10(R,M,T,F,A, dtB10_1, dtB10_2, dtB10_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: dtB10_1, dtB10_2, dtB10_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: sine, ocos 

      sine = A*sin(2.d0*pi*F*T)
      ocos = 2.d0*A*pi*F*cos(2.d0*pi*F*T)

      dtB10_1 = 0.d0
      dtB10_2 = 0.d0
      dtB10_3 = -4.d0/3.d0*M*ocos*R/(R+2.d0*M*sine)**2*dsqrt(3.d0)*dsqrt(pi)

  end subroutine sinSchdtBlm_10


  subroutine sinSchdtBlm_11(R,M,T,F,A, dtB11_1, dtB11_2, dtB11_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: dtB11_1, dtB11_2, dtB11_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: sine, ocos 

      sine = A*sin(2.d0*pi*F*T)
      ocos = 2.d0*A*pi*F*cos(2.d0*pi*F*T)

      dtB11_1 = 2.d0/3.d0*dsqrt(pi)*M*ocos*R/(R+2.d0*M*sine)**2*dsqrt(6.d0)
      dtB11_2 = dcmplx(0.d0,-2.d0/3.d0)*dsqrt(pi)*M*ocos*R/(R+2.d0*M*sine)**2*dsqrt(6.d0)
      dtB11_3 = 0.d0

  end subroutine sinSchdtBlm_11


  subroutine sinSchdtClm_2m2(R,M,T,F,A,&
             dtC2m2_11, dtC2m2_12, dtC2m2_13, dtC2m2_22, dtC2m2_23, dtC2m2_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: dtC2m2_11, dtC2m2_12, dtC2m2_13, dtC2m2_22, dtC2m2_23, dtC2m2_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: sine, ocos 

      sine = A*sin(2.d0*pi*F*T)
      ocos = 2.d0*A*pi*F*cos(2.d0*pi*F*T)

      dtC2m2_11 = 2.d0/15.d0*dsqrt(30.d0)/R*dsqrt(pi)*M*ocos
      dtC2m2_12 = dcmplx(0.d0,2.d0/15.d0)*dsqrt(pi)*M*ocos/R*dsqrt(30.d0)
      dtC2m2_13 = 0.d0
      dtC2m2_22 = -2.d0/15.d0*dsqrt(30.d0)/R*dsqrt(pi)*M*ocos
      dtC2m2_23 = 0.d0
      dtC2m2_33 = 0.d0

  end subroutine sinSchdtClm_2m2


  subroutine sinSchdtClm_2m1(R,M,T,F,A,&
             dtC2m1_11, dtC2m1_12, dtC2m1_13, dtC2m1_22, dtC2m1_23, dtC2m1_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: dtC2m1_11, dtC2m1_12, dtC2m1_13, dtC2m1_22, dtC2m1_23, dtC2m1_33 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: sine, ocos 

      sine = A*sin(2.d0*pi*F*T)
      ocos = 2.d0*A*pi*F*cos(2.d0*pi*F*T)

      dtC2m1_11 = 0.d0
      dtC2m1_12 = 0.d0
      dtC2m1_13 = 2.d0/15.d0*dsqrt(30.d0)/R*dsqrt(pi)*M*ocos
      dtC2m1_22 = 0.d0
      dtC2m1_23 = dcmplx(0.d0,2.d0/15.d0)*dsqrt(pi)*M*ocos/R*dsqrt(30.d0)
      dtC2m1_33 = 0.d0

  end subroutine sinSchdtClm_2m1


  subroutine sinSchdtClm_20(R,M,T,F,A,&
             dtC20_11, dtC20_12, dtC20_13, dtC20_22, dtC20_23, dtC20_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: dtC20_11, dtC20_12, dtC20_13, dtC20_22, dtC20_23, dtC20_33 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: sine, ocos 

      sine = A*sin(2.d0*pi*F*T)
      ocos = 2.d0*A*pi*F*cos(2.d0*pi*F*T)

      dtC20_11 = -4.d0/15.d0*dsqrt(pi)*M*ocos*dsqrt(5.d0)/R
      dtC20_12 = 0.d0
      dtC20_13 = 0.d0
      dtC20_22 = -4.d0/15.d0*dsqrt(pi)*M*ocos*dsqrt(5.d0)/R
      dtC20_23 = 0.d0
      dtC20_33 = 8.d0/15.d0*dsqrt(pi)*M*ocos*dsqrt(5.d0)/R

  end subroutine sinSchdtClm_20


  subroutine sinSchdtClm_21(R,M,T,F,A,&
             dtC21_11, dtC21_12, dtC21_13, dtC21_22, dtC21_23, dtC21_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: dtC21_11, dtC21_12, dtC21_13, dtC21_22, dtC21_23, dtC21_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: sine, ocos 

      sine = A*sin(2.d0*pi*F*T)
      ocos = 2.d0*A*pi*F*cos(2.d0*pi*F*T)

      dtC21_11 = 0.d0
      dtC21_12 = 0.d0
      dtC21_13 = -2.d0/15.d0*dsqrt(30.d0)/R*dsqrt(pi)*M*ocos
      dtC21_22 = 0.d0
      dtC21_23 = dcmplx(0.d0,2.d0/15.d0)*dsqrt(pi)*M*ocos/R*dsqrt(30.d0)
      dtC21_33 = 0.d0

  end subroutine sinSchdtClm_21


  subroutine sinSchdtClm_22(R,M,T,F,A,&
             dtC22_11, dtC22_12, dtC22_13, dtC22_22, dtC22_23, dtC22_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: dtC22_11, dtC22_12, dtC22_13, dtC22_22, dtC22_23, dtC22_33 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: sine, ocos 

      sine = A*sin(2.d0*pi*F*T)
      ocos = 2.d0*A*pi*F*cos(2.d0*pi*F*T)

      dtC22_11 = 2.d0/15.d0*dsqrt(30.d0)/R*dsqrt(pi)*M*ocos
      dtC22_12 = dcmplx(0.d0,-2.d0/15.d0)*dsqrt(pi)*M*ocos/R*dsqrt(30.d0)
      dtC22_13 = 0.d0
      dtC22_22 = -2.d0/15.d0*dsqrt(30.d0)/R*dsqrt(pi)*M*ocos
      dtC22_23 = 0.d0
      dtC22_33 = 0.d0

  end subroutine sinSchdtClm_22


!radial derivatives of spherical harmonic coefficients
  subroutine sinSchdrAlm_00(R,M,T,F,A, drA00)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: drA00 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: sine 

      sine = A*sin(2.d0*pi*F*T)

      drA00 = 2*dsqrt(pi)*M*sine/dsqrt(R+2*M*sine)**3/dsqrt(R) 

  end subroutine sinSchdrAlm_00  


  subroutine sinSchdrClm_00(R,M,T,F,A,&
             drC00_11, drC00_12, drC00_13, drC00_22, drC00_23, drC00_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: drC00_11, drC00_12, drC00_13, drC00_22, drC00_23, drC00_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: sine 

      sine = A*sin(2.d0*pi*F*T)


      drC00_11 = -4.d0/3.d0*dsqrt(pi)*M*sine/R**2
      drC00_12 = 0.d0
      drC00_13 = 0.d0
      drC00_22 = -4.d0/3.d0*dsqrt(pi)*M*sine/R**2
      drC00_23 = 0.d0
      drC00_33 = -4.d0/3.d0*dsqrt(pi)*M*sine/R**2

  end subroutine sinSchdrClm_00


  subroutine sinSchdrBlm_1m1(R,M,T,F,A, drB1m1_1, drB1m1_2, drB1m1_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: drB1m1_1, drB1m1_2, drB1m1_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: sine 

      sine = A*sin(2.d0*pi*F*T)

      drB1m1_1 = 2.d0/3.d0*M*sine/(R+2*M*sine)**2*dsqrt(6.d0)*dsqrt(pi)
      drB1m1_2 = dcmplx(0.d0,2.d0/3.d0)*M*sine/(R+2*M*sine)**2&
                *dsqrt(6.d0)*dsqrt(pi)
      drB1m1_3 = 0.d0

  end subroutine sinSchdrBlm_1m1


  subroutine sinSchdrBlm_10(R,M,T,F,A, drB10_1, drB10_2, drB10_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: drB10_1, drB10_2, drB10_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: sine 

      sine = A*sin(2.d0*pi*F*T)


      drB10_1 = 0.d0
      drB10_2 = 0.d0
      drB10_3 = 4.d0/3.d0*M*sine/(R+2*M*sine)**2*dsqrt(3.d0)*dsqrt(pi)

  end subroutine sinSchdrBlm_10


  subroutine sinSchdrBlm_11(R,M,T,F,A, drB11_1, drB11_2, drB11_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: drB11_1, drB11_2, drB11_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: sine 

      sine = A*sin(2.d0*pi*F*T)

      drB11_1 = -2.d0/3.d0*M*sine/(R+2*M*sine)**2*dsqrt(6.d0)*dsqrt(pi)
      drB11_2 = dcmplx(0.d0,2.d0/3.d0)*M*sine/(R+2*M*sine)**2*dsqrt(6.d0)*dsqrt(pi)
      drB11_3 = 0.d0

  end subroutine sinSchdrBlm_11


  subroutine sinSchdrClm_2m2(R,M,T,F,A,&
             drC2m2_11, drC2m2_12, drC2m2_13, drC2m2_22, drC2m2_23, drC2m2_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: drC2m2_11, drC2m2_12, drC2m2_13,&
                      drC2m2_22, drC2m2_23, drC2m2_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: sine 

      sine = A*sin(2.d0*pi*F*T)

      drC2m2_11 = -2.d0/15.d0*dsqrt(30.d0)/R**2*dsqrt(pi)*M*sine
      drC2m2_12 = dcmplx(0.d0,-2.d0/15.d0)*dsqrt(30.d0)/R**2*dsqrt(pi)*M*sine
      drC2m2_13 = 0.d0
      drC2m2_22 = 2.d0/15.d0*dsqrt(30.d0)/R**2*dsqrt(pi)*M*sine
      drC2m2_23 = 0.d0
      drC2m2_33 = 0.d0

  end subroutine sinSchdrClm_2m2


  subroutine sinSchdrClm_2m1(R,M,T,F,A,&
             drC2m1_11, drC2m1_12, drC2m1_13, drC2m1_22, drC2m1_23, drC2m1_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: drC2m1_11, drC2m1_12, drC2m1_13, &
                      drC2m1_22, drC2m1_23, drC2m1_33 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: sine 

      sine = A*sin(2.d0*pi*F*T)

      drC2m1_11 = 0.d0
      drC2m1_12 = 0.d0
      drC2m1_13 = -2.d0/15.d0*dsqrt(30.d0)/R**2*dsqrt(pi)*M*sine
      drC2m1_22 = 0.d0
      drC2m1_23 = dcmplx(0.d0,-2.d0/15.d0)*dsqrt(30.d0)/R**2*dsqrt(pi)*M*sine
      drC2m1_33 = 0.d0

  end subroutine sinSchdrClm_2m1


  subroutine sinSchdrClm_20(R,M,T,F,A,&
             drC20_11, drC20_12, drC20_13, drC20_22, drC20_23, drC20_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: drC20_11, drC20_12, drC20_13,&
                      drC20_22, drC20_23, drC20_33 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: sine 

      sine = A*sin(2.d0*pi*F*T)

      drC20_11 = 4.d0/15.d0*dsqrt(pi)*M*sine*dsqrt(5.d0)/R**2
      drC20_12 = 0.d0
      drC20_13 = 0.d0
      drC20_22 = 4.d0/15.d0*dsqrt(pi)*M*sine*dsqrt(5.d0)/R**2
      drC20_23 = 0.d0
      drC20_33 = -8.d0/15.d0*dsqrt(pi)*M*sine*dsqrt(5.d0)/R**2

  end subroutine sinSchdrClm_20


  subroutine sinSchdrClm_21(R,M,T,F,A,&
             drC21_11, drC21_12, drC21_13, drC21_22, drC21_23, drC21_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: drC21_11, drC21_12, drC21_13,&
                      drC21_22, drC21_23, drC21_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: sine 

      sine = A*sin(2.d0*pi*F*T)

      drC21_11 = 0.d0
      drC21_12 = 0.d0
      drC21_13 = 2.d0/15.d0*dsqrt(30.d0)/R**2*dsqrt(pi)*M*sine
      drC21_22 = 0.d0
      drC21_23 = dcmplx(0.d0,-2.d0/15.d0)*dsqrt(30.d0)/R**2*dsqrt(pi)*M*sine
      drC21_33 = 0.d0

  end subroutine sinSchdrClm_21


  subroutine sinSchdrClm_22(R,M,T,F,A,&
             drC22_11, drC22_12, drC22_13, drC22_22, drC22_23, drC22_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: drC22_11, drC22_12, drC22_13,&
                      drC22_22, drC22_23, drC22_33 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: sine 

      sine = A*sin(2.d0*pi*F*T)

      drC22_11 = -2.d0/15.d0*dsqrt(30.d0)/R**2*dsqrt(pi)*M*sine
      drC22_12 = dcmplx(0.d0,2.d0/15.d0)*dsqrt(30.d0)/R**2*dsqrt(pi)*M*sine
      drC22_13 = 0.d0
      drC22_22 = 2.d0/15.d0*dsqrt(30.d0)/R**2*dsqrt(pi)*M*sine
      drC22_23 = 0.d0
      drC22_33 = 0.d0

  end subroutine sinSchdrClm_22


end module NullSHRE_modsinSchClm

