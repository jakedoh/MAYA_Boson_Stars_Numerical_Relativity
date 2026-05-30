! vim: syntax=fortran
#include "cctk.h"

module NullSHRE_modvibSchClm

  ! gives the numerical spherical harmonic coefficients Clm for the IEF Schwarzschild metric 

  use cctk
  implicit none

contains

  subroutine vibSchAlm_00(R,M,T,F,A, A00)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: A00 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: V, Vt 

      V = A*sin(2.d0*pi*F*T)
      Vt = 2*A*cos(2*pi*F*T)*pi*T

      A00 = 2/sqrt(R+2*M)*sqrt(R)*sqrt(pi)

  end subroutine vibSchAlm_00  


  subroutine vibSchClm_00(R,M,T,F,A,&
             C00_11, C00_12, C00_13, C00_22, C00_23, C00_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: C00_11, C00_12, C00_13, C00_22, C00_23, C00_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: V, Vt 

      V = A*sin(2.d0*pi*F*T)
      Vt = 2*A*cos(2*pi*F*t)*pi*F

      C00_11 = 2.D0/3.D0*sqrt(pi)*(2*M+3*R)/R
      C00_12 = 0
      C00_13 = 0
      C00_22 = 2.D0/3.D0*sqrt(pi)*(2*M+3*R)/R
      C00_23 = 0
      C00_33 = 2.D0/3.D0*sqrt(pi)*(2*M+3*R)/R

  end subroutine vibSchClm_00


  subroutine vibSchBlm_1m1(R,M,T,F,A, B1m1_1, B1m1_2, B1m1_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: B1m1_1, B1m1_2, B1m1_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: V, Vt 

      V = A*sin(2.d0*pi*F*T)
      Vt = 2*A*cos(2*pi*F*t)*pi*F

      B1m1_1 = -2.D0/3.D0*sqrt(pi)*M/(R+2*M)*sqrt(6.D0)
      B1m1_2 = dcmplx(0.D0,-2.D0/3.D0)*sqrt(pi)*M/(R+2*M)*sqrt(6.D0)
      B1m1_3 = 0

  end subroutine vibSchBlm_1m1


  subroutine vibSchBlm_10(R,M,T,F,A, B10_1, B10_2, B10_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: B10_1, B10_2, B10_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: V, Vt 

      V = A*sin(2.d0*pi*F*T)
      Vt = 2*A*cos(2*pi*F*t)*pi*F

      B10_1 = 0
      B10_2 = 0
      B10_3 = -4.D0/3.D0*M/(R+2*M)*sqrt(3.D0)*sqrt(pi)

  end subroutine vibSchBlm_10


  subroutine vibSchBlm_11(R,M,T,F,A, B11_1, B11_2, B11_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: B11_1, B11_2, B11_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: V, Vt 

      V = A*sin(2.d0*pi*F*T)
      Vt = 2*A*cos(2*pi*F*t)*pi*F

      B11_1 = 2.D0/3.D0*sqrt(pi)*M/(R+2*M)*sqrt(6.D0)
      B11_2 = dcmplx(0.D0,-2.D0/3.D0)*sqrt(pi)*M/(R+2*M)*sqrt(6.D0)
      B11_3 = 0

  end subroutine vibSchBlm_11


  subroutine vibSchClm_2m2(R,M,T,F,A,&
             C2m2_11, C2m2_12, C2m2_13, C2m2_22, C2m2_23, C2m2_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: C2m2_11, C2m2_12, C2m2_13, C2m2_22, C2m2_23, C2m2_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: V, Vt 

      V = A*sin(2.d0*pi*F*T)
      Vt = 2*A*cos(2*pi*F*t)*pi*F

      C2m2_11 = 2.D0/15.D0*sqrt(30.D0)/R*sqrt(pi)*M
      C2m2_12 = dcmplx(0.D0,2.D0/15.D0)*sqrt(pi)*M/R*sqrt(30.D0)
      C2m2_13 = 0
      C2m2_22 = -2.D0/15.D0*sqrt(30.D0)/R*sqrt(pi)*M
      C2m2_23 = 0
      C2m2_33 = 0

  end subroutine vibSchClm_2m2


  subroutine vibSchClm_2m1(R,M,T,F,A,&
             C2m1_11, C2m1_12, C2m1_13, C2m1_22, C2m1_23, C2m1_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: C2m1_11, C2m1_12, C2m1_13, C2m1_22, C2m1_23, C2m1_33 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: V, Vt 

      V = A*sin(2.d0*pi*F*T)
      Vt = 2*A*cos(2*pi*F*t)*pi*F

      C2m1_11 = 0
      C2m1_12 = 0
      C2m1_13 = 2.D0/15.D0*sqrt(30.D0)/R*sqrt(pi)*M
      C2m1_22 = 0
      C2m1_23 = dcmplx(0.D0,2.D0/15.D0)*sqrt(pi)*M/R*sqrt(30.D0)
      C2m1_33 = 0

  end subroutine vibSchClm_2m1


  subroutine vibSchClm_20(R,M,T,F,A,&
             C20_11, C20_12, C20_13, C20_22, C20_23, C20_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: C20_11, C20_12, C20_13, C20_22, C20_23, C20_33 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: V, Vt 

      V = A*sin(2.d0*pi*F*T)
      Vt = 2*A*cos(2*pi*f*t)*pi*f

      C20_11 = -4.D0/15.D0*sqrt(pi)*M*sqrt(5.D0)/R
      C20_12 = 0
      C20_13 = 0
      C20_22 = -4.D0/15.D0*sqrt(pi)*M*sqrt(5.D0)/R
      C20_23 = 0
      C20_33 = 8.D0/15.D0*sqrt(pi)*M*sqrt(5.D0)/R

  end subroutine vibSchClm_20


  subroutine vibSchClm_21(R,M,T,F,A,&
             C21_11, C21_12, C21_13, C21_22, C21_23, C21_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: C21_11, C21_12, C21_13, C21_22, C21_23, C21_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: V, Vt 

      V = A*sin(2.d0*pi*F*T)
      Vt = 2*A*cos(2*pi*f*t)*pi*f

      C21_11 = 0
      C21_12 = 0
      C21_13 = -2.D0/15.D0*sqrt(30.D0)/R*sqrt(pi)*M
      C21_22 = 0
      C21_23 = dcmplx(0.D0,2.D0/15.D0)*sqrt(pi)*M/R*sqrt(30.D0)
      C21_33 = 0

  end subroutine vibSchClm_21


  subroutine vibSchClm_22(R,M,T,F,A,&
             C22_11, C22_12, C22_13, C22_22, C22_23, C22_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: C22_11, C22_12, C22_13, C22_22, C22_23, C22_33 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: V, Vt 

      V = A*sin(2.d0*pi*F*T)
      Vt = 2*A*cos(2*pi*f*t)*pi*f

      C22_11 = 2.D0/15.D0*sqrt(30.D0)/R*sqrt(pi)*M
      C22_12 = dcmplx(0.D0,-2.D0/15.D0)*sqrt(pi)*M/R*sqrt(30.D0)
      C22_13 = 0
      C22_22 = -2.D0/15.D0*sqrt(30.D0)/R*sqrt(pi)*M
      C22_23 = 0
      C22_33 = 0

  end subroutine vibSchClm_22

!time derivatives of spherical harmonic coefficients

  subroutine vibSchdtAlm_00(R,M,T,F,A, dtA00)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: dtA00 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: V, Vt, Vtt 

      V = A*sin(2.d0*pi*F*T)
      Vt = 2*A*cos(2*pi*f*t)*pi*f
      Vtt = -4*A*sin(2*pi*f*t)*pi**2*f**2

      dtA00 = 0 !2/sqrt(R+2*M)*sqrt(R)*sqrt(pi)

  end subroutine vibSchdtAlm_00  


  subroutine vibSchdtClm_00(R,M,T,F,A,&
             dtC00_11, dtC00_12, dtC00_13, dtC00_22, dtC00_23, dtC00_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: dtC00_11, dtC00_12, dtC00_13, dtC00_22, dtC00_23, dtC00_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: V, Vt, Vtt 

      V = A*sin(2.d0*pi*F*T)
      Vt = 2*A*cos(2*pi*f*t)*pi*f
      Vtt = -4*A*sin(2*pi*f*t)*pi**2*f**2

      dtC00_11 = 0 !2.D0/3.D0*sqrt(pi)*(2*M+3*R)/R
      dtC00_12 = 0
      dtC00_13 = 0
      dtC00_22 = 0 !2.D0/3.D0*sqrt(pi)*(2*M+3*R)/R
      dtC00_23 = 0
      dtC00_33 = 0 !2.D0/3.D0*sqrt(pi)*(2*M+3*R)/R

  end subroutine vibSchdtClm_00


  subroutine vibSchdtBlm_1m1(R,M,T,F,A, dtB1m1_1, dtB1m1_2, dtB1m1_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: dtB1m1_1, dtB1m1_2, dtB1m1_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: V, Vt, Vtt 

      V = A*sin(2.d0*pi*F*T)
      Vt = 2*A*cos(2*pi*f*t)*pi*f
      Vtt = -4*A*sin(2*pi*f*t)*pi**2*f**2

      dtB1m1_1 = 0 !-2.D0/3.D0*sqrt(pi)*M/(R+2*M)*sqrt(6.D0)
      dtB1m1_2 = 0 !dcmplx(0.D0,-2.D0/3.D0)*sqrt(pi)*M/(R+2*M)*sqrt(6.D0)
      dtB1m1_3 = 0

  end subroutine vibSchdtBlm_1m1


  subroutine vibSchdtBlm_10(R,M,T,F,A, dtB10_1, dtB10_2, dtB10_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: dtB10_1, dtB10_2, dtB10_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: V, Vt 

      V = A*sin(2.d0*pi*F*T)
      Vt = 2*A*cos(2*pi*f*t)*pi*f

      dtB10_1 = 0
      dtB10_2 = 0
      dtB10_3 = 0 !-4.D0/3.D0*M/(R+2*M)*sqrt(3.D0)*sqrt(pi)

  end subroutine vibSchdtBlm_10


  subroutine vibSchdtBlm_11(R,M,T,F,A, dtB11_1, dtB11_2, dtB11_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: dtB11_1, dtB11_2, dtB11_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: V, Vt, Vtt 

      V = A*sin(2.d0*pi*F*T)
      Vt = 2*A*cos(2*pi*f*t)*pi*f
      Vtt = -4*A*sin(2*pi*f*t)*pi**2*f**2

      dtB11_1 = 0 !2.D0/3.D0*sqrt(pi)*M/(R+2*M)*sqrt(6.D0)
      dtB11_2 = 0 !dcmplx(0.D0,-2.D0/3.D0)*sqrt(pi)*M/(R+2*M)*sqrt(6.D0)
      dtB11_3 = 0

  end subroutine vibSchdtBlm_11


  subroutine vibSchdtClm_2m2(R,M,T,F,A,&
             dtC2m2_11, dtC2m2_12, dtC2m2_13, dtC2m2_22, dtC2m2_23, dtC2m2_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: dtC2m2_11, dtC2m2_12, dtC2m2_13, dtC2m2_22, dtC2m2_23, dtC2m2_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: V, Vt, Vtt 

      V = A*sin(2.d0*pi*F*T)
      Vt = 2*A*cos(2*pi*f*t)*pi*f
      Vtt = -4*A*sin(2*pi*f*t)*pi**2*f**2

      dtC2m2_11 = 0 !2.D0/15.D0*sqrt(30.D0)/R*sqrt(pi)*M
      dtC2m2_12 = 0 !dcmplx(0.D0,2.D0/15.D0)*sqrt(pi)*M/R*sqrt(30.D0)
      dtC2m2_13 = 0
      dtC2m2_22 = 0 !-2.D0/15.D0*sqrt(30.D0)/R*sqrt(pi)*M
      dtC2m2_23 = 0
      dtC2m2_33 = 0

  end subroutine vibSchdtClm_2m2


  subroutine vibSchdtClm_2m1(R,M,T,F,A,&
             dtC2m1_11, dtC2m1_12, dtC2m1_13, dtC2m1_22, dtC2m1_23, dtC2m1_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: dtC2m1_11, dtC2m1_12, dtC2m1_13, dtC2m1_22, dtC2m1_23, dtC2m1_33 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: V, Vt, Vtt 

      V = A*sin(2.d0*pi*F*T)
      Vt = 2*A*cos(2*pi*f*t)*pi*f
      Vtt = -4*A*sin(2*pi*f*t)*pi**2*f**2

      dtC2m1_11 = 0
      dtC2m1_12 = 0
      dtC2m1_13 = 0 !2.D0/15.D0*sqrt(30.D0)/R*sqrt(pi)*M
      dtC2m1_22 = 0
      dtC2m1_23 = 0 !dcmplx(0.D0,2.D0/15.D0)*sqrt(pi)*M/R*sqrt(30.D0)
      dtC2m1_33 = 0

  end subroutine vibSchdtClm_2m1


  subroutine vibSchdtClm_20(R,M,T,F,A,&
             dtC20_11, dtC20_12, dtC20_13, dtC20_22, dtC20_23, dtC20_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: dtC20_11, dtC20_12, dtC20_13, dtC20_22, dtC20_23, dtC20_33 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: V, Vt, Vtt 

      V = A*sin(2.d0*pi*F*T)
      Vt = 2*A*cos(2*pi*f*t)*pi*f
      Vtt = -4*A*sin(2*pi*f*t)*pi**2*f**2

      dtC20_11 = 0 !-4.D0/15.D0*sqrt(pi)*M*sqrt(5.D0)/R
      dtC20_12 = 0
      dtC20_13 = 0
      dtC20_22 = 0 !-4.D0/15.D0*sqrt(pi)*M*sqrt(5.D0)/R
      dtC20_23 = 0
      dtC20_33 = 0 !8.D0/15.D0*sqrt(pi)*M*sqrt(5.D0)/R

  end subroutine vibSchdtClm_20


  subroutine vibSchdtClm_21(R,M,T,F,A,&
             dtC21_11, dtC21_12, dtC21_13, dtC21_22, dtC21_23, dtC21_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: dtC21_11, dtC21_12, dtC21_13, dtC21_22, dtC21_23, dtC21_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: V, Vt, Vtt 

      V = A*sin(2.d0*pi*F*T)
      Vt = 2*A*cos(2*pi*f*t)*pi*f
      Vtt = -4*A*sin(2*pi*f*t)*pi**2*f**2

      dtC21_11 = 0
      dtC21_12 = 0
      dtC21_13 = 0 !-2.D0/15.D0*sqrt(30.D0)/R*sqrt(pi)*M
      dtC21_22 = 0
      dtC21_23 = 0 !dcmplx(0.D0,2.D0/15.D0)*sqrt(pi)*M/R*sqrt(30.D0)
      dtC21_33 = 0

  end subroutine vibSchdtClm_21


  subroutine vibSchdtClm_22(R,M,T,F,A,&
             dtC22_11, dtC22_12, dtC22_13, dtC22_22, dtC22_23, dtC22_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: dtC22_11, dtC22_12, dtC22_13, dtC22_22, dtC22_23, dtC22_33 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: V, Vt, Vtt 

      V = A*sin(2.d0*pi*F*T)
      Vt = 2*A*cos(2*pi*f*t)*pi*f
      Vtt = -4*A*sin(2*pi*f*t)*pi**2*f**2

      dtC22_11 = 0 !2.D0/15.D0*sqrt(30.D0)/R*sqrt(pi)*M
      dtC22_12 = 0 !dcmplx(0.D0,-2.D0/15.D0)*sqrt(pi)*M/R*sqrt(30.D0)
      dtC22_13 = 0
      dtC22_22 = 0 !-2.D0/15.D0*sqrt(30.D0)/R*sqrt(pi)*M
      dtC22_23 = 0
      dtC22_33 = 0

  end subroutine vibSchdtClm_22


!radial derivatives of spherical harmonic coefficients
  subroutine vibSchdrAlm_00(R,M,T,F,A, drA00)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: drA00 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: V, Vt 

      V = A*sin(2.d0*pi*F*T)
      Vt = 2*A*cos(2*pi*f*t)*pi*f

      drA00 = 2*sqrt(pi)*M/sqrt(R+2*M)**3/sqrt(R) 

  end subroutine vibSchdrAlm_00  


  subroutine vibSchdrClm_00(R,M,T,F,A,&
             drC00_11, drC00_12, drC00_13, drC00_22, drC00_23, drC00_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: drC00_11, drC00_12, drC00_13, drC00_22, drC00_23, drC00_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: V, Vt 

      V = A*sin(2.d0*pi*F*T)
      Vt = 2*A*cos(2*pi*f*t)*pi*f


      drC00_11 = -4.D0/3.D0*sqrt(pi)*M/R**2
      drC00_12 = 0
      drC00_13 = 0
      drC00_22 = -4.D0/3.D0*sqrt(pi)*M/R**2
      drC00_23 = 0
      drC00_33 = -4.D0/3.D0*sqrt(pi)*M/R**2

  end subroutine vibSchdrClm_00


  subroutine vibSchdrBlm_1m1(R,M,T,F,A, drB1m1_1, drB1m1_2, drB1m1_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: drB1m1_1, drB1m1_2, drB1m1_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: V, Vt 

      V = A*sin(2.d0*pi*F*T)
      Vt = 2*A*cos(2*pi*f*t)*pi*f

      drB1m1_1 = 2.D0/3.D0*M/(R+2*M)**2*sqrt(6.D0)*sqrt(pi)
      drB1m1_2 = dcmplx(0.D0,2.D0/3.D0)*M/(R+2*M)**2*sqrt(6.D0)*sqrt(pi)
      drB1m1_3 = 0

  end subroutine vibSchdrBlm_1m1


  subroutine vibSchdrBlm_10(R,M,T,F,A, drB10_1, drB10_2, drB10_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: drB10_1, drB10_2, drB10_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: V, Vt 

      V = A*sin(2.d0*pi*F*T)
      Vt = 2*A*cos(2*pi*f*t)*pi*f


      drB10_1 = 0
      drB10_2 = 0
      drB10_3 = 4.D0/3.D0*M/(R+2*M)**2*sqrt(3.D0)*sqrt(pi)

  end subroutine vibSchdrBlm_10


  subroutine vibSchdrBlm_11(R,M,T,F,A, drB11_1, drB11_2, drB11_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: drB11_1, drB11_2, drB11_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: V, Vt 

      V = A*sin(2.d0*pi*F*T)
      Vt = 2*A*cos(2*pi*f*t)*pi*f

      drB11_1 = -2.D0/3.D0*M/(R+2*M)**2*sqrt(6.D0)*sqrt(pi)
      drB11_2 = dcmplx(0.D0,2.D0/3.D0)*M/(R+2*M)**2*sqrt(6.D0)*sqrt(pi)
      drB11_3 = 0

  end subroutine vibSchdrBlm_11


  subroutine vibSchdrClm_2m2(R,M,T,F,A,&
             drC2m2_11, drC2m2_12, drC2m2_13, drC2m2_22, drC2m2_23, drC2m2_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: drC2m2_11, drC2m2_12, drC2m2_13,&
                      drC2m2_22, drC2m2_23, drC2m2_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: V, Vt 

      V = A*sin(2.d0*pi*F*T)
      Vt = 2*A*cos(2*pi*f*t)*pi*f

      drC2m2_11 = -2.D0/15.D0*sqrt(30.D0)/R**2*sqrt(pi)*M
      drC2m2_12 = dcmplx(0.D0,-2.D0/15.D0)*sqrt(30.D0)/R**2*sqrt(pi)*M
      drC2m2_13 = 0
      drC2m2_22 = 2.D0/15.D0*sqrt(30.D0)/R**2*sqrt(pi)*M
      drC2m2_23 = 0
      drC2m2_33 = 0

  end subroutine vibSchdrClm_2m2


  subroutine vibSchdrClm_2m1(R,M,T,F,A,&
             drC2m1_11, drC2m1_12, drC2m1_13, drC2m1_22, drC2m1_23, drC2m1_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: drC2m1_11, drC2m1_12, drC2m1_13, &
                      drC2m1_22, drC2m1_23, drC2m1_33 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: V, Vt 

      V = A*sin(2.d0*pi*F*T)
      Vt = 2*A*cos(2*pi*f*t)*pi*f

      drC2m1_11 = 0
      drC2m1_12 = 0
      drC2m1_13 = -2.D0/15.D0*sqrt(30.D0)/R**2*sqrt(pi)*M
      drC2m1_22 = 0
      drC2m1_23 = dcmplx(0.D0,-2.D0/15.D0)*sqrt(30.D0)/R**2*sqrt(pi)*M
      drC2m1_33 = 0

  end subroutine vibSchdrClm_2m1


  subroutine vibSchdrClm_20(R,M,T,F,A,&
             drC20_11, drC20_12, drC20_13, drC20_22, drC20_23, drC20_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: drC20_11, drC20_12, drC20_13,&
                      drC20_22, drC20_23, drC20_33 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: V, Vt 

      V = A*sin(2.d0*pi*F*T)
      Vt = 2*A*cos(2*pi*f*t)*pi*f

      drC20_11 = 4.D0/15.D0*sqrt(pi)*M*sqrt(5.D0)/R**2
      drC20_12 = 0
      drC20_13 = 0
      drC20_22 = 4.D0/15.D0*sqrt(pi)*M*sqrt(5.D0)/R**2
      drC20_23 = 0
      drC20_33 = -8.D0/15.D0*sqrt(pi)*M*sqrt(5.D0)/R**2

  end subroutine vibSchdrClm_20


  subroutine vibSchdrClm_21(R,M,T,F,A,&
             drC21_11, drC21_12, drC21_13, drC21_22, drC21_23, drC21_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: drC21_11, drC21_12, drC21_13,&
                      drC21_22, drC21_23, drC21_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: V, Vt 

      V = A*sin(2.d0*pi*F*T)
      Vt = 2*A*cos(2*pi*f*t)*pi*f

      drC21_11 = 0
      drC21_12 = 0
      drC21_13 = 2.D0/15.D0*sqrt(30.D0)/R**2*sqrt(pi)*M
      drC21_22 = 0
      drC21_23 = dcmplx(0.D0,-2.D0/15.D0)*sqrt(30.D0)/R**2*sqrt(pi)*M
      drC21_33 = 0

  end subroutine vibSchdrClm_21


  subroutine vibSchdrClm_22(R,M,T,F,A,&
             drC22_11, drC22_12, drC22_13, drC22_22, drC22_23, drC22_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M, T, F, A
      CCTK_COMPLEX :: drC22_11, drC22_12, drC22_13,&
                      drC22_22, drC22_23, drC22_33 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
      CCTK_REAL :: V, Vt 

      V = A*sin(2.d0*pi*F*T)
      Vt = 2*A*cos(2*pi*f*t)*pi*f

      drC22_11 = -2.D0/15.D0*sqrt(30.D0)/R**2*sqrt(pi)*M
      drC22_12 = dcmplx(0.D0,2.D0/15.D0)*sqrt(30.D0)/R**2*sqrt(pi)*M
      drC22_13 = 0
      drC22_22 = 2.D0/15.D0*sqrt(30.D0)/R**2*sqrt(pi)*M
      drC22_23 = 0
      drC22_33 = 0

  end subroutine vibSchdrClm_22


end module NullSHRE_modvibSchClm

