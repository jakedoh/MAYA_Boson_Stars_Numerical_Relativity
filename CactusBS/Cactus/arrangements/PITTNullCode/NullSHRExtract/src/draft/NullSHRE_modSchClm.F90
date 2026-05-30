! vim: syntax=fortran
#include "cctk.h"

module NullSHRE_modSchClm

  ! gives the numerical spherical harmonic coefficients Clm for the IEF Schwarzschild metric 

contains

  subroutine SchAlm_00(R,M, A00)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: A00 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      A00 = 2/sqrt(R+2*M)*sqrt(R)*sqrt(pi)

  end subroutine SchAlm_00  

  subroutine SchAlm_1m1(R,M, A1m1)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: A1m1 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      A1m1 = 0.d0 

  end subroutine SchAlm_1m1  

  subroutine SchAlm_10(R,M, A10)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: A10 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      A10 = 0.d0

  end subroutine SchAlm_10  

  subroutine SchAlm_11(R,M, A11)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: A11 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      A11 = 0.d0

  end subroutine SchAlm_11  

  subroutine SchAlm_2m2(R,M, A2m2)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: A2m2 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      A2m2 = 0.d0 

  end subroutine SchAlm_2m2  

  subroutine SchAlm_2m1(R,M, A2m1)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: A2m1 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      A2m1 = 0.d0 

  end subroutine SchAlm_2m1  

  subroutine SchAlm_20(R,M, A20)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: A20 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      A20 = 0.d0

  end subroutine SchAlm_20  

  subroutine SchAlm_21(R,M, A21)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: A21 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      A21 = 0.d0

  end subroutine SchAlm_21  

  subroutine SchAlm_22(R,M, A22)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: A22 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      A22 = 0.d0

  end subroutine SchAlm_22  

  subroutine SchBlm_00(R,M, B00_1, B00_2, B00_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: B00_1, B00_2, B00_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      B00_1 = 0.d0 
      B00_2 = 0.d0
      B00_3 = 0.d0

  end subroutine SchBlm_00

  subroutine SchBlm_1m1(R,M, B1m1_1, B1m1_2, B1m1_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: B1m1_1, B1m1_2, B1m1_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      B1m1_1 = -2.D0/3.D0*sqrt(pi)*M/(R+2*M)*sqrt(6.D0)
      B1m1_2 = cmplx(0.D0,-2.D0/3.D0)*sqrt(pi)*M/(R+2*M)*sqrt(6.D0)
      B1m1_3 = 0.d0

  end subroutine SchBlm_1m1


  subroutine SchBlm_10(R,M, B10_1, B10_2, B10_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: B10_1, B10_2, B10_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      B10_1 = 0.d0
      B10_2 = 0.d0
      B10_3 = -4.D0/3.D0*M/(R+2*M)*sqrt(3.D0)*sqrt(pi)

  end subroutine SchBlm_10


  subroutine SchBlm_11(R,M, B11_1, B11_2, B11_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: B11_1, B11_2, B11_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      B11_1 = 2.D0/3.D0*sqrt(pi)*M/(R+2*M)*sqrt(6.D0)
      B11_2 = cmplx(0.D0,-2.D0/3.D0)*sqrt(pi)*M/(R+2*M)*sqrt(6.D0)
      B11_3 = 0.d0

  end subroutine SchBlm_11

  subroutine SchBlm_2m2(R,M, B2m2_1, B2m2_2, B2m2_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: B2m2_1, B2m2_2, B2m2_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      B2m2_1 = 0.d0 
      B2m2_2 = 0.d0
      B2m2_3 = 0.d0

  end subroutine SchBlm_2m2

  subroutine SchBlm_2m1(R,M, B2m1_1, B2m1_2, B2m1_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: B2m1_1, B2m1_2, B2m1_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      B2m1_1 = 0.d0 
      B2m1_2 = 0.d0
      B2m1_3 = 0.d0

  end subroutine SchBlm_2m1

  subroutine SchBlm_20(R,M, B20_1, B20_2, B20_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: B20_1, B20_2, B20_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      B20_1 = 0.d0 
      B20_2 = 0.d0
      B20_3 = 0.d0

  end subroutine SchBlm_20

  subroutine SchBlm_21(R,M, B21_1, B21_2, B21_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: B21_1, B21_2, B21_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      B21_1 = 0.d0 
      B21_2 = 0.d0
      B21_3 = 0.d0

  end subroutine SchBlm_21

  subroutine SchBlm_22(R,M, B22_1, B22_2, B22_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: B22_1, B22_2, B22_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      B22_1 = 0.d0 
      B22_2 = 0.d0
      B22_3 = 0.d0

  end subroutine SchBlm_22

  subroutine SchClm_00(R,M,&
             C00_11, C00_12, C00_13, C00_22, C00_23, C00_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: C00_11, C00_12, C00_13, C00_22, C00_23, C00_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      C00_11 = 2.D0/3.D0*sqrt(pi)*(2*M+3*R)/R
      C00_12 = 0.d0
      C00_13 = 0.d0
      C00_22 = 2.D0/3.D0*sqrt(pi)*(2*M+3*R)/R
      C00_23 = 0.d0
      C00_33 = 2.D0/3.D0*sqrt(pi)*(2*M+3*R)/R

  end subroutine SchClm_00

  subroutine SchClm_1m1(R,M,&
             C1m1_11, C1m1_12, C1m1_13, C1m1_22, C1m1_23, C1m1_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: C1m1_11, C1m1_12, C1m1_13, C1m1_22, C1m1_23, C1m1_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      C1m1_11 = 0.d0 
      C1m1_12 = 0.d0
      C1m1_13 = 0.d0
      C1m1_22 = 0.d0 
      C1m1_23 = 0.d0
      C1m1_33 = 0.d0 

  end subroutine SchClm_1m1

  subroutine SchClm_10(R,M,&
             C10_11, C10_12, C10_13, C10_22, C10_23, C10_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: C10_11, C10_12, C10_13, C10_22, C10_23, C10_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      C10_11 = 0.d0 
      C10_12 = 0.d0
      C10_13 = 0.d0
      C10_22 = 0.d0 
      C10_23 = 0.d0
      C10_33 = 0.d0 

  end subroutine SchClm_10

  subroutine SchClm_11(R,M,&
             C11_11, C11_12, C11_13, C11_22, C11_23, C11_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: C11_11, C11_12, C11_13, C11_22, C11_23, C11_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      C11_11 = 0.d0 
      C11_12 = 0.d0
      C11_13 = 0.d0
      C11_22 = 0.d0 
      C11_23 = 0.d0
      C11_33 = 0.d0 

  end subroutine SchClm_11


  subroutine SchClm_2m2(R,M,&
             C2m2_11, C2m2_12, C2m2_13, C2m2_22, C2m2_23, C2m2_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: C2m2_11, C2m2_12, C2m2_13, C2m2_22, C2m2_23, C2m2_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      C2m2_11 = 2.D0/15.D0*sqrt(30.D0)/R*sqrt(pi)*M
      C2m2_12 = cmplx(0.D0,2.D0/15.D0)*sqrt(pi)*M/R*sqrt(30.D0)
      C2m2_13 = 0.d0
      C2m2_22 = -2.D0/15.D0*sqrt(30.D0)/R*sqrt(pi)*M
      C2m2_23 = 0.d0
      C2m2_33 = 0.d0

  end subroutine SchClm_2m2


  subroutine SchClm_2m1(R,M,&
             C2m1_11, C2m1_12, C2m1_13, C2m1_22, C2m1_23, C2m1_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: C2m1_11, C2m1_12, C2m1_13, C2m1_22, C2m1_23, C2m1_33 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      C2m1_11 = 0.d0
      C2m1_12 = 0.d0
      C2m1_13 = 2.D0/15.D0*sqrt(30.D0)/R*sqrt(pi)*M
      C2m1_22 = 0.d0
      C2m1_23 = cmplx(0.D0,2.D0/15.D0)*sqrt(pi)*M/R*sqrt(30.D0)
      C2m1_33 = 0.d0

  end subroutine SchClm_2m1


  subroutine SchClm_20(R,M,&
             C20_11, C20_12, C20_13, C20_22, C20_23, C20_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: C20_11, C20_12, C20_13, C20_22, C20_23, C20_33 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      C20_11 = -4.D0/15.D0*sqrt(pi)*M*sqrt(5.D0)/R
      C20_12 = 0.d0
      C20_13 = 0.d0
      C20_22 = -4.D0/15.D0*sqrt(pi)*M*sqrt(5.D0)/R
      C20_23 = 0.d0
      C20_33 = 8.D0/15.D0*sqrt(pi)*M*sqrt(5.D0)/R

  end subroutine SchClm_20


  subroutine SchClm_21(R,M,&
             C21_11, C21_12, C21_13, C21_22, C21_23, C21_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: C21_11, C21_12, C21_13, C21_22, C21_23, C21_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      C21_11 = 0.d0
      C21_12 = 0.d0
      C21_13 = -2.D0/15.D0*sqrt(30.D0)/R*sqrt(pi)*M
      C21_22 = 0.d0
      C21_23 = cmplx(0.D0,2.D0/15.D0)*sqrt(pi)*M/R*sqrt(30.D0)
      C21_33 = 0.d0

  end subroutine SchClm_21


  subroutine SchClm_22(R,M,&
             C22_11, C22_12, C22_13, C22_22, C22_23, C22_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: C22_11, C22_12, C22_13, C22_22, C22_23, C22_33 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      C22_11 = 2.D0/15.D0*sqrt(30.D0)/R*sqrt(pi)*M
      C22_12 = cmplx(0.D0,-2.D0/15.D0)*sqrt(pi)*M/R*sqrt(30.D0)
      C22_13 = 0.d0
      C22_22 = -2.D0/15.D0*sqrt(30.D0)/R*sqrt(pi)*M
      C22_23 = 0.d0
      C22_33 = 0.d0

  end subroutine SchClm_22

!radial derivatives of spherical harmonic coefficients
  subroutine SchdrAlm_00(R,M, drA00)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: drA00 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      drA00 = 2*sqrt(pi)*M/sqrt(R+2*M)**3/sqrt(R) 

  end subroutine SchdrAlm_00  

  subroutine SchdrAlm_1m1(R,M, drA1m1)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: drA1m1 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      drA1m1 = 0.d0 

  end subroutine SchdrAlm_1m1  

  subroutine SchdrAlm_10(R,M, drA10)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: drA10
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      drA10 = 0.d0 

  end subroutine SchdrAlm_10  

  subroutine SchdrAlm_11(R,M, drA11)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: drA11
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      drA11 = 0.d0 

  end subroutine SchdrAlm_11  

  subroutine SchdrAlm_2m2(R,M, drA2m2)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: drA2m2
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      drA2m2 = 0.d0 

  end subroutine SchdrAlm_2m2  

  subroutine SchdrAlm_2m1(R,M, drA2m1)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: drA2m1
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      drA2m1 = 0.d0 

  end subroutine SchdrAlm_2m1  

  subroutine SchdrAlm_20(R,M, drA20)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: drA20
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      drA20 = 0.d0 

  end subroutine SchdrAlm_20  

  subroutine SchdrAlm_21(R,M, drA21)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: drA21
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      drA21 = 0.d0 

  end subroutine SchdrAlm_21  

  subroutine SchdrAlm_22(R,M, drA22)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: drA22
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      drA22 = 0.d0 

  end subroutine SchdrAlm_22  


  subroutine SchdrBlm_00(R,M, drB00_1, drB00_2, drB00_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: drB00_1, drB00_2, drB00_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      drB00_1 = 0.d0 
      drB00_2 = 0.d0
      drB00_3 = 0.d0

  end subroutine SchdrBlm_00


  subroutine SchdrBlm_1m1(R,M, drB1m1_1, drB1m1_2, drB1m1_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: drB1m1_1, drB1m1_2, drB1m1_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      drB1m1_1 = 2.D0/3.D0*M/(R+2*M)**2*sqrt(6.D0)*sqrt(pi)
      drB1m1_2 = cmplx(0.D0,2.D0/3.D0)*M/(R+2*M)**2*sqrt(6.D0)*sqrt(pi)
      drB1m1_3 = 0

  end subroutine SchdrBlm_1m1


  subroutine SchdrBlm_10(R,M, drB10_1, drB10_2, drB10_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: drB10_1, drB10_2, drB10_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0


      drB10_1 = 0
      drB10_2 = 0
      drB10_3 = 4.D0/3.D0*M/(R+2*M)**2*sqrt(3.D0)*sqrt(pi)

  end subroutine SchdrBlm_10


  subroutine SchdrBlm_11(R,M, drB11_1, drB11_2, drB11_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: drB11_1, drB11_2, drB11_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      drB11_1 = -2.D0/3.D0*M/(R+2*M)**2*sqrt(6.D0)*sqrt(pi)
      drB11_2 = cmplx(0.D0,2.D0/3.D0)*M/(R+2*M)**2*sqrt(6.D0)*sqrt(pi)
      drB11_3 = 0

  end subroutine SchdrBlm_11


  subroutine SchdrBlm_2m2(R,M, drB2m2_1, drB2m2_2, drB2m2_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: drB2m2_1, drB2m2_2, drB2m2_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      drB2m2_1 = 0.d0 
      drB2m2_2 = 0.d0
      drB2m2_3 = 0.d0

  end subroutine SchdrBlm_2m2

  subroutine SchdrBlm_2m1(R,M, drB2m1_1, drB2m1_2, drB2m1_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: drB2m1_1, drB2m1_2, drB2m1_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      drB2m1_1 = 0.d0 
      drB2m1_2 = 0.d0
      drB2m1_3 = 0.d0

  end subroutine SchdrBlm_2m1


  subroutine SchdrBlm_20(R,M, drB20_1, drB20_2, drB20_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: drB20_1, drB20_2, drB20_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      drB20_1 = 0.d0 
      drB20_2 = 0.d0
      drB20_3 = 0.d0

  end subroutine SchdrBlm_20


  subroutine SchdrBlm_21(R,M, drB21_1, drB21_2, drB21_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: drB21_1, drB21_2, drB21_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      drB21_1 = 0.d0 
      drB21_2 = 0.d0
      drB21_3 = 0.d0

  end subroutine SchdrBlm_21

  subroutine SchdrBlm_22(R,M, drB22_1, drB22_2, drB22_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: drB22_1, drB22_2, drB22_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      drB22_1 = 0.d0 
      drB22_2 = 0.d0
      drB22_3 = 0.d0

  end subroutine SchdrBlm_22


  subroutine SchdrClm_00(R,M,&
             drC00_11, drC00_12, drC00_13, drC00_22, drC00_23, drC00_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: drC00_11, drC00_12, drC00_13, drC00_22, drC00_23, drC00_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0


      drC00_11 = -4.D0/3.D0*sqrt(pi)*M/R**2
      drC00_12 = 0.d0
      drC00_13 = 0.d0
      drC00_22 = -4.D0/3.D0*sqrt(pi)*M/R**2
      drC00_23 = 0.d0
      drC00_33 = -4.D0/3.D0*sqrt(pi)*M/R**2

  end subroutine SchdrClm_00


  subroutine SchdrClm_1m1(R,M,&
             drC1m1_11, drC1m1_12, drC1m1_13, drC1m1_22, drC1m1_23, drC1m1_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: drC1m1_11, drC1m1_12, drC1m1_13, drC1m1_22, drC1m1_23, drC1m1_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0


      drC1m1_11 = 0.d0 
      drC1m1_12 = 0.d0
      drC1m1_13 = 0.d0
      drC1m1_22 = 0.d0 
      drC1m1_23 = 0.d0
      drC1m1_33 = 0.d0 

  end subroutine SchdrClm_1m1


  subroutine SchdrClm_10(R,M,&
             drC10_11, drC10_12, drC10_13, drC10_22, drC10_23, drC10_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: drC10_11, drC10_12, drC10_13, drC10_22, drC10_23, drC10_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0


      drC10_11 = 0.d0 
      drC10_12 = 0.d0
      drC10_13 = 0.d0
      drC10_22 = 0.d0 
      drC10_23 = 0.d0
      drC10_33 = 0.d0 

  end subroutine SchdrClm_10


  subroutine SchdrClm_11(R,M,&
             drC11_11, drC11_12, drC11_13, drC11_22, drC11_23, drC11_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: drC11_11, drC11_12, drC11_13, drC11_22, drC11_23, drC11_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0


      drC11_11 = 0.d0 
      drC11_12 = 0.d0
      drC11_13 = 0.d0
      drC11_22 = 0.d0 
      drC11_23 = 0.d0
      drC11_33 = 0.d0 

  end subroutine SchdrClm_11


  subroutine SchdrClm_2m2(R,M,&
             drC2m2_11, drC2m2_12, drC2m2_13, drC2m2_22, drC2m2_23, drC2m2_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: drC2m2_11, drC2m2_12, drC2m2_13,&
                      drC2m2_22, drC2m2_23, drC2m2_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      drC2m2_11 = -2.D0/15.D0*sqrt(30.D0)/R**2*sqrt(pi)*M
      drC2m2_12 = cmplx(0.D0,-2.D0/15.D0)*sqrt(30.D0)/R**2*sqrt(pi)*M
      drC2m2_13 = 0
      drC2m2_22 = 2.D0/15.D0*sqrt(30.D0)/R**2*sqrt(pi)*M
      drC2m2_23 = 0
      drC2m2_33 = 0

  end subroutine SchdrClm_2m2


  subroutine SchdrClm_2m1(R,M,&
             drC2m1_11, drC2m1_12, drC2m1_13, drC2m1_22, drC2m1_23, drC2m1_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: drC2m1_11, drC2m1_12, drC2m1_13, &
                      drC2m1_22, drC2m1_23, drC2m1_33 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      drC2m1_11 = 0
      drC2m1_12 = 0
      drC2m1_13 = -2.D0/15.D0*sqrt(30.D0)/R**2*sqrt(pi)*M
      drC2m1_22 = 0
      drC2m1_23 = cmplx(0.D0,-2.D0/15.D0)*sqrt(30.D0)/R**2*sqrt(pi)*M
      drC2m1_33 = 0

  end subroutine SchdrClm_2m1


  subroutine SchdrClm_20(R,M,&
             drC20_11, drC20_12, drC20_13, drC20_22, drC20_23, drC20_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: drC20_11, drC20_12, drC20_13,&
                      drC20_22, drC20_23, drC20_33 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      drC20_11 = 4.D0/15.D0*sqrt(pi)*M*sqrt(5.D0)/R**2
      drC20_12 = 0
      drC20_13 = 0
      drC20_22 = 4.D0/15.D0*sqrt(pi)*M*sqrt(5.D0)/R**2
      drC20_23 = 0
      drC20_33 = -8.D0/15.D0*sqrt(pi)*M*sqrt(5.D0)/R**2

  end subroutine SchdrClm_20


  subroutine SchdrClm_21(R,M,&
             drC21_11, drC21_12, drC21_13, drC21_22, drC21_23, drC21_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: drC21_11, drC21_12, drC21_13,&
                      drC21_22, drC21_23, drC21_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      drC21_11 = 0
      drC21_12 = 0
      drC21_13 = 2.D0/15.D0*sqrt(30.D0)/R**2*sqrt(pi)*M
      drC21_22 = 0
      drC21_23 = cmplx(0.D0,-2.D0/15.D0)*sqrt(30.D0)/R**2*sqrt(pi)*M
      drC21_33 = 0

  end subroutine SchdrClm_21


  subroutine SchdrClm_22(R,M,&
             drC22_11, drC22_12, drC22_13, drC22_22, drC22_23, drC22_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: drC22_11, drC22_12, drC22_13,&
                      drC22_22, drC22_23, drC22_33 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      drC22_11 = -2.D0/15.D0*sqrt(30.D0)/R**2*sqrt(pi)*M
      drC22_12 = cmplx(0.D0,2.D0/15.D0)*sqrt(30.D0)/R**2*sqrt(pi)*M
      drC22_13 = 0
      drC22_22 = 2.D0/15.D0*sqrt(30.D0)/R**2*sqrt(pi)*M
      drC22_23 = 0
      drC22_33 = 0

  end subroutine SchdrClm_22

end module NullSHRE_modSchClm

