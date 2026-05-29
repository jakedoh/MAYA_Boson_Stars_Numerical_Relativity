! vim: syntax=fortran
#include "cctk.h"

module NullSHRE_modSchClm

  use cctk
  implicit none

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


  subroutine SchClm_00(R,M,&
             C00_11, C00_12, C00_13, C00_22, C00_23, C00_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: C00_11, C00_12, C00_13, C00_22, C00_23, C00_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      C00_11 = 2.D0/3.D0*sqrt(pi)*(2*M+3*R)/R
      C00_12 = 0
      C00_13 = 0
      C00_22 = 2.D0/3.D0*sqrt(pi)*(2*M+3*R)/R
      C00_23 = 0
      C00_33 = 2.D0/3.D0*sqrt(pi)*(2*M+3*R)/R

  end subroutine SchClm_00


  subroutine SchBlm_1m1(R,M, B1m1_1, B1m1_2, B1m1_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: B1m1_1, B1m1_2, B1m1_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      B1m1_1 = -2.D0/3.D0*sqrt(pi)*M/(R+2*M)*sqrt(6.D0)
      B1m1_2 = dcmplx(0.D0,-2.D0/3.D0)*sqrt(pi)*M/(R+2*M)*sqrt(6.D0)
      B1m1_3 = 0

  end subroutine SchBlm_1m1


  subroutine SchBlm_10(R,M, B10_1, B10_2, B10_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: B10_1, B10_2, B10_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      B10_1 = 0
      B10_2 = 0
      B10_3 = -4.D0/3.D0*M/(R+2*M)*sqrt(3.D0)*sqrt(pi)

  end subroutine SchBlm_10


  subroutine SchBlm_11(R,M, B11_1, B11_2, B11_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: B11_1, B11_2, B11_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      B11_1 = 2.D0/3.D0*sqrt(pi)*M/(R+2*M)*sqrt(6.D0)
      B11_2 = dcmplx(0.D0,-2.D0/3.D0)*sqrt(pi)*M/(R+2*M)*sqrt(6.D0)
      B11_3 = 0

  end subroutine SchBlm_11


  subroutine SchClm_2m2(R,M,&
             C2m2_11, C2m2_12, C2m2_13, C2m2_22, C2m2_23, C2m2_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: C2m2_11, C2m2_12, C2m2_13, C2m2_22, C2m2_23, C2m2_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      C2m2_11 = 2.D0/15.D0*sqrt(30.D0)/R*sqrt(pi)*M
      C2m2_12 = dcmplx(0.D0,2.D0/15.D0)*sqrt(pi)*M/R*sqrt(30.D0)
      C2m2_13 = 0
      C2m2_22 = -2.D0/15.D0*sqrt(30.D0)/R*sqrt(pi)*M
      C2m2_23 = 0
      C2m2_33 = 0

  end subroutine SchClm_2m2


  subroutine SchClm_2m1(R,M,&
             C2m1_11, C2m1_12, C2m1_13, C2m1_22, C2m1_23, C2m1_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: C2m1_11, C2m1_12, C2m1_13, C2m1_22, C2m1_23, C2m1_33 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      C2m1_11 = 0
      C2m1_12 = 0
      C2m1_13 = 2.D0/15.D0*sqrt(30.D0)/R*sqrt(pi)*M
      C2m1_22 = 0
      C2m1_23 = dcmplx(0.D0,2.D0/15.D0)*sqrt(pi)*M/R*sqrt(30.D0)
      C2m1_33 = 0

  end subroutine SchClm_2m1


  subroutine SchClm_20(R,M,&
             C20_11, C20_12, C20_13, C20_22, C20_23, C20_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: C20_11, C20_12, C20_13, C20_22, C20_23, C20_33 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      C20_11 = -4.D0/15.D0*sqrt(pi)*M*sqrt(5.D0)/R
      C20_12 = 0
      C20_13 = 0
      C20_22 = -4.D0/15.D0*sqrt(pi)*M*sqrt(5.D0)/R
      C20_23 = 0
      C20_33 = 8.D0/15.D0*sqrt(pi)*M*sqrt(5.D0)/R

  end subroutine SchClm_20


  subroutine SchClm_21(R,M,&
             C21_11, C21_12, C21_13, C21_22, C21_23, C21_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: C21_11, C21_12, C21_13, C21_22, C21_23, C21_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      C21_11 = 0
      C21_12 = 0
      C21_13 = -2.D0/15.D0*sqrt(30.D0)/R*sqrt(pi)*M
      C21_22 = 0
      C21_23 = dcmplx(0.D0,2.D0/15.D0)*sqrt(pi)*M/R*sqrt(30.D0)
      C21_33 = 0

  end subroutine SchClm_21


  subroutine SchClm_22(R,M,&
             C22_11, C22_12, C22_13, C22_22, C22_23, C22_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: C22_11, C22_12, C22_13, C22_22, C22_23, C22_33 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      C22_11 = 2.D0/15.D0*sqrt(30.D0)/R*sqrt(pi)*M
      C22_12 = dcmplx(0.D0,-2.D0/15.D0)*sqrt(pi)*M/R*sqrt(30.D0)
      C22_13 = 0
      C22_22 = -2.D0/15.D0*sqrt(30.D0)/R*sqrt(pi)*M
      C22_23 = 0
      C22_33 = 0

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


  subroutine SchdrClm_00(R,M,&
             drC00_11, drC00_12, drC00_13, drC00_22, drC00_23, drC00_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: drC00_11, drC00_12, drC00_13, drC00_22, drC00_23, drC00_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0


      drC00_11 = -4.D0/3.D0*sqrt(pi)*M/R**2
      drC00_12 = 0
      drC00_13 = 0
      drC00_22 = -4.D0/3.D0*sqrt(pi)*M/R**2
      drC00_23 = 0
      drC00_33 = -4.D0/3.D0*sqrt(pi)*M/R**2

  end subroutine SchdrClm_00


  subroutine SchdrBlm_1m1(R,M, drB1m1_1, drB1m1_2, drB1m1_3)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: drB1m1_1, drB1m1_2, drB1m1_3 
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      drB1m1_1 = 2.D0/3.D0*M/(R+2*M)**2*sqrt(6.D0)*sqrt(pi)
      drB1m1_2 = dcmplx(0.D0,2.D0/3.D0)*M/(R+2*M)**2*sqrt(6.D0)*sqrt(pi)
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
      drB11_2 = dcmplx(0.D0,2.D0/3.D0)*M/(R+2*M)**2*sqrt(6.D0)*sqrt(pi)
      drB11_3 = 0

  end subroutine SchdrBlm_11


  subroutine SchdrClm_2m2(R,M,&
             drC2m2_11, drC2m2_12, drC2m2_13, drC2m2_22, drC2m2_23, drC2m2_33)

      implicit none

      CCTK_REAL, intent(in) :: R, M
      CCTK_COMPLEX :: drC2m2_11, drC2m2_12, drC2m2_13,&
                      drC2m2_22, drC2m2_23, drC2m2_33
      CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0

      drC2m2_11 = -2.D0/15.D0*sqrt(30.D0)/R**2*sqrt(pi)*M
      drC2m2_12 = dcmplx(0.D0,-2.D0/15.D0)*sqrt(30.D0)/R**2*sqrt(pi)*M
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
      drC2m1_23 = dcmplx(0.D0,-2.D0/15.D0)*sqrt(30.D0)/R**2*sqrt(pi)*M
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
      drC21_23 = dcmplx(0.D0,-2.D0/15.D0)*sqrt(30.D0)/R**2*sqrt(pi)*M
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
      drC22_12 = dcmplx(0.D0,2.D0/15.D0)*sqrt(30.D0)/R**2*sqrt(pi)*M
      drC22_13 = 0
      drC22_22 = 2.D0/15.D0*sqrt(30.D0)/R**2*sqrt(pi)*M
      drC22_23 = 0
      drC22_33 = 0

  end subroutine SchdrClm_22

end module NullSHRE_modSchClm

