! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"


subroutine NullSHRE_SchReadClm(CCTK_ARGUMENTS)

  use cctk
  use NullGrid_Vars
  use NullSHRE_modSchClm
  use NullSHRE_modsinSchClm
  use NullSHRE_modvibSchClm
  implicit none

    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    DECLARE_CCTK_FUNCTIONS
    CCTK_INT :: n, l, m

    CCTK_COMPLEX, dimension((l_max+1)*(l_max+1)) :: C1, C2, C3, C4, C5, &
                                                  C6, C7, C8, C9, C10
    CCTK_COMPLEX, dimension((l_max+1)*(l_max+1)) :: Ct1, Ct2, Ct3, Ct4, Ct5, &
                                                  Ct6, Ct7, Ct8, Ct9, Ct10
    CCTK_COMPLEX, dimension((l_max+1)*(l_max+1)) :: Cr1, Cr2, Cr3, Cr4, Cr5, &
                                                  Cr6, Cr7, Cr8, Cr9, Cr10

#define R cr
#define M mass 
#define T cctk_time
#define F fcoef 
#define A Afact 
 
    do l = 0, l_max
      do m = -l, l

         n = l*l + l + m + 1

         C1(n) = 0.d0; C2(n) = 0.d0
         C3(n) = 0.d0; C4(n) = 0.d0 
         C5(n) = 0.d0; C6(n) = 0.d0
         C7(n) = 0.d0; C8(n) = 0.d0
         C9(n) = 0.d0; C10(n) = 0.d0

         Ct1(n) = 0.d0; Ct2(n) = 0.d0
         Ct3(n) = 0.d0; Ct4(n) = 0.d0 
         Ct5(n) = 0.d0; Ct6(n) = 0.d0
         Ct7(n) = 0.d0; Ct8(n) = 0.d0
         Ct9(n) = 0.d0; Ct10(n) = 0.d0

         Cr1(n) = 0.d0; Cr2(n) = 0.d0
         Cr3(n) = 0.d0; Cr4(n) = 0.d0 
         Cr5(n) = 0.d0; Cr6(n) = 0.d0
         Cr7(n) = 0.d0; Cr8(n) = 0.d0
         Cr9(n) = 0.d0; Cr10(n) = 0.d0

      end do
    end do

    if (CCTK_EQUALS(SchIEF_time,"static")) then
 
       call SchAlm_00(R,M, C10(1))
       call SchBlm_1m1(R,M, C7(2), C8(2), C9(2))
       call SchBlm_10(R,M, C7(3), C8(3), C9(3))
       call SchBlm_11(R,M, C7(4), C8(4), C9(4))
       call SchClm_00(R,M, C1(1), C2(1), C3(1), C4(1), C5(1), C6(1))
       call SchClm_2m2(R,M, C1(5), C2(5), C3(5), C4(5), C5(5), C6(5))
       call SchClm_2m1(R,M, C1(6), C2(6), C3(6), C4(6), C5(6), C6(6))
       call SchClm_20(R,M, C1(7), C2(7), C3(7), C4(7), C5(7), C6(7))
       call SchClm_21(R,M, C1(8), C2(8), C3(8), C4(8), C5(8), C6(8))
       call SchClm_22(R,M, C1(9), C2(9), C3(9), C4(9), C5(9), C6(9))

       call SchdrAlm_00(R,M, Cr10(1))
       call SchdrBlm_1m1(R,M, Cr7(2), Cr8(2), Cr9(2))
       call SchdrBlm_10(R,M, Cr7(3), Cr8(3), Cr9(3))
       call SchdrBlm_11(R,M, Cr7(4), Cr8(4), Cr9(4))
       call SchdrClm_00(R,M, Cr1(1), Cr2(1), Cr3(1), Cr4(1), Cr5(1), Cr6(1))
       call SchdrClm_2m2(R,M, Cr1(5), Cr2(5), Cr3(5), Cr4(5), Cr5(5), Cr6(5))
       call SchdrClm_2m1(R,M, Cr1(6), Cr2(6), Cr3(6), Cr4(6), Cr5(6), Cr6(6))
       call SchdrClm_20(R,M, Cr1(7), Cr2(7), Cr3(7), Cr4(7), Cr5(7), Cr6(7))
       call SchdrClm_21(R,M, Cr1(8), Cr2(8), Cr3(8), Cr4(8), Cr5(8), Cr6(8))
       call SchdrClm_22(R,M, Cr1(9), Cr2(9), Cr3(9), Cr4(9), Cr5(9), Cr6(9))

    else if (CCTK_EQUALS(SchIEF_time, "sine_t")) then
 
       call sinSchAlm_00(R,M,T,F,A,C10(1))
       call sinSchBlm_1m1(R,M,T,F,A,C7(2),C8(2),C9(2))
       call sinSchBlm_10(R,M,T,F,A,C7(3),C8(3),C9(3))
       call sinSchBlm_11(R,M,T,F,A,C7(4),C8(4),C9(4))
       call sinSchClm_00(R,M,T,F,A,C1(1),C2(1),C3(1),C4(1),C5(1),C6(1))
       call sinSchClm_2m2(R,M,T,F,A,C1(5),C2(5),C3(5),C4(5),C5(5),C6(5))
       call sinSchClm_2m1(R,M,T,F,A,C1(6),C2(6),C3(6),C4(6),C5(6),C6(6))
       call sinSchClm_20(R,M,T,F,A,C1(7),C2(7),C3(7),C4(7),C5(7),C6(7))
       call sinSchClm_21(R,M,T,F,A,C1(8),C2(8),C3(8),C4(8),C5(8),C6(8))
       call sinSchClm_22(R,M,T,F,A,C1(9),C2(9),C3(9),C4(9),C5(9),C6(9))

       call sinSchdrAlm_00(R,M,T,F,A,Cr10(1))
       call sinSchdrBlm_1m1(R,M,T,F,A,Cr7(2),Cr8(2),Cr9(2))
       call sinSchdrBlm_10(R,M,T,F,A,Cr7(3),Cr8(3),Cr9(3))
       call sinSchdrBlm_11(R,M,T,F,A,Cr7(4),Cr8(4),Cr9(4))
       call sinSchdrClm_00(R,M,T,F,A,Cr1(1),Cr2(1),Cr3(1),Cr4(1),Cr5(1),Cr6(1))
       call sinSchdrClm_2m2(R,M,T,F,A,Cr1(5),Cr2(5),Cr3(5),Cr4(5),Cr5(5),Cr6(5))
       call sinSchdrClm_2m1(R,M,T,F,A,Cr1(6),Cr2(6),Cr3(6),Cr4(6),Cr5(6),Cr6(6))
       call sinSchdrClm_20(R,M,T,F,A,Cr1(7),Cr2(7),Cr3(7),Cr4(7),Cr5(7),Cr6(7))
       call sinSchdrClm_21(R,M,T,F,A,Cr1(8),Cr2(8),Cr3(8),Cr4(8),Cr5(8),Cr6(8))
       call sinSchdrClm_22(R,M,T,F,A,Cr1(9),Cr2(9),Cr3(9),Cr4(9),Cr5(9),Cr6(9))

       call sinSchdtAlm_00(R,M,T,F,A,Ct10(1))
       call sinSchdtBlm_1m1(R,M,T,F,A,Ct7(2),Ct8(2),Ct9(2))
       call sinSchdtBlm_10(R,M,T,F,A,Ct7(3),Ct8(3),Ct9(3))
       call sinSchdtBlm_11(R,M,T,F,A,Ct7(4),Ct8(4),Ct9(4))
       call sinSchdtClm_00(R,M,T,F,A,Ct1(1),Ct2(1),Ct3(1),Ct4(1),Ct5(1),Ct6(1))
       call sinSchdtClm_2m2(R,M,T,F,A,Ct1(5),Ct2(5),Ct3(5),Ct4(5),Ct5(5),Ct6(5))
       call sinSchdtClm_2m1(R,M,T,F,A,Ct1(6),Ct2(6),Ct3(6),Ct4(6),Ct5(6),Ct6(6))
       call sinSchdtClm_20(R,M,T,F,A,Ct1(7),Ct2(7),Ct3(7),Ct4(7),Ct5(7),Ct6(7))
       call sinSchdtClm_21(R,M,T,F,A,Ct1(8),Ct2(8),Ct3(8),Ct4(8),Ct5(8),Ct6(8))
       call sinSchdtClm_22(R,M,T,F,A,Ct1(9),Ct2(9),Ct3(9),Ct4(9),Ct5(9),Ct6(9))

    else if (CCTK_EQUALS(SchIEF_time, "vibe_t")) then
 
       call SchAlm_00(R,M, C10(1))
       call SchBlm_1m1(R,M, C7(2), C8(2), C9(2))
       call SchBlm_10(R,M, C7(3), C8(3), C9(3))
       call SchBlm_11(R,M, C7(4), C8(4), C9(4))
       call SchClm_00(R,M, C1(1), C2(1), C3(1), C4(1), C5(1), C6(1))
       call SchClm_2m2(R,M, C1(5), C2(5), C3(5), C4(5), C5(5), C6(5))
       call SchClm_2m1(R,M, C1(6), C2(6), C3(6), C4(6), C5(6), C6(6))
       call SchClm_20(R,M, C1(7), C2(7), C3(7), C4(7), C5(7), C6(7))
       call SchClm_21(R,M, C1(8), C2(8), C3(8), C4(8), C5(8), C6(8))
       call SchClm_22(R,M, C1(9), C2(9), C3(9), C4(9), C5(9), C6(9))
!       call vibSchAlm_00(R,M,T,F,A,C10(1))
!       call vibSchBlm_1m1(R,M,T,F,A,C7(2),C8(2),C9(2))
!       call vibSchBlm_10(R,M,T,F,A,C7(3),C8(3),C9(3))
!       call vibSchBlm_11(R,M,T,F,A,C7(4),C8(4),C9(4))
!       call vibSchClm_00(R,M,T,F,A,C1(1),C2(1),C3(1),C4(1),C5(1),C6(1))
!       call vibSchClm_2m2(R,M,T,F,A,C1(5),C2(5),C3(5),C4(5),C5(5),C6(5))
!       call vibSchClm_2m1(R,M,T,F,A,C1(6),C2(6),C3(6),C4(6),C5(6),C6(6))
!       call vibSchClm_20(R,M,T,F,A,C1(7),C2(7),C3(7),C4(7),C5(7),C6(7))
!       call vibSchClm_21(R,M,T,F,A,C1(8),C2(8),C3(8),C4(8),C5(8),C6(8))
!       call vibSchClm_22(R,M,T,F,A,C1(9),C2(9),C3(9),C4(9),C5(9),C6(9))

       call SchdrAlm_00(R,M, Cr10(1))
       call SchdrBlm_1m1(R,M, Cr7(2), Cr8(2), Cr9(2))
       call SchdrBlm_10(R,M, Cr7(3), Cr8(3), Cr9(3))
       call SchdrBlm_11(R,M, Cr7(4), Cr8(4), Cr9(4))
       call SchdrClm_00(R,M, Cr1(1), Cr2(1), Cr3(1), Cr4(1), Cr5(1), Cr6(1))
       call SchdrClm_2m2(R,M, Cr1(5), Cr2(5), Cr3(5), Cr4(5), Cr5(5), Cr6(5))
       call SchdrClm_2m1(R,M, Cr1(6), Cr2(6), Cr3(6), Cr4(6), Cr5(6), Cr6(6))
       call SchdrClm_20(R,M, Cr1(7), Cr2(7), Cr3(7), Cr4(7), Cr5(7), Cr6(7))
       call SchdrClm_21(R,M, Cr1(8), Cr2(8), Cr3(8), Cr4(8), Cr5(8), Cr6(8))
       call SchdrClm_22(R,M, Cr1(9), Cr2(9), Cr3(9), Cr4(9), Cr5(9), Cr6(9))
!       call vibSchdrAlm_00(R,M,T,F,A,Cr10(1))
!       call vibSchdrBlm_1m1(R,M,T,F,A,Cr7(2),Cr8(2),Cr9(2))
!       call vibSchdrBlm_10(R,M,T,F,A,Cr7(3),Cr8(3),Cr9(3))
!       call vibSchdrBlm_11(R,M,T,F,A,Cr7(4),Cr8(4),Cr9(4))
!       call vibSchdrClm_00(R,M,T,F,A,Cr1(1),Cr2(1),Cr3(1),Cr4(1),Cr5(1),Cr6(1))
!       call vibSchdrClm_2m2(R,M,T,F,A,Cr1(5),Cr2(5),Cr3(5),Cr4(5),Cr5(5),Cr6(5))
!       call vibSchdrClm_2m1(R,M,T,F,A,Cr1(6),Cr2(6),Cr3(6),Cr4(6),Cr5(6),Cr6(6))
!       call vibSchdrClm_20(R,M,T,F,A,Cr1(7),Cr2(7),Cr3(7),Cr4(7),Cr5(7),Cr6(7))
!       call vibSchdrClm_21(R,M,T,F,A,Cr1(8),Cr2(8),Cr3(8),Cr4(8),Cr5(8),Cr6(8))
!       call vibSchdrClm_22(R,M,T,F,A,Cr1(9),Cr2(9),Cr3(9),Cr4(9),Cr5(9),Cr6(9))

       Ct10(1) = 0.d0
       Ct7(2) = 0.d0; Ct8(2) =0.d0; Ct9(2) = 0.d0
       Ct7(3) = 0.d0; Ct8(3) =0.d0; Ct9(3) = 0.d0
       Ct7(4) = 0.d0; Ct8(4) =0.d0; Ct9(4) = 0.d0
       Ct1(1) = 0.d0; Ct2(1) = 0.d0; Ct3(1) = 0.d0; Ct4(1) = 0.d0; Ct5(1) = 0.d0; Ct6(1) = 0.d0
       Ct1(5) = 0.d0; Ct2(5) = 0.d0; Ct3(5) = 0.d0; Ct4(5) = 0.d0; Ct5(5) = 0.d0; Ct6(5) = 0.d0
       Ct1(6) = 0.d0; Ct2(6) = 0.d0; Ct3(6) = 0.d0; Ct4(6) = 0.d0; Ct5(6) = 0.d0; Ct6(6) = 0.d0
       Ct1(7) = 0.d0; Ct2(7) = 0.d0; Ct3(7) = 0.d0; Ct4(7) = 0.d0; Ct5(7) = 0.d0; Ct6(7) = 0.d0
       Ct1(8) = 0.d0; Ct2(8) = 0.d0; Ct3(8) = 0.d0; Ct4(8) = 0.d0; Ct5(8) = 0.d0; Ct6(8) = 0.d0
       Ct1(9) = 0.d0; Ct2(9) = 0.d0; Ct3(9) = 0.d0; Ct4(9) = 0.d0; Ct5(9) = 0.d0; Ct6(9) = 0.d0
!       call vibSchdtAlm_00(R,M,T,F,A, Ct10(1))
!       call vibSchdtBlm_1m1(R,M,T,F,A,Ct7(2),Ct8(2),Ct9(2))
!       call vibSchdtBlm_10(R,M,T,F,A,Ct7(3),Ct8(3),Ct9(3))
!       call vibSchdtBlm_11(R,M,T,F,A,Ct7(4),Ct8(4),Ct9(4))
!       call vibSchdtClm_00(R,M,T,F,A,Ct1(1),Ct2(1),Ct3(1),Ct4(1),Ct5(1),Ct6(1))
!       call vibSchdtClm_2m2(R,M,T,F,A,Ct1(5),Ct2(5),Ct3(5),Ct4(5),Ct5(5),Ct6(5))
!       call vibSchdtClm_2m1(R,M,T,F,A,Ct1(6),Ct2(6),Ct3(6),Ct4(6),Ct5(6),Ct6(6))
!       call vibSchdtClm_20(R,M,T,F,A,Ct1(7),Ct2(7),Ct3(7),Ct4(7),Ct5(7),Ct6(7))
!       call vibSchdtClm_21(R,M,T,F,A,Ct1(8),Ct2(8),Ct3(8),Ct4(8),Ct5(8),Ct6(8))
!       call vibSchdtClm_22(R,M,T,F,A,Ct1(9),Ct2(9),Ct3(9),Ct4(9),Ct5(9),Ct6(9))

    else
      call CCTK_WARN(0, "Coefficients for this metric not implemented")
    end if


    do l = 0, l_max
      do m = -l, l

         n = l*l + l + m + 1

         RTC(1,n) = dble(C1(n)); RTC(2,n) = dble(C2(n))
         RTC(3,n) = dble(C3(n)); RTC(4,n) = dble(C4(n))
         RTC(5,n) = dble(C5(n)); RTC(6,n) = dble(C6(n))
         RTC(7,n) = dble(C7(n)); RTC(8,n) = dble(C8(n))
         RTC(9,n) = dble(C9(n)); RTC(10,n) = dble(C10(n))

         ITC(1,n) = dimag(C1(n)); ITC(2,n) = dimag(C2(n))
         ITC(3,n) = dimag(C3(n)); ITC(4,n) = dimag(C4(n))
         ITC(5,n) = dimag(C5(n)); ITC(6,n) = dimag(C6(n))
         ITC(7,n) = dimag(C7(n)); ITC(8,n) = dimag(C8(n))
         ITC(9,n) = dimag(C9(n)); ITC(10,n) = dimag(C10(n))

         RTCt(1,n) = dble(Ct1(n)); RTCt(2,n) = dble(Ct2(n))
         RTCt(3,n) = dble(Ct3(n)); RTCt(4,n) = dble(Ct4(n))
         RTCt(5,n) = dble(Ct5(n)); RTCt(6,n) = dble(Ct6(n))
         RTCt(7,n) = dble(Ct7(n)); RTCt(8,n) = dble(Ct8(n))
         RTCt(9,n) = dble(Ct9(n)); RTCt(10,n) = dble(Ct10(n))

         ITCt(1,n) = dimag(Ct1(n)); ITCt(2,n) = dimag(Ct2(n))
         ITCt(3,n) = dimag(Ct3(n)); ITCt(4,n) = dimag(Ct4(n))
         ITCt(5,n) = dimag(Ct5(n)); ITCt(6,n) = dimag(Ct6(n))
         ITCt(7,n) = dimag(Ct7(n)); ITCt(8,n) = dimag(Ct8(n))
         ITCt(9,n) = dimag(Ct9(n)); ITCt(10,n) = dimag(Ct10(n))

         RTCr(1,n) = dble(Cr1(n)); RTCr(2,n) = dble(Cr2(n))
         RTCr(3,n) = dble(Cr3(n)); RTCr(4,n) = dble(Cr4(n))
         RTCr(5,n) = dble(Cr5(n)); RTCr(6,n) = dble(Cr6(n))
         RTCr(7,n) = dble(Cr7(n)); RTCr(8,n) = dble(Cr8(n))
         RTCr(9,n) = dble(Cr9(n)); RTCr(10,n) = dble(Cr10(n))

         ITCr(1,n) = dimag(Cr1(n)); ITCr(2,n) = dimag(Cr2(n))
         ITCr(3,n) = dimag(Cr3(n)); ITCr(4,n) = dimag(Cr4(n))
         ITCr(5,n) = dimag(Cr5(n)); ITCr(6,n) = dimag(Cr6(n))
         ITCr(7,n) = dimag(Cr7(n)); ITCr(8,n) = dimag(Cr8(n))
         ITCr(9,n) = dimag(Cr9(n)); ITCr(10,n) = dimag(Cr10(n))

      end do
    end do


 end subroutine NullSHRE_SchReadClm

