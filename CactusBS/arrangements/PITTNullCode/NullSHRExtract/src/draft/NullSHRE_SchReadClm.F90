! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"


subroutine NullSHRE_SchReadClm(CCTK_ARGUMENTS)

  use cctk
  use NullGrid_Vars
  use NullSHRE_modSchClm
  implicit none

    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    DECLARE_CCTK_FUNCTIONS
    CCTK_INT :: ll, n, l, m, k

    CCTK_COMPLEX, dimension((lmax+1)*(lmax+1)) :: &
                                   Clm_1, Clm_2, Clm_3, Clm_4, Clm_5, &
                                   Clm_6, Clm_7, Clm_8, Clm_9, Clm_10
    CCTK_COMPLEX, dimension((lmax+1)*(lmax+1)) :: &
                                   Ctlm_1, Ctlm_2, Ctlm_3, Ctlm_4, Ctlm_5, &
                                   Ctlm_6, Ctlm_7, Ctlm_8, Ctlm_9, Ctlm_10
    CCTK_COMPLEX, dimension((lmax+1)*(lmax+1)) :: &
                                   Crlm_1, Crlm_2, Crlm_3, Crlm_4, Crlm_5, &
                                   Crlm_6, Crlm_7, Crlm_8, Crlm_9, Crlm_10

#define R cr
#define M mass 

    Clm_1 = 0.d0; Clm_2 = 0.d0; Clm_3 = 0.d0; Clm_4 = 0.d0; Clm_5 = 0.d0
    Clm_6 = 0.d0; Clm_7 = 0.d0; Clm_8 = 0.d0; Clm_9 = 0.d0; Clm_10 = 0.d0
    Ctlm_1 = 0.d0; Ctlm_2 = 0.d0; Ctlm_3 = 0.d0; Ctlm_4 = 0.d0; Ctlm_5 = 0.d0
    Ctlm_6 = 0.d0; Ctlm_7 = 0.d0; Ctlm_8 = 0.d0; Ctlm_9 = 0.d0; Ctlm_10 = 0.d0
    Crlm_1 = 0.d0; Crlm_2 = 0.d0; Crlm_3 = 0.d0; Crlm_4 = 0.d0; Crlm_5 = 0.d0
    Crlm_6 = 0.d0; Crlm_7 = 0.d0; Crlm_8 = 0.d0; Crlm_9 = 0.d0; Crlm_10 = 0.d0
 
    reSCt(1,:) = real(Ctlm_1(:)) 
    reSCt(2,:) = real(Ctlm_2(:))
    reSCt(3,:) = real(Ctlm_3(:)) 
    reSCt(4,:) = real(Ctlm_4(:))
    reSCt(5,:) = real(Ctlm_5(:)) 
    reSCt(6,:) = real(Ctlm_6(:))
    reSCt(7,:) = real(Ctlm_7(:)) 
    reSCt(8,:) = real(Ctlm_8(:))
    reSCt(9,:) = real(Ctlm_9(:)) 
    reSCt(10,:) = real(Ctlm_10(:))
 
    imSCt(1,:) = imag(Ctlm_1(:)) 
    imSCt(2,:) = imag(Ctlm_2(:))
    imSCt(3,:) = imag(Ctlm_3(:)) 
    imSCt(4,:) = imag(Ctlm_4(:))
    imSCt(5,:) = imag(Ctlm_5(:)) 
    imSCt(6,:) = imag(Ctlm_6(:))
    imSCt(7,:) = imag(Ctlm_7(:)) 
    imSCt(8,:) = imag(Ctlm_8(:))
    imSCt(9,:) = imag(Ctlm_9(:)) 
    imSCt(10,:) = imag(Ctlm_10(:))


    call SchClm_00 (R,M, Clm_1(1),Clm_2(1),Clm_3(1),Clm_4(1),Clm_5(1),Clm_6(1)) 
    call SchClm_1m1(R,M, Clm_1(2),Clm_2(2),Clm_3(2),Clm_4(2),Clm_5(2),Clm_6(2)) 
    call SchClm_10 (R,M, Clm_1(3),Clm_2(3),Clm_3(3),Clm_4(3),Clm_5(3),Clm_6(3)) 
    call SchClm_11 (R,M, Clm_1(4),Clm_2(4),Clm_3(4),Clm_4(4),Clm_5(4),Clm_6(4)) 
    call SchClm_2m2(R,M, Clm_1(5),Clm_2(5),Clm_3(5),Clm_4(5),Clm_5(5),Clm_6(5)) 
    call SchClm_2m1(R,M, Clm_1(6),Clm_2(6),Clm_3(6),Clm_4(6),Clm_5(6),Clm_6(6)) 
    call SchClm_20 (R,M, Clm_1(7),Clm_2(7),Clm_3(7),Clm_4(7),Clm_5(7),Clm_6(7)) 
    call SchClm_21 (R,M, Clm_1(8),Clm_2(8),Clm_3(8),Clm_4(8),Clm_5(8),Clm_6(8)) 
    call SchClm_22 (R,M, Clm_1(9),Clm_2(9),Clm_3(9),Clm_4(9),Clm_5(9),Clm_6(9)) 

    call SchBlm_00 (R,M, Clm_7(1),Clm_8(1),Clm_9(1))
    call SchBlm_1m1(R,M, Clm_7(2),Clm_8(2),Clm_9(2))
    call SchBlm_10 (R,M, Clm_7(3),Clm_8(3),Clm_9(3))
    call SchBlm_11 (R,M, Clm_7(4),Clm_8(4),Clm_9(4))
    call SchBlm_2m2(R,M, Clm_7(5),Clm_8(5),Clm_9(5))
    call SchBlm_2m1(R,M, Clm_7(6),Clm_8(6),Clm_9(6))
    call SchBlm_20 (R,M, Clm_7(7),Clm_8(7),Clm_9(7))
    call SchBlm_21 (R,M, Clm_7(8),Clm_8(8),Clm_9(8))
    call SchBlm_22 (R,M, Clm_7(9),Clm_8(9),Clm_9(9))

    call SchAlm_00 (R,M, Clm_10(1))
    call SchAlm_1m1(R,M, Clm_10(2))
    call SchAlm_10 (R,M, Clm_10(3))
    call SchAlm_11 (R,M, Clm_10(4))
    call SchAlm_2m2(R,M, Clm_10(5))
    call SchAlm_2m1(R,M, Clm_10(6))
    call SchAlm_20 (R,M, Clm_10(7))
    call SchAlm_21 (R,M, Clm_10(8))
    call SchAlm_22 (R,M, Clm_10(9))

    call SchdrClm_00 (R,M, Crlm_1(1),Crlm_2(1),Crlm_3(1),Crlm_4(1),Crlm_5(1),Crlm_6(1)) 
    call SchdrClm_1m1(R,M, Crlm_1(2),Crlm_2(2),Crlm_3(2),Crlm_4(2),Crlm_5(2),Crlm_6(2)) 
    call SchdrClm_10 (R,M, Crlm_1(3),Crlm_2(3),Crlm_3(3),Crlm_4(3),Crlm_5(3),Crlm_6(3)) 
    call SchdrClm_11 (R,M, Crlm_1(4),Crlm_2(4),Crlm_3(4),Crlm_4(4),Crlm_5(4),Crlm_6(4)) 
    call SchdrClm_2m2(R,M, Crlm_1(5),Crlm_2(5),Crlm_3(5),Crlm_4(5),Crlm_5(5),Crlm_6(5)) 
    call SchdrClm_2m1(R,M, Crlm_1(6),Crlm_2(6),Crlm_3(6),Crlm_4(6),Crlm_5(6),Crlm_6(6)) 
    call SchdrClm_20 (R,M, Crlm_1(7),Crlm_2(7),Crlm_3(7),Crlm_4(7),Crlm_5(7),Crlm_6(7)) 
    call SchdrClm_21 (R,M, Crlm_1(8),Crlm_2(8),Crlm_3(8),Crlm_4(8),Crlm_5(8),Crlm_6(8)) 
    call SchdrClm_22 (R,M, Crlm_1(9),Crlm_2(9),Crlm_3(9),Crlm_4(9),Crlm_5(9),Crlm_6(9)) 

    call SchdrBlm_00 (R,M, Crlm_7(1),Crlm_8(1),Crlm_9(1))
    call SchdrBlm_1m1(R,M, Crlm_7(2),Crlm_8(2),Crlm_9(2))
    call SchdrBlm_10 (R,M, Crlm_7(3),Crlm_8(3),Crlm_9(3))
    call SchdrBlm_11 (R,M, Crlm_7(4),Crlm_8(4),Crlm_9(4))
    call SchdrBlm_2m2(R,M, Crlm_7(5),Crlm_8(5),Crlm_9(5))
    call SchdrBlm_2m1(R,M, Crlm_7(6),Crlm_8(6),Crlm_9(6))
    call SchdrBlm_20 (R,M, Crlm_7(7),Crlm_8(7),Crlm_9(7))
    call SchdrBlm_21 (R,M, Crlm_7(8),Crlm_8(8),Crlm_9(8))
    call SchdrBlm_22 (R,M, Crlm_7(9),Crlm_8(9),Crlm_9(9))

    call SchdrAlm_00 (R,M, Crlm_10(1))
    call SchdrAlm_1m1(R,M, Crlm_10(2))
    call SchdrAlm_10 (R,M, Crlm_10(3))
    call SchdrAlm_11 (R,M, Crlm_10(4))
    call SchdrAlm_2m2(R,M, Crlm_10(5))
    call SchdrAlm_2m1(R,M, Crlm_10(6))
    call SchdrAlm_20 (R,M, Crlm_10(7))
    call SchdrAlm_21 (R,M, Crlm_10(8))
    call SchdrAlm_22 (R,M, Crlm_10(9))


    do l = 0, lmax
      do m = -l, l

         n = l*l + l + m + 1

           reSC(1,n) = real(Clm_1(n)) 
           reSC(2,n) = real(Clm_2(n))
           reSC(3,n) = real(Clm_3(n)) 
           reSC(4,n) = real(Clm_4(n))
           reSC(5,n) = real(Clm_5(n)) 
           reSC(6,n) = real(Clm_6(n))
           reSC(7,n) = real(Clm_7(n)) 
           reSC(8,n) = real(Clm_8(n))
           reSC(9,n) = real(Clm_9(n)) 
           reSC(10,n) = real(Clm_10(n))

           imSC(1,n) = imag(Clm_1(n)) 
           imSC(2,n) = imag(Clm_2(n))
           imSC(3,n) = imag(Clm_3(n)) 
           imSC(4,n) = imag(Clm_4(n))
           imSC(5,n) = imag(Clm_5(n)) 
           imSC(6,n) = imag(Clm_6(n))
           imSC(7,n) = imag(Clm_7(n)) 
           imSC(8,n) = imag(Clm_8(n))
           imSC(9,n) = imag(Clm_9(n)) 
           imSC(10,n) = imag(Clm_10(n))

           reSCr(1,n) = real(Crlm_1(n)) 
           reSCr(2,n) = real(Crlm_2(n))
           reSCr(3,n) = real(Crlm_3(n)) 
           reSCr(4,n) = real(Crlm_4(n))
           reSCr(5,n) = real(Crlm_5(n)) 
           reSCr(6,n) = real(Crlm_6(n))
           reSCr(7,n) = real(Crlm_7(n)) 
           reSCr(8,n) = real(Crlm_8(n))
           reSCr(9,n) = real(Crlm_9(n)) 
           reSCr(10,n) = real(Crlm_10(n))

           imSCr(1,n) = imag(Crlm_1(n)) 
           imSCr(2,n) = imag(Crlm_2(n))
           imSCr(3,n) = imag(Crlm_3(n)) 
           imSCr(4,n) = imag(Crlm_4(n))
           imSCr(5,n) = imag(Crlm_5(n)) 
           imSCr(6,n) = imag(Crlm_6(n))
           imSCr(7,n) = imag(Crlm_7(n)) 
           imSCr(8,n) = imag(Crlm_8(n))
           imSCr(9,n) = imag(Crlm_9(n)) 
           imSCr(10,n) = imag(Crlm_10(n))

      end do
    end do

 end subroutine NullSHRE_SchReadClm

