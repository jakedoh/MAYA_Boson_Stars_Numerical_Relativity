 /*@@
   @file      Evolve.F77
   @date      June 13 2000
   @author    Gabrielle Allen
   @desc 
              Evolution routines for the 1d wave equation solver
   @enddesc 
 @@*/
  
#include "cctk.h" 
#include "cctk_Arguments.h"


 /*@@
   @routine    WaveToy1DF77_Evolve
   @date       June 13 2000
   @author     Gabrielle Allen
   @desc 
               Evolution for the 1d wave equation
   @enddesc 
   @calls      
   @calledby   
   @history 
 
   @endhistory 

@@*/

      subroutine WaveToy1DF77_Evolve(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS

      INTEGER   i
      integer   istart, iend
      CCTK_REAL dt,dx
      CCTK_REAL fac1,fac2

      dx     = CCTK_DELTA_SPACE(1)
      dt     = CCTK_DELTA_TIME

      istart = 2
      iend   = cctk_lsh(1)-1

      fac2   = (dt*dt)/(dx*dx)
      fac1   = 2.0d0*(1.0d0 - fac2)

      do i = istart, iend
         phi(i) = fac1*phi_p(i)-phi_p_p(i)+fac2*(phi_p(i+1)+phi_p(i-1))
      end do

      return
      end 








