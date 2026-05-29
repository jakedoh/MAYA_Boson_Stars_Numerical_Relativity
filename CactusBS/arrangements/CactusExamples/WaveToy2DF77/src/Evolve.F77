 /*@@
   @file      Evolve.F77
   @date      July 2 2000
   @author    Gabrielle Allen
   @desc 
              Evolution routines for the 2d wave equation solver
   @enddesc 
 @@*/
  
#include "cctk.h" 
#include "cctk_Arguments.h"


 /*@@
   @routine    WaveToy2DF77_Evolve
   @date       July 2 2000
   @author     Gabrielle Allen
   @desc 
               Evolution for the 2d wave equation
   @enddesc 
   @calls      
   @calledby   
   @history 
 
   @endhistory 

@@*/

      subroutine WaveToy2DF77_Evolve(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS

      INTEGER   i,j
      CCTK_REAL dt,dx,dy
      CCTK_REAL fac1,fac2,fac3

      dx     = CCTK_DELTA_SPACE(1)
      dy     = CCTK_DELTA_SPACE(2)
      dt     = CCTK_DELTA_TIME

      fac3   = (dt*dt)/(dy*dy)
      fac2   = (dt*dt)/(dx*dx)
      fac1   = 2.0d0*(1.0d0 - fac2 - fac3)

      do j = 2, cctk_lsh(2)-1
        do i = 2, cctk_lsh(1)-1
          phi(i,j) = fac1*phi_p(i,j)-phi_p_p(i,j)+
     &                 fac2*(phi_p(i+1,j)+phi_p(i-1,j))
     &                +fac3*(phi_p(i,j+1)+phi_p(i,j-1))
        end do
      end do

      return
      end 








