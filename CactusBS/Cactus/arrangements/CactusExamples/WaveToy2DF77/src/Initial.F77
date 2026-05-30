 /*@@
   @file      Initial.F77
   @date      June 13 2000
   @author    Gabrielle Allen
   @desc 
              Initial data for the 2d wave equation solver
   @enddesc 
 @@*/
  
#include "cctk.h" 
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"


 /*@@
   @routine    WaveToy2DF77_Initial
   @date       June 13 2000
   @author     Gabrielle Allen
   @desc 
               Initial data for the 2d wave equation. Assumes first
	       time derivative is zero.
   @enddesc 
   @calls      
   @calledby   
   @history 
 
   @endhistory 

@@*/

      subroutine WaveToy2DF77_Initial(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS 

      INTEGER   i,j,ierr
      CCTK_REAL dx,dy,dt	
      CCTK_REAL rad
      CCTK_REAL zero
      integer   sw(2)

      zero = 0.0
      sw(1) = 1
      sw(2) = 1

      dx = CCTK_DELTA_SPACE(1)
      dy = CCTK_DELTA_SPACE(2)
      dt = CCTK_DELTA_TIME

c     For now only have zero time derivative
      do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)
       	  timederiv(i,j) = 0.0
        end do
      end do               

      do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)
	   rad = sqrt((x(i,j)-xcenter)**2+(y(i,j)-ycenter)**2)
           phi(i,j) = amplitude*exp(-(rad)**2/sigma**2)  
        end do
      end do

c      call BndScalarVN(ierr,cctkGH, sw, zero, "wavetoy2df77::phi")

      do j = 2, cctk_lsh(2)-1
        do i = 2, cctk_lsh(1)-1
  	   phi_temp(i,j) = dt*timederiv(i,j) + 0.5*phi(i,j)+0.5*dt**2/dx**2*
     &	                (phi(i+1,j)-2.0*phi(i,j)+phi(i-1,j))
     &                  +0.5*dt**2/dy**2*(phi(i,j+1)-2.0*phi(i,j)+
     &                  phi(i,j-1))
        end do
      end do
 
      call CCTK_SyncGroup(ierr,cctkGH,"wavetoy2df77::dummy");
	
      do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)
  	   phi_p(i,j) = phi_temp(i,j)
        end do
      end do


      return
      end 


