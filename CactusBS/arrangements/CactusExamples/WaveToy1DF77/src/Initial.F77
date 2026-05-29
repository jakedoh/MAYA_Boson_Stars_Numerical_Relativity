 /*@@
   @file      Initial.F77
   @date      June 13 2000
   @author    Gabrielle Allen
   @desc 
              Initial data for the 1d wave equation solver
   @enddesc 
 @@*/
  
#include "cctk.h" 
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"


 /*@@
   @routine    WaveToy1DF77_Initial
   @date       June 13 2000
   @author     Gabrielle Allen
   @desc 
               Initial data for the 1d wave equation. Assumes first
	       time derivative is zero.
   @enddesc 
   @calls      
   @calledby   
   @history 
 
   @endhistory 

@@*/

      subroutine WaveToy1DF77_Initial(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS 

      INTEGER   i,ierr
      CCTK_REAL dx,dt
      CCTK_REAL zero
      integer   sw(1)

      sw(1) = 1
      zero = 0.0d0

      dx = CCTK_DELTA_SPACE(1)
      dt = CCTK_DELTA_TIME

c     For now only have zero time derivative
      do i = 1, cctk_lsh(1)
	timederiv(i) = 0.0d0
      end do               
    
      do i = 1, cctk_lsh(1)
         phi(i) = amplitude*exp(-(sqrt(x(i)**2)-center)**2/sigma**2)
      end do

c     This is horrible, use a temporary grid variable
      phi_temp(1) = 0.0d0
      do i = 2, cctk_lsh(1)-1
	 phi_temp(i) = dt*timederiv(i) + phi(i) + 0.5d0*dt**2/dx**2*
     &	            (phi(i+1)-2.0d0*phi(i)+phi(i-1))
      end do
      phi_temp(cctk_lsh(1)) = 0.0d0

      call CCTK_SyncGroup(ierr,cctkGH,"wavetoy1df77::dummy")

      do i = 1, cctk_lsh(1)
         phi_p(i) = phi_temp(i)
      end do

      return
      end 


