 /*@@
   @file      InitialData.F77
   @date      
   @author    Tom Goodale
   @desc 
              Initial data for the 3D Wave Equation
   @enddesc 
 @@*/

#include "cctk.h" 
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"


 /*@@
   @routine    IDScalarWave_InitialData
   @date       
   @author     Tom Goodale
   @desc 
               Set up initial data for the wave equation
   @enddesc 
   @calls      
   @calledby   
   @history 
 
   @endhistory 

@@*/

      subroutine IDScalarWave_InitialData(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

      CCTK_REAL  one
      parameter  (one = 1)
      INTEGER    i,j,k
      CCTK_REAL  dt,omega, pi
      CCTK_REAL  xp,yp,zp,rp

      pi = 4*atan(one)

      dt = CCTK_DELTA_TIME

      omega = sqrt(kx**2+ky**2+kz**2)

      if (CCTK_EQUALS(initial_data,"plane")) then 

         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)

                  phi(i,j,k)     = amplitude*cos(kx*x(i,j,k)+ky*y(i,j,k)
     &                             +kz*z(i,j,k)+omega*cctk_time)
                  phi_p(i,j,k) = amplitude*cos(kx*x(i,j,k)+ky*y(i,j,k)
     &                             +kz*z(i,j,k)+omega*(cctk_time-dt))
                  
               end do
            end do
         end do

      else if (CCTK_EQUALS(initial_data,"gaussian")) then 

         do k=1, cctk_lsh(3)
            do j=1, cctk_lsh(2)
               do i=1, cctk_lsh(1)
 
	          xp=x(i,j,k)
	          yp=y(i,j,k)
                  zp=z(i,j,k)

	          rp = sqrt(xp*xp+yp*yp+zp*zp)

                  phi(i,j,k) = amplitude*exp(-(rp-radius)**2/sigma**2)

		  if (rp .eq. 0.0) then
                    phi_p(i,j,k) = amplitude*(1.0 - 2.0*dt**2/sigma**2)*exp(-dt**2/sigma**2)
	          else
                    phi_p(i,j,k) = amplitude/2.0*(rp-dt)/rp*
     &                 exp( - ( (rp-radius-dt)/sigma)**2 )
     &                 + amplitude/2.0*(rp+dt)/rp*
     &                 exp( - ( (rp-radius+dt)/sigma)**2 )
                  endif

               end do
            end do
         end do

      else if (CCTK_EQUALS(initial_data, "box")) then

c        Use kx,ky,kz as number of modes in each direction.

         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)

                  phi(i,j,k) = amplitude*sin(kx*(x(i,j,k)-0.5)*pi)*
     &                 sin(ky*(y(i,j,k)-0.5)*pi)*
     &                 sin(kz*(z(i,j,k)-0.5)*pi)*
     &                 cos(omega*cctk_time*pi)
                  
                  phi_p(i,j,k)= amplitude*sin(kx*(x(i,j,k)-0.5)*pi)*
     &                 sin(ky*(y(i,j,k)-0.5)*pi)*
     &                 sin(kz*(z(i,j,k)-0.5)*pi)*
     &                 cos(omega*(cctk_time-dt)*pi)

               end do
            end do
         end do

      else if (CCTK_EQUALS(initial_data, "none")) then
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)
                  
                  phi(i,j,k)   = 0.0d0
                  phi_p(i,j,k) = 0.0d0

               end do
            end do
         end do 
         
      end if

      return
      end 


