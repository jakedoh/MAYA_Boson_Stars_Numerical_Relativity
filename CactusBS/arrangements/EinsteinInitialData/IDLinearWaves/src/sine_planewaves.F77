c     $Header$

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"


      subroutine IDLinearWaves_SinePlaneWaves (CCTK_ARGUMENTS)

      implicit none
      
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      
      CCTK_REAL pi
      CCTK_REAL the,phi
      CCTK_REAL st,ct,sp,cp
      CCTK_REAL kx,ky,kz,w
      CCTK_REAL b,bt
      
      INTEGER i,j,k
      
      CHARACTER*200 infoline
      
      pi = 3.14159265358979d0

c     determine the direction for the plane wave to travel
c     and convert it from degrees to radians
      the = pi*wavetheta/180
      phi = pi*wavephi/180

c     precalc
      st = sin(the)
      ct = cos(the)
      sp = sin(phi)
      cp = cos(phi)
      kx = 2*pi*st*cp/wavelength
      ky = 2*pi*st*sp/wavelength
      kz = 2*pi*ct/wavelength
      w = sqrt(kx*kx+ky*ky+kz*kz)

c     too lazy
      if (wavetheta .ne. 90) then
         call CCTK_WARN (0, "wavetheta.ne.90 is not implemented")
      end if

c  *************** plane waves ********************
      do k=1,cctk_lsh(3)
         do j=1,cctk_lsh(2)
            do i=1,cctk_lsh(1)           
               
               b = amplitude * sin(kx*x(i,j,k)+ky*y(i,j,k)+kz*z(i,j,k)-w*cctk_time)
               bt = - amplitude * w * cos(kx*x(i,j,k)+ky*y(i,j,k)+kz*z(i,j,k)-w*cctk_time)
               
c     the metric functions
               gxx(i,j,k) = sp*sp * (1+b) + cp*cp * (1)
               gxy(i,j,k) = cp*sp * (1-(1+b))
               gyy(i,j,k) = cp*cp * (1+b) + sp*sp * (1)
               gxz(i,j,k) = 0
               gyz(i,j,k) = 0
               gzz(i,j,k) = 1-b
        
c     and the extrinsic curvature
               
               kxx(i,j,k) = sp*sp * (bt/2) + cp*cp * (0)
               kxy(i,j,k) = cp*sp * (-bt/2)
               kyy(i,j,k) = cp*cp * (bt/2) + sp*sp * (0)
               kxz(i,j,k) = 0
               kyz(i,j,k) = 0
               kzz(i,j,k) = -bt/2
               
          enddo 
        enddo
      enddo
      
c     initialize the conformal factor
      if (CCTK_EQUALS(metric_type, "static conformal")) then 
         
         conformal_state = 1
         
         if (CCTK_EQUALS(conformal_storage, "factor+derivs")) then
            
            conformal_state = 2
            
         else if (CCTK_EQUALS(conformal_storage, "factor+derivs+2nd derivs")) then
            conformal_state = 3
            
         end if
         
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)                
                  
                  psi(i,j,k)   = 1d0
                  
                  if(conformal_state .gt. 1) then
                     psix(i,j,k)  = 0d0
                     psiy(i,j,k)  = 0d0
                     psiz(i,j,k)  = 0d0    
                  endif
                  
                  if(conformal_state .gt. 2) then
                     psixy(i,j,k) = 0d0
                     psixz(i,j,k) = 0d0
                     psiyz(i,j,k) = 0d0  
                     psixx(i,j,k) = 0d0
                     psiyy(i,j,k) = 0d0
                     psizz(i,j,k) = 0d0
                  endif
                  
               end do
            end do
         end do
         
      end if
      
      end
