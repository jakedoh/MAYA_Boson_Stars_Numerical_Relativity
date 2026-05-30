c     Written 2002-05-22 by Erik Schnetter <schnetter@uni-tuebingen.de>
c     $Header$

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
      
      subroutine IDLinearWaves_StandWaves (CCTK_ARGUMENTS)
      
      implicit none
      
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      
      CCTK_REAL pi
      parameter (pi = 3.14159265358979d0)
      CCTK_REAL aa, kk
      CCTK_REAL zz, tt
      CCTK_REAL bb, bbdot
      integer i,j,k
      
      aa = amplitude
      kk = 2*pi / wavelength
      
      tt = cctk_time
      
      do k=1,cctk_lsh(3)
         do j=1,cctk_lsh(2)
            do i=1,cctk_lsh(1)
               
               zz = z(i,j,k)
               
               bb = aa * sin(kk*zz) * cos(kk*tt)
               bbdot = -kk*aa * sin(kk*zz) * sin(kk*tt)
               
               gxx(i,j,k) = 1 + bb
               gxy(i,j,k) = 0
               gxz(i,j,k) = 0
               gyy(i,j,k) = 1 - bb
               gyz(i,j,k) = 0
               gzz(i,j,k) = 1
               
               kxx(i,j,k) = bbdot / (-2)
               kxy(i,j,k) = 0
               kxz(i,j,k) = 0
               kyy(i,j,k) = - bbdot / (-2)
               kyz(i,j,k) = 0
               kzz(i,j,k) = 0
               
            end do 
         end do
      end do
      
c     initialise the conformal factor
      
      if (CCTK_EQUALS(metric_type, "static conformal")) then 
         
         conformal_state = 1
         
         if (CCTK_EQUALS(conformal_storage,"factor+derivs")) then
            
            conformal_state = 2
            
         else if (CCTK_EQUALS(conformal_storage,"factor+derivs+2nd derivs")) then
            conformal_state = 3
            
         end if
         
         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)		 
                  
                  psi(i,j,k) = 1
                  
                  if (conformal_state .gt. 1) then
                     psix(i,j,k) = 0
                     psiy(i,j,k) = 0
                     psiz(i,j,k) = 0
                  end if
                  
                  if (conformal_state .gt. 2) then
                     psixy(i,j,k) = 0
                     psixz(i,j,k) = 0
                     psiyz(i,j,k) = 0
                     psixx(i,j,k) = 0
                     psiyy(i,j,k) = 0
                     psizz(i,j,k) = 0
                  end if
                  
               end do
            end do
         end do
         
      end if
      
      end
