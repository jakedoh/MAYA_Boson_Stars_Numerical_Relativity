 /*@@
   @file      Initial.F
   @date      May 2002
   @author    Gabrielle Allen
   @desc 
      Constraint violating initial data
   @enddesc 
   @version $Header$
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    IDConstraintViolate_Initial
   @date       May 2002
   @author     Gabrielle Allen
   @desc 
      Set up initial data which violates the constraints
   @enddesc 
   @calls     
@@*/

      subroutine IDConstraintViolate_Initial(CCTK_ARGUMENTS)

      implicit none

      integer i,j,k
      CCTK_REAL xp,yp,zp,rp

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS

      do k=1, cctk_lsh(3)
         do j=1, cctk_lsh(2)
            do i=1, cctk_lsh(1)

               xp=x(i,j,k)
               yp=y(i,j,k)
               zp=z(i,j,k)

               rp = sqrt(xp*xp+yp*yp+zp*zp)
               
               gxx(i,j,k) = 1.0d0+amplitude*exp(-(rp-radius)**2/sigma**2)
               gyy(i,j,k) = 1.0d0+amplitude*exp(-(rp-radius)**2/sigma**2)
               gzz(i,j,k) = 1.0d0+amplitude*exp(-(rp-radius)**2/sigma**2)

               gxy(i,j,k) = 0.0d0
               gxz(i,j,k) = 0.0d0
               gyz(i,j,k) = 0.0d0
                  
               kxx(i,j,k) = 0.0d0
               kxy(i,j,k) = 0.0d0
               kxz(i,j,k) = 0.0d0
               kyy(i,j,k) = 0.0d0
               kyz(i,j,k) = 0.0d0
               kzz(i,j,k) = 0.0d0

            end do
         end do
      end do
      
      end

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
