 /*@@
   @file      WaveToy.F77
   @date      
   @author    Tom Goodale
   @desc
              Evolution routines for the wave equation solver
   @enddesc 
   @version $Header$
 @@*/
  
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"


 /*@@
   @routine    WaveToyF77_Evolution
   @date       
   @author     Tom Goodale
   @desc 
               Evolution for the wave equation
   @enddesc 
   @calls      CCTK_SyncGroup, wavetoy_boundaries
   @calledby   
   @history 
 
   @endhistory 

@@*/

      subroutine WaveToyF77_Evolution(CCTK_ARGUMENTS)

      implicit none

c     Declare variables in argument list
      DECLARE_CCTK_ARGUMENTS

      CCTK_INT, parameter :: izero = 0
      integer, parameter :: ik = kind(izero)

      INTEGER   i,j,k
      INTEGER   istart, jstart, kstart, iend, jend, kend
      CCTK_REAL dx,dy,dz,dt
      CCTK_REAL dx2,dy2,dz2,dt2
      CCTK_REAL dx2i,dy2i,dz2i

      CCTK_REAL factor

c     Set up shorthands
c     -----------------
      dx = CCTK_DELTA_SPACE(1)
      dy = CCTK_DELTA_SPACE(2)
      dz = CCTK_DELTA_SPACE(3)
      dt = CCTK_DELTA_TIME

      dx2 = dx*dx
      dy2 = dy*dy
      dz2 = dz*dz
      dt2 = dt*dt

      dx2i = 1.0/dx2
      dy2i = 1.0/dy2
      dz2i = 1.0/dz2

      istart = 2
      jstart = 2
      kstart = 2
      
      iend = cctk_lsh(1)-1
      jend = cctk_lsh(2)-1
      kend = cctk_lsh(3)-1

      factor = 2*(1 - (dt2)*(dx2i + dy2i + dz2i))

c     Do the evolution
c     ----------------
      do k = kstart, kend
         do j = jstart, jend
            do i = istart, iend
               
               phi(i,j,k) = factor*phi_p(i,j,k) -
     &              phi_p_p(i,j,k) + (dt2) * 
     &              ((phi_p(i+1,j,k)+phi_p(i-1,j,k))*dx2i
     &              +(phi_p(i,j+1,k)+phi_p(i,j-1,k))*dy2i
     &              +(phi_p(i,j,k+1)+phi_p(i,j,k-1))*dz2i)

            end do
         end do
      end do
         
      return
      end 


 /*@@
   @routine    WaveToyF77_Boundaries
   @date       
   @author     Tom Goodale
   @desc 
               Boundary conditions for the wave equation
   @enddesc 
   @history 
 
   @endhistory 

@@*/

      subroutine WaveToyF77_Boundaries(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

      CCTK_INT, parameter :: izero = 0
      integer, parameter :: ik = kind(izero)

c     Local declarations
      CCTK_INT ierr

      CHARACTER*100 boundary
      INTEGER length

      ierr = 0

c     The "bound" parameter needs to be converted into a Fortran string.
      call CCTK_FortranString(length,bound,boundary)

c     Apply the outer boundary conditions
c     -----------------------------------
c     Note: In each of the following calls to Boundary_SelectVarForBC,
c     default arguments are used, so an invalid table handle of -1 can
c     be passed
      if (CCTK_EQUALS(bound,"flat") .or. CCTK_EQUALS(bound,"static") .or.
     $    CCTK_EQUALS(bound,"radiation") .or. CCTK_EQUALS(bound,"robin") .or.
     $    CCTK_EQUALS(bound,"none") ) then
         ierr = Boundary_SelectVarForBC(cctkGH, int(CCTK_ALL_FACES,ik), 1_ik, -1_ik,
     $        "wavetoy::phi", boundary)
      else if (CCTK_EQUALS(bound,"zero")) then
         ierr = Boundary_SelectVarForBC(cctkGH, int(CCTK_ALL_FACES,ik), 1_ik, -1_ik,
     $        "wavetoy::phi", "Scalar")
      end if

      if (ierr .lt. 0)  then
        call CCTK_WARN(0,"WaveToyF77_Boundaries: Error selecting boundary condition")
      end if

      return
      end 
