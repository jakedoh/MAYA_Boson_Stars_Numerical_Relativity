 /*@@
   @file      Whisky_TVDReconstruct.F90
   @date      Sat Jan 26 02:11:44 2002
   @author    Luca Baiotti
   @desc 
   The TVD reconstruction routine.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "SpaceMask.h"

/* #define WHISKY_UNSTABLE_TVD */

 /*@@
   @routine    tvdreconstruct
   @date       Sat Jan 26 02:12:12 2002
   @author     Luca Baiotti
   @desc 
   Performs slope limited TVD reconstruction on the given input GF
   @enddesc 
   @calls     
   @calledby   
   @history 
   Follows (in philosophy) old code by Ian Hawke
   @endhistory 

@@*/

subroutine tvdreconstruct(nx, ny, nz, xoffset, yoffset, zoffset, &
     orig, bextp, bextm, trivial_rp, space_mask, excision_descriptors)
  
  USE TVD_Scalars
  
  implicit none
  
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: i, j, k, xoffset, yoffset, zoffset, nx, ny, nz
  CCTK_REAL, dimension(nx, ny, nz) :: orig, bextp, bextm
  CCTK_REAL :: dupw, dloc, delta, ratio, hdelta
  logical, dimension(nx,ny,nz) :: trivial_rp
  CCTK_INT, dimension(nx,ny,nz) :: space_mask
  CCTK_INT :: excision_bits, excision_mask
  CCTK_INT, dimension(3) :: excision_descriptors


  bextp = 0.d0
  bextm = 0.d0

!!$ Initially all Riemann problems are NON-trivial

  trivial_rp = .false.
  do k = whisky_stencil, nz-whisky_stencil+1
    do j = whisky_stencil, ny-whisky_stencil+1
      do i = whisky_stencil, nx-whisky_stencil+1
          dupw = orig(i, j, k) - orig(i-xoffset, j-yoffset, k-zoffset)
          dloc = orig(i+xoffset, j+yoffset, k+zoffset) - orig(i, j, k)
!!$        For minmod set everything here.
!!$        Otherwise call the limiter function.
        
          if (MINMOD) then
            if (dupw*dloc < 0.d0) then
              delta=0.d0
            else if (abs(dupw)<abs(dloc)) then
              delta=dupw
            else
              delta=dloc
            end if
!!$        This is an alternative equivalent implementation 
!!$          of vanLeeer MC slopelimiter
          else if (MC2) then
            if (dupw*dloc < 0.d0) then
              delta=0.d0
            else 
              delta=sign(min(2.d0*abs(dupw),2.d0*abs(dloc),&
                   0.5d0*(abs(dupw)+abs(dloc))),dupw+dloc)
            end if
          else
            delta = 0.5d0*(dupw + dloc)
            if (abs(dupw) < myfloor ) dupw = myfloor*sign(1.d0, dupw) 
            if (abs(dloc) < myfloor ) dloc = myfloor*sign(1.d0, dloc) 
            ratio = dupw / dloc
            call slopelimiter(ratio, delta)
          end if
          hdelta = 0.5d0 * delta 
          bextm(i, j, k) = orig(i, j, k) - hdelta
          bextp(i, j, k) = orig(i, j, k) + hdelta
      end do
    end do
  end do

end subroutine tvdreconstruct

