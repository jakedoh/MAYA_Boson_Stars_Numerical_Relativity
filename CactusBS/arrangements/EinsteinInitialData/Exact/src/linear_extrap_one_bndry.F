c this subroutine linearly extrapolates one Cactus variable
c on one boundary of the Cactus grid box
C $Header$

#include "cctk.h"
#include "cctk_Arguments.h"
#include "Exact.inc"

c #define-ing the symbol EXACT_NO_F90 will turn this file into a no-op
#ifndef EXACT_NO_F90
      subroutine Exact__linear_extrap_one_bndry(CCTK_ARGUMENTS,var)

      implicit none

      DECLARE_CCTK_ARGUMENTS_Exact__linear_extrap_one_bndry

      integer nx,ny,nz

      CCTK_REAL var(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))

C     Grid parameters.

      nx = cctk_lsh(1)
      ny = cctk_lsh(2)
      nz = cctk_lsh(3)

C     Linear extrapolation from the interiors to the boundaries.
C     Does not support octant or quadrant.

C     6 faces.

      if (cctk_bbox(1).eq.1 .and. nx.ge.4 ) then 
         var(1,:,:) = 2.d0 * var(2,:,:) - var(3,:,:)
      endif
      if (cctk_bbox(3).eq.1 .and. ny.ge.4 ) then 
         var(:,1,:) = 2.d0 * var(:,2,:) - var(:,3,:)
      endif
      if (cctk_bbox(5).eq.1 .and. nz.ge.4 ) then 
         var(:,:,1) = 2.d0 * var(:,:,2) - var(:,:,3)
      endif
      if (cctk_bbox(2).eq.1 .and. nx.ge.4 ) then 
         var(nx,:,:) = 2.d0 * var(nx-1,:,:) - var(nx-2,:,:)
      endif
      if (cctk_bbox(4).eq.1 .and. ny.ge.4 ) then 
         var(:,ny,:) = 2.d0 * var(:,ny-1,:) - var(:,ny-2,:)
      endif
      if (cctk_bbox(6).eq.1 .and. nz.ge.4 ) then 
         var(:,:,nz) = 2.d0 * var(:,:,nz-1) - var(:,:,nz-2)
      endif

C     12 edges.
C     4 round face x=min.
      if ( cctk_bbox(1).eq.1 .and. nx.ge.4 
     $     .and. cctk_bbox(3).eq.1 .and. ny.ge.4 ) then 
         var(1,1,:) = 2.d0 * var(2,2,:) - var(3,3,:)
      end if

      if ( cctk_bbox(1).eq.1 .and. nx.ge.4 
     $     .and. cctk_bbox(4).eq.1 .and. ny.ge.4 ) then 
         var(1,ny,:) = 2.d0 * var(2,ny-1,:) - var(2,ny-2,:)
      end if

      if ( cctk_bbox(1).eq.1 .and. nx.ge.4 
     $     .and. cctk_bbox(5).eq.1 .and. nz.ge.4 ) then 
         var(1,:,1) = 2.d0 * var(2,:,2) - var(3,:,3)
      end if

      if ( cctk_bbox(1).eq.1 .and. nx.ge.4 
     $     .and. cctk_bbox(6).eq.1 .and. nz.ge.4 ) then 
         var(1,:,nz) = 2.d0 * var(2,:,nz-1) - var(3,:,nz-2)
      end if

C     4 around face x=max.
      if ( cctk_bbox(2).eq.1 .and. nx.ge.4 
     $     .and. cctk_bbox(3).eq.1 .and. ny.ge.4 ) then 
         var(nx,1,:) = 2.d0 * var(nx-1,2,:) - var(nx-2,3,:)
      end if

      if ( cctk_bbox(2).eq.1 .and. nx.ge.4 
     $     .and. cctk_bbox(4).eq.1 .and. ny.ge.4 ) then 
         var(nx,ny,:) = 2.d0 * var(nx-1,ny-1,:) - var(nx-2,ny-2,:)
      end if

      if ( cctk_bbox(2).eq.1 .and. nx.ge.4 
     $     .and. cctk_bbox(5).eq.1 .and. nz.ge.4 ) then 
         var(nx,:,1) = 2.d0 * var(nx-1,:,2) - var(nx-2,:,3)
      end if

      if ( cctk_bbox(2).eq.1 .and. nx.ge.4 
     $     .and. cctk_bbox(6).eq.1 .and. nz.ge.4 ) then 
         var(nx,:,nz) = 2.d0 * var(nx-1,:,nz-1) - var(nx-2,:,nz-2)
      end if

C     Remaining 2 in y=min.
      if ( cctk_bbox(3).eq.1 .and. ny.ge.4
     $     .and. cctk_bbox(5).eq.1 .and. nz.ge.4 ) then
         var(:,1,1) = 2.d0 * var(:,2,2) - var(:,3,3)
      end if

      if ( cctk_bbox(3).eq.1 .and. ny.ge.4
     $     .and. cctk_bbox(6).eq.1 .and. nz.ge.4 ) then 
         var(:,1,nz) = 2.d0 * var(:,2,nz-1) - var(:,2,nz-2)
      end if

C     Remaining 2 in y=ymax.      
      if ( cctk_bbox(4).eq.1 .and. ny.ge.4
     $     .and. cctk_bbox(5).eq.1 .and. nz.ge.4 ) then
         var(:,ny,1) = 2.d0 * var(:,ny-1,2) - var(:,ny-2,3)  
      end if

      if ( cctk_bbox(4).eq.1 .and. ny.ge.4
     $     .and. cctk_bbox(6).eq.1 .and. nz.ge.4 ) then 
         var(:,ny,nz) = 2.d0 * var(:,ny-1,nz-1) - var(:,ny-2,nz-2)
      end if

C     8 corners.

      if (nx.ge.4 .and. ny.ge.4 .and. nz.ge.4) then
         if (cctk_bbox(1).eq.1 .and. cctk_bbox(3).eq.1 .and. cctk_bbox(5).eq.1) then
            var(1,1,1) = 2.d0*var(2,2,2) - var(3,3,3)
         end if
         if (cctk_bbox(1).eq.1 .and. cctk_bbox(3).eq.1 .and. cctk_bbox(6).eq.1) then
            var(1,1,nz) = 2.d0*var(2,2,nz-1) - var(3,3,nz-2)
         end if
         if (cctk_bbox(1).eq.1 .and. cctk_bbox(4).eq.1 .and. cctk_bbox(5).eq.1) then
            var(1,ny,1) = 2.d0*var(2,ny-1,2) - var(3,ny-2,3)
         end if
         if (cctk_bbox(1).eq.1 .and. cctk_bbox(4).eq.1 .and. cctk_bbox(6).eq.1) then
            var(1,ny,nz) = 2.d0*var(2,ny-1,nz-1) - var(3,ny-2,nz-2)
         end if
        if (cctk_bbox(2).eq.1 .and. cctk_bbox(3).eq.1 .and. cctk_bbox(5).eq.1) then
            var(nx,1,1) = 2.d0*var(nx-1,2,2) - var(nx-2,3,3)
         end if
         if (cctk_bbox(2).eq.1 .and. cctk_bbox(3).eq.1 .and. cctk_bbox(6).eq.1) then
            var(nx,1,nz) = 2.d0*var(nx-1,2,nz-1) - var(nx-2,3,nz-2)
         end if
         if (cctk_bbox(2).eq.1 .and. cctk_bbox(4).eq.1 .and. cctk_bbox(5).eq.1) then
            var(nx,ny,1) = 2.d0*var(nx-1,ny-1,2) - var(nx-2,ny-2,3)
         end if
         if (cctk_bbox(2).eq.1 .and. cctk_bbox(4).eq.1 .and. cctk_bbox(6).eq.1) then
            var(nx,ny,nz) = 2.d0*var(nx-1,ny-1,nz-1) - var(nx-2,ny-2,nz-2)
         end if
      end if


      return
      end
#endif
