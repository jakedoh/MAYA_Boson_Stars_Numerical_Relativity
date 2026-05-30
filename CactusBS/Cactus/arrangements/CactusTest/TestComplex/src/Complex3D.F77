c     $Header$

#include "cctk.h"
#include "cctk_Arguments.h"

      subroutine InitComplex3D(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS

      integer grouplbnd(3), grouplsh(3)
      integer i0,j0,k0
      integer nx,ny,nz
      integer i,j,k
      integer ierr

      i0 = cctk_lbnd(1)
      j0 = cctk_lbnd(2)
      k0 = cctk_lbnd(3)
      nx = cctk_lsh(1)
      ny = cctk_lsh(2)
      nz = cctk_lsh(3)

      do k=1,nz
        do j=1,ny
          do i=1, nx

            comp3d_p(i,j,k) = dcmplx(sqrt(dble(i0+i)**2+
     &                       	          dble(j0+j)**2+
     &	                     	          dble(k0+k)**2
     &                                   )
     &                              ,
     &                               dble(i0+i)**2+
     &                               dble(j0+j)**2+
     &	                     	     dble(k0+k)**2
     &                              )

	    comp3d(i,j,k) = comp3d_p(i,j,k)

          end do
        end do
      end do

      call CCTK_GrouplbndGN (ierr, cctkGH, 3, grouplbnd, "TestComplex::carrays3d")
      if (ierr .ne. 0) call CCTK_ERROR ("Group not found")
      i0 = grouplbnd(1)-1
      j0 = grouplbnd(2)-1
      k0 = grouplbnd(3)-1

      call CCTK_GrouplshGN (ierr, cctkGH, 3, grouplsh, "TestComplex::carrays3d")
      if (ierr .ne. 0) call CCTK_ERROR ("Group not found")
      nx = grouplsh(1)
      ny = grouplsh(2)
      nz = grouplsh(3)

      do k=1,nz
         do j=1,ny
            do i=1,nx
               acomp3d(i,j,k) = dcmplx(i0+i+j0+j,j0+j+k0+k)
            end do
         end do
      end do

      call CCTK_GrouplbndGN (ierr, cctkGH, 2, grouplbnd, "TestComplex::carrays2d")
      if (ierr .ne. 0) call CCTK_ERROR ("Group not found")
      i0 = grouplbnd(1)-1
      j0 = grouplbnd(2)-1

      call CCTK_GrouplshGN (ierr, cctkGH, 2, grouplsh, "TestComplex::carrays2d")
      if (ierr .ne. 0) call CCTK_ERROR ("Group not found")
      nx = grouplsh(1)
      ny = grouplsh(2)

      do j=1,ny
         do i=1,nx
            acomp2d(i,j) = dcmplx(i0+i,j0+j)
         end do
      end do

      call CCTK_GrouplbndGN (ierr, cctkGH, 1, grouplbnd, "TestComplex::carrays1d")
      if (ierr .ne. 0) call CCTK_ERROR ("Group not found")
      i0 = grouplbnd(1)-1

      call CCTK_GrouplshGN (ierr, cctkGH, 1, grouplsh, "TestComplex::carrays1d")
      if (ierr .ne. 0) call CCTK_ERROR ("Group not found")
      nx = grouplsh(1)

      do i=1,nx
         acomp1d(i) = dcmplx(1,i0+i)
      end do

      scomp = dcmplx(1,2)

      end


      subroutine EvolveComplex3D(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS

      integer nx,ny,nz
      integer i,j,k

      nx = cctk_lsh(1)
      ny = cctk_lsh(2)
      nz = cctk_lsh(3)

c     provide boundary data
      comp3d = comp3d_p

      do k=2,nz-1
        do j=2,ny-1
          do i=2,nx-1

            comp3d(i,j,k) = 1.0/6.0*(
     &          comp3d_p(i-1,j,k)+comp3d_p(i+1,j,k)+
     &          comp3d_p(i,j-1,k)+comp3d_p(i,j+1,k)+
     &          comp3d_p(i,j,k-1)+comp3d_p(i,j,k+1)
     &      )

          end do
        end do
      end do

      end
