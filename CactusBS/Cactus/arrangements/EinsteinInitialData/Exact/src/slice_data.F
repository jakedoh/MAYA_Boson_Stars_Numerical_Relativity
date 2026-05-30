C     Extract Cauchy data from the slice x^A(x^i) stored in slicex,
C     slicey, slicez, slicet, and calculate dx^A/dt.
C $Header$

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "Exact.inc"

      subroutine Exact__slice_data(CCTK_ARGUMENTS)
      implicit none
      DECLARE_CCTK_ARGUMENTS_Exact__slice_data

c #define-ing the symbol EXACT_NO_F90 will turn this subroutine into a no-op
#ifndef EXACT_NO_F90
      integer i, j, k, l, m, n, p, q
      integer nx,ny,nz
      integer ierr

      CCTK_REAL s1d(4,3), nd(4), nu(4), norm, gd(4,4), gu(4,4), g3(3,3),
     $   gd_p(4,4), gd_m(4,4), gd1d(4,4,4), s2d(4,3,3), k3(3,3), 
     $   ex_eps, dx, dy, dz, exact_psi
      parameter (ex_eps=1.d-6)

C     Grid parameters.

      nx = cctk_lsh(1)
      ny = cctk_lsh(2)
      nz = cctk_lsh(3)

      dx = CCTK_DELTA_SPACE(1)
      dy = CCTK_DELTA_SPACE(2)
      dz = CCTK_DELTA_SPACE(3)

C     Sum over interior points on the slice.

      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1

C     Calculate first derivatives of slice coordinates.

               s1d(1,1) = 0.5d0 * (slicex(i+1,j,k) - slicex(i-1,j,k))/dx
               s1d(1,2) = 0.5d0 * (slicex(i,j+1,k) - slicex(i,j-1,k))/dy
               s1d(1,3) = 0.5d0 * (slicex(i,j,k+1) - slicex(i,j,k-1))/dz

               s1d(2,1) = 0.5d0 * (slicey(i+1,j,k) - slicey(i-1,j,k))/dx
               s1d(2,2) = 0.5d0 * (slicey(i,j+1,k) - slicey(i,j-1,k))/dy
               s1d(2,3) = 0.5d0 * (slicey(i,j,k+1) - slicey(i,j,k-1))/dz

               s1d(3,1) = 0.5d0 * (slicez(i+1,j,k) - slicez(i-1,j,k))/dx
               s1d(3,2) = 0.5d0 * (slicez(i,j+1,k) - slicez(i,j-1,k))/dy
               s1d(3,3) = 0.5d0 * (slicez(i,j,k+1) - slicez(i,j,k-1))/dz

               s1d(4,1) = 0.5d0 * (slicet(i+1,j,k) - slicet(i-1,j,k))/dx
               s1d(4,2) = 0.5d0 * (slicet(i,j+1,k) - slicet(i,j-1,k))/dy
               s1d(4,3) = 0.5d0 * (slicet(i,j,k+1) - slicet(i,j,k-1))/dz

C     Preset the conformal factor
               if (conformal_state .ne. 0) then
                 exact_psi=1.0d0
               else
                 exact_psi=0.0d0
               end if
               
C     Now we need the exact solution metric in the preferred coordinates x^A.

               call Exact__metric(
     $              decoded_exact_model,
     $              slicex(i,j,k), slicey(i,j,k), slicez(i,j,k), 
     $              slicet(i,j,k),
     $              gd(4,4), gd(1,4), gd(2,4), gd(3,4),
     $              gd(1,1), gd(2,2), gd(3,3), gd(1,2), gd(2,3), gd(1,3),
     $              gu(4,4), gu(1,4), gu(2,4), gu(3,4),
     $              gu(1,1), gu(2,2), gu(3,3), gu(1,2), gu(2,3), gu(1,3),
     $              exact_psi)

C     Calculate n^A and dx^A/dt

#include "include/slice_normal.inc"

C     Calculate g_ij = X^A,i X^B,j g_AB, at first in a temp.

               g3 = 0.d0

               do p=1,3
                  do q=p,3
                     do m=1,4
                        do n=1,4
                           g3(p,q) = g3(p,q) 
     $                          + s1d(m,p) * s1d(n,q) * gd(m,n)
                        end do
                     end do
                  end do
               end do

               gxx(i,j,k) = g3(1,1)
               gxy(i,j,k) = g3(1,2)
               gxz(i,j,k) = g3(1,3)
               gyy(i,j,k) = g3(2,2)
               gyz(i,j,k) = g3(2,3)
               gzz(i,j,k) = g3(3,3)

C     Calculate x^A,ij. Do this by hand, and loop over first index 
C     with an editor (hint for proofreading).

               s2d(1,1,1) = (slicex(i+1,j,k) + slicex(i-1,j,k) 
     $              - 2.d0 * slicex(i,j,k)) / dx**2 
               s2d(1,2,2) = (slicex(i,j+1,k) + slicex(i,j-1,k) 
     $              - 2.d0 * slicex(i,j,k)) / dy**2 
               s2d(1,3,3) = (slicex(i,j,k+1) + slicex(i,j,k-1) 
     $              - 2.d0 * slicex(i,j,k)) / dz**2 
               s2d(1,1,2) = (slicex(i+1,j+1,k) - slicex(i-1,j+1,k)
     $              - slicex(i+1,j-1,k) + slicex(i-1,j-1,k)) / (4.*dx*dy)
               s2d(1,1,3) = (slicex(i+1,j,k+1) - slicex(i-1,j,k+1)
     $              - slicex(i+1,j,k-1) + slicex(i-1,j,k-1)) / (4.*dy*dz)
               s2d(1,2,3) = (slicex(i,j+1,k+1) - slicex(i,j-1,k+1)
     $              - slicex(i,j+1,k-1) + slicex(i,j-1,k-1)) / (4.*dx*dz)
               s2d(1,2,1) = s2d(1,1,2)
               s2d(1,3,2) = s2d(1,2,3)
               s2d(1,3,1) = s2d(1,1,3)

               s2d(2,1,1) = (slicey(i+1,j,k) + slicey(i-1,j,k) 
     $              - 2.d0 * slicey(i,j,k)) / dx**2 
               s2d(2,2,2) = (slicey(i,j+1,k) + slicey(i,j-1,k) 
     $              - 2.d0 * slicey(i,j,k)) / dy**2 
               s2d(2,3,3) = (slicey(i,j,k+1) + slicey(i,j,k-1) 
     $              - 2.d0 * slicey(i,j,k)) / dz**2 
               s2d(2,1,2) = (slicey(i+1,j+1,k) - slicey(i-1,j+1,k)
     $              - slicey(i+1,j-1,k) + slicey(i-1,j-1,k)) / (4.*dx*dy)
               s2d(2,1,3) = (slicey(i+1,j,k+1) - slicey(i-1,j,k+1)
     $              - slicey(i+1,j,k-1) + slicey(i-1,j,k-1)) / (4.*dy*dz)
               s2d(2,2,3) = (slicey(i,j+1,k+1) - slicey(i,j-1,k+1)
     $              - slicey(i,j+1,k-1) + slicey(i,j-1,k-1)) / (4.*dx*dz)
               s2d(2,2,1) = s2d(2,1,2)
               s2d(2,3,2) = s2d(2,2,3)
               s2d(2,3,1) = s2d(2,1,3)

               s2d(3,1,1) = (slicez(i+1,j,k) + slicez(i-1,j,k) 
     $              - 2.d0 * slicez(i,j,k)) / dx**2 
               s2d(3,2,2) = (slicez(i,j+1,k) + slicez(i,j-1,k) 
     $              - 2.d0 * slicez(i,j,k)) / dy**2 
               s2d(3,3,3) = (slicez(i,j,k+1) + slicez(i,j,k-1) 
     $              - 2.d0 * slicez(i,j,k)) / dz**2 
               s2d(3,1,2) = (slicez(i+1,j+1,k) - slicez(i-1,j+1,k)
     $              - slicez(i+1,j-1,k) + slicez(i-1,j-1,k)) / (4.*dx*dy)
               s2d(3,1,3) = (slicez(i+1,j,k+1) - slicez(i-1,j,k+1)
     $              - slicez(i+1,j,k-1) + slicez(i-1,j,k-1)) / (4.*dy*dz)
               s2d(3,2,3) = (slicez(i,j+1,k+1) - slicez(i,j-1,k+1)
     $              - slicez(i,j+1,k-1) + slicez(i,j-1,k-1)) / (4.*dx*dz)
               s2d(3,2,1) = s2d(3,1,2)
               s2d(3,3,2) = s2d(3,2,3)
               s2d(3,3,1) = s2d(3,1,3)

               s2d(4,1,1) = (slicet(i+1,j,k) + slicet(i-1,j,k) 
     $              - 2.d0 * slicet(i,j,k)) / dx**2 
               s2d(4,2,2) = (slicet(i,j+1,k) + slicet(i,j-1,k) 
     $              - 2.d0 * slicet(i,j,k)) / dy**2 
               s2d(4,3,3) = (slicet(i,j,k+1) + slicet(i,j,k-1) 
     $              - 2.d0 * slicet(i,j,k)) / dz**2 
               s2d(4,1,2) = (slicet(i+1,j+1,k) - slicet(i-1,j+1,k)
     $              - slicet(i+1,j-1,k) + slicet(i-1,j-1,k)) / (4.*dx*dy)
               s2d(4,1,3) = (slicet(i+1,j,k+1) - slicet(i-1,j,k+1)
     $              - slicet(i+1,j,k-1) + slicet(i-1,j,k-1)) / (4.*dy*dz)
               s2d(4,2,3) = (slicet(i,j+1,k+1) - slicet(i,j-1,k+1)
     $              - slicet(i,j+1,k-1) + slicet(i,j-1,k-1)) / (4.*dx*dz)
               s2d(4,2,1) = s2d(4,1,2)
               s2d(4,3,2) = s2d(4,2,3)
               s2d(4,3,1) = s2d(4,1,3)


C     Calculate g_AB,C. Need to sum explicitly over C. Do this with 
C     the editor.

C     Preset the conformal factor
               if (conformal_state .ne. 0) then
                 exact_psi=1.0d0
               else
                 exact_psi=0.0d0
               end if
               call Exact__metric(
     $              decoded_exact_model,
     $              slicex(i,j,k)+ex_eps, slicey(i,j,k), slicez(i,j,k), 
     $              slicet(i,j,k),
     $              gd_p(4,4), gd_p(1,4), gd_p(2,4), gd_p(3,4),
     $              gd_p(1,1), gd_p(2,2), gd_p(3,3), 
     $              gd_p(1,2), gd_p(2,3), gd_p(1,3),
     $              gu(4,4), gu(1,4), gu(2,4), gu(3,4),
     $              gu(1,1), gu(2,2), gu(3,3), 
     $              gu(1,2), gu(2,3), gu(1,3), exact_psi)

C     Preset the conformal factor
               if (conformal_state .ne. 0) then
                 exact_psi=1.0d0
               else
                 exact_psi=0.0d0
               end if
               call Exact__metric(
     $              decoded_exact_model,
     $              slicex(i,j,k)-ex_eps, slicey(i,j,k), slicez(i,j,k), 
     $              slicet(i,j,k),
     $              gd_m(4,4), gd_m(1,4), gd_m(2,4), gd_m(3,4),
     $              gd_m(1,1), gd_m(2,2), gd_m(3,3), 
     $              gd_m(1,2), gd_m(2,3), gd_m(1,3),
     $              gu(4,4), gu(1,4), gu(2,4), gu(3,4),
     $              gu(1,1), gu(2,2), gu(3,3), 
     $              gu(1,2), gu(2,3), gu(1,3), exact_psi)
               do m=1,4
                  do n=m,4
                     gd1d(m,n,1) = (gd_p(m,n) - gd_m(m,n)) / (2.*ex_eps)
                  end do
               end do

C     Preset the conformal factor
               if (conformal_state .ne. 0) then
                 exact_psi=1.0d0
               else
                 exact_psi=0.0d0
               end if
               call Exact__metric(
     $              decoded_exact_model,
     $              slicex(i,j,k), slicey(i,j,k)+ex_eps, slicez(i,j,k), 
     $              slicet(i,j,k),
     $              gd_p(4,4), gd_p(1,4), gd_p(2,4), gd_p(3,4),
     $              gd_p(1,1), gd_p(2,2), gd_p(3,3), 
     $              gd_p(1,2), gd_p(2,3), gd_p(1,3),
     $              gu(4,4), gu(1,4), gu(2,4), gu(3,4),
     $              gu(1,1), gu(2,2), gu(3,3), 
     $              gu(1,2), gu(2,3), gu(1,3), exact_psi)
               
C     Preset the conformal factor
               if (conformal_state .ne. 0) then
                 exact_psi=1.0d0
               else
                 exact_psi=0.0d0
               end if
               call Exact__metric(
     $              decoded_exact_model,
     $              slicex(i,j,k), slicey(i,j,k)-ex_eps, slicez(i,j,k), 
     $              slicet(i,j,k),
     $              gd_m(4,4), gd_m(1,4), gd_m(2,4), gd_m(3,4),
     $              gd_m(1,1), gd_m(2,2), gd_m(3,3), 
     $              gd_m(1,2), gd_m(2,3), gd_m(1,3),
     $              gu(4,4), gu(1,4), gu(2,4), gu(3,4),
     $              gu(1,1), gu(2,2), gu(3,3), 
     $              gu(1,2), gu(2,3), gu(1,3), exact_psi)
               do m=1,4
                  do n=m,4
                     gd1d(m,n,2) = (gd_p(m,n) - gd_m(m,n)) / (2.*ex_eps)
                  end do
               end do

               
C     Preset the conformal factor
               if (conformal_state .ne. 0) then
                 exact_psi=1.0d0
               else
                 exact_psi=0.0d0
               end if
               call Exact__metric(
     $              decoded_exact_model,
     $              slicex(i,j,k), slicey(i,j,k), slicez(i,j,k)+ex_eps, 
     $              slicet(i,j,k),
     $              gd_p(4,4), gd_p(1,4), gd_p(2,4), gd_p(3,4),
     $              gd_p(1,1), gd_p(2,2), gd_p(3,3), 
     $              gd_p(1,2), gd_p(2,3), gd_p(1,3),
     $              gu(4,4), gu(1,4), gu(2,4), gu(3,4),
     $              gu(1,1), gu(2,2), gu(3,3), 
     $              gu(1,2), gu(2,3), gu(1,3), exact_psi)
               
C     Preset the conformal factor
               if (conformal_state .ne. 0) then
                 exact_psi=1.0d0
               else
                 exact_psi=0.0d0
               end if
               call Exact__metric(
     $              decoded_exact_model,
     $              slicex(i,j,k), slicey(i,j,k), slicez(i,j,k)-ex_eps, 
     $              slicet(i,j,k),
     $              gd_m(4,4), gd_m(1,4), gd_m(2,4), gd_m(3,4),
     $              gd_m(1,1), gd_m(2,2), gd_m(3,3), 
     $              gd_m(1,2), gd_m(2,3), gd_m(1,3),
     $              gu(4,4), gu(1,4), gu(2,4), gu(3,4),
     $              gu(1,1), gu(2,2), gu(3,3), 
     $              gu(1,2), gu(2,3), gu(1,3), exact_psi)
               do m=1,4
                  do n=m,4
                     gd1d(m,n,3) = (gd_p(m,n) - gd_m(m,n)) / (2.*ex_eps)
                  end do
               end do

               
C     Preset the conformal factor
               if (conformal_state .ne. 0) then
                 exact_psi=1.0d0
               else
                 exact_psi=0.0d0
               end if
               call Exact__metric(
     $              decoded_exact_model,
     $              slicex(i,j,k), slicey(i,j,k), slicez(i,j,k), 
     $              slicet(i,j,k)+ex_eps,
     $              gd_p(4,4), gd_p(1,4), gd_p(2,4), gd_p(3,4),
     $              gd_p(1,1), gd_p(2,2), gd_p(3,3), 
     $              gd_p(1,2), gd_p(2,3), gd_p(1,3),
     $              gu(4,4), gu(1,4), gu(2,4), gu(3,4),
     $              gu(1,1), gu(2,2), gu(3,3), 
     $              gu(1,2), gu(2,3), gu(1,3), exact_psi)
               
C     Preset the conformal factor
               if (conformal_state .ne. 0) then
                 exact_psi=1.0d0
               else
                 exact_psi=0.0d0
               end if
               call Exact__metric(
     $              decoded_exact_model,
     $              slicex(i,j,k), slicey(i,j,k), slicez(i,j,k), 
     $              slicet(i,j,k)-ex_eps,
     $              gd_m(4,4), gd_m(1,4), gd_m(2,4), gd_m(3,4),
     $              gd_m(1,1), gd_m(2,2), gd_m(3,3), 
     $              gd_m(1,2), gd_m(2,3), gd_m(1,3),
     $              gu(4,4), gu(1,4), gu(2,4), gu(3,4),
     $              gu(1,1), gu(2,2), gu(3,3), 
     $              gu(1,2), gu(2,3), gu(1,3), exact_psi)
               do m=1,4
                  do n=m,4
                     gd1d(m,n,4) = (gd_p(m,n) - gd_m(m,n)) / (2.*ex_eps)
                  end do
               end do

C     Fill in remaining components of metric derivative.

               do l=1,4
                  do m=1,4
                     do n=1,m-1
                        gd1d(m,n,l) = gd1d(n,m,l)
                     end do
                  end do
               end do

C     Calculate K_ij.

               k3 = 0.d0
               do p=1,3
                  do q=p,3
                     do l=1,4
                        do m=1,4
                           k3(p,q) = k3(p,q) 
     $                          + s2d(l,p,q) * nu(m) * gd(l,m)
                           do n=1,4
                              k3(p,q) = k3(p,q) + 0.5d0 * (
     $                             s1d(l,p) * s1d(n,q) * nu(m)
     $                             + s1d(n,p) * s1d(m,q) * nu(l)
     $                             - s1d(l,p) * s1d(m,q) * nu(n) )
     $                             * gd1d(l,m,n)
                           end do
                        end do
                     end do
                  end do
               end do

               kxx(i,j,k) = - k3(1,1)
               kxy(i,j,k) = - k3(1,2)
               kxz(i,j,k) = - k3(1,3)
               kyy(i,j,k) = - k3(2,2)
               kyz(i,j,k) = - k3(2,3)
               kzz(i,j,k) = - k3(3,3)

            end do
         end do
      end do

C     Synchronize and bound slicetmp2, which contains dx^A/dt.

      call Exact__linear_extrap_one_bndry(CCTK_ARGUMENTS,slicetmp2x)
      call Exact__linear_extrap_one_bndry(CCTK_ARGUMENTS,slicetmp2y)
      call Exact__linear_extrap_one_bndry(CCTK_ARGUMENTS,slicetmp2z)
      call Exact__linear_extrap_one_bndry(CCTK_ARGUMENTS,slicetmp2t)
 
      call CCTK_SyncGroup(ierr,cctkGH,"Exact::Exact_slicetemp2")

C     Bound and synchronize the 3-metric and extrinsic curvature.

      call Exact__linear_extrap_one_bndry(CCTK_ARGUMENTS,gxx)
      call Exact__linear_extrap_one_bndry(CCTK_ARGUMENTS,gxy)
      call Exact__linear_extrap_one_bndry(CCTK_ARGUMENTS,gxz)
      call Exact__linear_extrap_one_bndry(CCTK_ARGUMENTS,gyy)
      call Exact__linear_extrap_one_bndry(CCTK_ARGUMENTS,gyz)
      call Exact__linear_extrap_one_bndry(CCTK_ARGUMENTS,gzz)

      call CCTK_SyncGroup(ierr,cctkGH,"admbase::metric")

      call Exact__linear_extrap_one_bndry(CCTK_ARGUMENTS,kxx)
      call Exact__linear_extrap_one_bndry(CCTK_ARGUMENTS,kxy)
      call Exact__linear_extrap_one_bndry(CCTK_ARGUMENTS,kxz)
      call Exact__linear_extrap_one_bndry(CCTK_ARGUMENTS,kyy)
      call Exact__linear_extrap_one_bndry(CCTK_ARGUMENTS,kyz)
      call Exact__linear_extrap_one_bndry(CCTK_ARGUMENTS,kzz)

      call CCTK_SyncGroup(ierr,cctkGH,"admbase::curv")
#endif

      return
      end
