! $Header$

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine CalcK (CCTK_ARGUMENTS)
  use cctk
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  CCTK_INT, parameter :: idummy=0
  CCTK_INT, parameter :: ik=kind(idummy)
  
  integer, parameter :: bndwidth = 1
  
  integer         :: len_boundary, len_options
  character(1000) :: str_boundary, str_options
  integer         :: options
  
  CCTK_REAL :: dt, dx(3)
  
  CCTK_REAL :: gama(3,3), gama_dot(3,3), dgama(3,3,3)
  CCTK_REAL :: alfa
  CCTK_REAL :: beta(3), dbeta(3,3)
  CCTK_REAL :: kk(3,3)
  
  integer :: imin(3), imax(3)
  integer :: i, j, k
  integer :: a, b, c
  integer :: ierr
  
  dt = CCTK_DELTA_TIME
  dx(:) = CCTK_DELTA_SPACE(:)
  
  call CCTK_FortranString &
       (len_boundary, extcurv_boundary, str_boundary)
  call CCTK_FortranString &
       (len_options, extcurv_boundary_options, str_options)
  call Util_TableCreateFromString (options, str_options)
  if (options<0) then
     call CCTK_WARN (0, "Parameter ""extcurv_boundary_options"" has an illegal syntax")
  end if
  
  imin(:) = 1+cctk_nghostzones(:)
  imax(:) = cctk_lsh(:)-cctk_nghostzones(:)
  where (cctk_bbox(1::2)/=0) imin(:) = 1+bndwidth
  where (cctk_bbox(2::2)/=0) imax(:) = cctk_lsh(:)-bndwidth
  
  ! Convert to physical metric, if necessary
  if (CCTK_EQUALS(metric_type, "static conformal")) then
     call ConfToPhysInPlace &
          (int(cctk_lsh(1),ik), int(cctk_lsh(2),ik), int(cctk_lsh(3),ik), &
          psi, gxx, gxy, gxz, gyy, gyz, gzz)
     call ConfToPhysInPlace &
          (int(cctk_lsh(1),ik), int(cctk_lsh(2),ik), int(cctk_lsh(3),ik), &
          psi, gxx_prev, gxy_prev, gxz_prev, gyy_prev, gyz_prev, gzz_prev)
#if 0
     call ConfToPhysInPlace &
          (int(cctk_lsh(1),ik), int(cctk_lsh(2),ik), int(cctk_lsh(3),ik), &
          psi, gxx_prev2, gxy_prev2, gxz_prev2, gyy_prev2, gyz_prev2, gzz_prev2)
#else
     call ConfToPhysInPlace &
          (int(cctk_lsh(1),ik), int(cctk_lsh(2),ik), int(cctk_lsh(3),ik), &
          psi, gxx_next, gxy_next, gxz_next, gyy_next, gyz_next, gzz_next)
#endif
  end if
  
  do k = imin(3), imax(3)
     do j = imin(2), imax(2)
        do i = imin(1), imax(1)
           
           gama(1,1) = gxx(i,j,k)
           gama(1,2) = gxy(i,j,k)
           gama(1,3) = gxz(i,j,k)
           gama(2,2) = gyy(i,j,k)
           gama(2,3) = gyz(i,j,k)
           gama(3,3) = gzz(i,j,k)
           gama(2,1) = gama(1,2)
           gama(3,1) = gama(1,3)
           gama(3,2) = gama(2,3)
           
#if 0
           gama_dot(1,1) = (-3*gxx(i,j,k) + 4*gxx_prev(i,j,k) - gxx_prev2(i,j,k)) / (2*dt)
           gama_dot(1,2) = (-3*gxy(i,j,k) + 4*gxy_prev(i,j,k) - gxy_prev2(i,j,k)) / (2*dt)
           gama_dot(1,3) = (-3*gxz(i,j,k) + 4*gxz_prev(i,j,k) - gxz_prev2(i,j,k)) / (2*dt)
           gama_dot(2,2) = (-3*gyy(i,j,k) + 4*gyy_prev(i,j,k) - gyy_prev2(i,j,k)) / (2*dt)
           gama_dot(2,3) = (-3*gyz(i,j,k) + 4*gyz_prev(i,j,k) - gyz_prev2(i,j,k)) / (2*dt)
           gama_dot(3,3) = (-3*gzz(i,j,k) + 4*gzz_prev(i,j,k) - gzz_prev2(i,j,k)) / (2*dt)
#else
           gama_dot(1,1) = (gxx_next(i,j,k) - gxx_prev(i,j,k)) / (2*dt)
           gama_dot(1,2) = (gxy_next(i,j,k) - gxy_prev(i,j,k)) / (2*dt)
           gama_dot(1,3) = (gxz_next(i,j,k) - gxz_prev(i,j,k)) / (2*dt)
           gama_dot(2,2) = (gyy_next(i,j,k) - gyy_prev(i,j,k)) / (2*dt)
           gama_dot(2,3) = (gyz_next(i,j,k) - gyz_prev(i,j,k)) / (2*dt)
           gama_dot(3,3) = (gzz_next(i,j,k) - gzz_prev(i,j,k)) / (2*dt)
#endif
           gama_dot(2,1) = gama_dot(1,2)
           gama_dot(3,1) = gama_dot(1,3)
           gama_dot(3,2) = gama_dot(2,3)
           
           dgama(1,1,1) = (gxx(i+1,j,k) - gxx(i-1,j,k)) / (2*dx(1))
           dgama(1,2,1) = (gxy(i+1,j,k) - gxy(i-1,j,k)) / (2*dx(1))
           dgama(1,3,1) = (gxz(i+1,j,k) - gxz(i-1,j,k)) / (2*dx(1))
           dgama(2,2,1) = (gyy(i+1,j,k) - gyy(i-1,j,k)) / (2*dx(1))
           dgama(2,3,1) = (gyz(i+1,j,k) - gyz(i-1,j,k)) / (2*dx(1))
           dgama(3,3,1) = (gzz(i+1,j,k) - gzz(i-1,j,k)) / (2*dx(1))
           dgama(1,1,2) = (gxx(i,j+1,k) - gxx(i,j-1,k)) / (2*dx(2))
           dgama(1,2,2) = (gxy(i,j+1,k) - gxy(i,j-1,k)) / (2*dx(2))
           dgama(1,3,2) = (gxz(i,j+1,k) - gxz(i,j-1,k)) / (2*dx(2))
           dgama(2,2,2) = (gyy(i,j+1,k) - gyy(i,j-1,k)) / (2*dx(2))
           dgama(2,3,2) = (gyz(i,j+1,k) - gyz(i,j-1,k)) / (2*dx(2))
           dgama(3,3,2) = (gzz(i,j+1,k) - gzz(i,j-1,k)) / (2*dx(2))
           dgama(1,1,3) = (gxx(i,j,k+1) - gxx(i,j,k-1)) / (2*dx(3))
           dgama(1,2,3) = (gxy(i,j,k+1) - gxy(i,j,k-1)) / (2*dx(3))
           dgama(1,3,3) = (gxz(i,j,k+1) - gxz(i,j,k-1)) / (2*dx(3))
           dgama(2,2,3) = (gyy(i,j,k+1) - gyy(i,j,k-1)) / (2*dx(3))
           dgama(2,3,3) = (gyz(i,j,k+1) - gyz(i,j,k-1)) / (2*dx(3))
           dgama(3,3,3) = (gzz(i,j,k+1) - gzz(i,j,k-1)) / (2*dx(3))
           dgama(2,1,:) = dgama(1,2,:)
           dgama(3,1,:) = dgama(1,3,:)
           dgama(3,2,:) = dgama(2,3,:)
           
           alfa = alp(i,j,k)
           
           beta(1) = betax(i,j,k)
           beta(2) = betay(i,j,k)
           beta(3) = betaz(i,j,k)
           
           dbeta(1,1) = (betax(i+1,j,k) - betax(i-1,j,k)) / (2*dx(1))
           dbeta(2,1) = (betay(i+1,j,k) - betay(i-1,j,k)) / (2*dx(1))
           dbeta(3,1) = (betaz(i+1,j,k) - betaz(i-1,j,k)) / (2*dx(1))
           dbeta(1,2) = (betax(i,j+1,k) - betax(i,j-1,k)) / (2*dx(2))
           dbeta(2,2) = (betay(i,j+1,k) - betay(i,j-1,k)) / (2*dx(2))
           dbeta(3,2) = (betaz(i,j+1,k) - betaz(i,j-1,k)) / (2*dx(2))
           dbeta(1,3) = (betax(i,j,k+1) - betax(i,j,k-1)) / (2*dx(3))
           dbeta(2,3) = (betay(i,j,k+1) - betay(i,j,k-1)) / (2*dx(3))
           dbeta(3,3) = (betaz(i,j,k+1) - betaz(i,j,k-1)) / (2*dx(3))
           
           ! d/dt gamma_ij = - 2 alpha K_ij
           !                 + gamma_kj d_i beta^k
           !                 + gamma_kj d_i beta^k
           !                 + beta^k d_k gamma_ij
           
           do a=1,3
              do b=1,3
                 kk(a,b) = - gama_dot(a,b)
                 do c=1,3
                    kk(a,b) = kk(a,b) + gama(c,b) * dbeta(c,a) &
                         &            + gama(a,c) * dbeta(c,b) &
                         &            + beta(c) * dgama(a,b,c)
                 end do
                 kk(a,b) = kk(a,b) / (2*alfa)
              end do
           end do
           
           kxx(i,j,k) = kk(1,1)
           kxy(i,j,k) = kk(1,2)
           kxz(i,j,k) = kk(1,3)
           kyy(i,j,k) = kk(2,2)
           kyz(i,j,k) = kk(2,3)
           kzz(i,j,k) = kk(3,3)
           
        end do
     end do
  end do
  
  call CartSymGN (ierr, cctkGH, "ADMBase::curv")
  if (ierr /= 0) call CCTK_WARN (0, "internal error")
  ierr = Boundary_SelectGroupForBC (cctkGH, int(CCTK_ALL_FACES,ik), &
       int(bndwidth,ik), int(options,ik), "ADMBase::curv", str_boundary)
  if (ierr /= 0) call CCTK_WARN (0, "internal error")
  
  ! Convert back to conformal metric, if necessary
  if (CCTK_EQUALS(metric_type, "static conformal")) then
     call PhysToConfInPlace &
          (int(cctk_lsh(1),ik), int(cctk_lsh(2),ik), int(cctk_lsh(3),ik), &
          psi, gxx, gxy, gxz, gyy, gyz, gzz)
     call PhysToConfInPlace &
          (int(cctk_lsh(1),ik), int(cctk_lsh(2),ik), int(cctk_lsh(3),ik), &
          psi, gxx_prev, gxy_prev, gxz_prev, gyy_prev, gyz_prev, gzz_prev)
#if 0
     call PhysToConfInPlace &
          (int(cctk_lsh(1),ik), int(cctk_lsh(2),ik), int(cctk_lsh(3),ik), &
          psi, gxx_prev2, gxy_prev2, gxz_prev2, gyy_prev2, gyz_prev2, gzz_prev2)
#else
     call PhysToConfInPlace &
          (int(cctk_lsh(1),ik), int(cctk_lsh(2),ik), int(cctk_lsh(3),ik), &
          psi, gxx_next, gxy_next, gxz_next, gyy_next, gyz_next, gzz_next)
#endif
  end if
  
end subroutine CalcK
