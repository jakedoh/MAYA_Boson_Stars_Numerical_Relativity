c This subroutine calculates the 4-metric and its inverse at an event,
c taking into account an optional Lorentz boost and an optional rotation.
c The model is first rotated and then boosted, such that the boost is
c applied to the rotated model.
c $Header$
c
c The coordinates are
c   Cx(a) = Cactus $x^a$
c   Mx(a) = Model  $X^a$
c The 4-metrics are
c   Cgdd(a,b) = Cactus $g_{ab}$            Cguu(a,b) = Cactus $g^{ab}$
c   Mgdd(a,b) = Model  $g_{ab}$            Mguu(a,b) = Model  $g^{ab}$
c
c For a definition of the Euler angles in the conventions used below, see
c    http://mathworld.wolfram.com/EulerAngles.html
c Another useful resource may be
c    http://en.wikipedia.org/wiki/Euler_angles
c although this uses (on 2006-11-29) different conventions.
c
c This file is copyright (c) 2003 by Jonathan Thornburg <jthorn@aei.mpg.de>.
c This file is covered by the GNU GPL license; see the files ../README
c and ../COPYING for details. 
c

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#include "param_defs.inc"

      subroutine Exact__metric(
     $     decoded_exact_model,
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz,
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $     gutt, gutx, guty, gutz,
     $     guxx, guyy, guzz, guxy, guyz, guxz,
     $     psi)

      implicit none
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS

c input arguments
      CCTK_INT decoded_exact_model
      CCTK_REAL x, y, z, t

c output arguments
      CCTK_REAL gdtt, gdtx, gdty, gdtz,
     $          gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $          gutt, gutx, guty, gutz,
     $          guxx, guyy, guzz, guxy, guyz, guxz,
     $          psi

c static local variables describing Lorentz transformation
      logical   firstcall
      data      firstcall /.true./
      CCTK_REAL gamma
      CCTK_REAL vv(3), nn(3)
      CCTK_REAL parallel(3,3), perp(3,3)
      CCTK_REAL Cx_par(3), Cx_perp(3)
      CCTK_REAL partial_Mx_wrt_Cx(0:3,0:3)
      CCTK_REAL partial_Cx_wrt_Mx(0:3,0:3)
      CCTK_REAL R(0:3,0:3)
      save      firstcall
      save      gamma
      save      vv
      save      parallel, perp
      save      partial_Mx_wrt_Cx
      save      partial_Cx_wrt_Mx
      save      R
c$omp threadprivate (firstcall, gamma, vv, parallel, perp, 
c$omp+               partial_Mx_wrt_Cx,partial_Cx_wrt_Mx, R)

c coordinates and 4-metric
      CCTK_REAL Cx(0:3)
      CCTK_REAL Cgdd(0:3,0:3), Cguu(0:3,0:3)
      CCTK_REAL Mx(0:3)
      CCTK_REAL Mgdd(0:3,0:3), Mguu(0:3,0:3)
      CCTK_REAL Nx(0:3)
      CCTK_REAL Ngdd(0:3,0:3), Nguu(0:3,0:3)

c misc temps
      CCTK_REAL vnorm, vnormsq
      CCTK_REAL delta_ij
      CCTK_REAL Cx_par_i, Cx_perp_i
      CCTK_REAL vdotCx
      CCTK_REAL Cgdd_ab, Cguu_ab
      CCTK_REAL cos_phi, sin_phi
      CCTK_REAL cos_theta, sin_theta
      CCTK_REAL cos_psi, sin_psi
      CCTK_REAL R_phi(0:3,0:3), R_theta(0:3,0:3), R_psi(0:3,0:3)
      character*1000 warn_buffer

c flags, array indices, etc
      logical Tmunu_flag
      integer i, j, k, l
      integer Ca, Cb, MA, MB

c constants
      integer n
      parameter (n = 3)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c optimized fast-path if no Lorentz boost and no rotation
c
      if (       (boost_vx .eq. 0.0)
     $     .and. (boost_vy .eq. 0.0)
     $     .and. (boost_vz .eq. 0.0)
     $     .and. (rotation_euler_phi   .eq. 0.0)
     $     .and. (rotation_euler_theta .eq. 0.0)
     $     .and. (rotation_euler_psi   .eq. 0.0)) then
         call Exact__metric_for_model(
     $        decoded_exact_model,
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz,
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz,
     $        guxx, guyy, guzz, guxy, guyz, guxz,
     $        psi,
     $        Tmunu_flag)
         return
      end if

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c the rest of this function is the Lorentz-boost case:
c - Lorentz-transform Cactus coordinates --> Model coordinates
c - compute Model 4-metric and inverse at Model coordinates
c - tensor-transform 4-metric and inverse from Model coordinates
c   --> Cactus coordinates
c
c All the equations used are given in ../doc/documentation.tex
c

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c compute Lorentz transformation information on first call
c
      if (firstcall) then

c boost velocity
         vv(1) = boost_vx
         vv(2) = boost_vy
         vv(3) = boost_vz

c Lorentz gamma factor, unit vector in direction of boost velocity
         vnormsq = 0.0
         do i = 1,n
            vnormsq = vnormsq + vv(i)*vv(i)
         end do
         gamma = 1.0 / sqrt(1.0 - vnormsq)
         vnorm = sqrt(vnormsq)
         if (       (boost_vx .eq. 0.0)
     $        .and. (boost_vy .eq. 0.0)
     $        .and. (boost_vz .eq. 0.0)) then
            nn(1) = 1
            nn(2) = 0
            nn(3) = 0
         else
            do i = 1,n
               nn(i) = vv(i) / vnorm
            end do
         end if

c projection operators parallel(*,*) and perp(*,*)
         do j = 1,n
            do i = 1,n
               parallel(i,j) = nn(i) * nn(j)
               if (i .eq. j) then
                  delta_ij = 1.0
               else
                  delta_ij = 0.0
               end if
               perp(i,j) = delta_ij - parallel(i,j)
            end do
         end do

c partial derivatives of Model coordinates with respect to Cactus coordinates
         partial_Mx_wrt_Cx(0,0) = gamma
         do i = 1,n
            partial_Mx_wrt_Cx(0,i) = -gamma*vv(i)
         end do
         do i = 1,n
            partial_Mx_wrt_Cx(i,0) = -gamma*vv(i)
            do j=1,n
               partial_Mx_wrt_Cx(i,j) = gamma*parallel(i,j) + perp(i,j)
            end do
         end do

c partial derivatives of Cactus coordinates with respect to Model coordinates
         partial_Cx_wrt_Mx(0,0) = gamma
         do i = 1,n
            partial_Cx_wrt_Mx(0,i) = + gamma*vv(i)
         end do
         do i = 1,n
            partial_Cx_wrt_Mx(i,0) = + gamma*vv(i)
            do j=1,n
               partial_Cx_wrt_Mx(i,j) = gamma*parallel(i,j) + perp(i,j)
            end do
         end do

c Sines and cosines of rotation angles
         cos_phi   = cos (rotation_euler_phi)
         sin_phi   = sin (rotation_euler_phi)
         cos_theta = cos (rotation_euler_theta)
         sin_theta = sin (rotation_euler_theta)
         cos_psi   = cos (rotation_euler_psi)
         sin_psi   = sin (rotation_euler_psi)

c Set up individual rotation matrices
         R_phi(0,0) = 1
         R_phi(0,1) = 0
         R_phi(0,2) = 0
         R_phi(0,3) = 0
         R_phi(1,0) = 0
         R_phi(1,1) = + cos_phi
         R_phi(1,2) = + sin_phi
         R_phi(1,3) = 0
         R_phi(2,0) = 0
         R_phi(2,1) = - sin_phi
         R_phi(2,2) = + cos_phi
         R_phi(2,3) = 0
         R_phi(3,0) = 0
         R_phi(3,1) = 0
         R_phi(3,2) = 0
         R_phi(3,3) = 1

         R_theta(0,0) = 1
         R_theta(0,1) = 0
         R_theta(0,2) = 0
         R_theta(0,3) = 0
         R_theta(1,0) = 0
         R_theta(1,1) = 1
         R_theta(1,2) = 0
         R_theta(1,3) = 0
         R_theta(2,0) = 0
         R_theta(2,1) = 0
         R_theta(2,2) = + cos_theta
         R_theta(2,3) = + sin_theta
         R_theta(3,0) = 0
         R_theta(3,1) = 0
         R_theta(3,2) = - sin_theta
         R_theta(3,3) = + cos_theta

         R_psi(0,0) = 1
         R_psi(0,1) = 0
         R_psi(0,2) = 0
         R_psi(0,3) = 0
         R_psi(1,0) = 0
         R_psi(1,1) = + cos_psi
         R_psi(1,2) = + sin_psi
         R_psi(1,3) = 0
         R_psi(2,0) = 0
         R_psi(2,1) = - sin_psi
         R_psi(2,2) = + cos_psi
         R_psi(2,3) = 0
         R_psi(3,0) = 0
         R_psi(3,1) = 0
         R_psi(3,2) = 0
         R_psi(3,3) = 1

c Combine individual rotation matrices
         do i = 0,n
            do j = 0,n
               R(i,j) = 0
               do k = 0,n
                  do l = 0,n
                     R(i,j) = R(i,j) + R_psi(i,k) * R_theta(k,l) * R_phi(l,j)
                  end do
               end do
            end do
         end do

c     Notes that help me (Erik Schnetter) think:
c     This considers a rotation with phi=0, theta=pi/2, psi=0.
c     Nx(i) = Nx(i) + R(j,i) * Mx(j)
c     Nx(1) = - Mx(3)
c     Nx(3) =   Mx(1)
c     Mgxx(1,1) =   Ngxx(3,3)
c     Mgxx(1,3) = - Ngxx(1,3)
c     Mgxx(3,3) =   Ngxx(1,1)
c     Mgxx(i,j) = R(i,k) R(j,l) Ngxx(k,l)  [correct]
c     Mgxx(i,j) = R(k,i) R(l,j) Ngxx(k,l)  [correct]
c     Mbetax(1) = - Nbetax(3)
c     Mbetax(3) =   Nbetax(1)
c     Mbetax(i) + R(i,j) Nbetax(j)   [wrong]
c     Mbetax(i) + R(j,i) Nbetax(j)   [correct]

         firstcall = .false.
      end if

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c compute flat-space components of Cx(*) parallel and perpendicular to vv(*)
c
      Cx(0) = t
      Cx(1) = x
      Cx(2) = y
      Cx(3) = z

      do i=1,n
         Cx_par_i  = 0.0
         Cx_perp_i = 0.0
         do j=1,n
            Cx_par_i  = Cx_par_i  + parallel(i,j)*Cx(j)
            Cx_perp_i = Cx_perp_i + perp    (i,j)*Cx(j)
         end do
         Cx_par (i) = Cx_par_i
         Cx_perp(i) = Cx_perp_i
      end do

c
c Lorentz-transform and rotate the Cactus coordinate
c to get the Model coordinates
c

c Boost
      vdotCx = 0.0
      do i = 1,n
         vdotCx = vdotCx + vv(i)*Cx(i)
      end do

      Mx(0) = gamma * (Cx(0) - vdotCx)
      do i=1,n
         Mx(i) = gamma * (Cx_par(i) - vv(i)*Cx(0)) + Cx_perp(i)
      end do

c Rotation
      do i=0,n
         Nx(i) = 0
         do j = 0,n
            Nx(i) = Nx(i) + R(i,j) * Mx(j)
         end do
      end do


c
c compute the Model 4-metric and inverse 4-metric at the Model coordinates
c
      call Exact__metric_for_model(
     $     decoded_exact_model,
     $     Nx(1), Nx(2), Nx(3), Nx(0),
     $     Ngdd(0,0), Ngdd(0,1), Ngdd(0,2), Ngdd(0,3),
     $     Ngdd(1,1), Ngdd(2,2), Ngdd(3,3),
     $     Ngdd(1,2), Ngdd(2,3), Ngdd(1,3),
     $     Nguu(0,0), Nguu(0,1), Nguu(0,2), Nguu(0,3),
     $     Nguu(1,1), Nguu(2,2), Nguu(3,3),
     $     Nguu(1,2), Nguu(2,3), Nguu(1,3),
     $     psi, Tmunu_flag)

      if (Tmunu_flag) then
         write (warn_buffer, '(a,i8,a,a)')
     $        'exact_model = ', decoded_exact_model,
     $        'sets the stress-energy tensor',
     $        ' ==> we cannot Lorentz-boost or rotate it!'
         call CCTK_WARN(0, warn_buffer)
      end if

c
c symmetrize the Model 4-metric and inverse 4-metric arrays
c (the Exact__metric_for_model() call only set the upper triangles)
c
      Ngdd(1,0) = Ngdd(0,1)
      Ngdd(2,0) = Ngdd(0,2)
      Ngdd(2,1) = Ngdd(1,2)
      Ngdd(3,0) = Ngdd(0,3)
      Ngdd(3,1) = Ngdd(1,3)
      Ngdd(3,2) = Ngdd(2,3)

      Nguu(1,0) = Nguu(0,1)
      Nguu(2,0) = Nguu(0,2)
      Nguu(2,1) = Nguu(1,2)
      Nguu(3,0) = Nguu(0,3)
      Nguu(3,1) = Nguu(1,3)
      Nguu(3,2) = Nguu(2,3)

c
c tensor-transorm (the upper triangle of) the 4-metric and inverse 4-metric
c

c Rotations
      do i = 0,n
         do j = 0,n
            Mgdd(i,j) = 0
            Mguu(i,j) = 0
            do k = 0,n
               do l = 0,n
                  Mgdd(i,j) = Mgdd(i,j) + R(k,i) * R(l,j) * Ngdd(k,l)
c     The inverse of R is also its transpose.  That means that the
c     transpose of the inverse, which you would use for g^ij, is just R
c     again.
                  Mguu(i,j) = Mguu(i,j) + R(k,i) * R(l,j) * Nguu(k,l)
               end do
            end do
         end do
      end do

c Boost
      do Ca = 0,n
         do Cb = Ca,n
            Cgdd_ab = 0.0
            do Ma = 0,n
               do Mb = 0,n
                  Cgdd_ab = Cgdd_ab
     $                 +  Mgdd(Ma,Mb)
     $                    * partial_Mx_wrt_Cx(Ma,Ca)
     $                    * partial_Mx_wrt_Cx(Mb,Cb)
               end do
            end do
            Cgdd(Ca,Cb) = Cgdd_ab
         end do
      end do

      do Ca = 0,n
         do Cb = Ca,n
            Cguu_ab = 0.0
            do Ma = 0,n
               do Mb = 0,n
                  Cguu_ab = Cguu_ab
     $                 +  Mguu(Ma,Mb)
     $                    * partial_Cx_wrt_Mx(Ca,Ma)
     $                    * partial_Cx_wrt_Mx(Cb,Mb)
               end do
            end do
            Cguu(Ca,Cb) = Cguu_ab
         end do
      end do

c
c unpack the Cactus-coordinates 4-metric and inverse 4-metric
c into the corresponding output arguments
c
      gdtt = Cgdd(0,0)
      gdtx = Cgdd(0,1)
      gdty = Cgdd(0,2)
      gdtz = Cgdd(0,3)
      gdxx = Cgdd(1,1)
      gdxy = Cgdd(1,2)
      gdxz = Cgdd(1,3)
      gdyy = Cgdd(2,2)
      gdyz = Cgdd(2,3)
      gdzz = Cgdd(3,3)

      gutt = Cguu(0,0)
      gutx = Cguu(0,1)
      guty = Cguu(0,2)
      gutz = Cguu(0,3)
      guxx = Cguu(1,1)
      guxy = Cguu(1,2)
      guxz = Cguu(1,3)
      guyy = Cguu(2,2)
      guyz = Cguu(2,3)
      guzz = Cguu(3,3)

      end
