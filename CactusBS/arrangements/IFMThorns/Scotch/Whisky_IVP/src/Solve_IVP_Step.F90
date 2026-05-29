 /*@@
   @file      Solve_IVP_Step.F90
   @date      Sun Jul  7 15:55:26 2002
   @author    Ian Hawke
   @desc 
   A single step of the IVP routine using the calls to BAM.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

subroutine Solve_IVP_Step(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

#include "CactusEinstein/ADMMacros/src/macro/ADM_Spacing_declare.h"

  CCTK_INT :: i, j, k, nx, ny, nz
  CCTK_REAL :: DivergenceV, TwoPi

#include "CactusEinstein/ADMMacros/src/macro/CHR2_declare.h"

!!$  Here we need to calculate the things that may change at every step.
!!$  See Cook eqn 32.
!!$  This would be 
!!$  Atilde(ij)
!!$  Source terms for V
!!$  Note that the IVP_rho term does not change because it is the
!!$  conformal rho

#include "CactusEinstein/ADMMacros/src/macro/ADM_Spacing.h"

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  TwoPi = 8.d0 * atan(1.d0)

!!$ Initialise some variables, because loop is not going over the
!!$ whole grid
  IVP_Axx=0.0d0
  IVP_Axy=0.0d0
  IVP_Axz=0.0d0
  IVP_Ayy=0.0d0
  IVP_Ayz=0.0d0
  IVP_Azz=0.0d0

  do k = 3, nz - 2
    do j = 3, ny - 2
      do i = 3, nx - 2

!!$        Start with Atilde
!!$        Atilde^ij = (Ltilde V)^ij + Mtilde^ij

        IVP_Axx(i,j,k) = IVP_Mtildexx(i,j,k)
        IVP_Axy(i,j,k) = IVP_Mtildexy(i,j,k)
        IVP_Axz(i,j,k) = IVP_Mtildexz(i,j,k)
        IVP_Ayy(i,j,k) = IVP_Mtildeyy(i,j,k)
        IVP_Ayz(i,j,k) = IVP_Mtildeyz(i,j,k)
        IVP_Azz(i,j,k) = IVP_Mtildezz(i,j,k)
        
!!$        Add in the longitudinal derivatives of V.
!!$                         ---i  j   ---j i   2      ij ---   l
!!$        (Ltilde V)^ij == \ /  V  + \ / V  - - gamma   \ /  V
!!$                          V         V       3          V l
!!$
!!$        Divergence term first                         ^^^^^^^
!!$
!!$        ---   l    l    +--+l   m   /          i \     /
!!$        \ /  V  = V   + |      V  = | sqrt(g) V  |    / sqrt(g)
!!$         V l       ,l   |   lm      \            /,i /

        DivergenceV = ( &
            ( sqrt(IVP_detg(i+1,j,k)) * IVP_Vx(i+1,j,k) - &
              sqrt(IVP_detg(i-1,j,k)) * IVP_Vx(i-1,j,k) ) * i2dx + &
            ( sqrt(IVP_detg(i,j+1,k)) * IVP_Vy(i,j+1,k) - &
              sqrt(IVP_detg(i,j-1,k)) * IVP_Vy(i,j-1,k) ) * i2dy + &
            ( sqrt(IVP_detg(i,j,k+1)) * IVP_Vz(i,j,k+1) - &
              sqrt(IVP_detg(i,j,k-1)) * IVP_Vz(i,j,k-1) ) * i2dz ) / &
            sqrt(IVP_detg(i,j,k))

        IVP_Axx(i,j,k) = IVP_Axx(i,j,k) - &
             2.d0 / 3.d0 * DivergenceV * &
             IVP_Psi(i,j,k)**(-4) * IVP_uxx(i,j,k)
        IVP_Axy(i,j,k) = IVP_Axy(i,j,k) - &
             2.d0 / 3.d0 * DivergenceV * &
             IVP_Psi(i,j,k)**(-4) * IVP_uxy(i,j,k)
        IVP_Axz(i,j,k) = IVP_Axz(i,j,k) - &
             2.d0 / 3.d0 * DivergenceV * &
             IVP_Psi(i,j,k)**(-4) * IVP_uxz(i,j,k)
        IVP_Ayy(i,j,k) = IVP_Ayy(i,j,k) - &
             2.d0 / 3.d0 * DivergenceV * &
             IVP_Psi(i,j,k)**(-4) * IVP_uyy(i,j,k)
        IVP_Ayz(i,j,k) = IVP_Ayz(i,j,k) - &
             2.d0 / 3.d0 * DivergenceV * &
             IVP_Psi(i,j,k)**(-4) * IVP_uyz(i,j,k)
        IVP_Azz(i,j,k) = IVP_Azz(i,j,k) - &
             2.d0 / 3.d0 * DivergenceV * &
             IVP_Psi(i,j,k)**(-4) * IVP_uzz(i,j,k)

!!$        Then covariant derivatives (all covariant)
!!$    ---i  j  ---j  i        ik/  j   +--+j   l \        ik/  i   +--+i   l \
!!$    \ /  V + \ /  V  = gamma  | V  + |      V  | + gamma  | V  + |      V  |
!!$     V        V               \  ,k  |   kl    /          \  ,k  |   kl    /
!!$
!!$                        k     ij          ik  j         jk  i
!!$                     = V gamma     + gamma   V   + gamma   V
!!$                                ,k             ,k            ,k

        IVP_Axx(i,j,k) = IVP_Axx(i,j,k) + &
             IVP_Vx(i,j,k) * (IVP_uxx(i+1,j,k) - IVP_uxx(i-1,j,k)) * i2dx + &
             IVP_Vy(i,j,k) * (IVP_uxx(i,j+1,k) - IVP_uxx(i,j-1,k)) * i2dy + &
             IVP_Vz(i,j,k) * (IVP_uxx(i,j,k+1) - IVP_uxx(i,j,k-1)) * i2dz + &
             IVP_uxx(i,j,k) * (IVP_Vx(i+1,j,k) - IVP_Vx(i-1,j,k)) * i2dx + &
             IVP_uxy(i,j,k) * (IVP_Vx(i,j+1,k) - IVP_Vx(i,j-1,k)) * i2dy + &
             IVP_uxz(i,j,k) * (IVP_Vx(i,j,k+1) - IVP_Vx(i,j,k-1)) * i2dz + &
             IVP_uxx(i,j,k) * (IVP_Vx(i+1,j,k) - IVP_Vx(i-1,j,k)) * i2dx + &
             IVP_uxy(i,j,k) * (IVP_Vx(i,j+1,k) - IVP_Vx(i,j-1,k)) * i2dy + &
             IVP_uxz(i,j,k) * (IVP_Vx(i,j,k+1) - IVP_Vx(i,j,k-1)) * i2dz 

        IVP_Axy(i,j,k) = IVP_Axy(i,j,k) + &
             IVP_Vx(i,j,k) * (IVP_uxy(i+1,j,k) - IVP_uxy(i-1,j,k)) * i2dx + &
             IVP_Vy(i,j,k) * (IVP_uxy(i,j+1,k) - IVP_uxy(i,j-1,k)) * i2dy + &
             IVP_Vz(i,j,k) * (IVP_uxy(i,j,k+1) - IVP_uxy(i,j,k-1)) * i2dz + &
             IVP_uxx(i,j,k) * (IVP_Vy(i+1,j,k) - IVP_Vy(i-1,j,k)) * i2dx + &
             IVP_uxy(i,j,k) * (IVP_Vy(i,j+1,k) - IVP_Vy(i,j-1,k)) * i2dy + &
             IVP_uxz(i,j,k) * (IVP_Vy(i,j,k+1) - IVP_Vy(i,j,k-1)) * i2dz + &
             IVP_uxy(i,j,k) * (IVP_Vx(i+1,j,k) - IVP_Vx(i-1,j,k)) * i2dx + &
             IVP_uyy(i,j,k) * (IVP_Vx(i,j+1,k) - IVP_Vx(i,j-1,k)) * i2dy + &
             IVP_uyz(i,j,k) * (IVP_Vx(i,j,k+1) - IVP_Vx(i,j,k-1)) * i2dz 

        IVP_Axz(i,j,k) = IVP_Axz(i,j,k) + &
             IVP_Vx(i,j,k) * (IVP_uxz(i+1,j,k) - IVP_uxz(i-1,j,k)) * i2dx + &
             IVP_Vy(i,j,k) * (IVP_uxz(i,j+1,k) - IVP_uxz(i,j-1,k)) * i2dy + &
             IVP_Vz(i,j,k) * (IVP_uxz(i,j,k+1) - IVP_uxz(i,j,k-1)) * i2dz + &
             IVP_uxx(i,j,k) * (IVP_Vz(i+1,j,k) - IVP_Vz(i-1,j,k)) * i2dx + &
             IVP_uxy(i,j,k) * (IVP_Vz(i,j+1,k) - IVP_Vz(i,j-1,k)) * i2dy + &
             IVP_uxz(i,j,k) * (IVP_Vz(i,j,k+1) - IVP_Vz(i,j,k-1)) * i2dz + &
             IVP_uxz(i,j,k) * (IVP_Vx(i+1,j,k) - IVP_Vx(i-1,j,k)) * i2dx + &
             IVP_uyz(i,j,k) * (IVP_Vx(i,j+1,k) - IVP_Vx(i,j-1,k)) * i2dy + &
             IVP_uzz(i,j,k) * (IVP_Vx(i,j,k+1) - IVP_Vx(i,j,k-1)) * i2dz 

        IVP_Ayy(i,j,k) = IVP_Ayy(i,j,k) + &
             IVP_Vx(i,j,k) * (IVP_uyy(i+1,j,k) - IVP_uyy(i-1,j,k)) * i2dx + &
             IVP_Vy(i,j,k) * (IVP_uyy(i,j+1,k) - IVP_uyy(i,j-1,k)) * i2dy + &
             IVP_Vz(i,j,k) * (IVP_uyy(i,j,k+1) - IVP_uyy(i,j,k-1)) * i2dz + &
             IVP_uxy(i,j,k) * (IVP_Vy(i+1,j,k) - IVP_Vy(i-1,j,k)) * i2dx + &
             IVP_uyy(i,j,k) * (IVP_Vy(i,j+1,k) - IVP_Vy(i,j-1,k)) * i2dy + &
             IVP_uyz(i,j,k) * (IVP_Vy(i,j,k+1) - IVP_Vy(i,j,k-1)) * i2dz + &
             IVP_uxy(i,j,k) * (IVP_Vy(i+1,j,k) - IVP_Vy(i-1,j,k)) * i2dx + &
             IVP_uyy(i,j,k) * (IVP_Vy(i,j+1,k) - IVP_Vy(i,j-1,k)) * i2dy + &
             IVP_uyz(i,j,k) * (IVP_Vy(i,j,k+1) - IVP_Vy(i,j,k-1)) * i2dz 

        IVP_Ayz(i,j,k) = IVP_Ayz(i,j,k) + &
             IVP_Vx(i,j,k) * (IVP_uyz(i+1,j,k) - IVP_uyz(i-1,j,k)) * i2dx + &
             IVP_Vy(i,j,k) * (IVP_uyz(i,j+1,k) - IVP_uyz(i,j-1,k)) * i2dy + &
             IVP_Vz(i,j,k) * (IVP_uyz(i,j,k+1) - IVP_uyz(i,j,k-1)) * i2dz + &
             IVP_uxy(i,j,k) * (IVP_Vz(i+1,j,k) - IVP_Vz(i-1,j,k)) * i2dx + &
             IVP_uyy(i,j,k) * (IVP_Vz(i,j+1,k) - IVP_Vz(i,j-1,k)) * i2dy + &
             IVP_uyz(i,j,k) * (IVP_Vz(i,j,k+1) - IVP_Vz(i,j,k-1)) * i2dz + &
             IVP_uxz(i,j,k) * (IVP_Vy(i+1,j,k) - IVP_Vy(i-1,j,k)) * i2dx + &
             IVP_uyy(i,j,k) * (IVP_Vy(i,j+1,k) - IVP_Vy(i,j-1,k)) * i2dy + &
             IVP_uzz(i,j,k) * (IVP_Vy(i,j,k+1) - IVP_Vy(i,j,k-1)) * i2dz 

        IVP_Azz(i,j,k) = IVP_Azz(i,j,k) + &
             IVP_Vx(i,j,k) * (IVP_uzz(i+1,j,k) - IVP_uzz(i-1,j,k)) * i2dx + &
             IVP_Vy(i,j,k) * (IVP_uzz(i,j+1,k) - IVP_uzz(i,j-1,k)) * i2dy + &
             IVP_Vz(i,j,k) * (IVP_uzz(i,j,k+1) - IVP_uzz(i,j,k-1)) * i2dz + &
             IVP_uxz(i,j,k) * (IVP_Vz(i+1,j,k) - IVP_Vz(i-1,j,k)) * i2dx + &
             IVP_uyz(i,j,k) * (IVP_Vz(i,j+1,k) - IVP_Vz(i,j-1,k)) * i2dy + &
             IVP_uzz(i,j,k) * (IVP_Vz(i,j,k+1) - IVP_Vz(i,j,k-1)) * i2dz + &
             IVP_uxz(i,j,k) * (IVP_Vz(i+1,j,k) - IVP_Vz(i-1,j,k)) * i2dx + &
             IVP_uyz(i,j,k) * (IVP_Vz(i,j+1,k) - IVP_Vz(i,j-1,k)) * i2dy + &
             IVP_uzz(i,j,k) * (IVP_Vz(i,j,k+1) - IVP_Vz(i,j,k-1)) * i2dz 

      end do
    end do
  end do

!!$  Now calculate all the source terms.

!!$  First, the term with psi^-7
!!$  This is effectively the contraction of Atilde with itself (1/8 of it)
!!$  note: g?? is the conformal metric and has to be.

  IVP_Psim7Source = &
              gxx**2    * IVP_Axx**2 + &
       4.d0 * gxx * gxy * IVP_Axx * IVP_Axy + &
       2.d0 * gxy**2    * IVP_Axy**2 + &
       2.d0 * gxx * gyy * IVP_Axy**2 + &
       4.d0 * gxx * gxz * IVP_Axx * IVP_Axz + &
       4.d0 * gxy * gxz * IVP_Axy * IVP_Axz + &
       4.d0 * gxx * gyz * IVP_Axy * IVP_Axz + &
       2.d0 * gxz**2    * IVP_Axz**2 + &
       2.d0 * gxx * gzz * IVP_Axz**2 + &
       2.d0 * gxy**2    * IVP_Axx * IVP_Ayy + &
       4.d0 * gxy * gyy * IVP_Axy * IVP_Ayy + &
       4.d0 * gxy * gyz * IVP_Axz * IVP_Ayy + &
              gyy**2    * IVP_Ayy**2 + &
       4.d0 * gxy * gxz * IVP_Axx * IVP_Ayz + &
       4.d0 * gxz * gyy * IVP_Axy * IVP_Ayz + &
       4.d0 * gxy * gyz * IVP_Axy * IVP_Ayz + &
       4.d0 * gxz * gyz * IVP_Axz * IVP_Ayz + &
       4.d0 * gxy * gzz * IVP_Axz * IVP_Ayz + &
       4.d0 * gyy * gyz * IVP_Ayy * IVP_Ayz + &
       2.d0 * gyz**2    * IVP_Ayz**2 + &
       2.d0 * gyy * gzz * IVP_Ayz**2 + &
       2.d0 * gxz**2    * IVP_Axx * IVP_Azz + &
       4.d0 * gxz * gyz * IVP_Axy * IVP_Azz + &
       4.d0 * gxz * gzz * IVP_Axz * IVP_Azz + &
       2.d0 * gyz**2    * IVP_Ayy * IVP_Azz + &
       4.d0 * gyz * gzz * IVP_Ayz * IVP_Azz + &
              gzz**2    * IVP_Azz**2

  IVP_Psim7Source = IVP_Psim7Source / 8.d0

!!$  Scalar Curvature contributions

  IVP_Mlinear = - IVP_RicciScalar / 8.d0

!!$  Zero out the source term

  IVP_Nsource = 0.d0

!!$  Set up the psi^5 term, which is trK^2

  IVP_Psi5Source = - IVP_trK**2 / 12.d0

!!$  Add in the matter contributions

  if (.not. (Whisky_IVP_Ignore_Matter == 1)) then
    IVP_Psim3Source = TwoPi * IVP_rho
    IVP_Vx_Source = -4.d0 * TwoPi * IVP_jx
    IVP_Vy_Source = -4.d0 * TwoPi * IVP_jy
    IVP_Vz_Source = -4.d0 * TwoPi * IVP_jz
  else
    IVP_Psim3Source = 0.0d0
    IVP_Vx_Source = 0.0d0
    IVP_Vy_Source = 0.0d0
    IVP_Vz_Source = 0.0d0
  endif

!!$  Finish setting the vector source terms

  do k = 3, nz - 2
    do j = 3, ny - 2
      do i = 3, nx - 2

!!$     First add the divergence of trK
!!$     Again note, that we always have here conformal_state=0, because
!!$      we want to have the macros of the conformal metric

        IVP_Vx_Source(i,j,k) = IVP_Vx_Source(i,j,k) - 2.d0 / 3.d0 * &
             IVP_Psi(i,j,k)**6 * &
             ( IVP_uxx(i,j,k)*(IVP_trK(i+1,j,k) - IVP_trK(i-1,j,k)) * i2dx + &
               IVP_uxy(i,j,k)*(IVP_trK(i,j+1,k) - IVP_trK(i,j-1,k)) * i2dy + &
               IVP_uxz(i,j,k)*(IVP_trK(i,j,k+1) - IVP_trK(i,j,k-1)) * i2dz )
        IVP_Vy_Source(i,j,k) = IVP_Vy_Source(i,j,k) - 2.d0 / 3.d0 * &
             IVP_Psi(i,j,k)**6 * &
             ( IVP_uxy(i,j,k)*(IVP_trK(i+1,j,k) - IVP_trK(i-1,j,k)) * i2dx + &
               IVP_uyy(i,j,k)*(IVP_trK(i,j+1,k) - IVP_trK(i,j-1,k)) * i2dy + &
               IVP_uyz(i,j,k)*(IVP_trK(i,j,k+1) - IVP_trK(i,j,k-1)) * i2dz )
        IVP_Vz_Source(i,j,k) = IVP_Vz_Source(i,j,k) - 2.d0 / 3.d0 * &
             IVP_Psi(i,j,k)**6 * &
             ( IVP_uxz(i,j,k)*(IVP_trK(i+1,j,k) - IVP_trK(i-1,j,k)) * i2dx + &
               IVP_uyz(i,j,k)*(IVP_trK(i,j+1,k) - IVP_trK(i,j-1,k)) * i2dy + &
               IVP_uzz(i,j,k)*(IVP_trK(i,j,k+1) - IVP_trK(i,j,k-1)) * i2dz )

!!$        Now calculate the divergence of Mtilde

!!$        As we are taking a contraction it is written 
!!$           
!!$        ---
!!$        \ /  M^(ij) = {sqrt(detg) M^(ij)},i / sqrt(detg) 
!!$         V i           + Gamma^j_(ik) M^(ik)

!!$        Derivatives first

        IVP_Vx_Source(i,j,k) = IVP_Vx_Source(i,j,k) + &
             ( (sqrt(IVP_detg(i+1,j,k)) * IVP_Mtildexx(i+1,j,k) - &
                sqrt(IVP_detg(i-1,j,k)) * IVP_Mtildexx(i-1,j,k) ) * i2dx + &
               (sqrt(IVP_detg(i,j+1,k)) * IVP_Mtildexy(i,j+1,k) - &
                sqrt(IVP_detg(i,j-1,k)) * IVP_Mtildexy(i,j-1,k) ) * i2dy + &
               (sqrt(IVP_detg(i,j,k+1)) * IVP_Mtildexz(i,j,k+1) - &
                sqrt(IVP_detg(i,j,k-1)) * IVP_Mtildexz(i,j,k-1) ) * i2dz ) / &
             sqrt(IVP_detg(i,j,k))

        IVP_Vy_Source(i,j,k) = IVP_Vy_Source(i,j,k) + &
             ( (sqrt(IVP_detg(i+1,j,k)) * IVP_Mtildexy(i+1,j,k) - &
                sqrt(IVP_detg(i-1,j,k)) * IVP_Mtildexy(i-1,j,k) ) * i2dx + &
               (sqrt(IVP_detg(i,j+1,k)) * IVP_Mtildeyy(i,j+1,k) - &
                sqrt(IVP_detg(i,j-1,k)) * IVP_Mtildeyy(i,j-1,k) ) * i2dy + &
               (sqrt(IVP_detg(i,j,k+1)) * IVP_Mtildeyz(i,j,k+1) - &
                sqrt(IVP_detg(i,j,k-1)) * IVP_Mtildeyz(i,j,k-1) ) * i2dz ) / &
             sqrt(IVP_detg(i,j,k))

        IVP_Vz_Source(i,j,k) = IVP_Vz_Source(i,j,k) + &
             ( (sqrt(IVP_detg(i+1,j,k)) * IVP_Mtildexz(i+1,j,k) - &
                sqrt(IVP_detg(i-1,j,k)) * IVP_Mtildexz(i-1,j,k) ) * i2dx + &
               (sqrt(IVP_detg(i,j+1,k)) * IVP_Mtildeyz(i,j+1,k) - &
                sqrt(IVP_detg(i,j-1,k)) * IVP_Mtildeyz(i,j-1,k) ) * i2dy + &
               (sqrt(IVP_detg(i,j,k+1)) * IVP_Mtildezz(i,j,k+1) - &
                sqrt(IVP_detg(i,j,k-1)) * IVP_Mtildezz(i,j,k-1) ) * i2dz ) / &
             sqrt(IVP_detg(i,j,k))

!!$        Then get the Christoffel symbols and compute the last part

#include "CactusEinstein/ADMMacros/src/macro/CHR2_guts.h"        

        IVP_Vx_Source(i,j,k) = IVP_Vx_Source(i,j,k) + &
                    CHR2_XXX * IVP_Mtildexx(i,j,k) + &
             2.d0 * CHR2_XXY * IVP_Mtildexy(i,j,k) + &
             2.d0 * CHR2_XXZ * IVP_Mtildexz(i,j,k) + &
                    CHR2_XYY * IVP_Mtildeyy(i,j,k) + &
             2.d0 * CHR2_XYZ * IVP_Mtildeyz(i,j,k) + &
                    CHR2_XZZ * IVP_Mtildezz(i,j,k) 

        IVP_Vy_Source(i,j,k) = IVP_Vy_Source(i,j,k) + &
                    CHR2_YXX * IVP_Mtildexx(i,j,k) + &
             2.d0 * CHR2_YXY * IVP_Mtildexy(i,j,k) + &
             2.d0 * CHR2_YXZ * IVP_Mtildexz(i,j,k) + &
                    CHR2_YYY * IVP_Mtildeyy(i,j,k) + &
             2.d0 * CHR2_YYZ * IVP_Mtildeyz(i,j,k) + &
                    CHR2_YZZ * IVP_Mtildezz(i,j,k) 

        IVP_Vz_Source(i,j,k) = IVP_Vz_Source(i,j,k) + &
                    CHR2_ZXX * IVP_Mtildexx(i,j,k) + &
             2.d0 * CHR2_ZXY * IVP_Mtildexy(i,j,k) + &
             2.d0 * CHR2_ZXZ * IVP_Mtildexz(i,j,k) + &
                    CHR2_ZYY * IVP_Mtildeyy(i,j,k) + &
             2.d0 * CHR2_ZYZ * IVP_Mtildeyz(i,j,k) + &
                    CHR2_ZZZ * IVP_Mtildezz(i,j,k) 

      end do
    end do
  end do

#include "CactusEinstein/ADMMacros/src/macro/CHR2_undefine.h"

!!$  Finally, call BAM
!!$  conformal_state=3
  call Call_BAM_Nonlin(cctkGH, Whisky_IVP_Tolerance, IVP_Check_NaNs)
  call Call_BAM_VecLap(cctkGH, Whisky_IVP_Tolerance, IVP_Check_NaNs)
!!$  conformal_state=0
end subroutine Solve_IVP_Step

