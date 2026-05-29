  /*@@
   @file      Whisky_IVP.F90
   @date      Sun Jul  7 15:06:31 2002
   @author    Ian Hawke
   @desc 
   Solve the IVP for a given matter source using BAM.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

!!$#define DEBUG_WHISKY_IVP

subroutine Whisky_IVP(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

#include "CactusEinstein/ADMMacros/src/macro/ADM_Spacing_declare.h"

  CCTK_REAL :: Max_Change_Psi
  CCTK_INT :: Reduction_Handle, ChangePsiIndex, ierr
  CCTK_INT :: i, j, k, nx, ny, nz
  CCTK_REAL :: DivergenceV

  CCTK_INT :: count
  CCTK_REAL :: j2, rhoold, rhonew, pressnew, IVP_K, IVP_Gamma, &
       epsnew, w_lorentznew, rhoenthalpy, temp, dtemp, &
       drhoenthalpy, dw2, df, f, sqrtdet, vlowx, vlowy, vlowz, old_psi

  CCTK_INT, dimension(3) :: sw

!!$ ************************ ATTENTION **************************
!!$ the conformal state is saved here, because the variable
!!$ 'conformal_state' is set to 0 to force the macros to calculate
!!$ the conformal values instead of the physical
  CCTK_INT old_conformal_state

  CCTK_REAL :: IVP_Atmosphere

#include "CactusEinstein/ADMMacros/src/macro/UPPERMET_declare.h"
#include "CactusEinstein/ADMMacros/src/macro/DETG_declare.h"
#include "CactusEinstein/ADMMacros/src/macro/TRK_declare.h"
#include "CactusEinstein/ADMMacros/src/macro/TRRICCI_declare.h"

#include "CactusEOS/EOS_Base/src/EOS_Base.inc"

!!$ save the value of conformal_state
  old_conformal_state = conformal_state;

  if (.not. (Whisky_IVP_Ignore_Matter .eq. 1)) then
!!$  Set up the fluid constants

    IVP_K = EOS_Pressure(whisky_polytrope_handle, 1.d0, 1.d0)
    IVP_Gamma = 1.d0 + IVP_K / &
         EOS_SpecificIntEnergy(whisky_polytrope_handle, 1.d0, 1.d0)

!!$  Work out the atmosphere value

    if (rho_abs_min > 0.d0) then
      IVP_Atmosphere = rho_abs_min
    else
      IVP_Atmosphere = whisky_rho_central * rho_rel_min
    end if
  
    if (initial_rho_abs_min > 0.d0) then
      IVP_Atmosphere = initial_rho_abs_min
    else if (initial_rho_rel_min > 0.d0) then
      IVP_Atmosphere = whisky_rho_central * initial_rho_rel_min
    end if
  
    if (initial_atmosphere_factor > 0.0) then
      IVP_Atmosphere = IVP_Atmosphere * initial_atmosphere_factor
    end if
  end if
  
!!$  Set other constants

#include "CactusEinstein/ADMMacros/src/macro/ADM_Spacing.h"

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

!!$  Set up initial data. Something has already done this, so we do not
!!$  need to.

!!$ Initialise some internal variables since the following loop is not
!!$ going over the whole grid
  IVP_trK = 0.0d0
  IVP_Axx = 0.0d0
  IVP_Axy = 0.0d0
  IVP_Axz = 0.0d0
  IVP_Ayy = 0.0d0
  IVP_Ayz = 0.0d0
  IVP_Azz = 0.0d0
  IVP_Mtildexx = 0.0d0
  IVP_Mtildexy = 0.0d0
  IVP_Mtildexz = 0.0d0
  IVP_Mtildeyy = 0.0d0
  IVP_Mtildeyz = 0.0d0
  IVP_Mtildezz = 0.0d0
  IVP_uxx = 1.0d0
  IVP_uxy = 0.0d0
  IVP_uxz = 0.0d0
  IVP_uyy = 1.0d0
  IVP_uyz = 0.0d0
  IVP_uzz = 1.0d0
  IVP_detg = 1.0d0
  IVP_RicciScalar = 0.0d0
  IVP_rho = IVP_Atmosphere
  IVP_jx = 0.0d0
  IVP_jy = 0.0d0
  IVP_jz = 0.0d0

!!$  Calculate trace of extrinsic curvature and Mtilde
!!$  Set up scalar curvature, inverse metric. 
!!$  We are allowed to fix this.

  do k = 2, nz-1
    do j = 2, ny-1
      do i = 2, nx-1

!!$     We want to have the physical tr K
        conformal_state = old_conformal_state;
#include "CactusEinstein/ADMMacros/src/macro/TRK_guts.h"
        IVP_trK(i,j,k) = TRK_TRK
#include "CactusEinstein/ADMMacros/src/macro/TRK_undefine.h"
!!$     The rest should be computed conformal
        conformal_state = 0;

!!$     set up internal psi
        if (old_conformal_state .ne. 0) then
          old_psi = Psi(i,j,k)
        else
          old_psi = 1.d0
        endif

!!$     Use the A terms as temporary variables for the trace free Kij
!!$     (eq 19): A_ij = k_ij - 1/3 * gamma_ij * k

        IVP_Axx(i,j,k) = kxx(i,j,k) - &
                         IVP_trK(i,j,k) * old_psi**4 * gxx(i,j,k) / 3.d0
        IVP_Axy(i,j,k) = kxy(i,j,k) - &
                         IVP_trK(i,j,k) * old_psi**4 * gxy(i,j,k) / 3.d0
        IVP_Axz(i,j,k) = kxz(i,j,k) - &
                         IVP_trK(i,j,k) * old_psi**4 * gxz(i,j,k) / 3.d0
        IVP_Ayy(i,j,k) = kyy(i,j,k) - &
                         IVP_trK(i,j,k) * old_psi**4 * gyy(i,j,k) / 3.d0
        IVP_Ayz(i,j,k) = kyz(i,j,k) - &
                         IVP_trK(i,j,k) * old_psi**4 * gyz(i,j,k) / 3.d0
        IVP_Azz(i,j,k) = kzz(i,j,k) - &
                         IVP_trK(i,j,k) * old_psi**4 * gzz(i,j,k) / 3.d0

!!$     We have precomputed the uppermetric in the macro for trK
!!$     So lets upper the indices of A
!!$     Mtilde^ij = Atilde^ij = gammatilde^ik * gammatilde^jl * A_kl
!!$                             Psi^2
#include "CactusEinstein/ADMMacros/src/macro/UPPERMET_guts.h"

        IVP_Mtildexx(i,j,k) = old_psi**2 * (&
             UPPERMET_UXX**2 * IVP_Axx(i,j,k) + &
             2.d0 * UPPERMET_UXX * UPPERMET_UXY * IVP_Axy(i,j,k) + &
             2.d0 * UPPERMET_UXX * UPPERMET_UXZ * IVP_Axz(i,j,k) + &
             UPPERMET_UXY**2 * IVP_Ayy(i,j,k) + &
             2.d0 * UPPERMET_UXY * UPPERMET_UXZ * IVP_Ayz(i,j,k) + &
             UPPERMET_UXZ**2 * IVP_Azz(i,j,k) )

        IVP_Mtildexy(i,j,k) = old_psi**2 * (&
             UPPERMET_UXX * UPPERMET_UXY * IVP_Axx(i,j,k) + &
             UPPERMET_UXY**2 * IVP_Axy(i,j,k) + &
             UPPERMET_UXX * UPPERMET_UYY * IVP_Axy(i,j,k) + &
             UPPERMET_UXY * UPPERMET_UXZ * IVP_Axz(i,j,k) + &
             UPPERMET_UXX * UPPERMET_UYZ * IVP_Axz(i,j,k) + &
             UPPERMET_UXY * UPPERMET_UYY * IVP_Ayy(i,j,k) + &
             UPPERMET_UXZ * UPPERMET_UYY * IVP_Ayz(i,j,k) + &
             UPPERMET_UXY * UPPERMET_UYZ * IVP_Ayz(i,j,k) + &
             UPPERMET_UXZ * UPPERMET_UYZ * IVP_Azz(i,j,k))

        IVP_Mtildexz(i,j,k) = old_psi**2 * (&
             UPPERMET_UXX * UPPERMET_UXZ * IVP_Axx(i,j,k) + &
             UPPERMET_UXY * UPPERMET_UXZ * IVP_Axy(i,j,k) + &
             UPPERMET_UXX * UPPERMET_UYZ * IVP_Axy(i,j,k) + &
             UPPERMET_UXZ**2 * kxz(i,j,k) + &
             UPPERMET_UXX * UPPERMET_UZZ * IVP_Axz(i,j,k) + &
             UPPERMET_UXY * UPPERMET_UYZ * IVP_Ayy(i,j,k) + &
             UPPERMET_UXZ * UPPERMET_UYZ * IVP_Ayz(i,j,k) + &
             UPPERMET_UXY * UPPERMET_UZZ * IVP_Ayz(i,j,k) + &
             UPPERMET_UXZ * UPPERMET_UZZ * IVP_Azz(i,j,k))

        IVP_Mtildeyy(i,j,k) = old_psi**2 * (&
             UPPERMET_UXY**2 * IVP_Axx(i,j,k) + &
             2.d0 * UPPERMET_UXY * UPPERMET_UYY * IVP_Axy(i,j,k) + &
             2.d0 * UPPERMET_UXY * UPPERMET_UYZ * IVP_Axz(i,j,k) + &
             UPPERMET_UYY**2 * IVP_Ayy(i,j,k) + &
             2.d0 * UPPERMET_UYY * UPPERMET_UYZ * IVP_Ayz(i,j,k) + &
             UPPERMET_UYZ**2 * IVP_Azz(i,j,k))

        IVP_Mtildeyz(i,j,k) = old_psi**2 * (&
             UPPERMET_UXY * UPPERMET_UXZ * IVP_Axx(i,j,k) + &
             UPPERMET_UXZ * UPPERMET_UYY * IVP_Axy(i,j,k) + &
             UPPERMET_UXY * UPPERMET_UYZ * IVP_Axy(i,j,k) + &
             UPPERMET_UXZ * UPPERMET_UYZ * IVP_Axz(i,j,k) + &
             UPPERMET_UXY * UPPERMET_UZZ * IVP_Axz(i,j,k) + &
             UPPERMET_UYY * UPPERMET_UYZ * IVP_Ayy(i,j,k) + &
             UPPERMET_UYZ**2 * IVP_Ayz(i,j,k) + &
             UPPERMET_UYY * UPPERMET_UZZ * IVP_Ayz(i,j,k) + &
             UPPERMET_UYZ * UPPERMET_UZZ * IVP_Azz(i,j,k))

        IVP_Mtildezz(i,j,k) = old_psi**2 * (&
             UPPERMET_UXZ**2 * IVP_Axx(i,j,k) + &
             2.d0 * UPPERMET_UXZ * UPPERMET_UYZ * IVP_Axy(i,j,k) + &
             2.d0 * UPPERMET_UXZ * UPPERMET_UZZ * IVP_Axz(i,j,k) + &
             UPPERMET_UYZ**2 * IVP_Ayy(i,j,k) + &
             2.d0 * UPPERMET_UYZ * UPPERMET_UZZ * IVP_Ayz(i,j,k) + &
             UPPERMET_UZZ**2*IVP_Azz(i,j,k))

!!$     Dump the upper conformal metric and determinant into grid functions
!!$     These have been required for earlier calculations so will be around

        IVP_uxx(i,j,k) = UPPERMET_UXX
        IVP_uxy(i,j,k) = UPPERMET_UXY
        IVP_uxz(i,j,k) = UPPERMET_UXZ
        IVP_uyy(i,j,k) = UPPERMET_UYY
        IVP_uyz(i,j,k) = UPPERMET_UYZ
        IVP_uzz(i,j,k) = UPPERMET_UZZ
#include "CactusEinstein/ADMMacros/src/macro/UPPERMET_undefine.h"

#include "CactusEinstein/ADMMacros/src/macro/DETG_guts.h"
        IVP_detg(i,j,k) = DETG_DETG

#include "CactusEinstein/ADMMacros/src/macro/TRRICCI_guts.h"

!!$    The conformal RicciScalar is required inside the solve, but does
!!$    not change, so only compute this once - requires another GF
       IVP_RicciScalar(i,j,k) = TRRICCI_TRRICCI
!!$        IVP_RicciScalar(i,j,k) = TRRICCI_TRRICCI / IVP_uxx(i,j,k)

!!$     Set up the ADM matter density and momentum current. 
!!$     This only contains the Whisky terms
!!$     (and the conformal factor to get the conformal values, which
!!$      are fixed)

        if (.not. (Whisky_IVP_Ignore_Matter .eq. 1)) then
          IVP_rho(i,j,k)= old_psi**8 * &
                          (( rho(i,j,k) * (1.d0 + eps(i,j,k)) + &
                             press(i,j,k) ) * w_lorentz(i,j,k)**2 - &
                           press(i,j,k))
          IVP_jx(i,j,k) = old_psi**10 * &
                          (( rho(i,j,k) * (1.d0 + eps(i,j,k)) + &
                            press(i,j,k) ) * w_lorentz(i,j,k)**2 * &
                          velx(i,j,k))
          IVP_jy(i,j,k) = old_psi**10 * &
                          (( rho(i,j,k) * (1.d0 + eps(i,j,k)) + &
                            press(i,j,k) ) * w_lorentz(i,j,k)**2 * &
                          vely(i,j,k))
          IVP_jz(i,j,k) = old_psi**10 * &
                          (( rho(i,j,k) * (1.d0 + eps(i,j,k)) + &
                            press(i,j,k) ) * w_lorentz(i,j,k)**2 * &
                          velz(i,j,k))
        end if
      end do
    end do
  end do

#include "CactusEinstein/ADMMacros/src/macro/TRRICCI_undefine.h"

!!$  So, before solving we have to calculate the following only once.
!!$  See Cook, eqn 32.
!!$  trK
!!$  Mtilde(ij)
!!$  j(i)
!!$  rho_ADM
!!$  Rtilde
  
!!$  These should be set up by now. From these we can calculate the
!!$  terms in the power series for the psi equation,

!!$  ---2
!!$  \ / psi - psi Rtilde / 8 - psi^5 trK^2 / 12 + 
!!$   V
!!$
!!$          psi^-7 Atilde^(ij) Atilde_(ij) / 8 = -2 pi G psi^5 rho_ADM

!!$  As we are using the BAM_Power call, the only term which is not
!!$  fixed is that containing Atilde. The rest can be calculated here
!!$  for later use.
 
!!$  Set initial guess

  if (old_conformal_state .ne. 0) then
    IVP_Psi = Psi
  else
    IVP_Psi = 1.d0
  endif
  IVP_Vx = 0.d0
  IVP_Vy = 0.d0
  IVP_Vz = 0.d0

  IVP_Change_Psi = 1.d10
  Max_Change_Psi = 1.d10

!!$ just to be sure, should be 0 anyway
  conformal_state = 0;

  call CCTK_ReductionHandle(Reduction_Handle, "maximum")
  call CCTK_VarIndex(ChangePsiIndex, "Whisky_IVP::IVP_Change_Psi")

  if (IVP_Debug_Output .eq. 1) then
    call Call_IVP_Debug_Output(cctkGH, 0)
  endif

  do while(Max_Change_Psi > Whisky_IVP_Tolerance)

!!$ Single step solve, working out the change in Psi

    IVP_Change_Psi = IVP_Psi
    call Solve_IVP_Step(CCTK_PASS_FTOF)
    IVP_Change_Psi = abs(IVP_Psi - IVP_Change_Psi)

!!$ Find out the maximum global change in Psi

    call CCTK_Reduce(ierr, cctkGH, -1, Reduction_Handle,&
         1, CCTK_VARIABLE_REAL, Max_Change_Psi, 1, ChangePsiIndex)

    if (ierr .ne. 0) then
      call CCTK_WARN(0, "Failed to reduce the change in Psi")
    end if

  end do
  
  if (IVP_Debug_Output .eq. 1) then
    call Call_IVP_Debug_Output(cctkGH, 1)
  endif


!!$  Use the new conformal metric and extrinsic curvature to recover
!!$  the spacetime variables

!!$  First bring the metric to the physical form

  gxx = IVP_Psi**4 * gxx
  gxy = IVP_Psi**4 * gxy
  gxz = IVP_Psi**4 * gxz
  gyy = IVP_Psi**4 * gyy
  gyz = IVP_Psi**4 * gyz
  gzz = IVP_Psi**4 * gzz

!!$  Do the extrinsic curvature

  do k = 2, nz - 1
    do j = 2, ny - 1
      do i = 2, nx - 1

!!$     Start with Atilde
!!$     Atilde^ij = (Ltilde V)^ij + Mtilde^ij

        IVP_Axx(i,j,k) = IVP_Mtildexx(i,j,k)
        IVP_Axy(i,j,k) = IVP_Mtildexy(i,j,k)
        IVP_Axz(i,j,k) = IVP_Mtildexz(i,j,k)
        IVP_Ayy(i,j,k) = IVP_Mtildeyy(i,j,k)
        IVP_Ayz(i,j,k) = IVP_Mtildeyz(i,j,k)
        IVP_Azz(i,j,k) = IVP_Mtildezz(i,j,k)
        
!!$     Add in the longitudinal derivatives of V.
!!$                      ---i  j   ---j i   2      ij ---   l
!!$     (Ltilde V)^ij == \ /  V  + \ / V  - - gamma   \ /  V
!!$                       V         V       3          V l
!!$
!!$     Divergence term first
!!$
!!$     ---   l    l    +--+l   m   /          i \     /
!!$     \ /  V  = V   + |      V  = | sqrt(g) V  |    / sqrt(g)
!!$      V l       ,l   |   lm      \            /,i /
  
        DivergenceV = ( &
           ( sqrt(IVP_detg(i+1,j,k)) * IVP_Vx(i+1,j,k) - &
             sqrt(IVP_detg(i-1,j,k)) * IVP_Vx(i-1,j,k) ) * i2dx + &
           ( sqrt(IVP_detg(i,j+1,k)) * IVP_Vy(i,j+1,k) - &
             sqrt(IVP_detg(i,j-1,k)) * IVP_Vy(i,j-1,k) ) * i2dy + &
           ( sqrt(IVP_detg(i,j,k+1)) * IVP_Vz(i,j,k+1) - &
             sqrt(IVP_detg(i,j,k-1)) * IVP_Vz(i,j,k-1) ) * i2dz ) / &
             sqrt(IVP_detg(i,j,k))

        IVP_Axx(i,j,k) = IVP_Axx(i,j,k) - &
             2.d0 / 3.d0 * DivergenceV * IVP_uxx(i,j,k)
        IVP_Axy(i,j,k) = IVP_Axy(i,j,k) - &
             2.d0 / 3.d0 * DivergenceV * IVP_uxy(i,j,k)
        IVP_Axz(i,j,k) = IVP_Axz(i,j,k) - &
             2.d0 / 3.d0 * DivergenceV * IVP_uxz(i,j,k)
        IVP_Ayy(i,j,k) = IVP_Ayy(i,j,k) - &
             2.d0 / 3.d0 * DivergenceV * IVP_uyy(i,j,k)
        IVP_Ayz(i,j,k) = IVP_Ayz(i,j,k) - &
             2.d0 / 3.d0 * DivergenceV * IVP_uyz(i,j,k)
        IVP_Azz(i,j,k) = IVP_Azz(i,j,k) - &
             2.d0 / 3.d0 * DivergenceV * IVP_uzz(i,j,k)

!!$        Then covariant derivatives
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

!!$     Finally, reset the extrinsic curvature
!!$     It is simplest to reconstruct K with indices up 
!!$     in the Mtilde array as temporary storage, and then
!!$     covert to the lower indices

        IVP_Mtildexx(i,j,k) = IVP_Psi(i,j,k)**(-10) * IVP_Axx(i,j,k) + &
             IVP_Psi(i,j,k)**(-4) * IVP_uxx(i,j,k) * IVP_trK(i,j,k) / 3.d0
        IVP_Mtildexy(i,j,k) = IVP_Psi(i,j,k)**(-10) * IVP_Axy(i,j,k) + &
             IVP_Psi(i,j,k)**(-4) * IVP_uxy(i,j,k) * IVP_trK(i,j,k) / 3.d0
        IVP_Mtildexz(i,j,k) = IVP_Psi(i,j,k)**(-10) * IVP_Axz(i,j,k) + &
             IVP_Psi(i,j,k)**(-4) * IVP_uxz(i,j,k) * IVP_trK(i,j,k) / 3.d0
        IVP_Mtildeyy(i,j,k) = IVP_Psi(i,j,k)**(-10) * IVP_Ayy(i,j,k) + &
             IVP_Psi(i,j,k)**(-4) * IVP_uyy(i,j,k) * IVP_trK(i,j,k) / 3.d0
        IVP_Mtildeyz(i,j,k) = IVP_Psi(i,j,k)**(-10) * IVP_Ayz(i,j,k) + &
             IVP_Psi(i,j,k)**(-4) * IVP_uyz(i,j,k) * IVP_trK(i,j,k) / 3.d0
        IVP_Mtildezz(i,j,k) = IVP_Psi(i,j,k)**(-10) * IVP_Azz(i,j,k) + &
             IVP_Psi(i,j,k)**(-4) * IVP_uzz(i,j,k) * IVP_trK(i,j,k) / 3.d0
  
        kxx(i,j,k) = &
             gxx(i,j,k)**2 * IVP_Mtildexx(i,j,k) + &
             2.d0 * gxx(i,j,k) * gxy(i,j,k) * IVP_Mtildexy(i,j,k) + &
             2.d0 * gxx(i,j,k) * gxz(i,j,k) * IVP_Mtildexz(i,j,k) + &
             gxy(i,j,k)**2 * IVP_Mtildeyy(i,j,k) + &
             2.d0 * gxy(i,j,k) * gxz(i,j,k) * IVP_Mtildeyz(i,j,k) + &
             gxz(i,j,k)**2 * IVP_Mtildezz(i,j,k) 

        kxy(i,j,k) = &
             gxx(i,j,k) * gxy(i,j,k) * IVP_Mtildexx(i,j,k) + &
             gxy(i,j,k)**2 * IVP_Mtildexy(i,j,k) + &
             gxx(i,j,k) * gyy(i,j,k) * IVP_Mtildexy(i,j,k) + &
             gxy(i,j,k) * gxz(i,j,k) * IVP_Mtildexz(i,j,k) + &
             gxx(i,j,k) * gyz(i,j,k) * IVP_Mtildexz(i,j,k) + &
             gxy(i,j,k) * gyy(i,j,k) * IVP_Mtildeyy(i,j,k) + &
             gxz(i,j,k) * gyy(i,j,k) * IVP_Mtildeyz(i,j,k) + &
             gxy(i,j,k) * gyz(i,j,k) * IVP_Mtildeyz(i,j,k) + &
             gxz(i,j,k) * gyz(i,j,k) * IVP_Mtildezz(i,j,k)

        kxz(i,j,k) = &
             gxx(i,j,k) * gxz(i,j,k) * IVP_Mtildexx(i,j,k) + &
             gxy(i,j,k) * gxz(i,j,k) * IVP_Mtildexy(i,j,k) + &
             gxx(i,j,k) * gyz(i,j,k) * IVP_Mtildexy(i,j,k) + &
             gxz(i,j,k)**2 * IVP_Mtildexz(i,j,k) + &
             gxx(i,j,k) * gzz(i,j,k) * IVP_Mtildexz(i,j,k) + &
             gxy(i,j,k) * gyz(i,j,k) * IVP_Mtildeyy(i,j,k) + &
             gxz(i,j,k) * gyz(i,j,k) * IVP_Mtildeyz(i,j,k) + &
             gxy(i,j,k) * gzz(i,j,k) * IVP_Mtildeyz(i,j,k) + &
             gxz(i,j,k) * gzz(i,j,k) * IVP_Mtildezz(i,j,k)

        kyy(i,j,k) = &
             gxy(i,j,k)**2 * IVP_Mtildexx(i,j,k) + &
             2.d0 * gxy(i,j,k) * gyy(i,j,k) * IVP_Mtildexy(i,j,k) + &
             2.d0 * gxy(i,j,k) * gyz(i,j,k) * IVP_Mtildexz(i,j,k) + &
             gyy(i,j,k)**2 * IVP_Mtildeyy(i,j,k) + &
             2.d0 * gyy(i,j,k) * gyz(i,j,k) * IVP_Mtildeyz(i,j,k) + &
             gyz(i,j,k)**2 * IVP_Mtildezz(i,j,k)

        kyz(i,j,k) = &
             gxy(i,j,k) * gxz(i,j,k) * IVP_Mtildexx(i,j,k) + &
             gxz(i,j,k) * gyy(i,j,k) * IVP_Mtildexy(i,j,k) + &
             gxy(i,j,k) * gyz(i,j,k) * IVP_Mtildexy(i,j,k) + &
             gxz(i,j,k) * gyz(i,j,k) * IVP_Mtildexz(i,j,k) + &
             gxy(i,j,k) * gzz(i,j,k) * IVP_Mtildexz(i,j,k) + &
             gyy(i,j,k) * gyz(i,j,k) * IVP_Mtildeyy(i,j,k) + &
             gyz(i,j,k)**2 * IVP_Mtildeyz(i,j,k) + &
             gyy(i,j,k) * gzz(i,j,k) * IVP_Mtildeyz(i,j,k) + &
             gyz(i,j,k) * gzz(i,j,k) * IVP_Mtildezz(i,j,k)

        kzz(i,j,k) = &
             gxz(i,j,k)**2 * IVP_Mtildexx(i,j,k) + &
             2.d0 * gxz(i,j,k) * gyz(i,j,k) * IVP_Mtildexy(i,j,k) + &
             2.d0 * gxz(i,j,k) * gzz(i,j,k) * IVP_Mtildexz(i,j,k) + &
             gyz(i,j,k)**2 * IVP_Mtildeyy(i,j,k) + &
             2.d0 * gyz(i,j,k) * gzz(i,j,k) * IVP_Mtildeyz(i,j,k) + &
             gzz(i,j,k)**2*IVP_Mtildezz(i,j,k)

      end do
    end do
  end do

!!$  bring the metric to the conformal form, if wanted

  if (old_conformal_state .ne. 0) then
    gxx = Psi**(-4) * gxx
    gxy = Psi**(-4) * gxy
    gxz = Psi**(-4) * gxz
    gyy = Psi**(-4) * gyy
    gyz = Psi**(-4) * gyz
    gzz = Psi**(-4) * gzz
!!$ below this, conformal_state is again what it was before
    conformal_state = old_conformal_state;
  endif

!!$  Boundary conditions for the spacetime
!!$  Symmetry boundaries first  

  call CartSymGN(ierr,cctkGH,"ADMBase::lapse")
  call CartSymGN(ierr,cctkGH,"ADMBase::metric")
  call CartSymGN(ierr,cctkGH,"ADMBase::curv")
  if (shift_state .ne. 0) then
    call CartSymGN(ierr,cctkGH,"ADMBase::shift")
  end if

!!$  What should we do for the outer boundaries (if anything?)

!!$  call BndFlatGN(ierr,cctkGH,sw,"ADMBase::lapse")
!!$  call BndFlatGN(ierr,cctkGH,sw,"ADMBase::metric")
!!$  call BndFlatGN(ierr,cctkGH,sw,"ADMBase::curv")
!!$  if (shift_state .ne. 0) then
!!$    call BndFlatGN(ierr,cctkGH,sw,"ADMBase::shift")
!!$  end if

!!$  Then reconstruct the fluid variables

!!$  As it was really the scaled matter sources that remained
!!$  fixed, we first have to rescale the matter terms

  if (.not. (Whisky_IVP_Ignore_Matter .eq. 1)) then
    IVP_rho = IVP_rho * IVP_Psi**(-8)
    IVP_jx = IVP_jx * IVP_Psi**(-10)
    IVP_jy = IVP_jy * IVP_Psi**(-10)
    IVP_jz = IVP_jz * IVP_Psi**(-10)

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
  
!!$       set up internal psi
          if (old_conformal_state .ne. 0) then
              old_psi = Psi(i,j,k)
          else
            old_psi = 1.d0
          endif

          j2 = old_psi**4 * gxx(i,j,k) * IVP_jx(i,j,k) * IVP_jx(i,j,k) + &
               old_psi**4 * gyy(i,j,k) * IVP_jy(i,j,k) * IVP_jy(i,j,k) + &
               old_psi**4 * gzz(i,j,k) * IVP_jz(i,j,k) * IVP_jz(i,j,k) + &
               2.d0 * old_psi**4 * gxy(i,j,k) * IVP_jx(i,j,k) * IVP_jy(i,j,k) +&
               2.d0 * old_psi**4 * gxz(i,j,k) * IVP_jx(i,j,k) * IVP_jz(i,j,k) +&
               2.d0 * old_psi**4 * gyz(i,j,k) * IVP_jy(i,j,k) * IVP_jz(i,j,k)

          rhoold = rho(i,j,k)
          rhonew = rhoold
          pressnew = IVP_K * rhonew**IVP_Gamma
          epsnew = IVP_K * rhonew**(IVP_Gamma - 1.d0) / (IVP_Gamma - 1.d0)

          w_lorentznew = 1.d0 / sqrt(1.d0 - j2 / &
               (IVP_rho(i,j,k) + pressnew)**2)
          rhoenthalpy = rhonew * (1.d0 + epsnew) + pressnew

          f = rhoenthalpy * w_lorentznew - pressnew - IVP_rho(i,j,k)

          count = 0

          do while ( ((abs(rhonew - rhoold)/abs(rhonew) &
               .gt. IVP_rel_c2p_tol) .and. &
               (abs(rhonew - rhoold) .gt. IVP_abs_c2p_tol)) .or. &
               (count .lt. IVP_c2p_countmin) )
            count = count + 1
          
            if (count .gt. 100) then
              call CCTK_WARN(0, "Failed to converged IVP C2P")
            end if

            temp = j2 / (IVP_rho(i,j,k) + pressnew)**2
            dtemp = -2.d0 * j2 * IVP_K * IVP_Gamma * &
                 rhonew**(IVP_Gamma - 1.d0) / &
                 (IVP_rho(i,j,k) + pressnew)**3

            drhoenthalpy = 1.d0 + IVP_K * IVP_Gamma**2 * &
                 rhonew**(IVP_Gamma - 1.d0) / (IVP_Gamma - 1.d0)
!!$            dw2 = w_lorentznew**2 * dtemp / (1.d0 - temp)
            dw2 = dtemp / (1.d0 - temp)**2

            df = drhoenthalpy * w_lorentznew**2 + dw2 * rhoenthalpy

            rhoold = rhonew
            rhonew = rhoold - f / df
 
            pressnew = IVP_K * rhonew**IVP_Gamma
!!$            epsnew = IVP_K * pressnew * rhonew**(IVP_Gamma - 1.d0) &
            epsnew = IVP_K * rhonew**(IVP_Gamma - 1.d0) / &
                           (IVP_Gamma - 1.d0)
            rhoenthalpy = rhonew * (1.d0 + epsnew) + pressnew
            w_lorentznew = 1.d0 / sqrt(1.d0 - j2 / &
                                  (IVP_rho(i,j,k) + pressnew)**2)
          
!!$            f = rhoenthalpy * w_lorentznew**2 - IVP_rho(i,j,k)
            f = rhoenthalpy * w_lorentznew**2 - pressnew - IVP_rho(i,j,k)

          end do

          rho(i,j,k) = rhonew
          press(i,j,k) = pressnew
          eps(i,j,k) = epsnew
          velx(i,j,k) = IVP_jx(i,j,k) / rhoenthalpy / w_lorentznew**2
          vely(i,j,k) = IVP_jy(i,j,k) / rhoenthalpy / w_lorentznew**2
          velz(i,j,k) = IVP_jz(i,j,k) / rhoenthalpy / w_lorentznew**2
          w_lorentz(i,j,k) = w_lorentznew

!!$          Check that we are not below atmosphere

          if (rho(i,j,k) .le. IVP_Atmosphere) then

            rho(i,j,k) = IVP_Atmosphere
            press(i,j,k) = IVP_K * IVP_Atmosphere**IVP_Gamma
            eps(i,j,k) = IVP_K * IVP_Atmosphere**(IVP_Gamma - 1.d0) &
                               / (IVP_Gamma - 1.d0)
            velx(i,j,k) = 0.d0
            vely(i,j,k) = 0.d0
            velz(i,j,k) = 0.d0
            w_lorentz(i,j,k) = 1.d0

          end if

#ifdef DEBUG_WHISKY_IVP

          if ( (w_lorentz(i,j,k) - 1.d0 / sqrt(1.d0 - &
               (old_psi**4 * gxx(i,j,k) * velx(i,j,k) * velx(i,j,k) + &
                old_psi**4 * gyy(i,j,k) * vely(i,j,k) * vely(i,j,k) + &
                old_psi**4 * gzz(i,j,k) * velz(i,j,k) * velz(i,j,k) + &
                old_psi**4 * gxy(i,j,k) * velx(i,j,k) * vely(i,j,k) + &
                old_psi**4 * gxz(i,j,k) * velx(i,j,k) * velz(i,j,k) + &
                old_psi**4 * gyz(i,j,k) * vely(i,j,k) * velz(i,j,k))) ) / &
               w_lorentznew .gt. 1.d-10 ) then
            write(*,*) "We have not solved the C2P problem correctly.", &
                 w_lorentznew, 1.d0 / sqrt(1.d0 - &
               (old_psi**4 * gxx(i,j,k) * velx(i,j,k) * velx(i,j,k) + &
                old_psi**4 * gyy(i,j,k) * vely(i,j,k) * vely(i,j,k) + &
                old_psi**4 * gzz(i,j,k) * velz(i,j,k) * velz(i,j,k) + &
                old_psi**4 * gxy(i,j,k) * velx(i,j,k) * vely(i,j,k) + &
                old_psi**4 * gxz(i,j,k) * velx(i,j,k) * velz(i,j,k) + &
                old_psi**4 * gyz(i,j,k) * vely(i,j,k) * velz(i,j,k)))
            call CCTK_WARN(0, "Failure to solve IVP C2P")
          end if

#endif

#include "CactusEinstein/ADMMacros/src/macro/DETG_guts.h"

          sqrtdet = sqrt(DETG_DETG)

          vlowx = old_psi**4 * gxx(i,j,k) * velx(i,j,k) + &
                  old_psi**4 * gxy(i,j,k) * vely(i,j,k) + &
                  old_psi**4 * gxz(i,j,k) * velz(i,j,k)
          vlowy = old_psi**4 * gxy(i,j,k) * velx(i,j,k) + &
                  old_psi**4 * gyy(i,j,k) * vely(i,j,k) + &
                  old_psi**4 * gyz(i,j,k) * velz(i,j,k)
          vlowz = old_psi**4 * gxz(i,j,k) * velx(i,j,k) + &
                  old_psi**4 * gyz(i,j,k) * vely(i,j,k) + &
                  old_psi**4 * gzz(i,j,k) * velz(i,j,k)

          dens(i,j,k) = sqrtdet *  rho(i,j,k) * w_lorentz(i,j,k)
          sx(i,j,k)   = sqrtdet * (rho(i,j,k) * (1.d0 + eps(i,j,k)) + &
                          press(i,j,k)) * w_lorentz(i,j,k)**2 * vlowx
          sy(i,j,k)   = sqrtdet * (rho(i,j,k) * (1.d0 + eps(i,j,k)) + &
                          press(i,j,k)) * w_lorentz(i,j,k)**2 * vlowy
          sz(i,j,k)   = sqrtdet * (rho(i,j,k) * (1.d0 + eps(i,j,k)) + &
                          press(i,j,k)) * w_lorentz(i,j,k)**2 * vlowz
          !tau(i,j,k)  = sqrtdet * ( (rho(i,j,k) * (1.d0 + eps(i,j,k)) + &
          !     press(i,j,k)) * w_lorentz(i,j,k)**2 - press(i,j,k) ) - &
          !     dens(i,j,k)
          tau(i,j,k) = sqrtdetg * & 
                 ( rho(i,j,k) * (w_lorentz(i,j,k)-1.d0) * &
                   w_lorentz(i,j,k) + &
                   rho(i,j,k)*eps(i,j,k)*w_lorentz(i,j,k)**2 + & 
                   press(i,j,k)*(w_lorentz(i,j,k)**2-1.d0) )

        end do
      end do
    end do
  end if

#include "CactusEinstein/ADMMacros/src/macro/DETG_undefine.h"  

  sw = whisky_stencil

!!$  Symmetry boundaries first  

  call CartSymGN(ierr,cctkGH,"Whisky::Whisky_cons_scalar_var")
  call CartSymGN(ierr,cctkGH,"Whisky::Whisky_prim_scalar_var")
  call CartSymGN(ierr,cctkGH,"Whisky::Whisky_cons_vector_var")
  call CartSymGN(ierr,cctkGH,"Whisky::Whisky_prim_vector_var")

!!$  Set flat boundaries; what else could we do?

  call BndFlatGN(ierr,cctkGH,sw,"Whisky::Whisky_cons_scalar_var")
  call BndFlatGN(ierr,cctkGH,sw,"Whisky::Whisky_prim_scalar_var")
  call BndFlatGN(ierr,cctkGH,sw,"Whisky::Whisky_cons_vector_var")
  call BndFlatGN(ierr,cctkGH,sw,"Whisky::Whisky_prim_vector_var")

end subroutine Whisky_IVP
