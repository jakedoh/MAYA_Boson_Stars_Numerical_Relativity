 /*@@
   @file      Whisky_SetTmunu.F90
   @date      Wed Aug 19 18:53 2009
   @author    Tanja Bode
   @desc 
     Based on the include files, this function uses the 
     AddToTmunu mechanism to set the stress energy tensor.

     The calculation of the stress energy tensor.
     The version used here was worked out by Miguel Alcubierre. I
     think it was an extension of the routine from GR3D, written
     by Mark Miller.
     C version added by Ian Hawke.

     Lower components of the stress-energy tensor obtained from
     the hydro variables.  The components are given by:

     T      =  rho h  u   u    +  P  g
      mu nu            mu  nu         mu nu

     where rho is the energy density of the fluid, h the enthalpy
     and P the pressure.  The enthalpy is given in terms of the
     basic variables as:

     h  =  1  +  e  +  P/rho

     with e the internal energy (eps here).

     In the expresion for T_{mu,nu} we also have the four-velocity
     of the fluid given by (v_i is the 3-velocity field):

                                 i
     u  =  W ( - alpha +  v  beta  )
      0                    i

     u  =  W v
      i       i
                                                i  -1/2
     with W the Lorentz factor:   W = ( 1 -  v v  )
                                              i

     and where alpha and beta are the lapse and shift vector.

     Finally, the 4 metric is given by

                   2             i
     g   =  - alpha  + beta  beta
      00                   i

     g   =  beta
      0i        i


     g   =  gamma      (the spatial metric)
      ij        ij


     This has now been extended to add a magnetic field component 
     for MHD, taking our notation from Duez et al 2005

     T        += b^2 ( u  u  + 1/2 g     ) - b   b
      mu nu             mu nu       mu nu     mu  nu

     where b^{mu} = B^{mu}_(u) 
     ( the factor of 1/sqrt(4 pi) is absorbed into the B-field) 
     and

            mu         * nu mu
           B    = u   F        
            (u)    nu



   @enddesc 
 @@*/

#include "SpaceMask.h"
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#define Bvx(i,j,kk) Bvec(i,j,kk,1)
#define Bvy(i,j,kk) Bvec(i,j,kk,2)
#define Bvz(i,j,kk) Bvec(i,j,kk,3)

subroutine Whisky_SetTmunu (CCTK_ARGUMENTS)

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  CCTK_REAL vlow_x, vlow_y, vlow_z
  CCTK_REAL betaxlow, betaylow, betazlow, beta2
  CCTK_REAL ulow_t, ulow_x, ulow_y, ulow_z
  CCTK_REAL rhoenthalpy
  CCTK_INT :: type_bits, atmosphere

  integer i,j,kk

!     if (SpaceMask_CheckStateBitsF90(space_mask, i,j,k, \
!                                     atmosphere_field_descriptor, \
!                                     atmosphere_normal_descriptor)) then
  call SpaceMask_GetTypeBits(type_bits, "Hydro_Atmosphere")
  call SpaceMask_GetStateBits(atmosphere, "Hydro_Atmosphere", "in_atmosphere")

  !$OMP PARALLEL DO PRIVATE (vlow_x, vlow_y, vlow_z, betaxlow, betaylow, &
  !$OMP betazlow, beta2, ulow_t, ulow_x, ulow_y, ulow_z, rhoenthalpy, i,j,kk)
  do kk = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)

           vlow_x = (gxx(i,j,kk)*velx(i,j,kk) + gxy(i,j,kk)*vely(i,j,kk) &
             + gxz(i,j,kk)*velz(i,j,kk))
           vlow_y = (gxy(i,j,kk)*velx(i,j,kk) + gyy(i,j,kk)*vely(i,j,kk) &
             + gyz(i,j,kk)*velz(i,j,kk))
           vlow_z = (gxz(i,j,kk)*velx(i,j,kk) + gyz(i,j,kk)*vely(i,j,kk) &
             + gzz(i,j,kk)*velz(i,j,kk))

!          Calculate lower components and square of shift vector.

           betaxlow = (gxx(i,j,kk)*betax(i,j,kk) + gxy(i,j,kk)*betay(i,j,kk) &
                 + gxz(i,j,kk)*betaz(i,j,kk))
           betaylow = (gxy(i,j,kk)*betax(i,j,kk) + gyy(i,j,kk)*betay(i,j,kk) &
                 + gyz(i,j,kk)*betaz(i,j,kk))
           betazlow = (gxz(i,j,kk)*betax(i,j,kk) + gyz(i,j,kk)*betay(i,j,kk) &
                 + gzz(i,j,kk)*betaz(i,j,kk))

           beta2 = betax(i,j,kk)*betaxlow + betay(i,j,kk)*betaylow &
                 + betaz(i,j,kk)*betazlow 

!          Calculate the specific relativistic enthalpy times rho times the
!          square of the lorentz factor.

           rhoenthalpy = w_lorentz(i,j,kk)**2*&
                        (rho(i,j,kk)*(1.0D0 + eps(i,j,kk)) + press(i,j,kk))

!          Calculate lower components of 4-velocity (without the Lorent factor).

           ulow_t = (-alp(i,j,kk) + velx(i,j,kk)*betaxlow + vely(i,j,kk)*betaylow &
                   + velz(i,j,kk)*betazlow)

           ulow_x = vlow_x
           ulow_y = vlow_y
           ulow_z = vlow_z

!          Calculate Tmunu (the lower components!).

           if ( (whisky_atmosphere_is_vacuum .eq. 0) .or. &
              .not.( SpaceMask_CheckStateBitsF90(space_mask, i, j,kk, type_bits, atmosphere) .or. &
              (atmosphere_mask(i,j,kk) .ne. 0) )) then

              eTtt(i,j,kk) = eTtt(i,j,kk) + rhoenthalpy*ulow_t**2 &
                            + press(i,j,kk)*(beta2 - alp(i,j,kk)**2)
              eTtx(i,j,kk) = eTtx(i,j,kk) + rhoenthalpy*ulow_t*ulow_x &
                            + press(i,j,kk)*betaxlow
              eTty(i,j,kk) = eTty(i,j,kk) + rhoenthalpy*ulow_t*ulow_y &
                            + press(i,j,kk)*betaylow
              eTtz(i,j,kk) = eTtz(i,j,kk) + rhoenthalpy*ulow_t*ulow_z &
                            + press(i,j,kk)*betazlow

              eTxx(i,j,kk) = eTxx(i,j,kk) + rhoenthalpy*ulow_x**2 &
                            + press(i,j,kk)*gxx(i,j,kk)
              eTyy(i,j,kk) = eTyy(i,j,kk) + rhoenthalpy*ulow_y**2 &
                            + press(i,j,kk)*gyy(i,j,kk)
              eTzz(i,j,kk) = eTzz(i,j,kk) + rhoenthalpy*ulow_z**2 &
                            + press(i,j,kk)*gzz(i,j,kk)

              eTxy(i,j,kk) = eTxy(i,j,kk) + rhoenthalpy*ulow_x*ulow_y &
                            + press(i,j,kk)*gxy(i,j,kk)
              eTxz(i,j,kk) = eTxz(i,j,kk) + rhoenthalpy*ulow_x*ulow_z &
                            + press(i,j,kk)*gxz(i,j,kk)
              eTyz(i,j,kk) = eTyz(i,j,kk) + rhoenthalpy*ulow_y*ulow_z &
                            + press(i,j,kk)*gyz(i,j,kk)

          end if

        end do
     end do
  end do
  !$OMP END PARALLEL DO

end subroutine Whisky_SetTmunu

subroutine Whisky_MHD_SetTmunu (CCTK_ARGUMENTS)

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  CCTK_REAL vlow_x, vlow_y, vlow_z
  CCTK_REAL betaxlow, betaylow, betazlow, beta2
  CCTK_REAL ulow_t, ulow_x, ulow_y, ulow_z
  CCTK_REAL rhoenthalpy
  CCTK_REAL bu_0, bu_x, bu_y, bu_z, uB
  CCTK_REAL bulow_0, bulow_x, bulow_y, bulow_z, b2
  CCTK_INT :: type_bits, atmosphere

  integer i,j,kk

!     if (SpaceMask_CheckStateBitsF90(space_mask, i,j,k, \
!                                     atmosphere_field_descriptor, \
!                                     atmosphere_normal_descriptor)) then
  call SpaceMask_GetTypeBits(type_bits, "Hydro_Atmosphere")
  call SpaceMask_GetStateBits(atmosphere, "Hydro_Atmosphere", "in_atmosphere")

  !$OMP PARALLEL DO PRIVATE (vlow_x, vlow_y, vlow_z, betaxlow, betaylow, &
  !$OMP betazlow, beta2, ulow_t, ulow_x, ulow_y, ulow_z, rhoenthalpy, &
  !$OMP uB,bu_0,bu_x,bu_y,bu_z,bulow_0,bulow_x,bulow_y,bulow_z,b2, &
  !$OMP i,j,kk)
  do kk = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)

           vlow_x = (gxx(i,j,kk)*velx(i,j,kk) + gxy(i,j,kk)*vely(i,j,kk) &
             + gxz(i,j,kk)*velz(i,j,kk))
           vlow_y = (gxy(i,j,kk)*velx(i,j,kk) + gyy(i,j,kk)*vely(i,j,kk) &
             + gyz(i,j,kk)*velz(i,j,kk))
           vlow_z = (gxz(i,j,kk)*velx(i,j,kk) + gyz(i,j,kk)*vely(i,j,kk) &
             + gzz(i,j,kk)*velz(i,j,kk))

!          Calculate lower components and square of shift vector.

           betaxlow = (gxx(i,j,kk)*betax(i,j,kk) + gxy(i,j,kk)*betay(i,j,kk) &
                 + gxz(i,j,kk)*betaz(i,j,kk))
           betaylow = (gxy(i,j,kk)*betax(i,j,kk) + gyy(i,j,kk)*betay(i,j,kk) &
                 + gyz(i,j,kk)*betaz(i,j,kk))
           betazlow = (gxz(i,j,kk)*betax(i,j,kk) + gyz(i,j,kk)*betay(i,j,kk) &
                 + gzz(i,j,kk)*betaz(i,j,kk))

           beta2 = betax(i,j,kk)*betaxlow + betay(i,j,kk)*betaylow &
                 + betaz(i,j,kk)*betazlow 

!          Calculate the specific relativistic enthalpy times rho times the
!          square of the lorentz factor.

           rhoenthalpy = w_lorentz(i,j,kk)**2*&
                        (rho(i,j,kk)*(1.0D0 + eps(i,j,kk)) + press(i,j,kk))

!          Calculate lower components of 4-velocity (without the Lorent factor).

           ulow_t = (-alp(i,j,kk) + velx(i,j,kk)*betaxlow + vely(i,j,kk)*betaylow &
                   + velz(i,j,kk)*betazlow)

           ulow_x = vlow_x
           ulow_y = vlow_y
           ulow_z = vlow_z

!	   Calculate MHD components (little b^\mu = B^mu_(u) / sqrt(4pi)) )
           uB = ulow_x * Bvx(i,j,kk) + ulow_y * Bvy(i,j,kk) + ulow_z * Bvz(i,j,kk)
           bu_0 = uB / alp(i,j,kk)
           bu_x = Bvx(i,j,kk)/w_lorentz(i,j,kk) &
                     + uB * ( velx(i,j,kk) + betax(i,j,kk)/alp(i,j,kk) )
           bu_y = Bvy(i,j,kk)/w_lorentz(i,j,kk) &
                     + uB * ( vely(i,j,kk) + betay(i,j,kk)/alp(i,j,kk) )
           bu_z = Bvz(i,j,kk)/w_lorentz(i,j,kk) &
                     + uB * ( velz(i,j,kk) + betaz(i,j,kk)/alp(i,j,kk) )

           bulow_0 = (beta2 - alp(i,j,kk)**2)*bu_0 + betaxlow*bu_x &
                    + betaylow*bu_y + betazlow*bu_z
           bulow_x = betaxlow*bu_0 + gxx(i,j,kk)*bu_x + gxy(i,j,kk)*bu_y + gxz(i,j,kk)*bu_z
           bulow_y = betaylow*bu_0 + gxy(i,j,kk)*bu_x + gyy(i,j,kk)*bu_y + gyz(i,j,kk)*bu_z
           bulow_z = betazlow*bu_0 + gxz(i,j,kk)*bu_x + gyz(i,j,kk)*bu_y + gzz(i,j,kk)*bu_z

           b2 = bu_0*bulow_0 + bu_x*bulow_x + bu_y*bulow_y + bu_z*bulow_z
           
!          Calculate Tmunu (the lower components!).

           if ( (whisky_atmosphere_is_vacuum .eq. 0) .or. &
              .not.( SpaceMask_CheckStateBitsF90(space_mask, i, j,kk, type_bits, atmosphere) .or. &
              (atmosphere_mask(i,j,kk) .ne. 0) )) then

              eTtt(i,j,kk) = eTtt(i,j,kk) + rhoenthalpy*ulow_t**2 &
                            + press(i,j,kk)*(beta2 - alp(i,j,kk)**2) &
                            + b2*(ulow_t**2+(beta2-alp(i,j,kk)**2)/2.0d0 ) &
                            - bulow_0**2

              eTtx(i,j,kk) = eTtx(i,j,kk) + rhoenthalpy*ulow_t*ulow_x &
                            + press(i,j,kk)*betaxlow &
                            + b2*(ulow_t*ulow_x + betaxlow/2.0d0) &
                            - bulow_0*bulow_x

              eTty(i,j,kk) = eTty(i,j,kk) + rhoenthalpy*ulow_t*ulow_y &
                            + press(i,j,kk)*betaylow &
                            + b2*(ulow_t*ulow_y + betaylow/2.0d0) &
                            - bulow_0*bulow_y

              eTtz(i,j,kk) = eTtz(i,j,kk) + rhoenthalpy*ulow_t*ulow_z &
                            + press(i,j,kk)*betazlow &
                            + b2*(ulow_t*ulow_z + betazlow/2.0d0) &
                            - bulow_0*bulow_z

              eTxx(i,j,kk) = eTxx(i,j,kk) + rhoenthalpy*ulow_x**2 &
                            + press(i,j,kk)*gxx(i,j,kk) &
                            + b2*(ulow_x**2 + gxx(i,j,kk)/2.0d0) &
                            - bulow_x**2

              eTyy(i,j,kk) = eTyy(i,j,kk) + rhoenthalpy*ulow_y**2 &
                            + press(i,j,kk)*gyy(i,j,kk) &
                            + b2*(ulow_y**2 + gyy(i,j,kk)/2.0d0) &
                            - bulow_y**2 

              eTzz(i,j,kk) = eTzz(i,j,kk) + rhoenthalpy*ulow_z**2 &
                            + press(i,j,kk)*gzz(i,j,kk) &
                            + b2*(ulow_z**2 + gzz(i,j,kk)/2.0d0) &
                            - bulow_z**2

              eTxy(i,j,kk) = eTxy(i,j,kk) + rhoenthalpy*ulow_x*ulow_y &
                            + press(i,j,kk)*gxy(i,j,kk) &
                            + b2*(ulow_x*ulow_y + gxy(i,j,kk)/2.0d0) &
                            - bulow_x*bulow_y

              eTxz(i,j,kk) = eTxz(i,j,kk) + rhoenthalpy*ulow_x*ulow_z &
                            + press(i,j,kk)*gxz(i,j,kk) &
                            + b2*(ulow_x*ulow_z + gxz(i,j,kk)/2.0d0) &
                            - bulow_x*bulow_z

              eTyz(i,j,kk) = eTyz(i,j,kk) + rhoenthalpy*ulow_y*ulow_z &
                            + press(i,j,kk)*gyz(i,j,kk) &
                            + b2*(ulow_y*ulow_z + gyz(i,j,kk)/2.0d0) &
                            - bulow_y*bulow_z
           end if

        end do
     end do
  end do
  !$OMP END PARALLEL DO

end subroutine Whisky_MHD_SetTmunu
