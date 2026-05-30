 /*@@
   @file      Whisky_EMTensor.inc
   @date      Sat Jan 26 02:03:56 2002
   @author    Ian Hawke
   @desc 
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


   @enddesc 
 @@*/

#ifdef FCODE

/*       if (SpaceMask_CheckStateBitsF90(space_mask, i,j,k, \
                                       atmosphere_field_descriptor, \
                                       atmosphere_normal_descriptor)) then*/
       call SpaceMask_GetTypeBits(type_bits, "Hydro_Atmosphere")
       call SpaceMask_GetStateBits(atmosphere, "Hydro_Atmosphere", "in_atmosphere")
       velxlow = gxx(i,j,k)*velx(i,j,k) + gxy(i,j,k)*vely(i,j,k)
     .           + gxz(i,j,k)*velz(i,j,k)
       velylow = gxy(i,j,k)*velx(i,j,k) + gyy(i,j,k)*vely(i,j,k)
     .           + gyz(i,j,k)*velz(i,j,k)
       velzlow = gxz(i,j,k)*velx(i,j,k) + gyz(i,j,k)*vely(i,j,k)
     .           + gzz(i,j,k)*velz(i,j,k)

c        Calculate lower components and square of shift vector.

       betaxlow = gxx(i,j,k)*betax(i,j,k) + gxy(i,j,k)*betay(i,j,k)
     .              + gxz(i,j,k)*betaz(i,j,k)
       betaylow = gxy(i,j,k)*betax(i,j,k) + gyy(i,j,k)*betay(i,j,k)
     .              + gyz(i,j,k)*betaz(i,j,k)
       betazlow = gxz(i,j,k)*betax(i,j,k) + gyz(i,j,k)*betay(i,j,k)
     .              + gzz(i,j,k)*betaz(i,j,k)
      beta2 = betax(i,j,k)*betaxlow + betay(i,j,k)*betaylow + betaz(i,j,k)*betazlow 

c        Calculate the specific relativistic enthalpy times rho times the
c        square of the lorentz factor.

       rhoenthalpy = w_lorentz(i,j,k)**2*
     .       (rho(i,j,k)*(1.0D0 + eps(i,j,k)) + press(i,j,k))

c        Calculate lower components of 4-velocity (without the Lorent factor).

       utlow = (-alp(i,j,k) + velx(i,j,k)*betaxlow + vely(i,j,k)*betaylow
     .       + velz(i,j,k)*betazlow)

       uxlow = velxlow
       uylow = velylow
       uzlow = velzlow

c        Calculate Tmunu (the lower components!).

       if ( (whisky_atmosphere_is_vacuum .eq. 0) .or. 
     . .not.( SpaceMask_CheckStateBitsF90(space_mask, i, j, k, type_bits, atmosphere)  .or.
     .        (atmosphere_mask(i,j,k) .ne. 0) )) then

         Ttt = Ttt + rhoenthalpy*utlow**2
     .       + press(i,j,k)*(beta2 - alp(i,j,k)**2)

         Ttx = Ttx + rhoenthalpy*utlow*uxlow + press(i,j,k)*betaxlow
         Tty = Tty + rhoenthalpy*utlow*uylow + press(i,j,k)*betaylow
         Ttz = Ttz + rhoenthalpy*utlow*uzlow + press(i,j,k)*betazlow

         Txx = Txx + rhoenthalpy*uxlow**2 + press(i,j,k)*gxx(i,j,k)
         Tyy = Tyy + rhoenthalpy*uylow**2 + press(i,j,k)*gyy(i,j,k)
         Tzz = Tzz + rhoenthalpy*uzlow**2 + press(i,j,k)*gzz(i,j,k)

         Txy = Txy + rhoenthalpy*uxlow*uylow + press(i,j,k)*gxy(i,j,k)
         Txz = Txz + rhoenthalpy*uxlow*uzlow + press(i,j,k)*gxz(i,j,k)
         Tyz = Tyz + rhoenthalpy*uylow*uzlow + press(i,j,k)*gyz(i,j,k)

       end if
/*       end if*/

#endif

#ifdef CCODE

       type_bits = SpaceMask_GetTypeBits("Hydro_Atmosphere")
       atmosphere = SpaceMask_GetStateBits("Hydro_Atmosphere", "in_atmosphere")
       velxlow = gxx[ijk]*velx[ijk] + gxy[ijk]*vely[ijk] + 
                 gxz[ijk]*velz[ijk];
       velylow = gxy[ijk]*velx[ijk] + gyy[ijk]*vely[ijk] + 
                 gyz[ijk]*velz[ijk];
       velzlow = gxz[ijk]*velx[ijk] + gyz[ijk]*vely[ijk] + 
                 gzz[ijk]*velz[ijk];

/*
 *        Calculate lower components and square of shift vector.
 */

        betaxlow = gxx[ijk]*betax[ijk] + gxy[ijk]*betay[ijk] 
                  + gxz[ijk]*betaz[ijk];
        betaylow = gxy[ijk]*betax[ijk] + gyy[ijk]*betay[ijk]
                  + gyz[ijk]*betaz[ijk];
        betazlow = gxz[ijk]*betax[ijk] + gyz[ijk]*betay[ijk]
                  + gzz[ijk]*betaz[ijk];
        beta2 = betax[ijk]*betaxlow + betay[ijk]*betaylow
                + betaz[ijk]*betazlow ;

/*
 *        Calculate the specific relativistic enthalpy times rho times the
 *        square of the lorentz factor.
 */

       rhoenthalpy = w_lorentz[ijk]*w_lorentz[ijk]*
                      (rho[ijk]*(1.0 + eps[ijk]) + press[ijk]);

/*
 *        Calculate lower components of 4-velocity (without the Lorent factor).
 */
       utlow = (-alp[ijk] + velx[ijk]*betaxlow + vely[ijk]*betaylow
                          + velz[ijk]*betazlow);

       uxlow = velxlow;
       uylow = velylow;
       uzlow = velzlow;

/*
 *        Calculate Tmunu (the lower components!).
 */

       if ((!*whisky_atmosphere_is_vacuum) || 
	 !( (SpaceMask_CheckStateBits(space_mask, ijk, type_bits, atmosphere) ||
             atmosphere_mask[ijk] ) ) {

       Ttt = Ttt + rhoenthalpy*utlow*utlow +
             press[ijk]*(beta2 - alp[ijk]*alp[ijk]);

       Ttx = Ttx + rhoenthalpy*utlow*uxlow + press[ijk]*betaxlow;
       Tty = Tty + rhoenthalpy*utlow*uylow + press[ijk]*betaylow;
       Ttz = Ttz + rhoenthalpy*utlow*uzlow + press[ijk]*betazlow;

       Txx = Txx + rhoenthalpy*uxlow*uxlow + press[ijk]*gxx[ijk];
       Tyy = Tyy + rhoenthalpy*uylow*uylow + press[ijk]*gyy[ijk];
       Tzz = Tzz + rhoenthalpy*uzlow*uzlow + press[ijk]*gzz[ijk];

       Txy = Txy + rhoenthalpy*uxlow*uylow + press[ijk]*gxy[ijk];
       Txz = Txz + rhoenthalpy*uxlow*uzlow + press[ijk]*gxz[ijk];
       Tyz = Tyz + rhoenthalpy*uylow*uzlow + press[ijk]*gyz[ijk];

       }

#endif
