 /*@@
   @file      PPM.F90
   @date      Sun Feb 10 16:53:29 2002
   @author    Ian Hawke, Toni Font, Luca Baiotti, Frank Loeffler
   @desc 
   Routines to do PPM reconstruction.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "PPM_Utils.h"

module PPM_Simple1d_Routines

 CONTAINS

 /*@@
   @routine    SimplePPM_1d
   @date       Thu Feb 14 19:08:52 2002
   @author     Ian Hawke, Toni Font
   @desc 
   The simple PPM reconstruction routine that applies along
   each one dimensional slice.

   @enddesc 
   @calls     
   @calledby   
   @history 
   Written in frustration when IH couldn''t get Toni''s original code 
   to work.
   @endhistory 

@@*/

#define SpaceMask_CheckStateBitsF90_1D(mask,i,type_bits,state_bits) \
  (iand(mask((i)),(type_bits)).eq.(state_bits))


  subroutine SimplePPM_1d(handle,poly,&
       nx,dx,rho,velx,vely,velz,eps,press,rhominus,&
       velxminus,velyminus,velzminus,epsminus,rhoplus,velxplus,velyplus,&
       velzplus,epsplus,trivial_rp,space_mask, excision_descriptors,&
       gxx, gxy, gxz, gyy, gyz, gzz, beta, alp, w_lorentz, &
       dir, ni, nj, nrx, nry, nrz, ev_l, ev_r, xw)

    USE PPM_Scalars
    USE PPM_EigenvalueCalcs

    implicit none

    DECLARE_CCTK_PARAMETERS
    DECLARE_CCTK_FUNCTIONS

#ifdef _EOS_BASE_INC_
#undef _EOS_BASE_INC_
#endif
#include "EOS_Base.inc"

    CCTK_INT :: handle,poly,nx
    CCTK_REAL :: dx
    CCTK_REAL, dimension(nx) :: rho,velx,vely,velz,eps
    CCTK_REAL, dimension(nx) :: vx,vy,vz
    CCTK_REAL, dimension(nx) :: rhominus,velxminus,velyminus,velzminus,epsminus
    CCTK_REAL, dimension(nx) :: rhoplus,velxplus,velyplus,velzplus,epsplus
    CCTK_REAL, dimension(nx) :: rhominusl,velxminusl,velyminusl,velzminusl
    CCTK_REAL, dimension(nx) :: epsminusl
    CCTK_REAL, dimension(nx) :: rhoplusl,velxplusl,velyplusl,velzplusl,epsplusl
    CCTK_REAL, dimension(nx) :: rhominusr,velxminusr,velyminusr,velzminusr
    CCTK_REAL, dimension(nx) :: epsminusr
    CCTK_REAL, dimension(nx) :: rhoplusr,velxplusr,velyplusr,velzplusr,epsplusr

    CCTK_INT :: i,s
    CCTK_REAL, dimension(nx) :: drho,dvelx,dvely,dvelz,deps
    CCTK_REAL, dimension(nx) :: dmrho,dmvelx,dmvely,dmvelz,dmeps
    CCTK_REAL, dimension(nx) :: press,dpress,d2rho,tilde_flatten
    CCTK_REAL :: dpress2,dvel,w,flatten,eta,etatilde

    logical, dimension(nx) :: trivial_rp

    CCTK_INT, dimension(nx) :: space_mask
    CCTK_INT :: excision_bits, excision_mask
    CCTK_INT, dimension(3) :: excision_descriptors

    CCTK_REAL, dimension(nx) :: gxx, gxy, gxz, gyy, gyz, gzz, &
                              beta, alp, w_lorentz
    CCTK_INT :: dir, ni, nj, nrx, nry, nrz
    CCTK_REAL, dimension(nrx, nry, nrz) :: ev_l, ev_r, xw

    CCTK_REAL :: uxx, uxy, uxz, uyy, uyz, uzz, det
    CCTK_REAL, dimension(5) :: lam
    CCTK_REAL :: dupw, dloc, delta
    CCTK_REAL :: agxx, agxy, agxz, agyy, agyz, agzz
    CCTK_REAL, dimension(nx) :: xwind, l_ev_l, l_ev_r

    logical :: cond
    character(len=400) :: warnline


  !!$  Initially, all the Riemann problems will be trivial

  trivial_rp = .true.

  !!$  Loop to set the pressure and eps
  !!$  poly=0
  !!$  if (poly .eq. 0) then
  !!$    do i = 1, nx
  !!$      press(i) = EOS_Pressure(handle, rho(i), eps(i))
  !!$    end do
  !!$  else
  !!$    do i = 1, nx
  !!$      press(i) = EOS_Pressure(handle, rho(i), eps(i))
  !!$      eps(i) = EOS_SpecificIntEnergy(handle, rho(i), press(i))
  !!$    end do
  !!$  end if
  
  !!$ reconstruct w^i = w_lorentz*v^i to ensure slower than light speeds?
    if (reconstruct_Wv.ne.0) then
      vx = w_lorentz*velx
      vy = w_lorentz*vely
      vz = w_lorentz*velz
    else
      vx = velx
      vy = vely
      vz = velz
    end if

  !!$  Average slopes delta_m(a). See (1.7) of Colella and Woodward, p.178
  !!$  This is the expression for an even grid.
  
    do i = 2, nx - 1
      drho(i) = 0.5d0 * (rho(i+1) - rho(i-1))
      dvelx(i) = 0.5d0 * (vx(i+1) - vx(i-1))
      dvely(i) = 0.5d0 * (vy(i+1) - vy(i-1))
      dvelz(i) = 0.5d0 * (vz(i+1) - vz(i-1))
      dpress(i) = press(i+1) - press(i-1)
      d2rho(i) = (rho(i+1) - 2.d0 * rho(i) + rho(i-1))! / 6.d0 / dx / dx
      ! since we use d2rho only for the condition d2rho(i+1)*d2rhoi(i-1)<0 
      ! the denominator is not necessary
    end do
    if (poly .eq. 0) then
      do i = 2, nx - 1
        deps(i) = 0.5d0 * (eps(i+1) - eps(i-1))
      end do
    end if
  
  !!$  Steepened slope. See (1.8) of Colella and Woodward, p.178
  
    do i = 2, nx - 1
#define STEEP(x,dx,dmx)                                              \
    if ( (x(i+1) - x(i)) * (x(i) - x(i-1)) > 0.d0 ) then           &&\
      dmx(i) = sign(1.d0, dx(i)) *                                   \
           min(abs(dx(i)), 2.d0 * abs(x(i) - x(i-1)),                \
                           2.d0 * abs(x(i+1) - x(i)))              &&\
    else                                                           &&\
      dmx(i) = 0.d0                                                &&\
    end if
      STEEP(rho, drho, dmrho)
      STEEP(vx, dvelx, dmvelx)
      STEEP(vy, dvely, dmvely)
      STEEP(vz, dvelz, dmvelz)
    end do
    if (poly .eq. 0) then
      do i = 2, nx - 1
        STEEP(eps, deps, dmeps)
      end do
    end if
  
  !!$  Initial boundary states. See (1.9) of Colella and Woodward, p.178
  
    do i = 2, nx-2
      rhoplus(i) = 0.5d0 * (rho(i) + rho(i+1)) + &
           (dmrho(i) - dmrho(i+1)) / 6.d0
      rhominus(i+1) = rhoplus(i)
      velxplus(i) = 0.5d0 * (vx(i) + vx(i+1)) + &
           (dmvelx(i) - dmvelx(i+1)) / 6.d0
      velxminus(i+1) = velxplus(i)
      velyplus(i) = 0.5d0 * (vy(i) + vy(i+1)) + &
           (dmvely(i) - dmvely(i+1)) / 6.d0
      velyminus(i+1) = velyplus(i)
      velzplus(i) = 0.5d0 * (vz(i) + vz(i+1)) + &
           (dmvelz(i) - dmvelz(i+1)) / 6.d0
      velzminus(i+1) = velzplus(i)
    end do
    if (poly .eq. 0) then
      do i = 2, nx-2
        epsplus(i) = 0.5d0 * (eps(i) + eps(i+1)) + &
             (dmeps(i) - dmeps(i+1)) / 6.d0
        epsminus(i+1) = epsplus(i)
      end do
    end if
  
  !!$Discontinuity steepening. See (1.14-17) of C&W.
  !!$This is the detect routine which mat be activated with the ppm_detect parameter
  !!$Note that this part really also depends on the grid being even. 
  !!$Note also that we don''t have access to the gas constant gamma.
  !!$So this is just dropped from eq. (3.2) of C&W.
  !!$We can get around this by just rescaling the constant k0 (ppm_k0 here).
  
    if (ppm_detect .ne. 0) then
  
      do i = 3, nx - 2
        if ( (d2rho(i+1)*d2rho(i-1) < 0.d0).and.(abs(rho(i+1)-rho(i-1)) - &
             ppm_epsilon_shock * min(abs(rho(i+1)), abs(rho(i-1))) > 0.d0) ) then
          etatilde = (rho(i-2) - rho(i+2) + 4.d0 * drho(i)) / (drho(i) * 12.d0)
        else
          etatilde = 0.d0
        end if
        eta = max(0.d0, min(1.d0, ppm_eta1 * (etatilde - ppm_eta2)))
        if (ppm_k0 * abs(drho(i)) * min(press(i-1),press(i+1)) < &
             abs(dpress(i)) * min(rho(i-1), rho(i+1))) then
          eta = 0.d0
        end if
        if (eta > 0.d0) then
          trivial_rp(i-1) = .false.
          trivial_rp(i) = .false.
        end if
        rhominus(i) = rhominus(i) * (1.d0 - eta) + &
             (rho(i-1) + 0.5d0 * dmrho(i-1)) * eta
        rhoplus(i) = rhoplus(i) * (1.d0 - eta) + &
             (rho(i+1) - 0.5d0 * dmrho(i+1)) * eta
      end do
  
    end if
  
    !!$ mppm
#define D_UPW(x) (0.5d0 * (x(i) + x(i+1)))
#define LEFT1(x)  (13.d0*x(i+1)-5.d0*x(i+2)+x(i+3)+3.d0*x(i  ))/12.d0
#define RIGHT1(x) (13.d0*x(i  )-5.d0*x(i-1)+x(i-2)+3.d0*x(i+1))/12.d0
    if (ppm_mppm .gt. 0) then
      l_ev_l=0.d0
      l_ev_r=0.d0
      xwind=0.d0
      do i=3, nx - 3
        agxx = 0.5d0*( gxx(i) + gxx(i+1) )
        agxy = 0.5d0*( gxy(i) + gxy(i+1) )
        agxz = 0.5d0*( gxz(i) + gxz(i+1) )
        agyy = 0.5d0*( gyy(i) + gyy(i+1) )
        agyz = 0.5d0*( gyz(i) + gyz(i+1) )
        agzz = 0.5d0*( gzz(i) + gzz(i+1) )
        det = SPATIAL_DET(agxx,agxy,agxz,agyy,agyz,agzz)
        call UpperMetric (uxx, uxy, uxz, uyy, uyz, uzz, &
                          det, agxx, agxy, agxz, agyy, agyz, agzz)
  !!$ these always have to be vel[xyz], not v[xyz] since they are used to
  !!$ calculate the eigenvalues
        call eigenvalues_ppm(handle,&
                   D_UPW(rho), D_UPW(velx), D_UPW(vely), D_UPW(velz), &
                   D_UPW(eps), D_UPW(w_lorentz), lam, &
                   agxx, agxy, agxz, agyy, agyz, agzz, &
                   uxx, det, D_UPW(alp), D_UPW(beta))
        l_ev_l(i)=lam(1)
        l_ev_r(i)=lam(5)
        xwind(i) = (lam(1) + lam(5)) / (abs(lam(1)) + abs(lam(5)))
        xwind(i) = min(1.d0, max(-1.d0, xwind(i)))
#define LEFTPLUS(x,xplus)    xplus(i)   =       abs(xwind(i))  * LEFT1(x) + \
                                          (1.d0-abs(xwind(i))) * xplus(i)
#define LEFTMINUS(x,xminus)  xminus(i+1)=       abs(xwind(i))  * LEFT1(x) + \
                                          (1.d0-abs(xwind(i))) * xminus(i+1)
#define RIGHTPLUS(x,xplus)   xplus(i)   =       abs(xwind(i))  * RIGHT1(x) + \
                                          (1.d0-abs(xwind(i))) * xplus(i)
#define RIGHTMINUS(x,xminus) xminus(i+1)=       abs(xwind(i))  * RIGHT1(x) + \
                                          (1.d0-abs(xwind(i))) * xminus(i+1)
#define CHECK(x,xc) if (x(i+1) .gt. x(i)) then && xc=min(x(i+1),max(x(i),xc)) && else && xc=min(x(i),max(x(i+1),xc)) && endif
  !!$      xwind(i)=0.d0
        if (xwind(i) .lt. 0.0d0) then
          LEFTPLUS(rho, rhoplus)
          LEFTMINUS(rho, rhominus)
          LEFTPLUS(vx, velxplus)
          LEFTMINUS(vx, velxminus)
          LEFTPLUS(vy, velyplus)
          LEFTMINUS(vy, velyminus)
          LEFTPLUS(vz, velzplus)
          LEFTMINUS(vz, velzminus)
          if (poly .eq. 0) then
            LEFTPLUS(eps, epsplus)
            LEFTMINUS(eps, epsminus)
          end if
        else
          RIGHTPLUS(rho, rhoplus)
          RIGHTMINUS(rho, rhominus)
          RIGHTPLUS(vx, velxplus)
          RIGHTMINUS(vx, velxminus)
          RIGHTPLUS(vy, velyplus)
          RIGHTMINUS(vy, velyminus)
          RIGHTPLUS(vz, velzplus)
          RIGHTMINUS(vz, velzminus)
          if (poly .eq. 0) then
            RIGHTPLUS(eps, epsplus)
            RIGHTMINUS(eps, epsminus)
          end if
        end if
        CHECK(rho, rhoplus(i))
        CHECK(rho, rhominus(i+1))
        CHECK(vx, velxplus(i))
        CHECK(vx, velxminus(i+1))
        CHECK(vy, velyplus(i))
        CHECK(vy, velyminus(i+1))
        CHECK(vz, velzplus(i))
        CHECK(vz, velzminus(i+1))
        if (poly .eq. 0) then
          CHECK(eps, epsplus(i))
          CHECK(eps, epsminus(i+1))
        end if
  !!$      if ((dir .eq. 1) .and. (ni .eq. 4) .and. (nj .eq. 4)) then
  !!$        write (*,*) rhoplus(i), rhominus(i+1)
  !!$      end if
      end do
      !!$ mppm debug output
      if (mppm_debug_eigenvalues .gt. 0) then
        if (dir .eq. 1) then
          ev_l(:,ni,nj) = l_ev_l
          ev_r(:,ni,nj) = l_ev_r
          xw(:,ni,nj) = xwind
        else if (dir .eq. 2) then
          ev_l(ni,:,nj) = l_ev_l
          ev_r(ni,:,nj) = l_ev_r
          xw(ni,:,nj) = xwind
        else if (dir .eq. 3) then
          ev_l(ni,nj,:) = l_ev_l
          ev_r(ni,nj,:) = l_ev_r
          xw(ni,nj,:) = xwind
        else
          write (*,*) "flux direction not 1 to 3 ?"
        end if
      end if
    end if
  
  !!$  Zone flattening. See appendix of C&W, p. 197-8.
  
    do i = 3, nx - 2
      dpress2 = press(i+2) - press(i-2)
      dvel = vx(i+1) - vx(i-1)
      if ( (abs(dpress(i)) >  ppm_epsilon * min(press(i-1),press(i+1))) .and. &
           (dvel < 0.d0) ) then
        w = 1.d0
      else
        w = 0.d0
      end if
      if (abs(dpress2) < ppm_small) then
        tilde_flatten(i) = 1.d0
      else
        tilde_flatten(i) = max(0.d0, 1.d0 - w * max(0.d0, ppm_omega2 * &
             (dpress(i) / dpress2 - ppm_omega1)))
      end if
    end do
  
  
  
    if (PPM3) then !!$ Implement C&W, page 197, but with a workaround which allows to use stencil=3.
      do i = 3, nx - 2
        flatten = tilde_flatten(i)
        if (abs(1.d0 - flatten) > 0.d0) then
          trivial_rp(i-1) = .false.
          trivial_rp(i) = .false.
        end if
        rhoplus(i) = flatten * rhoplus(i) + (1.d0 - flatten) * rho(i)
        rhominus(i) = flatten * rhominus(i) + (1.d0 - flatten) * rho(i)
        velxplus(i) = flatten * velxplus(i) + (1.d0 - flatten) * vx(i)
        velxminus(i) = flatten * velxminus(i) + (1.d0 - flatten) * vx(i)
        velyplus(i) = flatten * velyplus(i) + (1.d0 - flatten) * vy(i)
        velyminus(i) = flatten * velyminus(i) + (1.d0 - flatten) * vy(i)
        velzplus(i) = flatten * velzplus(i) + (1.d0 - flatten) * vz(i)
        velzminus(i) = flatten * velzminus(i) + (1.d0 - flatten) * vz(i)
        if (poly .eq. 0) then
          epsplus(i) = flatten * epsplus(i) + (1.d0 - flatten) * eps(i)
          epsminus(i) = flatten * epsminus(i) + (1.d0 - flatten) * eps(i)
        end if
      end do
    else  !!$ Really implement C&W, page 197; which requires stencil 4.
      do i = 4, nx - 3
        s=sign(1.d0, -dpress(i))
        flatten = max(tilde_flatten(i), tilde_flatten(i+s))  
        if (abs(1.d0 - flatten) > 0.d0) then
          trivial_rp(i-1) = .false.
          trivial_rp(i) = .false.
        end if
        rhoplus(i) = flatten * rhoplus(i) + (1.d0 - flatten) * rho(i)
        rhominus(i) = flatten * rhominus(i) + (1.d0 - flatten) * rho(i)
        velxplus(i) = flatten * velxplus(i) + (1.d0 - flatten) * vx(i)
        velxminus(i) = flatten * velxminus(i) + (1.d0 - flatten) * vx(i)
        velyplus(i) = flatten * velyplus(i) + (1.d0 - flatten) * vy(i)
        velyminus(i) = flatten * velyminus(i) + (1.d0 - flatten) * vy(i)
        velzplus(i) = flatten * velzplus(i) + (1.d0 - flatten) * vz(i)
        velzminus(i) = flatten * velzminus(i) + (1.d0 - flatten) * vz(i)
        if (poly .eq. 0) then
          epsplus(i) = flatten * epsplus(i) + (1.d0 - flatten) * eps(i)
          epsminus(i) = flatten * epsminus(i) + (1.d0 - flatten) * eps(i)
        end if
      end do
    end if
  
  
  !!$ Monotonicity. See (1.10) of C&W.

  do i = whisky_stencil, nx - whisky_stencil + 1
#define MON(xminus,x,xplus)                                       \
    if (.not.( (xplus(i).eq.x(i)) .and. (x(i).eq.xminus(i)) )     \
        .and. ((xplus(i)-x(i))*(x(i)-xminus(i)) .le. 0.d0)) then&&\
      trivial_rp(i-1) = .false.                                 &&\
      trivial_rp(i) = .false.                                   &&\
      xminus(i) = x(i)                                          &&\
      xplus(i) = x(i)                                           &&\
    else if (6.d0 * (xplus(i) - xminus(i)) * (x(i) - 0.5d0 *      \
                (xplus(i) + xminus(i))) >                         \
                (xplus(i) - xminus(i))**2) then                 &&\
      xminus(i) = 3.d0 * x(i) - 2.d0 * xplus(i)                 &&\
      trivial_rp(i-1) = .false.                                 &&\
      trivial_rp(i) = .false.                                   &&\
    else if (6.d0 * (xplus(i) - xminus(i)) * (x(i) - 0.5d0 *      \
                (xplus(i) + xminus(i))) <                         \
               -(xplus(i) - xminus(i))**2) then                 &&\
      xplus(i) = 3.d0 * x(i) - 2.d0 * xminus(i)                 &&\
      trivial_rp(i-1) = .false.                                 &&\
      trivial_rp(i) = .false.                                   &&\
    end if                                                      &&\
    if (.not.( (xplus(i).eq.x(i)) .and. (x(i).eq.xminus(i)) ) ) then     &&\
      trivial_rp(i-1) = .false.                                 &&\
      trivial_rp(i) = .false.                                   &&\
    end if

      MON(rhominus,rho,rhoplus)
      MON(velxminus,vx,velxplus)
      MON(velyminus,vy,velyplus)
      MON(velzminus,vz,velzplus)
    end do
    if (poly .eq. 0) then
      do i = whisky_stencil, nx - whisky_stencil + 1
        MON(epsminus,eps,epsplus)
      end do
    end if
  
    if (check_for_trivial_rp .eq. 0) then
      trivial_rp = .false.
    end if
  
    !!$ excision
  
    !!$ recover primitive v^i from W*v^i at interfaces
    if (reconstruct_Wv.ne.0) then
      do i = whisky_stencil, nx - whisky_stencil + 1
        agxx = 0.5d0*( gxx(i) + gxx(i+1) )
        agxy = 0.5d0*( gxy(i) + gxy(i+1) )
        agxz = 0.5d0*( gxz(i) + gxz(i+1) )
        agyy = 0.5d0*( gyy(i) + gyy(i+1) )
        agyz = 0.5d0*( gyz(i) + gyz(i+1) )
        agzz = 0.5d0*( gzz(i) + gzz(i+1) )
        ! we re-use the variable w defined for zone flattening to hold the lorentz
        ! factor
        ! divide out the Loretnz factor for both the plus and minus quantities
        ! this should by construction ensure that any Lorentz factor calculated
        ! from them later on is physical (ie. > 1.d0)
        w = sqrt( 1.d0 + agxx*velxminus(i)*velxminus(i) + agyy*velyminus(i)*velyminus(i) &
                + agzz*velzminus(i)*velzminus(i) + 2.d0*agxy*velxminus(i)*velyminus(i) &
                + 2.d0*agxz*velxminus(i)*velzminus(i) + 2.d0*agyz*velyminus(i)*velzminus(i) )
        if (w.lt.1.d0 .or. w.ne.w) then ! .and. Whisky_CarpetWeights(i,j,k).eq.1.d0 (no CCTK_ARGUMENTS)
          write(warnline,'("SimplePPM_1d (minus): Unphysical Lorentz factor ",g15.6,&
            &" for data g = [",5(g15.6,","),g15.6,"] vel = [",&
            &2(g15.6,","),g15.6,"]")') &
            w, agxx,agxy,agxz,agyy,agyz,agzz,&
            velxminus(i),velyminus(i),velzminus(i)
          call CCTK_WARN(1,warnline)
        end if 
        velxminus(i) = velxminus(i)/w
        velyminus(i) = velyminus(i)/w
        velzminus(i) = velzminus(i)/w
        w = sqrt( 1.d0 + agxx*velxplus(i)*velxplus(i) + agyy*velyplus(i)*velyplus(i) &
                + agzz*velzplus(i)*velzplus(i) + 2.d0*agxy*velxplus(i)*velyplus(i) &
                + 2.d0*agxz*velxplus(i)*velzplus(i) + 2.d0*agyz*velyplus(i)*velzplus(i) )
        if (w.lt.1.d0 .or. w.ne.w) then ! .and. Whisky_CarpetWeights(i,j,k).eq.1.d0 (no CCTK_ARGUMENTS)
          write(warnline,'("SimplePPM_1d (plus): Unphysical Lorentz factor ",g15.6,&
            &" for data g = [",5(g15.6,","),g15.6,"] vel = [",&
            &2(g15.6,","),g15.6,"]")') &
            w, agxx,agxy,agxz,agyy,agyz,agzz,&
            velxplus(i),velyplus(i),velzplus(i)
          call CCTK_WARN(1,warnline)
        end if 
        velxplus(i) = velxplus(i)/w
        velyplus(i) = velyplus(i)/w
        velzplus(i) = velzplus(i)/w
      end do
    end if
  
    return
  
  end subroutine SimplePPM_1d

  /*@@
   @routine    SimplePPM_1dMHD
   @date       Thu Feb 14 19:08:52 2002
   @author     Ian Hawke, Toni Font, Tanja Bode
   @desc 
   Extended reconstruction for MHD:
   The simple PPM reconstruction routine that applies along
   each one dimensional slice.

   @enddesc 
   @calls     
   @calledby   
   @history 
   Written in frustration when IH couldn''t get Toni''s original code 
   to work.
   @endhistory 

  @@*/

  subroutine SimplePPM_1dMHD(handle,poly,&
       nx,dx,rho,velx,vely,velz,eps,press, &
       Bx, By, Bz, &
       rhominus,velxminus,velyminus,velzminus,epsminus,&
       Bxminus, Byminus, Bzminus,&
       rhoplus,velxplus,velyplus,velzplus,epsplus,&
       Bxplus, Byplus, Bzplus,&
       trivial_rp,space_mask, excision_descriptors,&
       gxx, gxy, gxz, gyy, gyz, gzz, beta, alp, w_lorentz, &
       dir, ni, nj, nrx, nry, nrz, ev_l, ev_r, xw)
  
    USE PPM_Scalars
    USE PPM_EigenvalueCalcs
   
    implicit none
  
    DECLARE_CCTK_PARAMETERS
    DECLARE_CCTK_FUNCTIONS
  
#ifdef _EOS_BASE_INC_
#undef _EOS_BASE_INC_
#endif
#include "EOS_Base.inc"

    CCTK_INT :: handle,poly,nx
    CCTK_REAL :: dx
    CCTK_REAL, dimension(nx) :: rho,velx,vely,velz,eps
    CCTK_REAL, dimension(nx) :: vx,vy,vz
    CCTK_REAL, dimension(nx) :: Bx,By,Bz
    CCTK_REAL, dimension(nx) :: Bxminus,Byminus,Bzminus
    CCTK_REAL, dimension(nx) :: Bxplus,Byplus,Bzplus
    CCTK_REAL, dimension(nx) :: rhominus,velxminus,velyminus,velzminus,epsminus
    CCTK_REAL, dimension(nx) :: rhoplus,velxplus,velyplus,velzplus,epsplus
    CCTK_REAL, dimension(nx) :: rhominusl,velxminusl,velyminusl,velzminusl
    CCTK_REAL, dimension(nx) :: epsminusl
    CCTK_REAL, dimension(nx) :: rhoplusl,velxplusl,velyplusl,velzplusl,epsplusl
    CCTK_REAL, dimension(nx) :: rhominusr,velxminusr,velyminusr,velzminusr
    CCTK_REAL, dimension(nx) :: epsminusr
    CCTK_REAL, dimension(nx) :: rhoplusr,velxplusr,velyplusr,velzplusr,epsplusr
  
    CCTK_INT :: i,s
    CCTK_REAL, dimension(nx) :: drho,dvelx,dvely,dvelz,deps
    CCTK_REAL, dimension(nx) :: dBx,dBy,dBz
    CCTK_REAL, dimension(nx) :: dmrho,dmvelx,dmvely,dmvelz,dmeps
    CCTK_REAL, dimension(nx) :: dmBx,dmBy,dmBz
    CCTK_REAL, dimension(nx) :: press,dpress,d2rho,tilde_flatten
    CCTK_REAL :: dpress2,dvel,w,flatten,eta,etatilde
  
    logical, dimension(nx) :: trivial_rp
  
    CCTK_INT, dimension(nx) :: space_mask
    CCTK_INT :: excision_bits, excision_mask
    CCTK_INT, dimension(3) :: excision_descriptors
  
    CCTK_REAL, dimension(nx) :: gxx, gxy, gxz, gyy, gyz, gzz, &
                                beta, alp, w_lorentz
    CCTK_INT :: dir, ni, nj, nrx, nry, nrz
    CCTK_REAL, dimension(nrx, nry, nrz) :: ev_l, ev_r, xw
  
    CCTK_REAL :: uxx, uxy, uxz, uyy, uyz, uzz, det
    CCTK_REAL, dimension(8) :: lam
    CCTK_REAL :: dupw, dloc, delta
    CCTK_REAL :: agxx, agxy, agxz, agyy, agyz, agzz
    CCTK_REAL, dimension(nx) :: xwind, l_ev_l, l_ev_r
  
    logical :: cond
    character(len=400) :: warnline
  
  
  !!$  Initially, all the Riemann problems will be trivial
  
  trivial_rp = .true.
  
  !!$  Loop to set the pressure and eps
  !!$  poly=0
  !!$  if (poly .eq. 0) then
  !!$    do i = 1, nx
  !!$      press(i) = EOS_Pressure(handle, rho(i), eps(i))
  !!$    end do
  !!$  else
  !!$    do i = 1, nx
  !!$      press(i) = EOS_Pressure(handle, rho(i), eps(i))
  !!$      eps(i) = EOS_SpecificIntEnergy(handle, rho(i), press(i))
  !!$    end do
  !!$  end if
  
  !!$ reconstruct w^i = w_lorentz*v^i to ensure slower than light speeds?
    if (reconstruct_Wv.ne.0) then
      vx = w_lorentz*velx
      vy = w_lorentz*vely
      vz = w_lorentz*velz
    else
      vx = velx
      vy = vely
      vz = velz
    end if
  
  !!$  Average slopes delta_m(a). See (1.7) of Colella and Woodward, p.178
  !!$  This is the expression for an even grid.
  
    do i = 2, nx - 1
      drho(i) = 0.5d0 * (rho(i+1) - rho(i-1))
      dvelx(i) = 0.5d0 * (vx(i+1) - vx(i-1))
      dvely(i) = 0.5d0 * (vy(i+1) - vy(i-1))
      dvelz(i) = 0.5d0 * (vz(i+1) - vz(i-1))
      dBx(i) = 0.5d0 * (Bx(i+1) - Bx(i-1))
      dBy(i) = 0.5d0 * (By(i+1) - By(i-1))
      dBz(i) = 0.5d0 * (Bz(i+1) - Bz(i-1))
      dpress(i) = press(i+1) - press(i-1)
      d2rho(i) = (rho(i+1) - 2.d0 * rho(i) + rho(i-1))! / 6.d0 / dx / dx
      ! since we use d2rho only for the condition d2rho(i+1)*d2rhoi(i-1)<0 
      ! the denominator is not necessary
    end do
    if (poly .eq. 0) then
      do i = 2, nx - 1
        deps(i) = 0.5d0 * (eps(i+1) - eps(i-1))
      end do
    end if
  
  !!$  Steepened slope. See (1.8) of Colella and Woodward, p.178
  
    do i = 2, nx - 1
#define STEEP(x,dx,dmx)                                              \
    if ( (x(i+1) - x(i)) * (x(i) - x(i-1)) > 0.d0 ) then           &&\
      dmx(i) = sign(1.d0, dx(i)) *                                   \
           min(abs(dx(i)), 2.d0 * abs(x(i) - x(i-1)),                \
                           2.d0 * abs(x(i+1) - x(i)))              &&\
    else                                                           &&\
      dmx(i) = 0.d0                                                &&\
    end if
      STEEP(rho, drho, dmrho)
      STEEP(vx, dvelx, dmvelx)
      STEEP(vy, dvely, dmvely)
      STEEP(vz, dvelz, dmvelz)
      STEEP(Bx, dBx, dmBx)
      STEEP(By, dBy, dmBy)
      STEEP(Bz, dBz, dmBz)
    end do
    if (poly .eq. 0) then
      do i = 2, nx - 1
        STEEP(eps, deps, dmeps)
      end do
    end if
  
  !!$  Initial boundary states. See (1.9) of Colella and Woodward, p.178
  
    do i = 2, nx-2
      rhoplus(i) = 0.5d0 * (rho(i) + rho(i+1)) + &
           (dmrho(i) - dmrho(i+1)) / 6.d0
      rhominus(i+1) = rhoplus(i)
      velxplus(i) = 0.5d0 * (vx(i) + vx(i+1)) + &
           (dmvelx(i) - dmvelx(i+1)) / 6.d0
      velxminus(i+1) = velxplus(i)
      velyplus(i) = 0.5d0 * (vy(i) + vy(i+1)) + &
           (dmvely(i) - dmvely(i+1)) / 6.d0
      velyminus(i+1) = velyplus(i)
      velzplus(i) = 0.5d0 * (vz(i) + vz(i+1)) + &
           (dmvelz(i) - dmvelz(i+1)) / 6.d0
      velzminus(i+1) = velzplus(i)
      Bxplus(i) = 0.5d0 * (Bx(i) + Bx(i+1)) + &
           (dmBx(i) - dmBx(i+1)) / 6.d0
      Bxminus(i+1) = Bxplus(i)
      Byplus(i) = 0.5d0 * (By(i) + By(i+1)) + &
           (dmBy(i) - dmBy(i+1)) / 6.d0
      Byminus(i+1) = Byplus(i)
      Bzplus(i) = 0.5d0 * (Bz(i) + Bz(i+1)) + &
           (dmBz(i) - dmBz(i+1)) / 6.d0
      Bzminus(i+1) = Bzplus(i)
    end do
    if (poly .eq. 0) then
      do i = 2, nx-2
        epsplus(i) = 0.5d0 * (eps(i) + eps(i+1)) + &
             (dmeps(i) - dmeps(i+1)) / 6.d0
        epsminus(i+1) = epsplus(i)
      end do
    end if
  
  !!$Discontinuity steepening. See (1.14-17) of C&W.
  !!$This is the detect routine which mat be activated with the ppm_detect parameter
  !!$Note that this part really also depends on the grid being even. 
  !!$Note also that we don''t have access to the gas constant gamma.
  !!$So this is just dropped from eq. (3.2) of C&W.
  !!$We can get around this by just rescaling the constant k0 (ppm_k0 here).
  
    if (ppm_detect .ne. 0) then
  
      do i = 3, nx - 2
        if ( (d2rho(i+1)*d2rho(i-1) < 0.d0).and.(abs(rho(i+1)-rho(i-1)) - &
             ppm_epsilon_shock * min(abs(rho(i+1)), abs(rho(i-1))) > 0.d0) ) then
          etatilde = (rho(i-2) - rho(i+2) + 4.d0 * drho(i)) / (drho(i) * 12.d0)
        else
          etatilde = 0.d0
        end if
        eta = max(0.d0, min(1.d0, ppm_eta1 * (etatilde - ppm_eta2)))
        if (ppm_k0 * abs(drho(i)) * min(press(i-1),press(i+1)) < &
             abs(dpress(i)) * min(rho(i-1), rho(i+1))) then
          eta = 0.d0
        end if
        if (eta > 0.d0) then
          trivial_rp(i-1) = .false.
          trivial_rp(i) = .false.
        end if
        rhominus(i) = rhominus(i) * (1.d0 - eta) + &
             (rho(i-1) + 0.5d0 * dmrho(i-1)) * eta
        rhoplus(i) = rhoplus(i) * (1.d0 - eta) + &
             (rho(i+1) - 0.5d0 * dmrho(i+1)) * eta
      end do
  
    end if
  
    !!$ mppm
#define D_UPW(x) (0.5d0 * (x(i) + x(i+1)))
#define LEFT1(x)  (13.d0*x(i+1)-5.d0*x(i+2)+x(i+3)+3.d0*x(i  ))/12.d0
#define RIGHT1(x) (13.d0*x(i  )-5.d0*x(i-1)+x(i-2)+3.d0*x(i+1))/12.d0
    if (ppm_mppm .gt. 0) then
      l_ev_l=0.d0
      l_ev_r=0.d0
      xwind=0.d0
      do i=3, nx - 3
        agxx = 0.5d0*( gxx(i) + gxx(i+1) )
        agxy = 0.5d0*( gxy(i) + gxy(i+1) )
        agxz = 0.5d0*( gxz(i) + gxz(i+1) )
        agyy = 0.5d0*( gyy(i) + gyy(i+1) )
        agyz = 0.5d0*( gyz(i) + gyz(i+1) )
        agzz = 0.5d0*( gzz(i) + gzz(i+1) )
        det = SPATIAL_DET(agxx,agxy,agxz,agyy,agyz,agzz)
        call UpperMetric (uxx, uxy, uxz, uyy, uyz, uzz, &
                          det, agxx, agxy, agxz, agyy, agyz, agzz)
  !!$ these always have to be vel[xyz], not v[xyz] since they are used to
  !!$ calculate the eigenvalues
        call eigenvalues_ppm(handle,&
                   D_UPW(rho), D_UPW(velx), D_UPW(vely), D_UPW(velz), &
                   D_UPW(eps), D_UPW(w_lorentz), lam, &
                   agxx, agxy, agxz, agyy, agyz, agzz, &
                   uxx, det, D_UPW(alp), D_UPW(beta))
        l_ev_l(i)=lam(1)
        l_ev_r(i)=lam(5)
        xwind(i) = (lam(1) + lam(5)) / (abs(lam(1)) + abs(lam(5)))
        xwind(i) = min(1.d0, max(-1.d0, xwind(i)))
#define LEFTPLUS(x,xplus)    xplus(i)   =       abs(xwind(i))  * LEFT1(x) + \
                                          (1.d0-abs(xwind(i))) * xplus(i)
#define LEFTMINUS(x,xminus)  xminus(i+1)=       abs(xwind(i))  * LEFT1(x) + \
                                          (1.d0-abs(xwind(i))) * xminus(i+1)
#define RIGHTPLUS(x,xplus)   xplus(i)   =       abs(xwind(i))  * RIGHT1(x) + \
                                          (1.d0-abs(xwind(i))) * xplus(i)
#define RIGHTMINUS(x,xminus) xminus(i+1)=       abs(xwind(i))  * RIGHT1(x) + \
                                          (1.d0-abs(xwind(i))) * xminus(i+1)
#define CHECK(x,xc) if (x(i+1) .gt. x(i)) then && xc=min(x(i+1),max(x(i),xc)) && else && xc=min(x(i),max(x(i+1),xc)) && endif
  !!$      xwind(i)=0.d0
        if (xwind(i) .lt. 0.0d0) then
          LEFTPLUS(rho, rhoplus)
          LEFTMINUS(rho, rhominus)
          LEFTPLUS(vx, velxplus)
          LEFTMINUS(vx, velxminus)
          LEFTPLUS(vy, velyplus)
          LEFTMINUS(vy, velyminus)
          LEFTPLUS(vz, velzplus)
          LEFTMINUS(vz, velzminus)
          if (poly .eq. 0) then
            LEFTPLUS(eps, epsplus)
            LEFTMINUS(eps, epsminus)
          end if
        else
          RIGHTPLUS(rho, rhoplus)
          RIGHTMINUS(rho, rhominus)
          RIGHTPLUS(vx, velxplus)
          RIGHTMINUS(vx, velxminus)
          RIGHTPLUS(vy, velyplus)
          RIGHTMINUS(vy, velyminus)
          RIGHTPLUS(vz, velzplus)
          RIGHTMINUS(vz, velzminus)
          if (poly .eq. 0) then
            RIGHTPLUS(eps, epsplus)
            RIGHTMINUS(eps, epsminus)
          end if
        end if
        CHECK(rho, rhoplus(i))
        CHECK(rho, rhominus(i+1))
        CHECK(vx, velxplus(i))
        CHECK(vx, velxminus(i+1))
        CHECK(vy, velyplus(i))
        CHECK(vy, velyminus(i+1))
        CHECK(vz, velzplus(i))
        CHECK(vz, velzminus(i+1))
        if (poly .eq. 0) then
          CHECK(eps, epsplus(i))
          CHECK(eps, epsminus(i+1))
        end if
  !!$      if ((dir .eq. 1) .and. (ni .eq. 4) .and. (nj .eq. 4)) then
  !!$        write (*,*) rhoplus(i), rhominus(i+1)
  !!$      end if
      end do
      !!$ mppm debug output
      if (mppm_debug_eigenvalues .gt. 0) then
        if (dir .eq. 1) then
          ev_l(:,ni,nj) = l_ev_l
          ev_r(:,ni,nj) = l_ev_r
          xw(:,ni,nj) = xwind
        else if (dir .eq. 2) then
          ev_l(ni,:,nj) = l_ev_l
          ev_r(ni,:,nj) = l_ev_r
          xw(ni,:,nj) = xwind
        else if (dir .eq. 3) then
          ev_l(ni,nj,:) = l_ev_l
          ev_r(ni,nj,:) = l_ev_r
          xw(ni,nj,:) = xwind
        else
          write (*,*) "flux direction not 1 to 3 ?"
        end if
      end if
    end if
  
  !!$  Zone flattening. See appendix of C&W, p. 197-8.
  
    do i = 3, nx - 2
      dpress2 = press(i+2) - press(i-2)
      dvel = vx(i+1) - vx(i-1)
      if ( (abs(dpress(i)) >  ppm_epsilon * min(press(i-1),press(i+1))) .and. &
           (dvel < 0.d0) ) then
        w = 1.d0
      else
        w = 0.d0
      end if
      if (abs(dpress2) < ppm_small) then
        tilde_flatten(i) = 1.d0
      else
        tilde_flatten(i) = max(0.d0, 1.d0 - w * max(0.d0, ppm_omega2 * &
             (dpress(i) / dpress2 - ppm_omega1)))
      end if
    end do
  
  
  
    if (PPM3) then !!$ Implement C&W, page 197, but with a workaround which allows to use stencil=3.
      do i = 3, nx - 2
        flatten = tilde_flatten(i)
        if (abs(1.d0 - flatten) > 0.d0) then
          trivial_rp(i-1) = .false.
          trivial_rp(i) = .false.
        end if
        rhoplus(i) = flatten * rhoplus(i) + (1.d0 - flatten) * rho(i)
        rhominus(i) = flatten * rhominus(i) + (1.d0 - flatten) * rho(i)
        velxplus(i) = flatten * velxplus(i) + (1.d0 - flatten) * vx(i)
        velxminus(i) = flatten * velxminus(i) + (1.d0 - flatten) * vx(i)
        velyplus(i) = flatten * velyplus(i) + (1.d0 - flatten) * vy(i)
        velyminus(i) = flatten * velyminus(i) + (1.d0 - flatten) * vy(i)
        velzplus(i) = flatten * velzplus(i) + (1.d0 - flatten) * vz(i)
        velzminus(i) = flatten * velzminus(i) + (1.d0 - flatten) * vz(i)
        Bxplus(i) = flatten * Bxplus(i) + (1.d0 - flatten) * Bx(i)
        Bxminus(i) = flatten * Bxminus(i) + (1.d0 - flatten) * Bx(i)
        Byplus(i) = flatten * Byplus(i) + (1.d0 - flatten) * By(i)
        Byminus(i) = flatten * Byminus(i) + (1.d0 - flatten) * By(i)
        Bzplus(i) = flatten * Bzplus(i) + (1.d0 - flatten) * Bz(i)
        Bzminus(i) = flatten * Bzminus(i) + (1.d0 - flatten) * Bz(i)
        if (poly .eq. 0) then
          epsplus(i) = flatten * epsplus(i) + (1.d0 - flatten) * eps(i)
          epsminus(i) = flatten * epsminus(i) + (1.d0 - flatten) * eps(i)
        end if
      end do
    else  !!$ Really implement C&W, page 197; which requires stencil 4.
      do i = 4, nx - 3
        s=sign(1.d0, -dpress(i))
        flatten = max(tilde_flatten(i), tilde_flatten(i+s))  
        if (abs(1.d0 - flatten) > 0.d0) then
          trivial_rp(i-1) = .false.
          trivial_rp(i) = .false.
        end if
        rhoplus(i) = flatten * rhoplus(i) + (1.d0 - flatten) * rho(i)
        rhominus(i) = flatten * rhominus(i) + (1.d0 - flatten) * rho(i)
        velxplus(i) = flatten * velxplus(i) + (1.d0 - flatten) * vx(i)
        velxminus(i) = flatten * velxminus(i) + (1.d0 - flatten) * vx(i)
        velyplus(i) = flatten * velyplus(i) + (1.d0 - flatten) * vy(i)
        velyminus(i) = flatten * velyminus(i) + (1.d0 - flatten) * vy(i)
        velzplus(i) = flatten * velzplus(i) + (1.d0 - flatten) * vz(i)
        velzminus(i) = flatten * velzminus(i) + (1.d0 - flatten) * vz(i)
        Bxplus(i) = flatten * Bxplus(i) + (1.d0 - flatten) * Bx(i)
        Bxminus(i) = flatten * Bxminus(i) + (1.d0 - flatten) * Bx(i)
        Byplus(i) = flatten * Byplus(i) + (1.d0 - flatten) * By(i)
        Byminus(i) = flatten * Byminus(i) + (1.d0 - flatten) * By(i)
        Bzplus(i) = flatten * Bzplus(i) + (1.d0 - flatten) * Bz(i)
        Bzminus(i) = flatten * Bzminus(i) + (1.d0 - flatten) * Bz(i)
        if (poly .eq. 0) then
          epsplus(i) = flatten * epsplus(i) + (1.d0 - flatten) * eps(i)
          epsminus(i) = flatten * epsminus(i) + (1.d0 - flatten) * eps(i)
        end if
      end do
    end if
  
  
  !!$ Monotonicity. See (1.10) of C&W.
  
  do i = whisky_stencil, nx - whisky_stencil + 1
#define MON(xminus,x,xplus)                                       \
    if (.not.( (xplus(i).eq.x(i)) .and. (x(i).eq.xminus(i)) )     \
        .and. ((xplus(i)-x(i))*(x(i)-xminus(i)) .le. 0.d0)) then&&\
      trivial_rp(i-1) = .false.                                 &&\
      trivial_rp(i) = .false.                                   &&\
      xminus(i) = x(i)                                          &&\
      xplus(i) = x(i)                                           &&\
    else if (6.d0 * (xplus(i) - xminus(i)) * (x(i) - 0.5d0 *      \
                (xplus(i) + xminus(i))) >                         \
                (xplus(i) - xminus(i))**2) then                 &&\
      xminus(i) = 3.d0 * x(i) - 2.d0 * xplus(i)                 &&\
      trivial_rp(i-1) = .false.                                 &&\
      trivial_rp(i) = .false.                                   &&\
    else if (6.d0 * (xplus(i) - xminus(i)) * (x(i) - 0.5d0 *      \
                (xplus(i) + xminus(i))) <                         \
               -(xplus(i) - xminus(i))**2) then                 &&\
      xplus(i) = 3.d0 * x(i) - 2.d0 * xminus(i)                 &&\
      trivial_rp(i-1) = .false.                                 &&\
      trivial_rp(i) = .false.                                   &&\
    end if                                                      &&\
    if (.not.( (xplus(i).eq.x(i)) .and. (x(i).eq.xminus(i)) ) ) then     &&\
      trivial_rp(i-1) = .false.                                 &&\
      trivial_rp(i) = .false.                                   &&\
    end if

      MON(rhominus,rho,rhoplus)
      MON(velxminus,vx,velxplus)
      MON(velyminus,vy,velyplus)
      MON(velzminus,vz,velzplus)
      MON(Bxminus,Bx,Bxplus)
      MON(Byminus,By,Byplus)
      MON(Bzminus,Bz,Bzplus)
    end do
    if (poly .eq. 0) then
      do i = whisky_stencil, nx - whisky_stencil + 1
        MON(epsminus,eps,epsplus)
      end do
    end if
  
    if (check_for_trivial_rp .eq. 0) then
      trivial_rp = .false.
    end if
  
    !!$ excision
  
    !!$ recover primitive v^i from W*v^i at interfaces
    if (reconstruct_Wv.ne.0) then
      do i = whisky_stencil, nx - whisky_stencil + 1
        agxx = 0.5d0*( gxx(i) + gxx(i+1) )
        agxy = 0.5d0*( gxy(i) + gxy(i+1) )
        agxz = 0.5d0*( gxz(i) + gxz(i+1) )
        agyy = 0.5d0*( gyy(i) + gyy(i+1) )
        agyz = 0.5d0*( gyz(i) + gyz(i+1) )
        agzz = 0.5d0*( gzz(i) + gzz(i+1) )
        ! we re-use the variable w defined for zone flattening to hold the lorentz
        ! factor
        ! divide out the Loretnz factor for both the plus and minus quantities
        ! this should by construction ensure that any Lorentz factor calculated
        ! from them later on is physical (ie. > 1.d0)
        w = sqrt( 1.d0 + agxx*velxminus(i)*velxminus(i) + agyy*velyminus(i)*velyminus(i) &
                + agzz*velzminus(i)*velzminus(i) + 2.d0*agxy*velxminus(i)*velyminus(i) &
                + 2.d0*agxz*velxminus(i)*velzminus(i) + 2.d0*agyz*velyminus(i)*velzminus(i) )
        if (w.lt.1.d0 .or. w.ne.w) then ! .and. Whisky_CarpetWeights(i,j,k).eq.1.d0 (no CCTK_ARGUMENTS)
          write(warnline,'("SimplePPM_1dMHD (minus): Unphysical Lorentz factor ",g15.6,&
            &" for data g = [",5(g15.6,","),g15.6,"] vel = [",&
            &2(g15.6,","),g15.6,"]")') &
            w, agxx,agxy,agxz,agyy,agyz,agzz,&
            velxminus(i),velyminus(i),velzminus(i)
          call CCTK_WARN(1,warnline)
        end if 
        velxminus(i) = velxminus(i)/w
        velyminus(i) = velyminus(i)/w
        velzminus(i) = velzminus(i)/w
        w = sqrt( 1.d0 + agxx*velxplus(i)*velxplus(i) + agyy*velyplus(i)*velyplus(i) &
                + agzz*velzplus(i)*velzplus(i) + 2.d0*agxy*velxplus(i)*velyplus(i) &
                + 2.d0*agxz*velxplus(i)*velzplus(i) + 2.d0*agyz*velyplus(i)*velzplus(i) )
        if (w.lt.1.d0 .or. w.ne.w) then ! .and. Whisky_CarpetWeights(i,j,k).eq.1.d0 (no CCTK_ARGUMENTS)
          write(warnline,'("SimplePPM_1dMHD (plus): Unphysical Lorentz factor ",g15.6,&
            &" for data g = [",5(g15.6,","),g15.6,"] vel = [",&
            &2(g15.6,","),g15.6,"]")') &
            w, agxx,agxy,agxz,agyy,agyz,agzz,&
            velxplus(i),velyplus(i),velzplus(i)
          call CCTK_WARN(1,warnline)
        end if 
        velxplus(i) = velxplus(i)/w
        velyplus(i) = velyplus(i)/w
        velzplus(i) = velzplus(i)/w
      end do
    end if
  
    return
  
  end subroutine SimplePPM_1dMHD
  
   /*@@
     @routine    SimplePPM_1dMHD_wDC
     @date       Thu Feb 14 19:08:52 2002
     @author     Ian Hawke, Toni Font, Tanja Bode
     @desc 
     Extended reconstruction for MHD:
     The simple PPM reconstruction routine that applies along
     each one dimensional slice.  With MHD and Divergence cleaning.
  
     @enddesc 
     @calls     
     @calledby   
     @history 
     Written in frustration when IH couldn''t get Toni''s original code 
     to work.
     @endhistory 
  
  @@*/

  subroutine SimplePPM_1dMHD_wDC(handle,poly,&
       nx,dx,rho,velx,vely,velz,eps,press, &
       Bx, By, Bz, psi, &
       rhominus,velxminus,velyminus,velzminus,epsminus,&
       Bxminus, Byminus, Bzminus, psiminus,&
       rhoplus,velxplus,velyplus,velzplus,epsplus,&
       Bxplus, Byplus, Bzplus, psiplus,&
       trivial_rp,space_mask, excision_descriptors,&
       gxx, gxy, gxz, gyy, gyz, gzz, beta, alp, w_lorentz, &
       dir, ni, nj, nrx, nry, nrz, ev_l, ev_r, xw)
  
    USE PPM_Scalars
    USE PPM_EigenvalueCalcs
  
    implicit none
  
    DECLARE_CCTK_PARAMETERS
    DECLARE_CCTK_FUNCTIONS
  
#ifdef _EOS_BASE_INC_
#undef _EOS_BASE_INC_
#endif
#include "EOS_Base.inc"

    CCTK_INT :: handle,poly,nx
    CCTK_REAL :: dx
    CCTK_REAL, dimension(nx) :: rho,velx,vely,velz,eps
    CCTK_REAL, dimension(nx) :: vx,vy,vz
    CCTK_REAL, dimension(nx) :: Bx,By,Bz
    CCTK_REAL, dimension(nx) :: Bxminus,Byminus,Bzminus
    CCTK_REAL, dimension(nx) :: Bxplus,Byplus,Bzplus
    CCTK_REAL, dimension(nx) :: rhominus,velxminus,velyminus,velzminus,epsminus
    CCTK_REAL, dimension(nx) :: rhoplus,velxplus,velyplus,velzplus,epsplus
    CCTK_REAL, dimension(nx) :: rhominusl,velxminusl,velyminusl,velzminusl
    CCTK_REAL, dimension(nx) :: epsminusl
    CCTK_REAL, dimension(nx) :: rhoplusl,velxplusl,velyplusl,velzplusl,epsplusl
    CCTK_REAL, dimension(nx) :: rhominusr,velxminusr,velyminusr,velzminusr
    CCTK_REAL, dimension(nx) :: epsminusr
    CCTK_REAL, dimension(nx) :: rhoplusr,velxplusr,velyplusr,velzplusr,epsplusr
    CCTK_REAL, dimension(nx) :: psi, psiplus, psiminus
  
    CCTK_INT :: i,s
    CCTK_REAL, dimension(nx) :: drho,dvelx,dvely,dvelz,deps
    CCTK_REAL, dimension(nx) :: dBx,dBy,dBz
    CCTK_REAL, dimension(nx) :: dmrho,dmvelx,dmvely,dmvelz,dmeps
    CCTK_REAL, dimension(nx) :: dmBx,dmBy,dmBz
    CCTK_REAL, dimension(nx) :: press,dpress,d2rho,tilde_flatten
    CCTK_REAL, dimension(nx) :: dpsi,dmpsi
    CCTK_REAL :: dpress2,dvel,w,flatten,eta,etatilde
  
    logical, dimension(nx) :: trivial_rp
  
    CCTK_INT, dimension(nx) :: space_mask
    CCTK_INT :: excision_bits, excision_mask
    CCTK_INT, dimension(3) :: excision_descriptors
  
    CCTK_REAL, dimension(nx) :: gxx, gxy, gxz, gyy, gyz, gzz, &
                                beta, alp, w_lorentz
    CCTK_INT :: dir, ni, nj, nrx, nry, nrz
    CCTK_REAL, dimension(nrx, nry, nrz) :: ev_l, ev_r, xw
  
    CCTK_REAL :: uxx, uxy, uxz, uyy, uyz, uzz, det
    CCTK_REAL, dimension(8) :: lam
    CCTK_REAL :: dupw, dloc, delta
    CCTK_REAL :: agxx, agxy, agxz, agyy, agyz, agzz
    CCTK_REAL, dimension(nx) :: xwind, l_ev_l, l_ev_r
  
    logical :: cond
    character(len=400) :: warnline
    CCTK_REAL :: B2, B2plus, B2minus
  
  
  !!$  Initially, all the Riemann problems will be trivial
  
  trivial_rp = .true.
  
  !!$  Loop to set the pressure and eps
  !!$  poly=0
  !!$  if (poly .eq. 0) then
  !!$    do i = 1, nx
  !!$      press(i) = EOS_Pressure(handle, rho(i), eps(i))
  !!$    end do
  !!$  else
  !!$    do i = 1, nx
  !!$      press(i) = EOS_Pressure(handle, rho(i), eps(i))
  !!$      eps(i) = EOS_SpecificIntEnergy(handle, rho(i), press(i))
  !!$    end do
  !!$  end if
  
  !!$ reconstruct w^i = w_lorentz*v^i to ensure slower than light speeds?
    if (reconstruct_Wv.ne.0) then
      vx = w_lorentz*velx
      vy = w_lorentz*vely
      vz = w_lorentz*velz
    else
      vx = velx
      vy = vely
      vz = velz
    end if
  
  !!$  Average slopes delta_m(a). See (1.7) of Colella and Woodward, p.178
  !!$  This is the expression for an even grid.
  
    do i = 2, nx - 1
      drho(i) = 0.5d0 * (rho(i+1) - rho(i-1))
      dvelx(i) = 0.5d0 * (vx(i+1) - vx(i-1))
      dvely(i) = 0.5d0 * (vy(i+1) - vy(i-1))
      dvelz(i) = 0.5d0 * (vz(i+1) - vz(i-1))
      dBx(i) = 0.5d0 * (Bx(i+1) - Bx(i-1))
      dBy(i) = 0.5d0 * (By(i+1) - By(i-1))
      dBz(i) = 0.5d0 * (Bz(i+1) - Bz(i-1))
      dpsi(i) = 0.5d0 * (psi(i+1) - psi(i-1))
      dpress(i) = press(i+1) - press(i-1)
      d2rho(i) = (rho(i+1) - 2.d0 * rho(i) + rho(i-1))! / 6.d0 / dx / dx
      ! since we use d2rho only for the condition d2rho(i+1)*d2rhoi(i-1)<0 
      ! the denominator is not necessary
    end do
    if (poly .eq. 0) then
      do i = 2, nx - 1
        deps(i) = 0.5d0 * (eps(i+1) - eps(i-1))
      end do
    end if
  
  !!$  Steepened slope. See (1.8) of Colella and Woodward, p.178
  
    do i = 2, nx - 1
#define STEEP(x,dx,dmx)                                              \
      if ( (x(i+1) - x(i)) * (x(i) - x(i-1)) > 0.d0 ) then           &&\
        dmx(i) = sign(1.d0, dx(i)) *                                   \
             min(abs(dx(i)), 2.d0 * abs(x(i) - x(i-1)),                \
                             2.d0 * abs(x(i+1) - x(i)))              &&\
      else                                                           &&\
        dmx(i) = 0.d0                                                &&\
      end if
      STEEP(rho, drho, dmrho)
      STEEP(vx, dvelx, dmvelx)
      STEEP(vy, dvely, dmvely)
      STEEP(vz, dvelz, dmvelz)
      STEEP(Bx, dBx, dmBx)
      STEEP(By, dBy, dmBy)
      STEEP(Bz, dBz, dmBz)
      STEEP(psi, dpsi, dmpsi)
    end do
    if (poly .eq. 0) then
      do i = 2, nx - 1
        STEEP(eps, deps, dmeps)
      end do
    end if
  
  !!$  Initial boundary states. See (1.9) of Colella and Woodward, p.178
  
    do i = 2, nx-2
      rhoplus(i) = 0.5d0 * (rho(i) + rho(i+1)) + &
           (dmrho(i) - dmrho(i+1)) / 6.d0
      rhominus(i+1) = rhoplus(i)
      velxplus(i) = 0.5d0 * (vx(i) + vx(i+1)) + &
           (dmvelx(i) - dmvelx(i+1)) / 6.d0
      velxminus(i+1) = velxplus(i)
      velyplus(i) = 0.5d0 * (vy(i) + vy(i+1)) + &
           (dmvely(i) - dmvely(i+1)) / 6.d0
      velyminus(i+1) = velyplus(i)
      velzplus(i) = 0.5d0 * (vz(i) + vz(i+1)) + &
           (dmvelz(i) - dmvelz(i+1)) / 6.d0
      velzminus(i+1) = velzplus(i)
      Bxplus(i) = 0.5d0 * (Bx(i) + Bx(i+1)) + &
           (dmBx(i) - dmBx(i+1)) / 6.d0
      Bxminus(i+1) = Bxplus(i)
      Byplus(i) = 0.5d0 * (By(i) + By(i+1)) + &
           (dmBy(i) - dmBy(i+1)) / 6.d0
      Byminus(i+1) = Byplus(i)
      Bzplus(i) = 0.5d0 * (Bz(i) + Bz(i+1)) + &
           (dmBz(i) - dmBz(i+1)) / 6.d0
      Bzminus(i+1) = Bzplus(i)
      psiplus(i) = 0.5d0 * (psi(i) + psi(i+1)) + &
           (dmpsi(i) - dmpsi(i+1)) / 6.d0
      psiminus(i+1) = psiplus(i)
    end do
    if (poly .eq. 0) then
      do i = 2, nx-2
        epsplus(i) = 0.5d0 * (eps(i) + eps(i+1)) + &
             (dmeps(i) - dmeps(i+1)) / 6.d0
        epsminus(i+1) = epsplus(i)
      end do
    end if
  
  !!$Discontinuity steepening. See (1.14-17) of C&W.
  !!$This is the detect routine which mat be activated with the ppm_detect parameter
  !!$Note that this part really also depends on the grid being even. 
  !!$Note also that we don''t have access to the gas constant gamma.
  !!$So this is just dropped from eq. (3.2) of C&W.
  !!$We can get around this by just rescaling the constant k0 (ppm_k0 here).
  
    if (ppm_detect .ne. 0) then
  
      do i = 3, nx - 2
        if ( (d2rho(i+1)*d2rho(i-1) < 0.d0).and.(abs(rho(i+1)-rho(i-1)) - &
             ppm_epsilon_shock * min(abs(rho(i+1)), abs(rho(i-1))) > 0.d0) ) then
          etatilde = (rho(i-2) - rho(i+2) + 4.d0 * drho(i)) / (drho(i) * 12.d0)
        else
          etatilde = 0.d0
        end if
        eta = max(0.d0, min(1.d0, ppm_eta1 * (etatilde - ppm_eta2)))
        if (ppm_k0 * abs(drho(i)) * min(press(i-1),press(i+1)) < &
             abs(dpress(i)) * min(rho(i-1), rho(i+1))) then
          eta = 0.d0
        end if
        if (eta > 0.d0) then
          trivial_rp(i-1) = .false.
          trivial_rp(i) = .false.
        end if
        rhominus(i) = rhominus(i) * (1.d0 - eta) + &
             (rho(i-1) + 0.5d0 * dmrho(i-1)) * eta
        rhoplus(i) = rhoplus(i) * (1.d0 - eta) + &
             (rho(i+1) - 0.5d0 * dmrho(i+1)) * eta
      end do
  
    end if
  
    !!$ mppm
#define D_UPW(x) (0.5d0 * (x(i) + x(i+1)))
#define LEFT1(x)  (13.d0*x(i+1)-5.d0*x(i+2)+x(i+3)+3.d0*x(i  ))/12.d0
#define RIGHT1(x) (13.d0*x(i  )-5.d0*x(i-1)+x(i-2)+3.d0*x(i+1))/12.d0
    if (ppm_mppm .gt. 0) then
      l_ev_l=0.d0
      l_ev_r=0.d0
      xwind=0.d0
      do i=3, nx - 3
        agxx = 0.5d0*( gxx(i) + gxx(i+1) )
        agxy = 0.5d0*( gxy(i) + gxy(i+1) )
        agxz = 0.5d0*( gxz(i) + gxz(i+1) )
        agyy = 0.5d0*( gyy(i) + gyy(i+1) )
        agyz = 0.5d0*( gyz(i) + gyz(i+1) )
        agzz = 0.5d0*( gzz(i) + gzz(i+1) )
        det=SPATIAL_DET(agxx,agxy,agxz,agyy,agyz,agzz)
        call UpperMetric (uxx, uxy, uxz, uyy, uyz, uzz, &
                          det, agxx, agxy, agxz, agyy, agyz, agzz)
  !!$ these always have to be vel[xyz], not v[xyz] since they are used to
  !!$ calculate the eigenvalues
        call eigenvalues_ppm(handle,&
                   D_UPW(rho), D_UPW(velx), D_UPW(vely), D_UPW(velz), &
                   D_UPW(eps), D_UPW(w_lorentz), lam, &
                   agxx, agxy, agxz, agyy, agyz, agzz, &
                   uxx, det, D_UPW(alp), D_UPW(beta))
        l_ev_l(i)=lam(1)
        l_ev_r(i)=lam(5)
        xwind(i) = (lam(1) + lam(5)) / (abs(lam(1)) + abs(lam(5)))
        xwind(i) = min(1.d0, max(-1.d0, xwind(i)))
#define LEFTPLUS(x,xplus)    xplus(i)   =       abs(xwind(i))  * LEFT1(x) + \
                                          (1.d0-abs(xwind(i))) * xplus(i)
#define LEFTMINUS(x,xminus)  xminus(i+1)=       abs(xwind(i))  * LEFT1(x) + \
                                          (1.d0-abs(xwind(i))) * xminus(i+1)
#define RIGHTPLUS(x,xplus)   xplus(i)   =       abs(xwind(i))  * RIGHT1(x) + \
                                          (1.d0-abs(xwind(i))) * xplus(i)
#define RIGHTMINUS(x,xminus) xminus(i+1)=       abs(xwind(i))  * RIGHT1(x) + \
                                          (1.d0-abs(xwind(i))) * xminus(i+1)
#define CHECK(x,xc) if (x(i+1) .gt. x(i)) then && xc=min(x(i+1),max(x(i),xc)) && else && xc=min(x(i),max(x(i+1),xc)) && endif
  !!$      xwind(i)=0.d0
        if (xwind(i) .lt. 0.0d0) then
          LEFTPLUS(rho, rhoplus)
          LEFTMINUS(rho, rhominus)
          LEFTPLUS(vx, velxplus)
          LEFTMINUS(vx, velxminus)
          LEFTPLUS(vy, velyplus)
          LEFTMINUS(vy, velyminus)
          LEFTPLUS(vz, velzplus)
          LEFTMINUS(vz, velzminus)
          if (poly .eq. 0) then
            LEFTPLUS(eps, epsplus)
            LEFTMINUS(eps, epsminus)
          end if
        else
          RIGHTPLUS(rho, rhoplus)
          RIGHTMINUS(rho, rhominus)
          RIGHTPLUS(vx, velxplus)
          RIGHTMINUS(vx, velxminus)
          RIGHTPLUS(vy, velyplus)
          RIGHTMINUS(vy, velyminus)
          RIGHTPLUS(vz, velzplus)
          RIGHTMINUS(vz, velzminus)
          if (poly .eq. 0) then
            RIGHTPLUS(eps, epsplus)
            RIGHTMINUS(eps, epsminus)
          end if
        end if
        CHECK(rho, rhoplus(i))
        CHECK(rho, rhominus(i+1))
        CHECK(vx, velxplus(i))
        CHECK(vx, velxminus(i+1))
        CHECK(vy, velyplus(i))
        CHECK(vy, velyminus(i+1))
        CHECK(vz, velzplus(i))
        CHECK(vz, velzminus(i+1))
        if (poly .eq. 0) then
          CHECK(eps, epsplus(i))
          CHECK(eps, epsminus(i+1))
        end if
  !!$      if ((dir .eq. 1) .and. (ni .eq. 4) .and. (nj .eq. 4)) then
  !!$        write (*,*) rhoplus(i), rhominus(i+1)
  !!$      end if
      end do
      !!$ mppm debug output
      if (mppm_debug_eigenvalues .gt. 0) then
        if (dir .eq. 1) then
          ev_l(:,ni,nj) = l_ev_l
          ev_r(:,ni,nj) = l_ev_r
          xw(:,ni,nj) = xwind
        else if (dir .eq. 2) then
          ev_l(ni,:,nj) = l_ev_l
          ev_r(ni,:,nj) = l_ev_r
          xw(ni,:,nj) = xwind
        else if (dir .eq. 3) then
          ev_l(ni,nj,:) = l_ev_l
          ev_r(ni,nj,:) = l_ev_r
          xw(ni,nj,:) = xwind
        else
          write (*,*) "flux direction not 1 to 3 ?"
        end if
      end if
    end if
  
  !!$  Zone flattening. See appendix of C&W, p. 197-8.
  
    do i = 3, nx - 2
      dpress2 = press(i+2) - press(i-2)
      dvel = vx(i+1) - vx(i-1)
      if ( (abs(dpress(i)) >  ppm_epsilon * min(press(i-1),press(i+1))) .and. &
           (dvel < 0.d0) ) then
        w = 1.d0
      else
        w = 0.d0
      end if
      if (abs(dpress2) < ppm_small) then
        tilde_flatten(i) = 1.d0
      else
        tilde_flatten(i) = max(0.d0, 1.d0 - w * max(0.d0, ppm_omega2 * &
             (dpress(i) / dpress2 - ppm_omega1)))
      end if
    end do
  
  
  
    if (PPM3) then !!$ Implement C&W, page 197, but with a workaround which allows to use stencil=3.
      do i = 3, nx - 2
        flatten = tilde_flatten(i)
        if (abs(1.d0 - flatten) > 0.d0) then
          trivial_rp(i-1) = .false.
          trivial_rp(i) = .false.
        end if
        rhoplus(i) = flatten * rhoplus(i) + (1.d0 - flatten) * rho(i)
        rhominus(i) = flatten * rhominus(i) + (1.d0 - flatten) * rho(i)
        velxplus(i) = flatten * velxplus(i) + (1.d0 - flatten) * vx(i)
        velxminus(i) = flatten * velxminus(i) + (1.d0 - flatten) * vx(i)
        velyplus(i) = flatten * velyplus(i) + (1.d0 - flatten) * vy(i)
        velyminus(i) = flatten * velyminus(i) + (1.d0 - flatten) * vy(i)
        velzplus(i) = flatten * velzplus(i) + (1.d0 - flatten) * vz(i)
        velzminus(i) = flatten * velzminus(i) + (1.d0 - flatten) * vz(i)
        Bxplus(i) = flatten * Bxplus(i) + (1.d0 - flatten) * Bx(i)
        Bxminus(i) = flatten * Bxminus(i) + (1.d0 - flatten) * Bx(i)
        Byplus(i) = flatten * Byplus(i) + (1.d0 - flatten) * By(i)
        Byminus(i) = flatten * Byminus(i) + (1.d0 - flatten) * By(i)
        Bzplus(i) = flatten * Bzplus(i) + (1.d0 - flatten) * Bz(i)
        Bzminus(i) = flatten * Bzminus(i) + (1.d0 - flatten) * Bz(i)
        psiplus(i) = flatten * psiplus(i) + (1.d0 - flatten) * psi(i)
        psiminus(i) = flatten * psiminus(i) + (1.d0 - flatten) * psi(i)
        if (poly .eq. 0) then
          epsplus(i) = flatten * epsplus(i) + (1.d0 - flatten) * eps(i)
          epsminus(i) = flatten * epsminus(i) + (1.d0 - flatten) * eps(i)
        end if
      end do
    else  !!$ Really implement C&W, page 197; which requires stencil 4.
      do i = 4, nx - 3
        s=sign(1.d0, -dpress(i))
        flatten = max(tilde_flatten(i), tilde_flatten(i+s))  
        if (abs(1.d0 - flatten) > 0.d0) then
          trivial_rp(i-1) = .false.
          trivial_rp(i) = .false.
        end if
        rhoplus(i) = flatten * rhoplus(i) + (1.d0 - flatten) * rho(i)
        rhominus(i) = flatten * rhominus(i) + (1.d0 - flatten) * rho(i)
        velxplus(i) = flatten * velxplus(i) + (1.d0 - flatten) * vx(i)
        velxminus(i) = flatten * velxminus(i) + (1.d0 - flatten) * vx(i)
        velyplus(i) = flatten * velyplus(i) + (1.d0 - flatten) * vy(i)
        velyminus(i) = flatten * velyminus(i) + (1.d0 - flatten) * vy(i)
        velzplus(i) = flatten * velzplus(i) + (1.d0 - flatten) * vz(i)
        velzminus(i) = flatten * velzminus(i) + (1.d0 - flatten) * vz(i)
        Bxplus(i) = flatten * Bxplus(i) + (1.d0 - flatten) * Bx(i)
        Bxminus(i) = flatten * Bxminus(i) + (1.d0 - flatten) * Bx(i)
        Byplus(i) = flatten * Byplus(i) + (1.d0 - flatten) * By(i)
        Byminus(i) = flatten * Byminus(i) + (1.d0 - flatten) * By(i)
        Bzplus(i) = flatten * Bzplus(i) + (1.d0 - flatten) * Bz(i)
        Bzminus(i) = flatten * Bzminus(i) + (1.d0 - flatten) * Bz(i)
        psiplus(i) = flatten * psiplus(i) + (1.d0 - flatten) * psi(i)
        psiminus(i) = flatten * psiminus(i) + (1.d0 - flatten) * psi(i)
        if (poly .eq. 0) then
          epsplus(i) = flatten * epsplus(i) + (1.d0 - flatten) * eps(i)
          epsminus(i) = flatten * epsminus(i) + (1.d0 - flatten) * eps(i)
        end if
      end do
    end if
  
  
  !!$ Monotonicity. See (1.10) of C&W.
  
  do i = whisky_stencil, nx - whisky_stencil + 1
#define MON(xminus,x,xplus)                                       \
    if (.not.( (xplus(i).eq.x(i)) .and. (x(i).eq.xminus(i)) )     \
        .and. ((xplus(i)-x(i))*(x(i)-xminus(i)) .le. 0.d0)) then&&\
      trivial_rp(i-1) = .false.                                 &&\
      trivial_rp(i) = .false.                                   &&\
      xminus(i) = x(i)                                          &&\
      xplus(i) = x(i)                                           &&\
    else if (6.d0 * (xplus(i) - xminus(i)) * (x(i) - 0.5d0 *      \
                (xplus(i) + xminus(i))) >                         \
                (xplus(i) - xminus(i))**2) then                 &&\
      xminus(i) = 3.d0 * x(i) - 2.d0 * xplus(i)                 &&\
      trivial_rp(i-1) = .false.                                 &&\
      trivial_rp(i) = .false.                                   &&\
    else if (6.d0 * (xplus(i) - xminus(i)) * (x(i) - 0.5d0 *      \
                (xplus(i) + xminus(i))) <                         \
               -(xplus(i) - xminus(i))**2) then                 &&\
      xplus(i) = 3.d0 * x(i) - 2.d0 * xminus(i)                 &&\
      trivial_rp(i-1) = .false.                                 &&\
      trivial_rp(i) = .false.                                   &&\
    end if                                                      &&\
    if (.not.( (xplus(i).eq.x(i)) .and. (x(i).eq.xminus(i)) ) ) then     &&\
      trivial_rp(i-1) = .false.                                 &&\
      trivial_rp(i) = .false.                                   &&\
    end if

      MON(rhominus,rho,rhoplus)
      MON(velxminus,vx,velxplus)
      MON(velyminus,vy,velyplus)
      MON(velzminus,vz,velzplus)
      MON(Bxminus,Bx,Bxplus)
      MON(Byminus,By,Byplus)
      MON(Bzminus,Bz,Bzplus)
      MON(psiminus,psi,psiplus)
    end do
    if (poly .eq. 0) then
      do i = whisky_stencil, nx - whisky_stencil + 1
        MON(epsminus,eps,epsplus)
      end do
    end if
  
    if (check_for_trivial_rp .eq. 0) then
      trivial_rp = .false.
    end if
  
    !!$ excision -- check Bfield here.
    B2 = gxx(i)*Bx(i)**2 + gyy(i)*By(i)**2 + gzz(i)*Bz(i)**2 &
             + 2.d0*gxy(i)*Bx(i)*By(i) + 2.d0*gxz(i)*Bx(i)*Bz(i) &
             + 2.d0*gyz(i)*By(i)*Bz(i)
    B2plus = gxx(i)*Bxplus(i)**2 + gyy(i)*Byplus(i)**2 + gzz(i)*Bzplus(i)**2 &
             + 2.d0*gxy(i)*Bxplus(i)*Byplus(i) + 2.d0*gxz(i)*Bxplus(i)*Bzplus(i) &
             + 2.d0*gyz(i)*Byplus(i)*Bzplus(i)
    B2minus = gxx(i)*Bxminus(i)**2 + gyy(i)*Byminus(i)**2 + gzz(i)*Bzminus(i)**2 &
             + 2.d0*gxy(i)*Bxminus(i)*Byminus(i) + 2.d0*gxz(i)*Bxminus(i)*Bzminus(i) &
             + 2.d0*gyz(i)*Byminus(i)*Bzminus(i)
    if ( B2plus .lt. 0.d0 ) then
       Bxplus(i) = Bx(i)
       Byplus(i) = By(i)
       Bzplus(i) = Bz(i)
    end if
    if ( B2minus .lt. 0.d0 ) then
       Bxminus(i) = Bx(i)
       Byminus(i) = By(i)
       Bzminus(i) = Bz(i)
    end if
    if ( B2 .lt. 0.d0 ) then
       call CCTK_WARN(1,'B2 is negative even at the cell center ... ')
    end if
  
    !!$ recover primitive v^i from W*v^i at interfaces
    if (reconstruct_Wv.ne.0) then
      do i = whisky_stencil, nx - whisky_stencil + 1
        agxx = 0.5d0*( gxx(i) + gxx(i+1) )
        agxy = 0.5d0*( gxy(i) + gxy(i+1) )
        agxz = 0.5d0*( gxz(i) + gxz(i+1) )
        agyy = 0.5d0*( gyy(i) + gyy(i+1) )
        agyz = 0.5d0*( gyz(i) + gyz(i+1) )
        agzz = 0.5d0*( gzz(i) + gzz(i+1) )
        ! we re-use the variable w defined for zone flattening to hold the lorentz
        ! factor
        ! divide out the Loretnz factor for both the plus and minus quantities
        ! this should by construction ensure that any Lorentz factor calculated
        ! from them later on is physical (ie. > 1.d0)
        w = sqrt( 1.d0 + agxx*velxminus(i)*velxminus(i) + agyy*velyminus(i)*velyminus(i) &
                + agzz*velzminus(i)*velzminus(i) + 2.d0*agxy*velxminus(i)*velyminus(i) &
                + 2.d0*agxz*velxminus(i)*velzminus(i) + 2.d0*agyz*velyminus(i)*velzminus(i) )
        if (w.lt.1.d0 .or. w.ne.w) then ! .and. Whisky_CarpetWeights(i,j,k).eq.1.d0 (no CCTK_ARGUMENTS)
          write(warnline,'("SimplePPM_1dMHD_wDC (minus): Unphysical Lorentz factor ",g15.6,&
            &" for data g = [",5(g15.6,","),g15.6,"] vel = [",&
            &2(g15.6,","),g15.6,"]")') &
            w, agxx,agxy,agxz,agyy,agyz,agzz,&
            velxminus(i),velyminus(i),velzminus(i)
          call CCTK_WARN(1,warnline)
        end if 
        velxminus(i) = velxminus(i)/w
        velyminus(i) = velyminus(i)/w
        velzminus(i) = velzminus(i)/w
        w = sqrt( 1.d0 + agxx*velxplus(i)*velxplus(i) + agyy*velyplus(i)*velyplus(i) &
                + agzz*velzplus(i)*velzplus(i) + 2.d0*agxy*velxplus(i)*velyplus(i) &
                + 2.d0*agxz*velxplus(i)*velzplus(i) + 2.d0*agyz*velyplus(i)*velzplus(i) )
        if (w.lt.1.d0 .or. w.ne.w) then ! .and. Whisky_CarpetWeights(i,j,k).eq.1.d0 (no CCTK_ARGUMENTS)
          write(warnline,'("SimplePPM_1dMHD_wDC (plus): Unphysical Lorentz factor ",g15.6,&
            &" for data g = [",5(g15.6,","),g15.6,"] vel = [",&
            &2(g15.6,","),g15.6,"]")') &
            w, agxx,agxy,agxz,agyy,agyz,agzz,&
            velxplus(i),velyplus(i),velzplus(i)
          call CCTK_WARN(1,warnline)
        end if 
        velxplus(i) = velxplus(i)/w
        velyplus(i) = velyplus(i)/w
        velzplus(i) = velzplus(i)/w
      end do
    end if
  
    return
  
  end subroutine SimplePPM_1dMHD_wDC

end module PPM_Simple1d_Routines

