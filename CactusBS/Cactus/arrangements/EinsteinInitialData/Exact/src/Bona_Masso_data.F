c This routine calculates Bona-Masso initial data, making use of the
c subroutine Exact__metric() to calculate the spacetime metric and its
c inverse.  Note that this use of the Bona-Masso variables is independent
c of how (or even if) we are evolving the Einstein equations -- here
c the Bona-Masso variables are "just" used as intermediate variables.

c $Header$

#include "cctk.h"
#include "cctk_Parameters.h"

      subroutine  Exact__Bona_Masso_data(
     $     decoded_exact_model,
     $     x, y, z, t,
     $     gxx, gyy, gzz, gxy, gyz, gxz,
     $     hxx, hyy, hzz, hxy, hyz, hxz,
     $     psi, psix, psiy, psiz,
     $     psixx, psiyy, psizz, psixy, psiyz, psixz,
     $     dxgxx, dxgyy, dxgzz, dxgxy, dxgyz, dxgxz,
     $     dygxx, dygyy, dygzz, dygxy, dygyz, dygxz,
     $     dzgxx, dzgyy, dzgzz, dzgxy, dzgyz, dzgxz,
     $     alp, dtalp, ax, ay, az,
     $     betax, betay, betaz, dtbetax, dtbetay, dtbetaz,
     $     bxx, bxy, bxz, byx, byy, byz, bzx, bzy, bzz)

      implicit none
      CCTK_INT, parameter :: izero = 0
      integer, parameter :: ik = kind(izero)
      CCTK_INT decoded_exact_model
      logical psi_on
      CCTK_REAL x, y, z, t,
     $       gxx, gyy, gzz, gxy, gyz, gxz,
     $       hxx, hyy, hzz, hxy, hyz, hxz,
     $       psi, psix, psiy, psiz,
     $       psixx, psiyy, psizz, psixy, psiyz, psixz,
     $       dxgxx, dxgyy, dxgzz, dxgxy, dxgyz, dxgxz,
     $       dygxx, dygyy, dygzz, dygxy, dygyz, dygxz,
     $       dzgxx, dzgyy, dzgzz, dzgxy, dzgyz, dzgxz,
     $       alp, dtalp, ax, ay, az,
     $       betax, betay, betaz, dtbetax, dtbetay, dtbetaz,
     $       bxx, bxy, bxz, byx, byy, byz, bzx, bzy, bzz

C     gxx is g_xx etc.
C     hxx is K_xx etc.
C     dxgyy is (/2) dg_yy / dx (sic!)
C     alp is N, betax is N^x etc.
C     bxy is (/2) dN^y / dx (sic and sic!)
C     ax is dN / dx / N (sic!)

      CCTK_REAL
     $       gdtt, gdtx, gdty, gdtz, 
     $       gutt, gutx, guty, gutz, 
     $       guxx, guyy, guzz, guxy, guyz, guxz,
     $       dgdtt, dgdtx, dgdty, dgdtz,
     $       dgdxx, dgdyy, dgdzz, dgdxy, dgdyz, dgdxz,
     $       dgutt, dgutx, dguty, dgutz,
     $       dguxx, dguyy, dguzz, dguxy, dguyz, dguxz,
     $       dpsi

C     Save, if we have to provide the conformal factor

      psi_on = psi .gt. 0.0d0

C     Get the spacetime metric and its inverse at the base point.

      call Exact__metric(
     $     decoded_exact_model,
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz, 
     $     gxx, gyy, gzz, gxy, gyz, gxz,
     $     gutt, gutx, guty, gutz, 
     $     guxx, guyy, guzz, guxy, guyz, guxz,
     $     psi)

C     Calculate lapse and shift from the upper metric.

      alp = 1.d0 / sqrt(-gutt)

      betax = - gutx / gutt
      betay = - guty / gutt
      betaz = - gutz / gutt

C     In order to calculate the derivatives of the lapse and shift from
C     the contravariant metric, use g^tt = -1/N^2 and g^i0 = N^i / N^2
C     Note that ax is dN/dx / N and that bxy is 1/2 dN^y / dx.

C     Calculate x-derivatives. 

      call Exact__metric_deriv(
     $     decoded_exact_model,
     $     1_ik,
     $     x, y, z, t,
     $     dgdtt, dgdtx, dgdty, dgdtz,
     $     dgdxx, dgdyy, dgdzz, dgdxy, dgdyz, dgdxz,
     $     dgutt, dgutx, dguty, dgutz,
     $     dguxx, dguyy, dguzz, dguxy, dguyz, dguxz,
     $     dpsi)

      dxgxx = 0.5d0 * dgdxx
      dxgyy = 0.5d0 * dgdyy
      dxgzz = 0.5d0 * dgdzz
      dxgxy = 0.5d0 * dgdxy
      dxgyz = 0.5d0 * dgdyz
      dxgxz = 0.5d0 * dgdxz

      ax = - 0.5d0 * dgutt / gutt

      bxx = 0.5d0 * (- dgutx * gutt + gutx * dgutt) / gutt**2
      bxy = 0.5d0 * (- dguty * gutt + guty * dgutt) / gutt**2
      bxz = 0.5d0 * (- dgutz * gutt + gutz * dgutt) / gutt**2

      if (psi_on) then
         psix = dpsi
      end if

C     Calculate y-derivatives.

      call Exact__metric_deriv(
     $     decoded_exact_model,
     $     2_ik,
     $     x, y, z, t,
     $     dgdtt, dgdtx, dgdty, dgdtz,
     $     dgdxx, dgdyy, dgdzz, dgdxy, dgdyz, dgdxz,
     $     dgutt, dgutx, dguty, dgutz,
     $     dguxx, dguyy, dguzz, dguxy, dguyz, dguxz,
     $     dpsi)

      dygxx = 0.5d0 * dgdxx
      dygyy = 0.5d0 * dgdyy
      dygzz = 0.5d0 * dgdzz
      dygxy = 0.5d0 * dgdxy
      dygyz = 0.5d0 * dgdyz
      dygxz = 0.5d0 * dgdxz

      ay = - 0.5d0 * dgutt / gutt

      byx = 0.5d0 * (- dgutx * gutt + gutx * dgutt) / gutt**2
      byy = 0.5d0 * (- dguty * gutt + guty * dgutt) / gutt**2
      byz = 0.5d0 * (- dgutz * gutt + gutz * dgutt) / gutt**2

      if (psi_on) then
         psiy = dpsi
      end if

C     Calculate z-derivatives.

      call Exact__metric_deriv(
     $     decoded_exact_model,
     $     3_ik,
     $     x, y, z, t,
     $     dgdtt, dgdtx, dgdty, dgdtz,
     $     dgdxx, dgdyy, dgdzz, dgdxy, dgdyz, dgdxz,
     $     dgutt, dgutx, dguty, dgutz,
     $     dguxx, dguyy, dguzz, dguxy, dguyz, dguxz,
     $     dpsi)

      dzgxx = 0.5d0 * dgdxx
      dzgyy = 0.5d0 * dgdyy
      dzgzz = 0.5d0 * dgdzz
      dzgxy = 0.5d0 * dgdxy
      dzgyz = 0.5d0 * dgdyz
      dzgxz = 0.5d0 * dgdxz

      az = - 0.5d0 * dgutt / gutt

      bzx = 0.5d0 * (- dgutx * gutt + gutx * dgutt) / gutt**2
      bzy = 0.5d0 * (- dguty * gutt + guty * dgutt) / gutt**2
      bzz = 0.5d0 * (- dgutz * gutt + gutz * dgutt) / gutt**2

      if (psi_on) then
         psiz = dpsi
      end if

C     Calculate t-derivatives, and extrinsic curvature.

      call Exact__metric_deriv(
     $     decoded_exact_model,
     $     0_ik,
     $     x, y, z, t,
     $     dgdtt, dgdtx, dgdty, dgdtz,
     $     dgdxx, dgdyy, dgdzz, dgdxy, dgdyz, dgdxz,
     $     dgutt, dgutx, dguty, dgutz,
     $     dguxx, dguyy, dguzz, dguxy, dguyz, dguxz,
     $     dpsi)

      hxx = - 0.5d0 * dgdxx / alp
     $     + (dxgxx * betax + dygxx * betay + dzgxx * betaz
     $     + 2.d0 * (bxx * gxx + bxy * gxy + bxz * gxz)) / alp

      hyy = - 0.5d0 * dgdyy / alp
     $     + (dxgyy * betax + dygyy * betay + dzgyy * betaz
     $     + 2.d0 * (byx * gxy + byy * gyy + byz * gyz)) / alp

      hzz = - 0.5d0 * dgdzz / alp
     $     + (dxgzz * betax + dygzz * betay + dzgzz * betaz
     $     + 2.d0 * (bzx * gxz + bzy * gyz + bzz * gzz)) / alp

      hxy = - 0.5d0 * dgdxy / alp
     $     + (dxgxy * betax + dygxy * betay + dzgxy * betaz
     $     + bxx * gxy + bxy * gyy + bxz * gyz
     $     + byx * gxx + byy * gxy + byz * gxz) / alp

      hyz = - 0.5d0 * dgdyz / alp
     $     + (dxgyz * betax + dygyz * betay + dzgyz * betaz
     $     + byx * gxz + byy * gyz + byz * gzz
     $     + bzx * gxy + bzy * gyy + bzz * gyz) / alp

      hxz = - 0.5d0 * dgdxz / alp
     $     + (dxgxz * betax + dygxz * betay + dzgxz * betaz
     $     + bxx * gxz + bxy * gyz + bxz * gzz
     $     + bzx * gxx + bzy * gxy + bzz * gxz) / alp

C     Calculate time derivatives of  lapse and shift

C     alp = 1.d0 / sqrt(-gutt)
      dtalp = 0.5d0 / sqrt(-gutt)**3 * dgutt

C     betax = - gutx / gutt
C     betay = - guty / gutt
C     betaz = - gutz / gutt
      dtbetax = (- dgutx * gutt + gutx * dgutt) / gutt**2
      dtbetay = (- dguty * gutt + guty * dgutt) / gutt**2
      dtbetaz = (- dgutz * gutt + gutz * dgutt) / gutt**2

C     Calculate second derivatives of the conformal factor

      if (psi_on) then
        call Exact__metric_deriv2(
     $        decoded_exact_model,
     $        1_ik, 1_ik,
     $        x, y, z, t,
     $        dgdtt, dgdtx, dgdty, dgdtz,
     $        dgdxx, dgdyy, dgdzz, dgdxy, dgdyz, dgdxz,
     $        dgutt, dgutx, dguty, dgutz,
     $        dguxx, dguyy, dguzz, dguxy, dguyz, dguxz,
     $        psixx)
        call Exact__metric_deriv2(
     $        decoded_exact_model,
     $        1_ik, 2_ik,
     $        x, y, z, t,
     $        dgdtt, dgdtx, dgdty, dgdtz,
     $        dgdxx, dgdyy, dgdzz, dgdxy, dgdyz, dgdxz,
     $        dgutt, dgutx, dguty, dgutz,
     $        dguxx, dguyy, dguzz, dguxy, dguyz, dguxz,
     $        psixy)
        call Exact__metric_deriv2(
     $        decoded_exact_model,
     $        1_ik, 3_ik,
     $        x, y, z, t,
     $        dgdtt, dgdtx, dgdty, dgdtz,
     $        dgdxx, dgdyy, dgdzz, dgdxy, dgdyz, dgdxz,
     $        dgutt, dgutx, dguty, dgutz,
     $        dguxx, dguyy, dguzz, dguxy, dguyz, dguxz,
     $        psixz)
        call Exact__metric_deriv2(
     $        decoded_exact_model,
     $        2_ik, 2_ik,
     $        x, y, z, t,
     $        dgdtt, dgdtx, dgdty, dgdtz,
     $        dgdxx, dgdyy, dgdzz, dgdxy, dgdyz, dgdxz,
     $        dgutt, dgutx, dguty, dgutz,
     $        dguxx, dguyy, dguzz, dguxy, dguyz, dguxz,
     $        psiyy)
        call Exact__metric_deriv2(
     $        decoded_exact_model,
     $        2_ik, 3_ik,
     $        x, y, z, t,
     $        dgdtt, dgdtx, dgdty, dgdtz,
     $        dgdxx, dgdyy, dgdzz, dgdxy, dgdyz, dgdxz,
     $        dgutt, dgutx, dguty, dgutz,
     $        dguxx, dguyy, dguzz, dguxy, dguyz, dguxz,
     $        psiyz)
        call Exact__metric_deriv2(
     $        decoded_exact_model,
     $        3_ik, 3_ik,
     $        x, y, z, t,
     $        dgdtt, dgdtx, dgdty, dgdtz,
     $        dgdxx, dgdyy, dgdzz, dgdxy, dgdyz, dgdxz,
     $        dgutt, dgutx, dguty, dgutz,
     $        dguxx, dguyy, dguzz, dguxy, dguyz, dguxz,
     $        psizz)
        
      end if
      return
      end



      subroutine Exact__metric_deriv(
     $     decoded_exact_model,
     $     dir,
     $     x, y, z, t,
     $     dgdtt, dgdtx, dgdty, dgdtz,
     $     dgdxx, dgdyy, dgdzz, dgdxy, dgdyz, dgdxz,
     $     dgutt, dgutx, dguty, dgutz,
     $     dguxx, dguyy, dguzz, dguxy, dguyz, dguxz,
     $     dpsi)
      
      implicit none
      DECLARE_CCTK_PARAMETERS
      
      CCTK_INT
     $     decoded_exact_model,
     $     dir
      CCTK_REAL
     $     x, y, z, t,
     $     dgdtt, dgdtx, dgdty, dgdtz,
     $     dgdxx, dgdyy, dgdzz, dgdxy, dgdyz, dgdxz,
     $     dgutt, dgutx, dguty, dgutz,
     $     dguxx, dguyy, dguzz, dguxy, dguyz, dguxz,
     $     dpsi,
     $     gdtt_p, gdtx_p, gdty_p, gdtz_p, 
     $     gdxx_p, gdyy_p, gdzz_p, gdxy_p, gdyz_p, gdxz_p,
     $     gutt_p, gutx_p, guty_p, gutz_p, 
     $     guxx_p, guyy_p, guzz_p, guxy_p, guyz_p, guxz_p,
     $     psi_p,
     $     gdtt_m, gdtx_m, gdty_m, gdtz_m, 
     $     gdxx_m, gdyy_m, gdzz_m, gdxy_m, gdyz_m, gdxz_m,
     $     gutt_m, gutx_m, guty_m, gutz_m, 
     $     guxx_m, guyy_m, guzz_m, guxy_m, guyz_m, guxz_m,
     $     psi_m,
     $     gdtt_p_p, gdtx_p_p, gdty_p_p, gdtz_p_p, 
     $     gdxx_p_p, gdyy_p_p, gdzz_p_p, gdxy_p_p, gdyz_p_p, gdxz_p_p,
     $     gutt_p_p, gutx_p_p, guty_p_p, gutz_p_p, 
     $     guxx_p_p, guyy_p_p, guzz_p_p, guxy_p_p, guyz_p_p, guxz_p_p,
     $     psi_p_p,
     $     gdtt_m_m, gdtx_m_m, gdty_m_m, gdtz_m_m, 
     $     gdxx_m_m, gdyy_m_m, gdzz_m_m, gdxy_m_m, gdyz_m_m, gdxz_m_m,
     $     gutt_m_m, gutx_m_m, guty_m_m, gutz_m_m, 
     $     guxx_m_m, guyy_m_m, guzz_m_m, guxy_m_m, guyz_m_m, guxz_m_m,
     $     psi_m_m,
     $     eps,
     $     dx, dy, dz, dt
      
      eps = exact_eps
      
      dx = 0
      dy = 0
      dz = 0
      dt = 0
      if (dir.eq.0) dt = eps
      if (dir.eq.1) dx = eps
      if (dir.eq.2) dy = eps
      if (dir.eq.3) dz = eps
      
      if (exact_order .eq. 2) then
         
         call Exact__metric(
     $        decoded_exact_model,
     $        x-dx, y-dy, z-dz, t-dt,
     $        gdtt_m, gdtx_m, gdty_m, gdtz_m,
     $        gdxx_m, gdyy_m, gdzz_m, gdxy_m, gdyz_m, gdxz_m,
     $        gutt_m, gutx_m, guty_m, gutz_m,
     $        guxx_m, guyy_m, guzz_m, guxy_m, guyz_m, guxz_m,
     $        psi_m)
         call Exact__metric(
     $        decoded_exact_model,
     $        x+dx, y+dy, z+dz, t+dt,
     $        gdtt_p, gdtx_p, gdty_p, gdtz_p,
     $        gdxx_p, gdyy_p, gdzz_p, gdxy_p, gdyz_p, gdxz_p,
     $        gutt_p, gutx_p, guty_p, gutz_p,
     $        guxx_p, guyy_p, guzz_p, guxy_p, guyz_p, guxz_p,
     $        psi_p)
         
         dgdtt = (gdtt_p - gdtt_m) / (2*eps)
         dgdtx = (gdtx_p - gdtx_m) / (2*eps)
         dgdty = (gdty_p - gdty_m) / (2*eps)
         dgdtz = (gdtz_p - gdtz_m) / (2*eps)
         dgdxx = (gdxx_p - gdxx_m) / (2*eps)
         dgdyy = (gdyy_p - gdyy_m) / (2*eps)
         dgdzz = (gdzz_p - gdzz_m) / (2*eps)
         dgdxy = (gdxy_p - gdxy_m) / (2*eps)
         dgdyz = (gdyz_p - gdyz_m) / (2*eps)
         dgdxz = (gdxz_p - gdxz_m) / (2*eps)
         dgutt = (gutt_p - gutt_m) / (2*eps)
         dgutx = (gutx_p - gutx_m) / (2*eps)
         dguty = (guty_p - guty_m) / (2*eps)
         dgutz = (gutz_p - gutz_m) / (2*eps)
         dguxx = (guxx_p - guxx_m) / (2*eps)
         dguyy = (guyy_p - guyy_m) / (2*eps)
         dguzz = (guzz_p - guzz_m) / (2*eps)
         dguxy = (guxy_p - guxy_m) / (2*eps)
         dguyz = (guyz_p - guyz_m) / (2*eps)
         dguxz = (guxz_p - guxz_m) / (2*eps)
         dpsi  = (psi_p  - psi_m ) / (2*eps)
         
      else if (exact_order .eq. 4) then
         
         call Exact__metric(
     $        decoded_exact_model,
     $        x-2*dx, y-2*dy, z-2*dz, t-2*dt,
     $        gdtt_m_m, gdtx_m_m, gdty_m_m, gdtz_m_m,
     $        gdxx_m_m, gdyy_m_m, gdzz_m_m, gdxy_m_m, gdyz_m_m, gdxz_m_m,
     $        gutt_m_m, gutx_m_m, guty_m_m, gutz_m_m,
     $        guxx_m_m, guyy_m_m, guzz_m_m, guxy_m_m, guyz_m_m, guxz_m_m,
     $        psi_m_m)
         call Exact__metric(
     $        decoded_exact_model,
     $        x-dx, y-dy, z-dz, t-dt,
     $        gdtt_m, gdtx_m, gdty_m, gdtz_m,
     $        gdxx_m, gdyy_m, gdzz_m, gdxy_m, gdyz_m, gdxz_m,
     $        gutt_m, gutx_m, guty_m, gutz_m,
     $        guxx_m, guyy_m, guzz_m, guxy_m, guyz_m, guxz_m,
     $        psi_m)
         call Exact__metric(
     $        decoded_exact_model,
     $        x+dx, y+dy, z+dz, t+dt,
     $        gdtt_p, gdtx_p, gdty_p, gdtz_p,
     $        gdxx_p, gdyy_p, gdzz_p, gdxy_p, gdyz_p, gdxz_p,
     $        gutt_p, gutx_p, guty_p, gutz_p,
     $        guxx_p, guyy_p, guzz_p, guxy_p, guyz_p, guxz_p,
     $        psi_p)
         call Exact__metric(
     $        decoded_exact_model,
     $        x+2*dx, y+2*dy, z+2*dz, t+2*dt,
     $        gdtt_p_p, gdtx_p_p, gdty_p_p, gdtz_p_p,
     $        gdxx_p_p, gdyy_p_p, gdzz_p_p, gdxy_p_p, gdyz_p_p, gdxz_p_p,
     $        gutt_p_p, gutx_p_p, guty_p_p, gutz_p_p,
     $        guxx_p_p, guyy_p_p, guzz_p_p, guxy_p_p, guyz_p_p, guxz_p_p,
     $        psi_p_p)
         
         dgdtt = (- gdtt_p_p + 8*gdtt_p - 8*gdtt_m + gdtt_m_m) / (12*eps)
         dgdtx = (- gdtx_p_p + 8*gdtx_p - 8*gdtx_m + gdtx_m_m) / (12*eps)
         dgdty = (- gdty_p_p + 8*gdty_p - 8*gdty_m + gdty_m_m) / (12*eps)
         dgdtz = (- gdtz_p_p + 8*gdtz_p - 8*gdtz_m + gdtz_m_m) / (12*eps)
         dgdxx = (- gdxx_p_p + 8*gdxx_p - 8*gdxx_m + gdxx_m_m) / (12*eps)
         dgdyy = (- gdyy_p_p + 8*gdyy_p - 8*gdyy_m + gdyy_m_m) / (12*eps)
         dgdzz = (- gdzz_p_p + 8*gdzz_p - 8*gdzz_m + gdzz_m_m) / (12*eps)
         dgdxy = (- gdxy_p_p + 8*gdxy_p - 8*gdxy_m + gdxy_m_m) / (12*eps)
         dgdyz = (- gdyz_p_p + 8*gdyz_p - 8*gdyz_m + gdyz_m_m) / (12*eps)
         dgdxz = (- gdxz_p_p + 8*gdxz_p - 8*gdxz_m + gdxz_m_m) / (12*eps)
         dgutt = (- gutt_p_p + 8*gutt_p - 8*gutt_m + gutt_m_m) / (12*eps)
         dgutx = (- gutx_p_p + 8*gutx_p - 8*gutx_m + gutx_m_m) / (12*eps)
         dguty = (- guty_p_p + 8*guty_p - 8*guty_m + guty_m_m) / (12*eps)
         dgutz = (- gutz_p_p + 8*gutz_p - 8*gutz_m + gutz_m_m) / (12*eps)
         dguxx = (- guxx_p_p + 8*guxx_p - 8*guxx_m + guxx_m_m) / (12*eps)
         dguyy = (- guyy_p_p + 8*guyy_p - 8*guyy_m + guyy_m_m) / (12*eps)
         dguzz = (- guzz_p_p + 8*guzz_p - 8*guzz_m + guzz_m_m) / (12*eps)
         dguxy = (- guxy_p_p + 8*guxy_p - 8*guxy_m + guxy_m_m) / (12*eps)
         dguyz = (- guyz_p_p + 8*guyz_p - 8*guyz_m + guyz_m_m) / (12*eps)
         dguxz = (- guxz_p_p + 8*guxz_p - 8*guxz_m + guxz_m_m) / (12*eps)
         dpsi  = (- psi_p_p  + 8*psi_p  - 8*psi_m  + psi_m_m ) / (12*eps)
         
      else
         call CCTK_WARN (CCTK_WARN_ABORT, "internal error")
      end if
      
      end



      subroutine Exact__metric_deriv2(
     $     decoded_exact_model,
     $     dir1, dir2,
     $     x, y, z, t,
     $     dgdtt, dgdtx, dgdty, dgdtz,
     $     dgdxx, dgdyy, dgdzz, dgdxy, dgdyz, dgdxz,
     $     dgutt, dgutx, dguty, dgutz,
     $     dguxx, dguyy, dguzz, dguxy, dguyz, dguxz,
     $     dpsi)
      
      implicit none
      DECLARE_CCTK_PARAMETERS
      
      CCTK_INT
     $     decoded_exact_model,
     $     dir1, dir2
      CCTK_REAL
     $     x, y, z, t,
     $     dgdtt, dgdtx, dgdty, dgdtz,
     $     dgdxx, dgdyy, dgdzz, dgdxy, dgdyz, dgdxz,
     $     dgutt, dgutx, dguty, dgutz,
     $     dguxx, dguyy, dguzz, dguxy, dguyz, dguxz,
     $     dpsi,
     $     gdtt_0, gdtx_0, gdty_0, gdtz_0, 
     $     gdxx_0, gdyy_0, gdzz_0, gdxy_0, gdyz_0, gdxz_0,
     $     gutt_0, gutx_0, guty_0, gutz_0, 
     $     guxx_0, guyy_0, guzz_0, guxy_0, guyz_0, guxz_0,
     $     psi_0,
     $     gdtt_p, gdtx_p, gdty_p, gdtz_p, 
     $     gdxx_p, gdyy_p, gdzz_p, gdxy_p, gdyz_p, gdxz_p,
     $     gutt_p, gutx_p, guty_p, gutz_p, 
     $     guxx_p, guyy_p, guzz_p, guxy_p, guyz_p, guxz_p,
     $     psi_p,
     $     gdtt_m, gdtx_m, gdty_m, gdtz_m, 
     $     gdxx_m, gdyy_m, gdzz_m, gdxy_m, gdyz_m, gdxz_m,
     $     gutt_m, gutx_m, guty_m, gutz_m, 
     $     guxx_m, guyy_m, guzz_m, guxy_m, guyz_m, guxz_m,
     $     psi_m,
     $     gdtt_p_p, gdtx_p_p, gdty_p_p, gdtz_p_p, 
     $     gdxx_p_p, gdyy_p_p, gdzz_p_p, gdxy_p_p, gdyz_p_p, gdxz_p_p,
     $     gutt_p_p, gutx_p_p, guty_p_p, gutz_p_p, 
     $     guxx_p_p, guyy_p_p, guzz_p_p, guxy_p_p, guyz_p_p, guxz_p_p,
     $     psi_p_p,
     $     gdtt_m_m, gdtx_m_m, gdty_m_m, gdtz_m_m, 
     $     gdxx_m_m, gdyy_m_m, gdzz_m_m, gdxy_m_m, gdyz_m_m, gdxz_m_m,
     $     gutt_m_m, gutx_m_m, guty_m_m, gutz_m_m, 
     $     guxx_m_m, guyy_m_m, guzz_m_m, guxy_m_m, guyz_m_m, guxz_m_m,
     $     psi_m_m,
     $     eps,
     $     dx, dy, dz, dt
      
      eps = exact_eps
      
      dx = 0
      dy = 0
      dz = 0
      dt = 0
      if (dir1.eq.0) dt = eps
      if (dir1.eq.1) dx = eps
      if (dir1.eq.2) dy = eps
      if (dir1.eq.3) dz = eps
      if (dir1.lt.0 .or. dir1.gt.3) then
         call CCTK_WARN (CCTK_WARN_ABORT, "internal error")
      end if
      if (dir2.lt.0 .or. dir2.gt.3) then
         call CCTK_WARN (CCTK_WARN_ABORT, "internal error")
      end if
      
      if (exact_order .eq. 2) then
         
         if (dir1 .eq. dir2) then
            
            call Exact__metric(
     $           decoded_exact_model,
     $           x-dx, y-dy, z-dz, t-dt,
     $           gdtt_m, gdtx_m, gdty_m, gdtz_m,
     $           gdxx_m, gdyy_m, gdzz_m, gdxy_m, gdyz_m, gdxz_m,
     $           gutt_m, gutx_m, guty_m, gutz_m,
     $           guxx_m, guyy_m, guzz_m, guxy_m, guyz_m, guxz_m,
     $           psi_m)
            call Exact__metric(
     $           decoded_exact_model,
     $           x, y, z, t,
     $           gdtt_0, gdtx_0, gdty_0, gdtz_0,
     $           gdxx_0, gdyy_0, gdzz_0, gdxy_0, gdyz_0, gdxz_0,
     $           gutt_0, gutx_0, guty_0, gutz_0,
     $           guxx_0, guyy_0, guzz_0, guxy_0, guyz_0, guxz_0,
     $           psi_0)
            call Exact__metric(
     $           decoded_exact_model,
     $           x+dx, y+dy, z+dz, t+dt,
     $           gdtt_p, gdtx_p, gdty_p, gdtz_p,
     $           gdxx_p, gdyy_p, gdzz_p, gdxy_p, gdyz_p, gdxz_p,
     $           gutt_p, gutx_p, guty_p, gutz_p,
     $           guxx_p, guyy_p, guzz_p, guxy_p, guyz_p, guxz_p,
     $           psi_p)
            
            dgdtt = (gdtt_m - 2*gdtt_0 + gdtt_p) / eps**2
            dgdtx = (gdtx_m - 2*gdtx_0 + gdtx_p) / eps**2
            dgdty = (gdty_m - 2*gdty_0 + gdty_p) / eps**2
            dgdtz = (gdtz_m - 2*gdtz_0 + gdtz_p) / eps**2
            dgdxx = (gdxx_m - 2*gdxx_0 + gdxx_p) / eps**2
            dgdyy = (gdyy_m - 2*gdyy_0 + gdyy_p) / eps**2
            dgdzz = (gdzz_m - 2*gdzz_0 + gdzz_p) / eps**2
            dgdxy = (gdxy_m - 2*gdxy_0 + gdxy_p) / eps**2
            dgdyz = (gdyz_m - 2*gdyz_0 + gdyz_p) / eps**2
            dgdxz = (gdxz_m - 2*gdxz_0 + gdxz_p) / eps**2
            dgutt = (gutt_m - 2*gutt_0 + gutt_p) / eps**2
            dgutx = (gutx_m - 2*gutx_0 + gutx_p) / eps**2
            dguty = (guty_m - 2*guty_0 + guty_p) / eps**2
            dgutz = (gutz_m - 2*gutz_0 + gutz_p) / eps**2
            dguxx = (guxx_m - 2*guxx_0 + guxx_p) / eps**2
            dguyy = (guyy_m - 2*guyy_0 + guyy_p) / eps**2
            dguzz = (guzz_m - 2*guzz_0 + guzz_p) / eps**2
            dguxy = (guxy_m - 2*guxy_0 + guxy_p) / eps**2
            dguyz = (guyz_m - 2*guyz_0 + guyz_p) / eps**2
            dguxz = (guxz_m - 2*guxz_0 + guxz_p) / eps**2
            dpsi  = (psi_m  - 2*psi_0  + psi_p ) / eps**2
            
         else
            
            call Exact__metric_deriv(
     $           decoded_exact_model,
     $           dir2,
     $           x-dx, y-dy, z-dz, t-dt,
     $           gdtt_m, gdtx_m, gdty_m, gdtz_m,
     $           gdxx_m, gdyy_m, gdzz_m, gdxy_m, gdyz_m, gdxz_m,
     $           gutt_m, gutx_m, guty_m, gutz_m,
     $           guxx_m, guyy_m, guzz_m, guxy_m, guyz_m, guxz_m,
     $           psi_m)
            call Exact__metric_deriv(
     $           decoded_exact_model,
     $           dir2,
     $           x+dx, y+dy, z+dz, t+dt,
     $           gdtt_p, gdtx_p, gdty_p, gdtz_p,
     $           gdxx_p, gdyy_p, gdzz_p, gdxy_p, gdyz_p, gdxz_p,
     $           gutt_p, gutx_p, guty_p, gutz_p,
     $           guxx_p, guyy_p, guzz_p, guxy_p, guyz_p, guxz_p,
     $           psi_p)
         
            dgdtt = (gdtt_p - gdtt_m) / (2*eps)
            dgdtx = (gdtx_p - gdtx_m) / (2*eps)
            dgdty = (gdty_p - gdty_m) / (2*eps)
            dgdtz = (gdtz_p - gdtz_m) / (2*eps)
            dgdxx = (gdxx_p - gdxx_m) / (2*eps)
            dgdyy = (gdyy_p - gdyy_m) / (2*eps)
            dgdzz = (gdzz_p - gdzz_m) / (2*eps)
            dgdxy = (gdxy_p - gdxy_m) / (2*eps)
            dgdyz = (gdyz_p - gdyz_m) / (2*eps)
            dgdxz = (gdxz_p - gdxz_m) / (2*eps)
            dgutt = (gutt_p - gutt_m) / (2*eps)
            dgutx = (gutx_p - gutx_m) / (2*eps)
            dguty = (guty_p - guty_m) / (2*eps)
            dgutz = (gutz_p - gutz_m) / (2*eps)
            dguxx = (guxx_p - guxx_m) / (2*eps)
            dguyy = (guyy_p - guyy_m) / (2*eps)
            dguzz = (guzz_p - guzz_m) / (2*eps)
            dguxy = (guxy_p - guxy_m) / (2*eps)
            dguyz = (guyz_p - guyz_m) / (2*eps)
            dguxz = (guxz_p - guxz_m) / (2*eps)
            dpsi  = (psi_p  - psi_m ) / (2*eps)
            
         end if
         
      else if (exact_order .eq. 4) then
         
         if (dir1 .eq. dir2) then
            
            call Exact__metric(
     $           decoded_exact_model,
     $           x-2*dx, y-2*dy, z-2*dz, t-2*dt,
     $           gdtt_m_m, gdtx_m_m, gdty_m_m, gdtz_m_m,
     $           gdxx_m_m, gdyy_m_m, gdzz_m_m, gdxy_m_m, gdyz_m_m, gdxz_m_m,
     $           gutt_m_m, gutx_m_m, guty_m_m, gutz_m_m,
     $           guxx_m_m, guyy_m_m, guzz_m_m, guxy_m_m, guyz_m_m, guxz_m_m,
     $           psi_m_m)
            call Exact__metric(
     $           decoded_exact_model,
     $           x-dx, y-dy, z-dz, t-dt,
     $           gdtt_m, gdtx_m, gdty_m, gdtz_m,
     $           gdxx_m, gdyy_m, gdzz_m, gdxy_m, gdyz_m, gdxz_m,
     $           gutt_m, gutx_m, guty_m, gutz_m,
     $           guxx_m, guyy_m, guzz_m, guxy_m, guyz_m, guxz_m,
     $           psi_m)
            call Exact__metric(
     $           decoded_exact_model,
     $           x, y, z, t,
     $           gdtt_0, gdtx_0, gdty_0, gdtz_0,
     $           gdxx_0, gdyy_0, gdzz_0, gdxy_0, gdyz_0, gdxz_0,
     $           gutt_0, gutx_0, guty_0, gutz_0,
     $           guxx_0, guyy_0, guzz_0, guxy_0, guyz_0, guxz_0,
     $           psi_0)
            call Exact__metric(
     $           decoded_exact_model,
     $           x+dx, y+dy, z+dz, t+dt,
     $           gdtt_p, gdtx_p, gdty_p, gdtz_p,
     $           gdxx_p, gdyy_p, gdzz_p, gdxy_p, gdyz_p, gdxz_p,
     $           gutt_p, gutx_p, guty_p, gutz_p,
     $           guxx_p, guyy_p, guzz_p, guxy_p, guyz_p, guxz_p,
     $           psi_p)
            call Exact__metric(
     $           decoded_exact_model,
     $           x+2*dx, y+2*dy, z+2*dz, t+2*dt,
     $           gdtt_p_p, gdtx_p_p, gdty_p_p, gdtz_p_p,
     $           gdxx_p_p, gdyy_p_p, gdzz_p_p, gdxy_p_p, gdyz_p_p, gdxz_p_p,
     $           gutt_p_p, gutx_p_p, guty_p_p, gutz_p_p,
     $           guxx_p_p, guyy_p_p, guzz_p_p, guxy_p_p, guyz_p_p, guxz_p_p,
     $           psi_p_p)
            
            dgdtt = (- gdtt_m_m - 16*gdtt_m + 30*gdtt_0 - 16*gdtt_p - gdtt_p_p) / (12*eps**2)
            dgdtx = (- gdtx_m_m - 16*gdtx_m + 30*gdtx_0 - 16*gdtx_p - gdtx_p_p) / (12*eps**2)
            dgdty = (- gdty_m_m - 16*gdty_m + 30*gdty_0 - 16*gdty_p - gdty_p_p) / (12*eps**2)
            dgdtz = (- gdtz_m_m - 16*gdtz_m + 30*gdtz_0 - 16*gdtz_p - gdtz_p_p) / (12*eps**2)
            dgdxx = (- gdxx_m_m - 16*gdxx_m + 30*gdxx_0 - 16*gdxx_p - gdxx_p_p) / (12*eps**2)
            dgdyy = (- gdyy_m_m - 16*gdyy_m + 30*gdyy_0 - 16*gdyy_p - gdyy_p_p) / (12*eps**2)
            dgdzz = (- gdzz_m_m - 16*gdzz_m + 30*gdzz_0 - 16*gdzz_p - gdzz_p_p) / (12*eps**2)
            dgdxy = (- gdxy_m_m - 16*gdxy_m + 30*gdxy_0 - 16*gdxy_p - gdxy_p_p) / (12*eps**2)
            dgdyz = (- gdyz_m_m - 16*gdyz_m + 30*gdyz_0 - 16*gdyz_p - gdyz_p_p) / (12*eps**2)
            dgdxz = (- gdxz_m_m - 16*gdxz_m + 30*gdxz_0 - 16*gdxz_p - gdxz_p_p) / (12*eps**2)
            dgutt = (- gutt_m_m - 16*gutt_m + 30*gutt_0 - 16*gutt_p - gutt_p_p) / (12*eps**2)
            dgutx = (- gutx_m_m - 16*gutx_m + 30*gutx_0 - 16*gutx_p - gutx_p_p) / (12*eps**2)
            dguty = (- guty_m_m - 16*guty_m + 30*guty_0 - 16*guty_p - guty_p_p) / (12*eps**2)
            dgutz = (- gutz_m_m - 16*gutz_m + 30*gutz_0 - 16*gutz_p - gutz_p_p) / (12*eps**2)
            dguxx = (- guxx_m_m - 16*guxx_m + 30*guxx_0 - 16*guxx_p - guxx_p_p) / (12*eps**2)
            dguyy = (- guyy_m_m - 16*guyy_m + 30*guyy_0 - 16*guyy_p - guyy_p_p) / (12*eps**2)
            dguzz = (- guzz_m_m - 16*guzz_m + 30*guzz_0 - 16*guzz_p - guzz_p_p) / (12*eps**2)
            dguxy = (- guxy_m_m - 16*guxy_m + 30*guxy_0 - 16*guxy_p - guxy_p_p) / (12*eps**2)
            dguyz = (- guyz_m_m - 16*guyz_m + 30*guyz_0 - 16*guyz_p - guyz_p_p) / (12*eps**2)
            dguxz = (- guxz_m_m - 16*guxz_m + 30*guxz_0 - 16*guxz_p - guxz_p_p) / (12*eps**2)
            dpsi  = (- psi_m_m  - 16*psi_m  + 30*psi_0  - 16*psi_p  - psi_p_p ) / (12*eps**2)
            
         else
            
            call Exact__metric_deriv(
     $           decoded_exact_model,
     $           dir2,
     $           x-2*dx, y-2*dy, z-2*dz, t-2*dt,
     $           gdtt_m_m, gdtx_m_m, gdty_m_m, gdtz_m_m,
     $           gdxx_m_m, gdyy_m_m, gdzz_m_m, gdxy_m_m, gdyz_m_m, gdxz_m_m,
     $           gutt_m_m, gutx_m_m, guty_m_m, gutz_m_m,
     $           guxx_m_m, guyy_m_m, guzz_m_m, guxy_m_m, guyz_m_m, guxz_m_m,
     $           psi_m_m)
            call Exact__metric_deriv(
     $           decoded_exact_model,
     $           dir2,
     $           x-dx, y-dy, z-dz, t-dt,
     $           gdtt_m, gdtx_m, gdty_m, gdtz_m,
     $           gdxx_m, gdyy_m, gdzz_m, gdxy_m, gdyz_m, gdxz_m,
     $           gutt_m, gutx_m, guty_m, gutz_m,
     $           guxx_m, guyy_m, guzz_m, guxy_m, guyz_m, guxz_m,
     $           psi_m)
            call Exact__metric_deriv(
     $           decoded_exact_model,
     $           dir2,
     $           x+dx, y+dy, z+dz, t+dt,
     $           gdtt_p, gdtx_p, gdty_p, gdtz_p,
     $           gdxx_p, gdyy_p, gdzz_p, gdxy_p, gdyz_p, gdxz_p,
     $           gutt_p, gutx_p, guty_p, gutz_p,
     $           guxx_p, guyy_p, guzz_p, guxy_p, guyz_p, guxz_p,
     $           psi_p)
            call Exact__metric_deriv(
     $           decoded_exact_model,
     $           dir2,
     $           x+2*dx, y+2*dy, z+2*dz, t+2*dt,
     $           gdtt_p_p, gdtx_p_p, gdty_p_p, gdtz_p_p,
     $           gdxx_p_p, gdyy_p_p, gdzz_p_p, gdxy_p_p, gdyz_p_p, gdxz_p_p,
     $           gutt_p_p, gutx_p_p, guty_p_p, gutz_p_p,
     $           guxx_p_p, guyy_p_p, guzz_p_p, guxy_p_p, guyz_p_p, guxz_p_p,
     $           psi_p_p)
            
            dgdtt = (- gdtt_p_p + 8*gdtt_p - 8*gdtt_m + gdtt_m_m) / (12*eps)
            dgdtx = (- gdtx_p_p + 8*gdtx_p - 8*gdtx_m + gdtx_m_m) / (12*eps)
            dgdty = (- gdty_p_p + 8*gdty_p - 8*gdty_m + gdty_m_m) / (12*eps)
            dgdtz = (- gdtz_p_p + 8*gdtz_p - 8*gdtz_m + gdtz_m_m) / (12*eps)
            dgdxx = (- gdxx_p_p + 8*gdxx_p - 8*gdxx_m + gdxx_m_m) / (12*eps)
            dgdyy = (- gdyy_p_p + 8*gdyy_p - 8*gdyy_m + gdyy_m_m) / (12*eps)
            dgdzz = (- gdzz_p_p + 8*gdzz_p - 8*gdzz_m + gdzz_m_m) / (12*eps)
            dgdxy = (- gdxy_p_p + 8*gdxy_p - 8*gdxy_m + gdxy_m_m) / (12*eps)
            dgdyz = (- gdyz_p_p + 8*gdyz_p - 8*gdyz_m + gdyz_m_m) / (12*eps)
            dgdxz = (- gdxz_p_p + 8*gdxz_p - 8*gdxz_m + gdxz_m_m) / (12*eps)
            dgutt = (- gutt_p_p + 8*gutt_p - 8*gutt_m + gutt_m_m) / (12*eps)
            dgutx = (- gutx_p_p + 8*gutx_p - 8*gutx_m + gutx_m_m) / (12*eps)
            dguty = (- guty_p_p + 8*guty_p - 8*guty_m + guty_m_m) / (12*eps)
            dgutz = (- gutz_p_p + 8*gutz_p - 8*gutz_m + gutz_m_m) / (12*eps)
            dguxx = (- guxx_p_p + 8*guxx_p - 8*guxx_m + guxx_m_m) / (12*eps)
            dguyy = (- guyy_p_p + 8*guyy_p - 8*guyy_m + guyy_m_m) / (12*eps)
            dguzz = (- guzz_p_p + 8*guzz_p - 8*guzz_m + guzz_m_m) / (12*eps)
            dguxy = (- guxy_p_p + 8*guxy_p - 8*guxy_m + guxy_m_m) / (12*eps)
            dguyz = (- guyz_p_p + 8*guyz_p - 8*guyz_m + guyz_m_m) / (12*eps)
            dguxz = (- guxz_p_p + 8*guxz_p - 8*guxz_m + guxz_m_m) / (12*eps)
            dpsi  = (- psi_p_p  + 8*psi_p  - 8*psi_m  + psi_m_m ) / (12*eps)
            
         end if
         
      else
         call CCTK_WARN (CCTK_WARN_ABORT, "internal error")
      end if
      
      end
