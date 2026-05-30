c This subroutine calculates the 4-metric and its inverse at an event,
c for a given model, by decoding  decoded_exact_model  and calling the
c appropriate subroutine for that model.
C $Header$

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#include "param_defs.inc"

      subroutine Exact__metric_for_model(
     $     decoded_exact_model,
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz,
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $     gutt, gutx, guty, gutz,
     $     guxx, guyy, guzz, guxy, guyz, guxz,
     $     psi, Tmunu_flag)

      implicit none
      DECLARE_CCTK_FUNCTIONS

c arguments
      CCTK_INT decoded_exact_model
      CCTK_REAL x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz,
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $     gutt, gutx, guty, gutz,
     $     guxx, guyy, guzz, guxy, guyz, guxz
      CCTK_REAL psi
      LOGICAL   Tmunu_flag

c local variables
      character*100 warn_buffer

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c Minkowski spacetime
c

      if     (decoded_exact_model .eq. EXACT__Minkowski) then
         call Exact__Minkowski(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz,
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz,
     $        guxx, guyy, guzz, guxy, guyz, guxz,
     $        psi, Tmunu_flag)

      elseif (decoded_exact_model .eq. EXACT__Minkowski_shift) then
         call Exact__Minkowski_shift(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz,
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz,
     $        guxx, guyy, guzz, guxy, guyz, guxz,
     $        psi, Tmunu_flag)

      elseif (decoded_exact_model .eq. EXACT__Minkowski_funny) then
         call Exact__Minkowski_funny(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz,
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz,
     $        guxx, guyy, guzz, guxy, guyz, guxz,
     $        psi, Tmunu_flag)

      elseif (decoded_exact_model .eq. EXACT__Minkowski_gauge_wave) then
         call Exact__Minkowski_gauge_wave(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz,
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz,
     $        guxx, guyy, guzz, guxy, guyz, guxz,
     $        psi, Tmunu_flag)

      elseif (decoded_exact_model .eq. EXACT__Minkowski_shifted_gauge_wave) then
         call Exact__Minkowski_shifted_gauge_wave(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz,
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz,
     $        guxx, guyy, guzz, guxy, guyz, guxz,
     $        psi, Tmunu_flag)

      elseif (decoded_exact_model .eq. EXACT__Minkowski_conf_wave) then
         call Exact__Minkowski_conf_wave(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz,
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz,
     $        guxx, guyy, guzz, guxy, guyz, guxz,
     $        psi, Tmunu_flag)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c black hole spacetimes
c

      elseif (decoded_exact_model .eq. EXACT__Schwarzschild_EF) then
         call Exact__Schwarzschild_EF(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz,
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz,
     $        guxx, guyy, guzz, guxy, guyz, guxz,
     $        psi, Tmunu_flag)

      elseif (decoded_exact_model .eq. EXACT__Schwarzschild_PG) then
         call Exact__Schwarzschild_PG(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz,
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz,
     $        guxx, guyy, guzz, guxy, guyz, guxz,
     $        psi, Tmunu_flag)

      elseif (decoded_exact_model .eq. EXACT__Schwarzschild_BL) then
         call Exact__Schwarzschild_BL(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz,
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz,
     $        guxx, guyy, guzz, guxy, guyz, guxz,
     $        psi, Tmunu_flag)

      elseif (decoded_exact_model .eq. EXACT__Schwarzschild_Novikov) then
         call Exact__Schwarzschild_Novikov(x,y,z,t,
     $        gdtt, gdtx, gdty, gdtz,
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz,
     $        guxx, guyy, guzz, guxy, guyz, guxz,
     $        psi, Tmunu_flag)

      elseif (decoded_exact_model .eq. EXACT__Kerr_BoyerLindquist) then
         call Exact__Kerr_BoyerLindquist(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz,
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz,
     $        guxx, guyy, guzz, guxy, guyz, guxz,
     $        psi, Tmunu_flag)

      elseif (decoded_exact_model .eq. EXACT__Kerr_KerrSchild) then
         call Exact__Kerr_KerrSchild(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz,
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz,
     $        guxx, guyy, guzz, guxy, guyz, guxz,
     $        psi, Tmunu_flag)

      elseif (decoded_exact_model .eq. EXACT__Kerr_KerrSchild_spherical) then
         call Exact__Kerr_KerrSchild_spherical(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz,
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz,
     $        guxx, guyy, guzz, guxy, guyz, guxz,
     $        psi, Tmunu_flag)

      elseif (decoded_exact_model .eq. EXACT__Schwarzschild_Lemaitre) then
         call Exact__Schwarzschild_Lemaitre(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz,
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz,
     $        guxx, guyy, guzz, guxy, guyz, guxz,
     $        psi, Tmunu_flag)

      elseif (decoded_exact_model .eq. EXACT__multi_BH) then
         call Exact__multi_BH(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz,
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz,
     $        guxx, guyy, guzz, guxy, guyz, guxz,
     $        psi, Tmunu_flag)

c
c not fully implemented yet -- see Nina Jansen for details
c
c     elseif (decoded_exact_model .eq. EXACT__Alvi) then
c        call Exact__Alvi(
c    $        x, y, z, t,
c    $        gdtt, gdtx, gdty, gdtz,
c    $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
c    $        gutt, gutx, guty, gutz,
c    $        guxx, guyy, guzz, guxy, guyz, guxz,
c    $        psi, Tmunu_flag)
c

      elseif (decoded_exact_model .eq. EXACT__Thorne_fakebinary) then
         call Exact__Thorne_fakebinary(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz,
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz,
     $        guxx, guyy, guzz, guxy, guyz, guxz,
     $        psi, Tmunu_flag)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c cosmological spacetimes
c

      elseif (decoded_exact_model .eq. EXACT__Lemaitre) then
         call Exact__Lemaitre(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz,
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz,
     $        guxx, guyy, guzz, guxy, guyz, guxz,
     $        psi, Tmunu_flag)

      elseif (decoded_exact_model .eq. EXACT__de_Sitter) then
         call Exact__de_Sitter(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz,
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz,
     $        guxx, guyy, guzz, guxy, guyz, guxz,
     $        psi, Tmunu_flag)

      elseif (decoded_exact_model .eq. EXACT__de_Sitter_Lambda) then
         call Exact__de_Sitter_Lambda(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz,
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz,
     $        guxx, guyy, guzz, guxy, guyz, guxz,
     $        psi, Tmunu_flag)

      elseif (decoded_exact_model .eq. EXACT__anti_de_Sitter_Lambda) then
         call Exact__anti_de_Sitter_Lambda(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz,
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz,
     $        guxx, guyy, guzz, guxy, guyz, guxz,
     $        psi, Tmunu_flag)

      elseif (decoded_exact_model .eq. EXACT__Bianchi_I) then
         call Exact__Bianchi_I(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz,
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz,
     $        guxx, guyy, guzz, guxy, guyz, guxz,
     $        psi, Tmunu_flag)

      elseif (decoded_exact_model .eq. EXACT__Goedel) then
         call Exact__Goedel(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz,
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz,
     $        guxx, guyy, guzz, guxy, guyz, guxz,
     $        psi, Tmunu_flag)

      elseif (decoded_exact_model .eq. EXACT__Bertotti) then
         call Exact__Bertotti(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz,
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz,
     $        guxx, guyy, guzz, guxy, guyz, guxz,
     $        psi, Tmunu_flag)

      elseif (decoded_exact_model .eq. EXACT__Kasner_like) then
         call Exact__Kasner_like(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz,
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz,
     $        guxx, guyy, guzz, guxy, guyz, guxz,
     $        psi, Tmunu_flag)

      elseif (decoded_exact_model .eq. EXACT__Kasner_axisymmetric) then
         call Exact__Kasner_axisymmetric(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz,
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz,
     $        guxx, guyy, guzz, guxy, guyz, guxz,
     $        psi, Tmunu_flag)

      elseif (decoded_exact_model .eq. EXACT__Kasner_generalized) then
         call Exact__Kasner_generalized(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz,
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz,
     $        guxx, guyy, guzz, guxy, guyz, guxz,
     $        psi, Tmunu_flag)

      elseif (decoded_exact_model .eq. EXACT__Gowdy_wave) then
         call Exact__Gowdy_wave(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz,
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz,
     $        guxx, guyy, guzz, guxy, guyz, guxz,
     $        psi, Tmunu_flag)

      elseif (decoded_exact_model .eq. EXACT__Milne) then
         call Exact__Milne(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz,
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz,
     $        guxx, guyy, guzz, guxy, guyz, guxz,
     $        psi, Tmunu_flag)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c miscellaneous spacetimes
c

      elseif (decoded_exact_model .eq. EXACT__boost_rotation_symmetric) then
         call Exact__boost_rotation_symmetric(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz,
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz,
     $        guxx, guyy, guzz, guxy, guyz, guxz,
     $        psi, Tmunu_flag)

      elseif (decoded_exact_model .eq. EXACT__bowl) then
         call Exact__bowl(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz,
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz,
     $        guxx, guyy, guzz, guxy, guyz, guxz,
     $        psi, Tmunu_flag)

      elseif (decoded_exact_model .eq. EXACT__constant_density_star) then
         call Exact__constant_density_star(
     $        x, y, z, t,
     $        gdtt, gdtx, gdty, gdtz,
     $        gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $        gutt, gutx, guty, gutz,
     $        guxx, guyy, guzz, guxy, guyz, guxz,
     $        psi, Tmunu_flag)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      else
         write (warn_buffer, '(a,i8)')
     $         'Unknown decoded_exact_model = ', decoded_exact_model
         call CCTK_WARN(0, warn_buffer)
      endif

      return
      end
