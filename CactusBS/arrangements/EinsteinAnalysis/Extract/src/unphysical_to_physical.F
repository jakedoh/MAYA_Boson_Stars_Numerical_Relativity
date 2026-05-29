
#include "cctk.h"

c     ==================================================================

      SUBROUTINE unphysical_to_physical(grr,grt,grp,gtt,gtp,gpp,dgtt,
     &     dgtp,dgpp,Psis,dPsis,Psi_power)

c     ------------------------------------------------------------------
c
c     Convert unphysical metric components and radial (eta) derivatives 
c     on the 2-sphere into physical quantities, using the conformal factor 
c     on the sphere
c
c     ------------------------------------------------------------------

      IMPLICIT NONE

c     Input variables
      INTEGER :: 
     &    Psi_power
      CCTK_REAL,INTENT(INOUT),DIMENSION (:,:) :: 
     &    grr,grt,grp,gtt,gtp,gpp,dgtt,dgtp,dgpp
      CCTK_REAL,INTENT(IN),DIMENSION (:,:) :: 
     &    Psis,dPsis

c     Output variables
c     WE ARE CHANGING THE INPUT VARIABLES !!!!

c     ------------------------------------------------------------------

      dgtt = Psis**4*dgtt + 4.0d0*Psis**3*dPsis*gtt
      dgtp = Psis**4*dgtp + 4.0d0*Psis**3*dPsis*gtp
      dgpp = Psis**4*dgpp + 4.0d0*Psis**3*dPsis*gpp

      grr  = Psis**4*grr
      grt  = Psis**4*grt
      grp  = Psis**4*grp
      gtt  = Psis**4*gtt
      gtp  = Psis**4*gtp
      gpp  = Psis**4*gpp


      END SUBROUTINE unphysical_to_physical

c     ==================================================================

      








