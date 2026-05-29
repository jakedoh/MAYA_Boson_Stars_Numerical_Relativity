! vim: syntax=fortran
! $Header$

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

  subroutine NullGrid_Setup(CCTK_ARGUMENTS)
    use NullGrid_Vars
    implicit none

    TARGET :: null_xb, null_xbh, stereo_q, stereo_p, stereo_pp,&
       zeta, guard_mask, EG_mask, null_rb, null_rbh

    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    gsh = null_gsh
    lsh = null_lsh

    nx = N_radial_pts
    dx = null_dx
    rwt = null_rwt        
    xbin = null_xin

    xb => null_xb
    xbh => null_xbh

    rb => null_rb
    rbh => null_rbh

    lbnd = null_lbnd
    ubnd = null_ubnd
    delta = null_delta
    qs => stereo_q
    ps => stereo_p
    pp => stereo_pp
    zz => zeta

    guard => guard_mask
    EG    => EG_mask

  end subroutine NullGrid_Setup

