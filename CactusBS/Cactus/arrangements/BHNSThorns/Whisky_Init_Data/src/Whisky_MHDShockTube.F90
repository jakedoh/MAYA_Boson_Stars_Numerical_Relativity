 /*@@
   @file      Whisky_MHDShockTube.F90
   @date      Sat Jan 26 02:53:25 2002
   @author    Tanja Bode
   @desc 
   Initial data of the shock tube type.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

 /*@@
   @routine    Whisky_MHDshocktube
   @date       
   @author     Tanja Bode
   @desc 
   Initial data for magnetized shock tubes. Either diagonal or parallel to
   a coordinate axis. Either Sods problem or the standard shock tube.
   @enddesc 
   @calls     
   @calledby   
   @history 
   Expansion and alteration of the test code from GRAstro_Hydro, 
   written by Mark Miller.
   @endhistory 

@@*/

subroutine Whisky_MHDshocktube(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

#include "EOS_Base.inc"
  
  CCTK_INT :: i, j, k, nx, ny, nz
  CCTK_REAL :: dist_from_plane, pcoord_x, pcoord_y, pcoord_z, det
  CCTK_REAL :: shiftx, shifty, shiftz
  CCTK_REAL :: rhol, rhor, velxl, velxr, velyl, velyr, &
       velzl, velzr, epsl, epsr
  CCTK_REAL :: pressl, pressr, ut,ux,uy,uz
  CCTK_REAL :: bvecxl, bvecyl, bveczl, &
               bvecxr, bvecyr, bveczr
  CCTK_REAL :: avecxl, avecyl, aveczl, &
               avecxr, avecyr, aveczr
  CCTK_REAL, dimension(4) :: ul, ur, Bvecl, Bvecr, upos, Bpos
  logical   :: set_Bfield, set_Afield
  character(len=2) :: shockdir
  CCTK_REAL :: invroot3 = 1./sqrt(3.)

  !! Bfield type
  if ( CCTK_EQUALS(initial_Bvec,"shocktube") ) then
      set_Bfield = .true.
  else
      set_Bfield = .false.
  end if

  if ( CCTK_EQUALS(initial_Avec,"shocktube") ) then
      set_Afield = .true.
      if (CCTK_EQUALS(shocktube_type,"sphere")) then
         call CCTK_WARN(0,'Current spherical shocktube with Avec evolution not set up')
      end if
  else
      set_Afield = .false.
  end if
 

  !! Key to shock direction 
  if (CCTK_EQUALS(shocktube_type,"xshock")) then
     if ( change_shock_direction.eq.0 ) then
        shockdir="+x"
     else
        shockdir="-x"
     endif
  else if (CCTK_EQUALS(shocktube_type,"yshock")) then
     if ( change_shock_direction.eq.0 ) then
        shockdir="+y"
     else
        shockdir="-y"
     endif
  else if (CCTK_EQUALS(shocktube_type,"zshock")) then
     if ( change_shock_direction.eq.0 ) then
        shockdir="+z"
     else
        shockdir="-z"
     endif
  else if (CCTK_EQUALS(shocktube_type,"diagshock")) then
     if ( change_shock_direction.eq.0 ) then
        shockdir="+d"
     else
        shockdir="-d"
     endif
  else if (CCTK_EQUALS(shocktube_type,"sphere")) then
     if ( change_shock_direction.eq.0 ) then
        shockdir="+r"
     else
        shockdir="-r"
     endif
  end if

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

 
  do i=1,nx
    do j=1,ny
      do k=1,nz
        !! d>0=Right, d<0=Left
        if (CCTK_EQUALS(shocktube_type,"diagshock")) then
          dist_from_plane = invroot3*( x(i,j,k) - shock_xpos + &
               y(i,j,k) - shock_ypos + z(i,j,k) - shock_zpos )
          pcoord_x = x(i,j,k) - shock_xpos - dist_from_plane/3.d0
          pcoord_y = y(i,j,k) - shock_ypos - dist_from_plane/3.d0
          pcoord_z = z(i,j,k) - shock_zpos - dist_from_plane/3.d0
        else if (CCTK_EQUALS(shocktube_type,"xshock")) then
          dist_from_plane = x(i,j,k) - shock_xpos
          pcoord_x = dist_from_plane
          pcoord_y = y(i,j,k) - shock_ypos
          pcoord_z = z(i,j,k) - shock_zpos
        else if (CCTK_EQUALS(shocktube_type,"yshock")) then
          dist_from_plane = y(i,j,k) - shock_ypos
          pcoord_x = dist_from_plane
          pcoord_y = z(i,j,k) - shock_zpos
          pcoord_z = x(i,j,k) - shock_xpos
        else if (CCTK_EQUALS(shocktube_type,"zshock")) then
          dist_from_plane = z(i,j,k) - shock_zpos
          pcoord_x = dist_from_plane
          pcoord_y = x(i,j,k) - shock_xpos
          pcoord_z = y(i,j,k) - shock_ypos
        else if (CCTK_EQUALS(shocktube_type,"sphere")) then
          dist_from_plane = sqrt((x(i,j,k)-shock_xpos)**2+&
                           (y(i,j,k)-shock_ypos)**2+&
                           (z(i,j,k)-shock_zpos)**2)-shock_radius
        end if
        if (CCTK_EQUALS(shock_case,"Simple")) then
          rhol = 10.d0
          rhor = 1.d0
          velxl = 0.d0
          velxr = 0.d0
          velyl = 0.d0
          velyr = 0.d0
          velzl = 0.d0
          velzr = 0.d0
          epsl = 2.d0
          epsr = 1.d-6
          pressl = EOS_Pressure(whisky_eos_handle,rhol,epsl)
          pressr = EOS_Pressure(whisky_eos_handle,rhor,epsr)
          bvecxl = 0.d0
          bvecyl = 0.d0
          bveczl = 0.d0
          bvecxr = 0.d0
          bvecyr = 0.d0
          bveczr = 0.d0
          avecxl = 0.d0
          avecyl = 0.d0
          aveczl = 0.d0
          avecxr = 0.d0
          avecyr = 0.d0
          aveczr = 0.d0
        else if (CCTK_EQUALS(shock_case,"Sod")) then
          rhol = 1.d0
          rhor = 0.125d0
          velxl = 0.d0
          velxr = 0.d0
          velyl = 0.d0
          velyr = 0.d0
          velzl = 0.d0
          velzr = 0.d0
          epsl = 1.5d0
          epsr = 1.2d0
          pressl = EOS_Pressure(whisky_eos_handle,rhol,epsl)
          pressr = EOS_Pressure(whisky_eos_handle,rhor,epsr)
          bvecxl = 0.d0
          bvecyl = 0.d0
          bveczl = 0.d0
          bvecxr = 0.d0
          bvecyr = 0.d0
          bveczr = 0.d0
          avecxl = 0.d0
          avecyl = 0.d0
          aveczl = 0.d0
          avecxr = 0.d0
          avecyr = 0.d0
          aveczr = 0.d0
!!$This line only for polytrope, k=1
!!$          epsr = 0.375d0
        else if (CCTK_EQUALS(shock_case,"Blast")) then
          rhol = 1.d0
          rhor = 1.d0
          velxl = 0.d0
          velxr = 0.d0
          velyl = 0.d0
          velyr = 0.d0
          velzl = 0.d0
          velzr = 0.d0
          epsl = 1500.d0
          epsr = 1.5d-2
          pressl = EOS_Pressure(whisky_eos_handle,rhol,epsl)
          pressr = EOS_Pressure(whisky_eos_handle,rhor,epsr)
          bvecxl = 0.d0
          bvecyl = 0.d0
          bveczl = 0.d0
          bvecxr = 0.d0
          bvecyr = 0.d0
          bveczr = 0.d0
          avecxl = 0.d0
          avecyl = 0.d0
          aveczl = 0.d0
          avecxr = 0.d0
          avecyr = 0.d0
          aveczr = 0.d0

        else if (CCTK_EQUALS(shock_case,"DuezMHDA") .or. CCTK_EQUALS(shock_case,"Komissarov1")) then !! Fast Shock ( v(shock)=0.2 )
          ux=25.0d0
          uy=0.0d0
          uz=0.0d0
          bvecxl = 20.0d0
          bvecyl = 25.02d0
          bveczl = 0.0d0
          pressl = 1.0d0
          rhol = 1.0d0

          ut = sqrt( 1. + ux*ux + uy*uy + uz*uz )
          velxl = ux/ut
          velyl = uy/ut
          velzl = uz/ut
          epsl = 3.0d0*pressl/rhol

          ux=1.091d0
          uy=0.3923d0
          uz=0.0d0
          bvecxr = 20.0d0
          bvecyr = 49.0d0
          bveczr = 0.0d0
          pressr = 367.5d0
          rhor = 25.48d0

          ut = sqrt( 1. + ux*ux + uy*uy + uz*uz )
          velxr = ux/ut
          velyr = uy/ut
          velzr = uz/ut
          epsr = 3.0d0*pressr/rhor
          
          avecxl = 0 
          avecyl = bveczl*dist_from_plane
          aveczl = bvecxl*pcoord_y - bvecyl*dist_from_plane
          avecxr = 0 
          avecyr = bveczr*dist_from_plane
          aveczr = bvecxr*pcoord_y - bvecyr*dist_from_plane

        else if (CCTK_EQUALS(shock_case,"DuezMHDB") .or. CCTK_EQUALS(shock_case,"Komissarov2") ) then !! Slow Shock ( v(shock) = 0.5 )
          ux=1.53d0
          uy=0.0d0
          uz=0.0d0
          bvecxl = 10.0d0
          bvecyl = 18.28d0
          bveczl = 0.0d0
          pressl = 10.0d0
          rhol = 1.0d0

          ut = sqrt( 1. + ux*ux + uy*uy + uz*uz )
          velxl = ux/ut
          velyl = uy/ut
          velzl = uz/ut
          epsl = 3.0d0*pressl/rhol

          ux=0.9571d0
          uy=-0.6822d0
          uz=0.0d0
          bvecxr = 10.0d0
          bvecyr = 14.49d0
          bveczr = 0.0d0
          pressr = 55.36d0
          rhor = 3.323d0

          ut = sqrt( 1. + ux*ux + uy*uy + uz*uz )
          velxr = ux/ut
          velyr = uy/ut
          velzr = uz/ut
          epsr = 3.0d0*pressr/rhor
          
          avecxl = 0 
          avecyl = bveczl*dist_from_plane
          aveczl = bvecxl*pcoord_y - bvecyl*dist_from_plane
          avecxr = 0 
          avecyr = bveczr*dist_from_plane
          aveczr = bvecxr*pcoord_y - bvecyr*dist_from_plane

        else if (CCTK_EQUALS(shock_case,"DuezMHDC") .or. CCTK_EQUALS(shock_case,"Komissarov3")) then !! Switch-off Fast Rarefaction
          ux=-2.0d0
          uy=0.0d0
          uz=0.0d0
          bvecxl = 2.0d0
          bvecyl = 0.0d0
          bveczl = 0.0d0
          pressl = 1.0d0
          rhol = 0.1d0

          ut = sqrt( 1. + ux*ux + uy*uy + uz*uz )
          velxl = ux/ut
          velyl = uy/ut
          velzl = uz/ut
          epsl = 3.0d0*pressl/rhol

          ux=-0.212d0
          uy=-0.590d0
          uz=0.0d0
          bvecxr = 2.0d0
          bvecyr = 4.710d0
          bveczr = 0.0d0
          pressr = 10.0d0
          rhor = 0.562d0

          ut = sqrt( 1. + ux*ux + uy*uy + uz*uz )
          velxr = ux/ut
          velyr = uy/ut
          velzr = uz/ut
          epsr = 3.0d0*pressr/rhor
          
          avecxl = 0 
          avecyl = bveczl*dist_from_plane
          aveczl = bvecxl*pcoord_y - bvecyl*dist_from_plane
          avecxr = 0 
          avecyr = bveczr*dist_from_plane
          aveczr = bvecxr*pcoord_y - bvecyr*dist_from_plane

        else if (CCTK_EQUALS(shock_case,"DuezMHDD") .or. CCTK_EQUALS(shock_case,"Komissarov4")) then  !! Switch-off Slow Rarefaction 
          ux=-0.765d0
          uy=-1.386d0
          uz=0.0d0
          bvecxl = 1.0d0
          bvecyl = 1.022d0
          bveczl = 0.0d0
          pressl = 0.1d0
          rhol = 0.00178d0

          ut = sqrt( 1. + ux*ux + uy*uy + uz*uz )
          velxl = ux/ut
          velyl = uy/ut
          velzl = uz/ut
          epsl = 3.0d0*pressl/rhol

          ux=0.0d0
          uy=0.0d0
          uz=0.0d0
          bvecxr = 1.0d0
          bvecyr = 0.0d0
          bveczr = 0.0d0
          pressr = 1.0d0
          rhor = 0.01d0

          ut = sqrt( 1. + ux*ux + uy*uy + uz*uz )
          velxr = ux/ut
          velyr = uy/ut
          velzr = uz/ut
          epsr = 3.0d0*pressr/rhor
          
          avecxl = 0 
          avecyl = bveczl*dist_from_plane
          aveczl = bvecxl*pcoord_y - bvecyl*dist_from_plane
          avecxr = 0 
          avecyr = bveczr*dist_from_plane
          aveczr = bvecxr*pcoord_y - bvecyr*dist_from_plane

        else if (CCTK_EQUALS(shock_case,"DuezMHDE") .or. CCTK_EQUALS(shock_case,"Komissarov7")) then  !! Shock Tube 1
          ux=0.0d0
          uy=0.0d0
          uz=0.0d0
          bvecxl = 1.0d0
          bvecyl = 0.0d0
          bveczl = 0.0d0
          pressl = 1000.0d0
          rhol = 1.0d0

          ut = sqrt( 1. + ux*ux + uy*uy + uz*uz )
          velxl = ux/ut
          velyl = uy/ut
          velzl = uz/ut
          epsl = 3.0d0*pressl/rhol

          ux=0.0d0
          uy=0.0d0
          uz=0.0d0
          bvecxr = 1.0d0
          bvecyr = 0.0d0
          bveczr = 0.0d0
          pressr = 1.0d0
          rhor = 0.1d0

          ut = sqrt( 1. + ux*ux + uy*uy + uz*uz )
          velxr = ux/ut
          velyr = uy/ut
          velzr = uz/ut
          epsr = 3.0d0*pressr/rhor
          
          avecxl = 0 
          avecyl = bveczl*dist_from_plane
          aveczl = bvecxl*pcoord_y - bvecyl*dist_from_plane
          avecxr = 0 
          avecyr = bveczr*dist_from_plane
          aveczr = bvecxr*pcoord_y - bvecyr*dist_from_plane

        else if (CCTK_EQUALS(shock_case,"DuezMHDF") .or. CCTK_EQUALS(shock_case,"Komissarov8")) then  !! Shock Tube 2
          ux=0.0d0
          uy=0.0d0
          uz=0.0d0
          bvecxl = 0.0d0
          bvecyl = 20.0d0
          bveczl = 0.0d0
          pressl = 30.0d0
          rhol = 1.0d0

          ut = sqrt( 1. + ux*ux + uy*uy + uz*uz )
          velxl = ux/ut
          velyl = uy/ut
          velzl = uz/ut
          epsl = 3.0d0*pressl/rhol

          ux=0.0d0
          uy=0.0d0
          uz=0.0d0
          bvecxr = 0.0d0
          bvecyr = 0.0d0
          bveczr = 0.0d0
          pressr = 1.0d0
          rhor = 0.1d0

          ut = sqrt( 1. + ux*ux + uy*uy + uz*uz )
          velxr = ux/ut
          velyr = uy/ut
          velzr = uz/ut
          epsr = 3.0d0*pressr/rhor
          
          avecxl = 0 
          avecyl = bveczl*dist_from_plane
          aveczl = bvecxl*pcoord_y - bvecyl*dist_from_plane
          avecxr = 0 
          avecyr = bveczr*dist_from_plane
          aveczr = bvecxr*pcoord_y - bvecyr*dist_from_plane

        else if (CCTK_EQUALS(shock_case,"DuezMHDG") .or. CCTK_EQUALS(shock_case,"Komissarov9")) then  !! Collision
          ux=5.0d0
          uy=0.0d0
          uz=0.0d0
          bvecxl = 10.0d0
          bvecyl = 10.0d0
          bveczl = 0.0d0
          pressl = 1.0d0
          rhol = 1.0d0

          ut = sqrt( 1. + ux*ux + uy*uy + uz*uz )
          velxl = ux/ut
          velyl = uy/ut
          velzl = uz/ut
          epsl = 3.0d0*pressl/rhol

          ux=-5.0d0
          uy=0.0d0
          uz=0.0d0
          bvecxr = 10.0d0
          bvecyr = -10.0d0
          bveczr = 0.0d0
          pressr = 1.0d0
          rhor = 1.0d0

          ut = sqrt( 1. + ux*ux + uy*uy + uz*uz )
          velxr = ux/ut
          velyr = uy/ut
          velzr = uz/ut
          epsr = 3.0d0*pressr/rhor
          
          avecxl = 0 
          avecyl = bveczl*dist_from_plane
          aveczl = bvecxl*pcoord_y - bvecyl*dist_from_plane
          avecxr = 0 
          avecyr = bveczr*dist_from_plane
          aveczr = bvecxr*pcoord_y - bvecyr*dist_from_plane

        else if (CCTK_EQUALS(shock_case,"DuezMHDH") .or. CCTK_EQUALS(shock_case,"Komissarov5")) then  !! Alfven wave ( v(wave)=0.626 )

          ul=(/ -1.0d0, 0.0d0, 0.0d0, 0.0d0 /) !! u^i
          Bvecl=(/ 0.0d0, 3.0d0, 3.0d0, 0.0d0 /)
          rhol = 1.0d0
          pressl = 1.0d0
          epsl = 3.0d0*pressl/rhol !! p=rho*eps*(gam-1)

          ur=(/ -1.0d0, 3.70d0, 5.76d0, 0.0d0 /) !! u^i
          ur(1)=-1.d0*sqrt( 1.0d0 + ur(2)*ur(2) + ur(3)*ur(3) + ur(4)*ur(4) ) !! u^0=-W
          Bvecr=(/ 0.0d0, 3.0d0, -6.857d0, 0.0d0 /)
          rhor = 1.0d0
          pressr = 1.0d0
          epsr = epsl

          call AlfvenWave( rhol, pressl, ul, Bvecl, ur, Bvecr, pcoord_x, 0.5d0, upos, Bpos ) 

          velxl = -1.d0*upos(2)/upos(1)
          velyl = -1.d0*upos(3)/upos(1)
          velzl = -1.d0*upos(4)/upos(1)
          bvecxl = Bpos(2)
          bvecyl = Bpos(3)
          bveczl = Bpos(4)

          velxr = velxl
          velyr = velyl
          velzr = velzl
          bvecxr = bvecxl
          bvecyr = bvecyl
          bveczr = bveczl

          avecxl = 0 
          avecyl = bveczl*dist_from_plane
          aveczl = bvecxl*pcoord_y - bvecyl*dist_from_plane
          avecxr = 0 
          avecyr = bveczr*dist_from_plane
          aveczr = bvecxr*pcoord_y - bvecyr*dist_from_plane

        else if (CCTK_EQUALS(shock_case,"Komissarov6")) then  !! Compound wave
          ux=0.0d0
          uy=0.0d0
          uz=0.0d0
          bvecxl = 3.0d0
          bvecyl = 3.0d0
          bveczl = 0.0d0
          pressl = 1.0d0
          rhol = 1.0d0

          ut = sqrt( 1. + ux*ux + uy*uy + uz*uz )
          velxl = ux/ut
          velyl = uy/ut
          velzl = uz/ut
          epsl = 3.0d0*pressl/rhol

          ux=3.70d0
          uy=5.76d0
          uz=0.0d0
          bvecxr = 3.0d0
          bvecyr = -6.857d0
          bveczr = 0.0d0
          pressr = 1.0d0
          rhor = 1.0d0

          ut = sqrt( 1. + ux*ux + uy*uy + uz*uz )
          velxr = ux/ut
          velyr = uy/ut
          velzr = uz/ut
          epsr = 3.0d0*pressr/rhor
          
          avecxl = 0 
          avecyl = bveczl*dist_from_plane
          aveczl = bvecxl*pcoord_y - bvecyl*dist_from_plane
          avecxr = 0 
          avecyr = bveczr*dist_from_plane
          aveczr = bvecxr*pcoord_y - bvecyr*dist_from_plane

        else if (CCTK_EQUALS(shock_case,"Balsara1")) then
          rhol = 1.0d0
          velxl = 0.0d0
          velyl = 0.0d0
          velzl = 0.0d0
          bvecxl=0.5d0
          bvecyl=1.0d0
          bveczl=0.0d0
          !!epsl = 1.0d0/rhol
          pressl=1.0d0
          epsl = EOS_SpecificIntEnergy( whisky_eos_handle, rhol, pressl );

          rhor = 0.125d0
          velxr = 0.0d0
          velyr = 0.0d0
          velzr = 0.0d0
          bvecxr=0.5d0
          bvecyr=-1.0d0
          bveczr=0.0d0
          pressr=0.1d0
          epsr = EOS_SpecificIntEnergy( whisky_eos_handle, rhol, pressl );

          avecxl = 0 
          avecyl = bveczl*dist_from_plane
          aveczl = bvecxl*pcoord_y - bvecyl*dist_from_plane
          avecxr = 0 
          avecyr = bveczr*dist_from_plane
          aveczr = bvecxr*pcoord_y - bvecyr*dist_from_plane

       else if (CCTK_EQUALS(shock_case,"Balsara2")) then
          rhol = 1.d0
          rhor = 1.d0
          velxl = 0.d0
          velxr = 0.d0
          velyl = 0.d0
          velyr = 0.d0
          velzl = 0.d0
          velzr = 0.d0
          bvecxl=5.0d0
          bvecxr=5.0d0
          bvecyl=6.d0
          bvecyr=0.7d0
          bveczl=6.d0
          bveczr=0.7d0
          epsl = 1.5d0*30.0d0/rhol
          epsr = 1.5d0*1.0d0/rhor

          avecxl = 0 
          avecyl = bveczl*dist_from_plane
          aveczl = bvecxl*pcoord_y - bvecyl*dist_from_plane
          avecxr = 0 
          avecyr = bveczr*dist_from_plane
          aveczr = bvecxr*pcoord_y - bvecyr*dist_from_plane

!!$ Test 3 (blast wave) of Balsara 2001  -- compare at t=0.4
       else if (CCTK_EQUALS(shock_case,"Balsara3")) then
          rhol = 1.d0
          rhor = 1.d0
          velxl = 0.d0
          velxr = 0.d0
          velyl = 0.d0
          velyr = 0.d0
          velzl = 0.d0
          velzr = 0.d0
          bvecxl=10.0d0
          bvecxr=10.0d0
          bvecyl=7.d0
          bvecyr=0.7d0
          bveczl=7.d0
          bveczr=0.7d0
          epsl = 1.5d0*1000.0d0/rhol
          epsr = 1.5d0*0.1d0/rhor
          pressl = EOS_Pressure(whisky_eos_handle,rhol,epsl)
          pressr = EOS_Pressure(whisky_eos_handle,rhor,epsr)

          avecxl = 0 
          avecyl = bveczl*dist_from_plane
          aveczl = bvecxl*pcoord_y - bvecyl*dist_from_plane
          avecxr = 0 
          avecyr = bveczr*dist_from_plane
          aveczr = bvecxr*pcoord_y - bvecyr*dist_from_plane

!!$ Test 4 (rel. version of Noh 1987) of Balsara 2001  -- compare at t=0.4
       else if (CCTK_EQUALS(shock_case,"Balsara4")) then
          rhol = 1.d0
          rhor = 1.d0
          velxl = 0.999d0
          velxr = -0.999d0
          velyl = 0.d0
          velyr = 0.d0
          velzl = 0.d0
          velzr = 0.d0
          bvecxl=10.0d0
          bvecxr=10.0d0
          bvecyl=7.d0
          bvecyr=-7.d0
          bveczl=7.d0
          bveczr=-7.d0
          epsl = 1.5d0*0.1d0/rhol
          epsr = 1.5d0*0.1d0/rhor
          pressl = EOS_Pressure(whisky_eos_handle,rhol,epsl)
          pressr = EOS_Pressure(whisky_eos_handle,rhor,epsr)

          avecxl = 0 
          avecyl = bveczl*dist_from_plane
          aveczl = bvecxl*pcoord_y - bvecyl*dist_from_plane
          avecxr = 0 
          avecyr = bveczr*dist_from_plane
          aveczr = bvecxr*pcoord_y - bvecyr*dist_from_plane

!!$ Test 5 (non-coplanar set of waves) of Balsara 2001  -- compare at t=0.55
       else if (CCTK_EQUALS(shock_case,"Balsara5")) then
          rhol = 1.08d0
          rhor = 1.d0
          velxl = 0.4d0
          velxr = -0.45d0
          velyl = 0.3d0
          velyr = -0.2d0
          velzl = 0.2d0
          velzr = 0.2d0
          bvecxl=2.0d0
          bvecxr=2.0d0
          bvecyl=0.3d0
          bvecyr=-0.7d0
          bveczl=0.3d0
          bveczr=0.5d0
          epsl = 1.5d0*0.95d0/rhol
          epsr = 1.5d0*1.0d0/rhor
          pressl = EOS_Pressure(whisky_eos_handle,rhol,epsl)
          pressr = EOS_Pressure(whisky_eos_handle,rhor,epsr)

          avecxl = 0 
          avecyl = bveczl*dist_from_plane
          aveczl = bvecxl*pcoord_y - bvecyl*dist_from_plane
          avecxr = 0 
          avecyr = bveczr*dist_from_plane
          aveczr = bvecxr*pcoord_y - bvecyr*dist_from_plane

!!$  "Generic Alfven Test of  Giacomazzo and Rezzolla J.Comp.Phys (2006)
       else if (CCTK_EQUALS(shock_case,"Alfven")) then
          rhol = 1.d0
          rhor = 0.9d0
          velxl = 0.d0
          velxr = 0.d0
          velyl = 0.3d0
          velyr = 0.d0
          velzl = 0.4d0
          velzr = 0.d0
          bvecxl=1.0d0
          bvecxr=1.0d0
          bvecyl=6.d0
          bvecyr=5.d0
          bveczl=2.d0
          bveczr=2.d0
          epsl = 1.5d0*5.d0/rhol
          epsr = 1.5d0*5.3d0/rhor
          pressl = EOS_Pressure(whisky_eos_handle,rhol,epsl)
          pressr = EOS_Pressure(whisky_eos_handle,rhor,epsr)

          avecxl = 0 
          avecyl = bveczl*dist_from_plane
          aveczl = bvecxl*pcoord_y - bvecyl*dist_from_plane
          avecxr = 0 
          avecyr = bveczr*dist_from_plane
          aveczr = bvecxr*pcoord_y - bvecyr*dist_from_plane

        else if (CCTK_EQUALS(shock_case,"Custom")) then
          rhol = left_rho
          rhor = right_rho
          velxl = left_vel(1)
          velxr = right_vel(1) 
          velyl = left_vel(2)
          velyr = right_vel(2)
          velzl = left_vel(3)
          velzr = right_vel(3)
          epsl = left_eps
          epsr = right_eps
          bvecxl = left_bvec(1)
          bvecyl = left_bvec(2)
          bveczl = left_bvec(3)
          bvecxr = right_bvec(1)
          bvecyr = right_bvec(2)
          bveczr = right_bvec(3)
          pressl = EOS_Pressure(whisky_eos_handle,rhol,epsl)
          pressr = EOS_Pressure(whisky_eos_handle,rhor,epsr)
          avecxl = 0 
          avecyl = bveczl*dist_from_plane
          aveczl = bvecxl*pcoord_y - bvecyl*dist_from_plane
          avecxr = 0 
          avecyr = bveczr*dist_from_plane
          aveczr = bvecxr*pcoord_y - bvecyr*dist_from_plane

        else
          call CCTK_WARN(0,"Shock case not recognized")
        end if

        !! This doesn't work for vectors! They need to be rotated!
        !! Assume velx and bvecx to be what we want perpendicular
        !! to the shock!
        if ( ((change_shock_direction==0).and.(dist_from_plane .lt. 0.0d0)).or.& 
             ((change_shock_direction==1).and.(dist_from_plane .gt. 0.0d0)) ) then
          rho(i,j,k) = rhol
          eps(i,j,k) = epsl
          press(i,j,k) = pressl
          velx(i,j,k) = velxl
          vely(i,j,k) = velyl
          velz(i,j,k) = velzl
          bvec(i,j,k,1) = bvecxl
          bvec(i,j,k,2) = bvecyl
          bvec(i,j,k,3) = bveczl
          if ( set_Afield ) then
             avec(i,j,k,1) = avecxl
             avec(i,j,k,2) = avecyl
             avec(i,j,k,3) = aveczl
          end if
        else
          rho(i,j,k) = rhor
          eps(i,j,k) = epsr
          press(i,j,k) = pressr
          velx(i,j,k) = velxr
          vely(i,j,k) = velyr
          velz(i,j,k) = velzr
          bvec(i,j,k,1) = bvecxr
          bvec(i,j,k,2) = bvecyr
          bvec(i,j,k,3) = bveczr
          if ( set_Afield ) then
             avec(i,j,k,1) = avecxr
             avec(i,j,k,2) = avecyr
             avec(i,j,k,3) = aveczr
          end if
        end if

        !! Need to rotate v^i, B^i
        !! Before, each change in case creates a completely different shock system.
        !! Whereas we want to make sure direction doesn't cause issues
        call RotateVector( velx(i,j,k), vely(i,j,k), velz(i,j,k), 1, shockdir, &
                           x(i,j,k), y(i,j,k), z(i,j,k) )
        call RotateVector( bvec(i,j,k,1), bvec(i,j,k,2), bvec(i,j,k,3), -1, shockdir, &
                           x(i,j,k), y(i,j,k), z(i,j,k) )
        if ( set_Afield ) then
           call RotateVector( avec(i,j,k,1), avec(i,j,k,2), avec(i,j,k,3), 1, shockdir, &
                           x(i,j,k), y(i,j,k), z(i,j,k) )
        end if

        shiftx = betax(i,j,k)
        shifty = betax(i,j,k)
        shiftz = betax(i,j,k)

        call SpatialDet(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),&
             gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),det)

        if (CCTK_EQUALS(whisky_eos_type,"Polytype")) then
          call Prim2ConMHDPoly(whisky_eos_handle,gxx(i,j,k),gxy(i,j,k),&
               gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
               det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
               tau(i,j,k),rho(i,j,k),&
               velx(i,j,k),vely(i,j,k),velz(i,j,k),&
               eps(i,j,k),press(i,j,k),w_lorentz(i,j,k),&
               bnx(i,j,k),bny(i,j,k),bnz(i,j,k),&
               bvec(i,j,k,1),bvec(i,j,k,2),bvec(i,j,k,3))
        else
          call Prim2ConMHDGen(whisky_eos_handle,gxx(i,j,k),gxy(i,j,k),&
               gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
               det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
               tau(i,j,k),rho(i,j,k),&
               velx(i,j,k),vely(i,j,k),velz(i,j,k),&
               eps(i,j,k),press(i,j,k),w_lorentz(i,j,k),&
               bnx(i,j,k),bny(i,j,k),bnz(i,j,k),&
               bvec(i,j,k,1),bvec(i,j,k,2),bvec(i,j,k,3))
        end if

    enddo
    enddo
  enddo

  if (clean_divergence.ne.0) then
    divclean_psi = 0.d0
  end if

  densrhs = 0.d0
  sxrhs = 0.d0
  syrhs = 0.d0
  szrhs = 0.d0
  taurhs = 0.d0
  if ( CCTK_EQUALS(evolution_method,"Scotch_MHD_Bvec") ) then
     bnxrhs = 0.d0
     bnyrhs = 0.d0
     bnzrhs = 0.d0
  end if
  if ( CCTK_EQUALS(evolution_method,"Scotch_MHD_Avec") ) then
     avecrhs = 0.d0
  end if

  return
  
end subroutine Whisky_MHDshocktube

 /*@@
   @routine    AlfvenWave
   @date       Fri Feb 10 2012
   @author     Tanja Bode 
   @desc 
   Create an Alfvenwave and return the position relative to the shock, assuming
   the shock is in the direction x and spacetime is Minkowski
   @enddesc 
   @calls     
   @calledby   
   @history 
   @endhistory 

@@*/

subroutine AlfvenWave( rho, press, ul, Bl, ur, Br, xpos, width, upos, Bpos )

   implicit none
   DECLARE_CCTK_PARAMETERS

   CCTK_REAL, intent(in) :: rho, press, xpos, width
   CCTK_REAL, dimension(4), intent(in) :: ul, Bl, ur, Br !! u^i, B^i
   CCTK_REAL, dimension(4), intent(out) :: upos, Bpos

   CCTK_REAL, dimension(4) :: Phi, bul, bupos
   CCTK_REAL, dimension(4) :: bulHF, ulHF, BlHF, uposHF, buposHF
   CCTK_REAL, dimension(4) :: bur, urHF, burHF
   CCTK_REAL, dimension(4,4) :: Lab2HF, HF2Lab
   CCTK_REAL :: a, ScB, ScE, vAlfven
   CCTK_REAL :: wl, uBl, bvecl2, bul2, uBr
   CCTK_REAL :: wHF
   CCTK_REAL :: chi, ay, az, cb, db, invd, thl, thetax, byzth, byc, bzc
   CCTK_REAL, dimension(3,3) :: aij
   CCTK_INT  :: i, j

   !! Hard-coded parameters
   CCTK_REAL :: pi=4.d0 * atan(1.0)
   CCTK_REAL :: gam=4./3., Aph !! Hard-coded EOS and phase-shift
   Aph=pi

   !! Shorthands and b^mu
   wl = -1.d0*ul(1) 
   uBl = ul(2)*Bl(2) + ul(3)*Bl(3) + ul(4)*Bl(4)
   uBr = ur(2)*Br(2) + ur(3)*Br(3) + ur(4)*Br(4)
   bvecl2 = Bl(2)*Bl(2) + Bl(3)*Bl(3) + Bl(4)*Bl(4)
   bul2 = ( bvecl2 + uBl*uBl ) / (wl*wl)
   do i=1,4
      bul(i) = ( Bl(i) - ul(i)*uBl )/wl
      bur(i) = ( Br(i) - ul(i)*uBr )/wl
   end do

   !! Calculate a= u.Phi, B=b.Phi, E=rho*h+b2 with vAlfven from a^2 E - B^2 = 0
   ScE = rho + gam*press/(gam-1.0d0) + bul2
   call solvequad( wl*wl*ScE - uBl*uBl,   2.d0*(wl*ScE*ul(2) + uBl*bul(2) ), &
         ScE*ul(2)*ul(2) - bul(2)*bul(2),   vAlfven )
   Phi = (/ -1.0d0*vAlfven, 1.0d0, 0.d0, 0.d0 /) !! Phi_a
   a = ul(1)*Phi(1) + ul(2)*Phi(2) + ul(3)*Phi(3) + ul(4)*Phi(4)
   ScB = sqrt( a*a*ScE )

   !! Boost to wave from (u', b')|L and compute B'|L with these
   wHF=1.d0/sqrt( 1.d0 - vAlfven*vAlfven )
   Lab2HF=0.d0
   Lab2HF(1,1)=wHF
   Lab2HF(2,1)=-vAlfven*wHF
   Lab2HF(1,2)=-vAlfven*wHF
   Lab2HF(2,2)=wHF
   Lab2HF(3,3)=1.0d0
   Lab2HF(4,4)=1.0d0
   bulHF=0.d0
   ulHF=0.d0
   burHF=0.d0
   urHF=0.d0
   do i=1,4
      do j=1,4
        bulHF(i) = bulHF(i) + Lab2HF(i,j)*bul(j)
        ulHF(i)  = ulHF(i)  + Lab2HF(i,j)* ul(j)
        burHF(i) = burHF(i) + Lab2HF(i,j)*bur(j)
        urHF(i)  = urHF(i)  + Lab2HF(i,j)* ur(j)
      enddo
   enddo
   BlHF=0.d0
   do i=1,4
      BlHF(i) = ulHF(i)*bulHF(1) - ulHF(1)*bulHF(i)
   enddo

   !! Start creating the u(x), b(x) in HF frame: ux and bx are constant across an Alfven wave
   uposHF=0.d0
   buposHF=0.d0 
   uposHF(2) = ulHF(2)
   buposHF(2) = bulHF(2)

   !! Parameters of the ellipse which the transverse components of b lie on  
   chi = uposHF(2)/buposHF(2)
   ay = ( ulHF(3) - chi*bulHF(3))/( ulHF(1) - chi*bulHF(1)) 
   az = ( ulHF(4) - chi*bulHF(4))/( ulHF(1) - chi*bulHF(1)) 
   cb = chi*bul2 / ( ulHF(1) - chi*bulHF(1))
   db = bul2 - buposHF(2)*buposHF(2)
   aij(1,1) = 1.0d0 - ay*ay
   aij(2,2) = 1.0d0 - az*az
   aij(3,3) = -1.d0*( cb*cb + db )
   aij(1,2) = -ay*az
   aij(1,3) = -cb*ay
   aij(2,3) = -cb*az
   aij(2,1) = aij(1,2)
   aij(3,1) = aij(1,3)
   aij(3,2) = aij(2,3)

!!!! CHECKED TO HERE !!!!!

   !! Center of ellipse
   invd = BlHF(2)*BlHF(2) / ( buposHF(2)*buposHF(2) - bul2*uposHF(2)*uposHF(2) )
   byc = cb*invd*ay
   bzc = cb*invd*az

   !! Find the rotation, theta(x), of the transverse components of b^i in HF frame
   thl = atan( (bulHF(4) - bzc) / (bulHF(3)-byc) )
   if ( xpos.lt.(-0.5d0*width) ) then
      thetax = thl
   else if ( xpos.gt.(0.5d0*width) ) then
      thetax = thl + Aph
   else
      thetax = thl + Aph*sin(pi*(xpos+0.5d0*width)/(2.d0*width))**2
   end if 
   byzth = sqrt( (db + cb*cb*invd )/( aij(1,1)*cos(thetax)**2 + &
                 2.d0*aij(1,2)*sin(thetax)*cos(thetax) + &
                 aij(2,2)*sin(thetax)**2 ) )

   !! Compute the rotated magnetic field in HF frame
   buposHF(3) = byc + byzth*cos(thetax) 
   buposHF(4) = bzc + byzth*sin(thetax) 

   !! Fill in rest of uHF(x) from [u] = chi*[b]
   !! Complete u by u.u=-1 and b by u.b=0
   uposHF(3) = ulHF(3) + chi*(buposHF(3)-bulHF(3))
   uposHF(4) = ulHF(4) + chi*(buposHF(4)-bulHF(4))
   uposHF(1) = -1.d0*sqrt( 1.0d0 + uposHF(2)*uposHF(2) + uposHF(3)*uposHF(3) + uposHF(4)*uposHF(4) )
   buposHF(1) = ( uposHF(2)*buposHF(2) + uposHF(3)*buposHF(3) + uposHF(4)*buposHF(4) ) / uposHF(1)

   !! Boost back to wave from (u', b') and compute B(x)
   HF2Lab=0.d0
   HF2Lab(1,1)=wHF
   HF2Lab(2,1)=vAlfven*wHF
   HF2Lab(1,2)=vAlfven*wHF
   HF2Lab(2,2)=wHF
   HF2Lab(3,3)=1.0d0
   HF2Lab(4,4)=1.0d0
   bupos=0.d0
   upos=0.d0
   do i=1,4
      do j=1,4
          bupos(i) = bupos(i) + HF2Lab(i,j)*buposHF(j)
          upos(i)  = upos(i)  + HF2Lab(i,j)* uposHF(j)
      enddo
   enddo
   Bpos=0.d0
   Bpos(2) = (-1.d0*upos(1))*bupos(2) + upos(2)*bupos(1)
   Bpos(3) = (-1.d0*upos(1))*bupos(3) + upos(3)*bupos(1)
   Bpos(4) = (-1.d0*upos(1))*bupos(4) + upos(4)*bupos(1)

   return

end subroutine AlfvenWave

subroutine solvequad( a, b, c, x )

  CCTK_REAL, intent(in) :: a,b,c
  CCTK_REAL, intent(out) :: x
  CCTK_REAL :: del, xplus, xminus

  del = b*b - 4.0d0*a*c

  if ( del.lt.0.d0 ) then
     x = 0
     return
  else if ( del .eq. 0.d0 ) then
     x = - 0.5d0*b/a
     return
  end if

  xplus = ( -b + sqrt(del) )/(2*a)
  xminus = ( -b - sqrt(del) )/(2*a)

  if ( xminus.lt.0.d0 .and. xplus.gt.0.d0 ) then
     x = xminus !! Prefer right-going waves
  else 
     x = xplus !! Prefer faster of two waves
  end if

  return 

end subroutine
