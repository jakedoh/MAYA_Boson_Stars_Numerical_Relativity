 /*@@
   @file      Whisky_ShockTube.F90
   @date      Sat Jan 26 02:53:25 2002
   @author    Ian Hawke
   @desc 
   Initial data of the shock tube type.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

 /*@@
   @routine    Whisky_shocktube
   @date       Sat Jan 26 02:53:49 2002
   @author     Ian Hawke
   @desc 
   Initial data for shock tubes. Either diagonal or parallel to
   a coordinate axis. Either Sods problem or the standard shock tube.
   @enddesc 
   @calls     
   @calledby   
   @history 
   Expansion and alteration of the test code from GRAstro_Hydro, 
   written by Mark Miller.
   @endhistory 

@@*/

subroutine Whisky_shocktube(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT :: i, j, k, nx, ny, nz
  CCTK_REAL :: direction, det
  CCTK_REAL :: rhol, rhor, velxl, velxr, velyl, velyr, &
       velzl, velzr, epsl, epsr

#include "EOS_Base.inc"
  !! Key to shock direction 
  character(len=2) :: shockdir
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
        if (CCTK_EQUALS(shocktube_type,"diagshock")) then
          direction = x(i,j,k) - shock_xpos + &
               y(i,j,k) - shock_ypos + z(i,j,k) - shock_zpos
        else if (CCTK_EQUALS(shocktube_type,"xshock")) then
          direction = x(i,j,k) - shock_xpos
        else if (CCTK_EQUALS(shocktube_type,"yshock")) then
          direction = y(i,j,k) - shock_ypos
        else if (CCTK_EQUALS(shocktube_type,"zshock")) then
          direction = z(i,j,k) - shock_zpos
        else if (CCTK_EQUALS(shocktube_type,"sphere")) then
          direction = sqrt((x(i,j,k)-shock_xpos)**2+&
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
        else
          call CCTK_WARN(0,"Shock case not recognized.  If you wanted MHD, change whisky_initial_data.")
        end if

        if ( ((change_shock_direction==0).and.(direction .lt. 0.0d0)).or.& 
             ((change_shock_direction==1).and.(direction .gt. 0.0d0)) ) then
          rho(i,j,k) = rhol
          velx(i,j,k) = velxl
          vely(i,j,k) = velyl
          velz(i,j,k) = velzl
          eps(i,j,k) = epsl
        else
          rho(i,j,k) = rhor
          velx(i,j,k) = velxr
          vely(i,j,k) = velyr
          velz(i,j,k) = velzr
          eps(i,j,k) = epsr
        end if

        !! Need to rotate v^i
        !! Before, each change in case creates a completely different shock system.
        !! Whereas we want to make sure direction doesn't cause issues
        call RotateVector( velx(i,j,k), vely(i,j,k), velz(i,j,k), 1, shockdir, &
                           x(i,j,k), y(i,j,k), z(i,j,k) )

        call SpatialDet(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),&
             gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),det)

        if (CCTK_EQUALS(whisky_eos_type,"Polytype")) then
          call Prim2ConPoly(whisky_eos_handle,gxx(i,j,k),gxy(i,j,k),&
               gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
               det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
               tau(i,j,k),rho(i,j,k),&
               velx(i,j,k),vely(i,j,k),velz(i,j,k),&
               eps(i,j,k),press(i,j,k),w_lorentz(i,j,k))
        else
          call Prim2ConGen(whisky_eos_handle,gxx(i,j,k),gxy(i,j,k),&
               gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
               det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
               tau(i,j,k),rho(i,j,k),&
               velx(i,j,k),vely(i,j,k),velz(i,j,k),&
               eps(i,j,k),press(i,j,k),w_lorentz(i,j,k))
        end if
    enddo
    enddo
  enddo

  densrhs = 0.d0
  sxrhs = 0.d0
  syrhs = 0.d0
  szrhs = 0.d0
  taurhs = 0.d0

  return
  
end subroutine Whisky_shocktube
