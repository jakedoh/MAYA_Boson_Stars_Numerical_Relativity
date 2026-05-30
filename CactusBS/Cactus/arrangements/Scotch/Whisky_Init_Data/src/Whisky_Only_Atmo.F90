 /*@@
   @file      Whisky_Only_Atmo.F90
   @date      Sat Jan 29 04:44:25 2002
   @author    Luca Baiotti
   @desc 
   Hydro initial data for vacuum.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

 /*@@
   @routine    Whisky_Only_Atmo
   @date       Sat Jan 29 04:53:43 2002
   @author     Luca Baiotti
   @desc 
   Hydro initial data for vacuum.
   @enddesc 
   @calls     
   @calledby   
   @history 

   @endhistory 

@@*/

subroutine Whisky_Only_Atmo(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  CCTK_INT i,j,k,nx,ny,nz
  CCTK_REAL det, atmo

  call CCTK_INFO("Setting only atmosphere as initial data.")

  
  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)
  
  if (rho_abs_min > 0.0) then
    atmo = rho_abs_min
    whisky_rho_min = rho_abs_min
  else 
    call CCTK_WARN(0,"For initial data 'only_atmo' the parameter 'rho_abs_min' must be set")
  end if
    
  if (initial_rho_abs_min > 0.0) then
    atmo = initial_rho_abs_min
  else if (initial_rho_rel_min > 0.0) then
    call CCTK_WARN(3,"For initial data 'only_atmo' the parameter 'initial_rho_rel_min' is ignored")
  end if

  if (initial_atmosphere_factor > 0.0) then
    atmo = atmo * initial_atmosphere_factor
  end if

  rho = atmo
  if (rho_abs_min < initial_rho_abs_min) then
    velx = atmosphere_vel(1)
    vely = atmosphere_vel(2)
    velz = atmosphere_vel(3)
  else
    velx = 0.d0
    vely = 0.d0
    velz = 0.d0
  end if
  if (attenuate_atmosphere .ne. 0) then
    velx = velx * (1.d0 - exp(-2.d0 * (r - 1.d0)**2))
    vely = vely * (1.d0 - exp(-2.d0 * (r - 1.d0)**2))
    velz = velz * (1.d0 - exp(-2.d0 * (r - 1.d0)**2))
  end if
  eps = whisky_eps_min

  do i=1,nx
    do j=1,ny
      do k=1,nz

!!$        if (x(i,j,k) < 0.0) then
!!$          velx(i,j,k) = -velx(i,j,k)
!!$          vely(i,j,k) = -vely(i,j,k)
!!$          velz(i,j,k) = -velz(i,j,k)
!!$        end if
        
        call SpatialDeterminant(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),&
             gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),det)
        call prim2con(whisky_eos_handle,gxx(i,j,k),gxy(i,j,k),&
             gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
             det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
             tau(i,j,k),rho(i,j,k),velx(i,j,k),vely(i,j,k),velz(i,j,k),&
             eps(i,j,k),press(i,j,k),w_lorentz(i,j,k))
      enddo
    enddo
  enddo

  densrhs = 0.d0
  sxrhs = 0.d0
  syrhs = 0.d0
  szrhs = 0.d0
  taurhs = 0.d0

  return
  
end subroutine Whisky_Only_Atmo
