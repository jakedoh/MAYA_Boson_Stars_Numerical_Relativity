#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#define Bvecx(i,j,k) Bvec(i,j,k,1)
#define Bvecy(i,j,k) Bvec(i,j,k,2)
#define Bvecz(i,j,k) Bvec(i,j,k,3)
#define Avecx(i,j,k) Avec(i,j,k,1)
#define Avecy(i,j,k) Avec(i,j,k,2)
#define Avecz(i,j,k) Avec(i,j,k,3)

#define SPATIAL_DET(gxx_,gxy_,gxz_,gyy_,gyz_,gzz_) \
                   (-(gxz_)**2*(gyy_) + 2*(gxy_)*(gxz_)*(gyz_) - (gxx_)*(gyz_)**2 - (gxy_)**2*(gzz_) \
                   + (gxx_)*(gyy_)*(gzz_))

subroutine Whisky_CylindricalExplosion(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL :: Bx, By, Bz
  CCTK_REAL :: rho_inner, rho_outer, press_inner, press_outer
  CCTK_REAL :: localrho, localp
  CCTK_REAL :: gam, direction, det
  CCTK_INT  :: nx, ny, nz, i, j, k
  logical   :: set_Afield, set_Bfield 

  if ( CCTK_EQUALS(initial_Bvec,"cylExplosion") ) then
      set_Bfield = .true.
  else
      set_Bfield = .false.
  end if

  if ( CCTK_EQUALS(initial_Avec,"cylExplosion") ) then
      set_Afield = .true.
  else
      set_Afield = .false.
  end if

  !!$ The magnetic field is constant across the entire domain
  !!$ Original tests varied Bx as 0.1,1.0 with By,Bz=0
  Bx = Bx_init
  By = By_init
  Bz = Bz_init

  !!$ Hard-code eos gamma
  gam = (4.d0/3.d0)

  !!$ rho and press defined inside and outside regions
  rho_outer   = 1.d-4
  press_outer = 3.d-5
  rho_inner   = 1.d-2
  press_inner = 1.d0 

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  do k=1,nz
    do j=1,ny
      do i=1,nx

        if ( CCTK_EQUALS(shocktube_type,"xshock")) then
           direction = sqrt( (x(i,j,k)-shock_xpos)**2 + (y(i,j,k)-shock_ypos)**2 )
        else if ( CCTK_EQUALS(shocktube_type,"yshock")) then
           direction = sqrt( (y(i,j,k)-shock_ypos)**2 + (z(i,j,k)-shock_zpos)**2 )
        else if ( CCTK_EQUALS(shocktube_type,"zshock")) then
           direction = sqrt( (x(i,j,k)-shock_xpos)**2 + (z(i,j,k)-shock_zpos)**2 )
        else
           call CCTK_WARN(0,"That shocktube type has not been implemented.")
        end if 
      
        if ( direction.gt.cylexp_radius_outer ) then
           localrho = rho_outer
           localp   = press_outer
        else if ( direction.lt.cylexp_radius_inner ) then
           localrho = rho_inner
           localp   = press_inner
        else !! Smooth exponential linking inner and outer regions
           localrho = exp( ( log(rho_inner)*(cylexp_radius_outer-direction) &
                          +  log(rho_outer)*(direction - cylexp_radius_inner) ) /&
                          ( cylexp_radius_outer - cylexp_radius_inner ) ) 
           localp   = exp( ( log(press_inner)*(cylexp_radius_outer-direction) &
                          +  log(press_outer)*(direction - cylexp_radius_inner) ) /&
                          ( cylexp_radius_outer - cylexp_radius_inner ) ) 
        end if                           

        rho(i,j,k)   = localrho
        press(i,j,k) = localp
        eps(i,j,k)   = localp/( (gam - 1.d0)*localrho )

        velx(i,j,k) = 0.d0
        vely(i,j,k) = 0.d0
        velz(i,j,k) = 0.d0

        if ( set_Bfield ) then
          Bvecx(i,j,k) = Bx
          Bvecy(i,j,k) = By
          Bvecz(i,j,k) = Bz
        end if
        if ( set_Afield ) then
          Avecx(i,j,k) = By*z(i,j,k)
          Avecy(i,j,k) = Bz*x(i,j,k)
          Avecz(i,j,k) = Bx*y(i,j,k)
        end if

        !!$ This should be run on top of Minkowski, but we'll allow some other metric
        det = SPATIAL_DET(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))

        if (CCTK_EQUALS(whisky_eos_type,"Polytype")) then
          call Prim2ConMHDPoly(whisky_eos_handle,gxx(i,j,k),gxy(i,j,k),&
               gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
               det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
               tau(i,j,k),rho(i,j,k),&
               velx(i,j,k),vely(i,j,k),velz(i,j,k),&
               eps(i,j,k),press(i,j,k),w_lorentz(i,j,k),&
               bnx(i,j,k),bny(i,j,k),bnz(i,j,k),&
               Bx,By,Bz)
        else
          call Prim2ConMHDGen(whisky_eos_handle,gxx(i,j,k),gxy(i,j,k),&
               gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
               det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
               tau(i,j,k),rho(i,j,k),&
               velx(i,j,k),vely(i,j,k),velz(i,j,k),&
               eps(i,j,k),press(i,j,k),w_lorentz(i,j,k),&
               bnx(i,j,k),bny(i,j,k),bnz(i,j,k),&
               Bx,By,Bz)
        end if


      end do
    end do
  end do

  

end subroutine Whisky_CylindricalExplosion

