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

subroutine Whisky_MagneticRotor(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT  :: nx, ny, nz, i, j, k
  CCTK_REAL :: xx, yy, rr, det, rfact
  logical   :: set_Afield, set_Bfield

#include "EOS_Base.inc"

  if ( CCTK_EQUALS(initial_Bvec,"magRotor") ) then
      set_Bfield = .true.
  else
      set_Bfield = .false.
  end if

  if ( CCTK_EQUALS(initial_Avec,"magRotor") ) then
      set_Afield = .true.
  else
      set_Afield = .false.
  end if

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  do k=1,nz
    do j=1,ny
      do i=1,nx

        xx = x(i,j,k)-rotor_xc
        yy = y(i,j,k)-rotor_yc
        rr = sqrt( xx**2 + yy**2 )

        if ( rr.le.rotor_r_rot ) then

           rho(i,j,k)   = rotor_rhoin
           press(i,j,k) = rotor_pressin
           velx(i,j,k)  = -1.d0*rotor_v_max/rotor_r_rot*yy
           vely(i,j,k)  =  1.d0*rotor_v_max/rotor_r_rot*xx
           velz(i,j,k)  =  0.d0 

        else if ( (rotor_use_smoothing.eq.1) .and. &
                ( (rr.gt.rotor_r_rot).and.(rr.le.((1.d0+rotor_rsmooth_rel)*rotor_r_rot))) ) then

           rfact = ( rr/rotor_r_rot - 1.d0 ) / rotor_rsmooth_rel
           rho(i,j,k) = rfact*rotor_rhoout + (1.d0-rfact)*rotor_rhoin
           press(i,j,k) = rfact*rotor_pressout + (1.d0-rfact)*rotor_pressin
           velx(i,j,k) = -1.d0*(1.d0-rfact)*rotor_v_max * yy/rr
           vely(i,j,k) =  1.d0*(1.d0-rfact)*rotor_v_max * xx/rr
           velz(i,j,k) = 0.d0

        else

           rho(i,j,k) = rotor_rhoout
           press(i,j,k) = rotor_pressout
           velx(i,j,k) = 0.d0
           vely(i,j,k) = 0.d0
           velz(i,j,k) = 0.d0

        end if

        eps(i,j,k)   = EOS_SpecificIntEnergy(whisky_eos_handle, rho(i,j,k), press(i,j,k) )

        if ( set_Bfield ) then
          Bvecx(i,j,k) = rotor_bx
          Bvecy(i,j,k) = rotor_by
          Bvecz(i,j,k) = rotor_bz
        end if
        if ( set_Afield ) then
          Avecx(i,j,k) = rotor_by*z(i,j,k)
          Avecy(i,j,k) = rotor_bz*x(i,j,k)
          Avecz(i,j,k) = rotor_bx*y(i,j,k)
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
               rotor_bx,rotor_by,rotor_bz)
        else
          call Prim2ConMHDGen(whisky_eos_handle,gxx(i,j,k),gxy(i,j,k),&
               gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
               det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
               tau(i,j,k),rho(i,j,k),&
               velx(i,j,k),vely(i,j,k),velz(i,j,k),&
               eps(i,j,k),press(i,j,k),w_lorentz(i,j,k),&
               bnx(i,j,k),bny(i,j,k),bnz(i,j,k),&
               rotor_bx,rotor_by,rotor_bz)
        end if


      end do
    end do
  end do

  

end subroutine Whisky_MagneticRotor

