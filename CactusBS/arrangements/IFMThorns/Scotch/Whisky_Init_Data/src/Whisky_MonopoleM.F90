 /*@@
   @file      Whisky_MonopoleM.F90
   @date      Wed Feb  8 2012 
   @author    Tanja Bode 
   @desc 
   Monopole initial data.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

 /*@@
   @routine    Whisky_Monopole
   @date       Wed Feb  8 2012 
   @author     Tanja Bode
   @desc 
   Initial data for magnetized tests, monopole.
   @enddesc 
   @calls     
   @calledby   
   @history 
   Expansion and alteration of the test code from GRAstro_Hydro, 
   written by Mark Miller.
   @endhistory 

@@*/

subroutine Whisky_Monopole(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

#include "EOS_Base.inc"
  
  CCTK_INT :: i, j, k, nx, ny, nz
  CCTK_REAL :: rpole2, det; 
  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)
  
  
  do i=1,nx
    do j=1,ny
      do k=1,nz

        !! Base system on top of which we add the monopole
        rho(i,j,k) = 1.0
        eps(i,j,k) = 0.1
        velx(i,j,k) = 0d0
        vely(i,j,k) = 0d0
        velz(i,j,k) = 0d0
        bvec(i,j,k,1) = 0d0
        bvec(i,j,k,2) = 0d0
        bvec(i,j,k,3) = 0d0

        rpole2 = (x(i,j,k)-monopole_pos(1))**2 + &
                 (y(i,j,k)-monopole_pos(2))**2 + &
                 (z(i,j,k)-monopole_pos(3))**2

        !! Add monopole as pure Bvecx.
        if ( CCTK_EQUALS(monopole_type, "Point" )) then
           if ( sqrt(rpole2).lt.(0.5d0*CCTK_DELTA_SPACE(1))) then
              bvec(i,j,k,3) = bvec(i,j,k,3) + monopole_val 
           end if
        else if (CCTK_EQUALS(monopole_type, "Gaussian" )) then
           if ( sqrt(rpole2).le.monopole_sig ) then
              bvec(i,j,k,3) = bvec(i,j,k,3) + &
                  monopole_val*(exp( -rpole2/(monopole_sig)**2 ) - exp(-1.d0))
           end if
        else
           call CCTK_WARN(0,"Unrecognized monopole_type")
        end if 

        !! Create Conservatives from these primitives
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
  bnxrhs = 0.d0
  bnyrhs = 0.d0
  bnzrhs = 0.d0

  return
  
end subroutine Whisky_Monopole
