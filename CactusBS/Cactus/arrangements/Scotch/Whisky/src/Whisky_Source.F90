 /*@@
   @file      Whisky_Source.F90
   @date      Sat Jan 26 02:03:56 2002
   @author    Ian Hawke
   @desc 
   The geometric source terms for the matter evolution
   @enddesc 
 @@*/

! Second order f.d.

#define DIFF_X_2(q) 0.5d0 * (q(i+1,j,k) - q(i-1,j,k)) * idx
#define DIFF_Y_2(q) 0.5d0 * (q(i,j+1,k) - q(i,j-1,k)) * idy
#define DIFF_Z_2(q) 0.5d0 * (q(i,j,k+1) - q(i,j,k-1)) * idz

! Fourth order f.d.

#define DIFF_X_4(q) (-q(i+2,j,k) + 8.d0 * q(i+1,j,k) - 8.d0 * q(i-1,j,k) + \
                      q(i-2,j,k)) / 12.d0 * idx
#define DIFF_Y_4(q) (-q(i,j+2,k) + 8.d0 * q(i,j+1,k) - 8.d0 * q(i,j-1,k) + \
                      q(i,j-2,k)) / 12.d0 * idy
#define DIFF_Z_4(q) (-q(i,j,k+2) + 8.d0 * q(i,j,k+1) - 8.d0 * q(i,j,k-1) + \
                      q(i,j,k-2)) / 12.d0 * idz

! B-field
#define Bvx(i,j,k) Bvec(i,j,k,1)
#define Bvy(i,j,k) Bvec(i,j,k,2)
#define Bvz(i,j,k) Bvec(i,j,k,3)

! A-field
#define Ax(i,j,k) Avec(i,j,k,1)
#define Ay(i,j,k) Avec(i,j,k,2)
#define Az(i,j,k) Avec(i,j,k,3)
#define Axrhs(i,j,k) Avecrhs(i,j,k,1)
#define Ayrhs(i,j,k) Avecrhs(i,j,k,2)
#define Azrhs(i,j,k) Avecrhs(i,j,k,3)

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "Whisky_Utils.h"

 /*@@
   @routine    SourceTerms
   @date       Sat Jan 26 02:04:21 2002
   @author     Ian Hawke
   @desc 
   Calculate the geometric source terms and add to the update GFs
   @enddesc 
   @calls     
   @calledby   
   @history 
   Minor alterations of routine from GR3D.
   @endhistory 

@@*/

subroutine SourceTerms(CCTK_ARGUMENTS)
      
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  CCTK_INT :: i, j, k, nx, ny, nz
  CCTK_REAL :: one, two, half
  CCTK_REAL :: t00, t0x, t0y, t0z, txx, txy, txz, tyy, tyz, tzz
  CCTK_REAL :: sqrtdet, det, uxx, uxy, uxz, uyy, uyz, uzz, rhoenthalpyW2
  CCTK_REAL :: shiftx, shifty, shiftz, velxshift, velyshift, velzshift 
  CCTK_REAL :: vlowx, vlowy, vlowz
  CCTK_REAL :: dx_betax, dx_betay, dx_betaz, dy_betax, dy_betay,&
       dy_betaz, dz_betax, dz_betay, dz_betaz
  CCTK_REAL :: dx_alp, dy_alp, dz_alp
  CCTK_REAL :: dx_psi, dy_psi, dz_psi, djBnj, gdxg, gdyg, gdzg
  CCTK_REAL :: div_ug_lx, div_ug_ly, div_ug_lz 
  CCTK_REAL :: tau_source, sx_source, sy_source, sz_source
  CCTK_REAL :: localgxx,localgxy,localgxz,localgyy,localgyz,localgzz
  CCTK_REAL :: localalp, localpress, localeps
  CCTK_REAL :: localvx, localvy, localvz, local_wlorentz
  CCTK_REAL :: localBnx, localBny, localBnz
  CCTK_REAL :: localBx, localBy, localBz, local_dcpsi
  CCTK_REAL :: localAx, localAy, localAz, localAphi
  CCTK_REAL :: upAx, upAy, upAz
  CCTK_REAL :: dx_Ax, dx_Ay, dx_Az, dy_Ax, dy_Ay, dy_Az
  CCTK_REAL :: dz_Ax, dz_Ay, dz_Az, dx_Aphi, dy_Aphi, dz_Aphi
  CCTK_REAL :: dx_gxx, dx_gxy, dx_gxz, dx_gyy, dx_gyz, dx_gzz
  CCTK_REAL :: dy_gxx, dy_gxy, dy_gxz, dy_gyy, dy_gyz, dy_gzz
  CCTK_REAL :: dz_gxx, dz_gxy, dz_gxz, dz_gyy, dz_gyz, dz_gzz
  CCTK_REAL :: dx, dy, dz, idx, idy, idz
  CCTK_REAL :: shiftshiftk, shiftkx, shiftky, shiftkz
  CCTK_REAL :: sumTK
  CCTK_REAL :: halfshiftdgx, halfshiftdgy, halfshiftdgz
  CCTK_REAL :: halfTdgx, halfTdgy, halfTdgz
  CCTK_INT :: whisky_local_spatial_order
  CCTK_REAL :: udotB, bx, by, bz, b2, w2, blowx, blowy, blowz
  CCTK_REAL :: rr

  logical, allocatable, dimension (:,:,:) :: force_spatial_second_order
  CHARACTER(len=400) :: warnline

  one = 1.0d0
  two = 2.0d0
  half = 0.5d0
  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)
  dx = CCTK_DELTA_SPACE(1)
  dy = CCTK_DELTA_SPACE(2)
  dz = CCTK_DELTA_SPACE(3)
  idx = 1.d0/dx
  idy = 1.d0/dy
  idz = 1.d0/dz
 
!!$  Initialize the update terms to be zero.
!!$  This will guarantee that no garbage in the boundaries is updated.

  densrhs = 0.d0
  sxrhs = 0.d0
  syrhs = 0.d0
  szrhs = 0.d0
  taurhs = 0.d0
  if (whisky_mhd_handle.eq.2) then
    Bnxrhs = 0.d0
    Bnyrhs = 0.d0
    Bnzrhs = 0.d0
  else if (whisky_mhd_handle.ge.3) then
    Avecrhs = 0.d0 
  end if

  if (evolve_tracer .ne. 0) then

    cons_tracerrhs = 0.d0

  end if
  
!!$  Set up the array for checking the order. We switch to second order
!!$  differencing at boundaries and near excision regions.
!!$  Copied straight from BSSN.

  allocate (force_spatial_second_order(nx,ny,nz))
  force_spatial_second_order = .FALSE.
  
  if (spatial_order > 2) then
    !$OMP PARALLEL DO PRIVATE (i,j,k)
    do k = 1 + whisky_stencil, nz - whisky_stencil
      do j = 1 + whisky_stencil, ny - whisky_stencil
        do i = 1 + whisky_stencil, nx - whisky_stencil
          if ((i < 3).or.(i > cctk_lsh(1) - 2).or. &
               (j < 3).or.(j > cctk_lsh(2) - 2).or. &
               (k < 3).or.(k > cctk_lsh(3) - 2) ) then
            force_spatial_second_order(i,j,k) = .TRUE.
          else if ( use_mask > 0 ) then
            if (minval(emask(i-2:i+2,j-2:j+2,k-2:k+2)) < 0.75d0) then
              force_spatial_second_order(i,j,k) = .TRUE.
            end if
          end if
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end if
  
  !$OMP PARALLEL DO PRIVATE (i, j, k, &
  !$OMP t00, t0x, t0y, t0z, txx, txy, txz, tyy, tyz, tzz, &
  !$OMP sqrtdet, det, uxx, uxy, uxz, uyy, uyz, uzz, rhoenthalpyW2, &
  !$OMP shiftx, shifty, shiftz, velxshift, velyshift, velzshift , &
  !$OMP vlowx, vlowy, vlowz, &
  !$OMP dx_betax, dx_betay, dx_betaz, dy_betax, dy_betay, &
  !$OMP dy_betaz, dz_betax, dz_betay, dz_betaz, &
  !$OMP dx_alp, dy_alp, dz_alp, &
  !$OMP tau_source, sx_source, sy_source, sz_source, &
  !$OMP localalp,localgxx,localgxy,localgxz,localgyy,localgyz,localgzz, &
  !$OMP localpress, localeps, localvx, localvy, localvz, local_wlorentz, &
  !$OMP localBnx, localBny, localBnz, localBx, localBy, localBz, local_dcpsi, &
  !$OMP localAx, localAy, localAz, localAphi, &
  !$OMP dx_Ax, dx_Ay, dx_Az, &
  !$OMP dy_Ax, dy_Ay, dy_Az, &
  !$OMP dz_Ax, dz_Ay, dz_Az, &
  !$OMP dx_Aphi, dy_Aphi, dz_Aphi, &
  !$OMP upAx, upAy, upAz, &
  !$OMP dx_gxx, dx_gxy, dx_gxz, dx_gyy, dx_gyz, dx_gzz, &
  !$OMP dy_gxx, dy_gxy, dy_gxz, dy_gyy, dy_gyz, dy_gzz, &
  !$OMP dz_gxx, dz_gxy, dz_gxz, dz_gyy, dz_gyz, dz_gzz, &
  !$OMP shiftshiftk, shiftkx, shiftky, shiftkz, &
  !$OMP sumTK, &
  !$OMP b2,w2,udotB,djbnj, bx,by,bz,blowx,blowy,blowz,warnline, &
  !$OMP gdxg, gdyg, gdzg, div_ug_lx,div_ug_ly,div_ug_lz,&
  !$OMP halfshiftdgx, halfshiftdgy, halfshiftdgz, &
  !$OMP halfTdgx, halfTdgy, halfTdgz, whisky_local_spatial_order)
  do k=1 + whisky_stencil,nz - whisky_stencil
    do j=1 + whisky_stencil,ny - whisky_stencil
      do i=1 + whisky_stencil,nx - whisky_stencil

        whisky_local_spatial_order = spatial_order
        if (force_spatial_second_order(i,j,k)) then
          whisky_local_spatial_order = 2
        end if
        
!!$        Set the metric terms.

        localalp = alp(i,j,k)
        localgxx = gxx(i,j,k)
        localgxy = gxy(i,j,k)
        localgxz = gxz(i,j,k)
        localgyy = gyy(i,j,k)
        localgyz = gyz(i,j,k)
        localgzz = gzz(i,j,k)

!!$     Set local hydro variables
        localpress = press(i,j,k)
        localeps   = eps(i,j,k)
        localvx    = velx(i,j,k)
        localvy    = vely(i,j,k)
        localvz    = velz(i,j,k)
        local_wlorentz = w_lorentz(i,j,k)

        if ( whisky_mhd_handle.gt.0 ) then
           localBx = Bvec(i,j,k,1)
           localBy = Bvec(i,j,k,2)
           localBz = Bvec(i,j,k,3)
           if ( whisky_mhd_handle.gt.1 ) then
              localBnx = Bnx(i,j,k)
              localBny = Bny(i,j,k)
              localBnz = Bnz(i,j,k)
           end if 
        end if

        det = SPATIAL_DET(localgxx,localgxy,localgxz,localgyy,localgyz,localgzz)
        sqrtdet = sqrt(det)
        call UpperMetric(uxx, uxy, uxz, uyy, uyz, uzz, det, localgxx,&
             localgxy, localgxz, localgyy, localgyz, localgzz)
        
!!$        All the matter variables (except velocity) always appear
!!$        together in this form

        rhoenthalpyW2 = (rho(i,j,k)*(one + localeps) + localpress)*&
             local_wlorentz**2
        
        shiftx = betax(i,j,k)
        shifty = betay(i,j,k)
        shiftz = betaz(i,j,k)

        if (whisky_local_spatial_order .eq. 2) then

          dx_betax = DIFF_X_2(betax)
          dx_betay = DIFF_X_2(betay)
          dx_betaz = DIFF_X_2(betaz)
            
          dy_betax = DIFF_Y_2(betax)
          dy_betay = DIFF_Y_2(betay)
          dy_betaz = DIFF_Y_2(betaz)
            
          dz_betax = DIFF_Z_2(betax)
          dz_betay = DIFF_Z_2(betay)
          dz_betaz = DIFF_Z_2(betaz)

        else

          dx_betax = DIFF_X_4(betax)
          dx_betay = DIFF_X_4(betay)
          dx_betaz = DIFF_X_4(betaz)
            
          dy_betax = DIFF_Y_4(betax)
          dy_betay = DIFF_Y_4(betay)
          dy_betaz = DIFF_Y_4(betaz)
            
          dz_betax = DIFF_Z_4(betax)
          dz_betay = DIFF_Z_4(betay)
          dz_betaz = DIFF_Z_4(betaz)

        end if
          
        velxshift = localvx - shiftx/localalp
        velyshift = localvy - shifty/localalp
        velzshift = localvz - shiftz/localalp
        vlowx = localvx*localgxx + localvy*localgxy +&
             localvz*localgxz
        vlowy = localvx*localgxy + localvy*localgyy +&
             localvz*localgyz
        vlowz = localvx*localgxz + localvy*localgyz +&
             localvz*localgzz

!!$        For a change, these are T^{ij}

        t00 = (rhoenthalpyW2 - localpress)/(localalp**2)
        t0x = rhoenthalpyW2*velxshift/localalp +&
             localpress*shiftx/(localalp**2)
        t0y = rhoenthalpyW2*velyshift/localalp +&
             localpress*shifty/(localalp**2)
        t0z = rhoenthalpyW2*velzshift/localalp +&
             localpress*shiftz/(localalp**2)
        txx = rhoenthalpyW2*velxshift*velxshift +&
             localpress*(uxx - shiftx*shiftx/(localalp**2))
        txy = rhoenthalpyW2*velxshift*velyshift +&
             localpress*(uxy - shiftx*shifty/(localalp**2))
        txz = rhoenthalpyW2*velxshift*velzshift +&
             localpress*(uxz - shiftx*shiftz/(localalp**2))
        tyy = rhoenthalpyW2*velyshift*velyshift +&
             localpress*(uyy - shifty*shifty/(localalp**2))
        tyz = rhoenthalpyW2*velyshift*velzshift +&
             localpress*(uyz - shifty*shiftz/(localalp**2))
        tzz = rhoenthalpyW2*velzshift*velzshift +&
             localpress*(uzz - shiftz*shiftz/(localalp**2))

        !! Add B-field contribution to Tmunu
        if ( whisky_mhd_handle.gt.0 ) then 

          w2 = local_wlorentz**2
          udotB = local_wlorentz * ( vlowx*localBx &
                    + vlowy*localBy + vlowz*localBz )
          b2 = ( localgxx*localBx**2 + localgyy*localBy**2 &
                 + localgzz*localBz**2 &
                 + 2.*localgxy*localBx*localBy &
                 + 2.*localgxz*localBx*localBz &
                 + 2.*localgyz*localBy*localBz &
                 + udotB**2 ) / w2
          bx = localBx/local_wlorentz + udotB*velxshift
          by = localBy/local_wlorentz + udotB*velyshift
          bz = localBz/local_wlorentz + udotB*velzshift
          ! b_mu is a four vector so its covariant components are not just
          ! gamma_ij b^i
          blowx = (localBx*localgxx + localBy*localgxy +&
               localBz*localgxz) / local_wlorentz + udotB *vlowx
          blowy = (localBx*localgxy + localBy*localgyy +&
               localBz*localgyz) / local_wlorentz + udotB *vlowy
          blowz = (localBx*localgxz + localBy*localgyz +&
               localBz*localgzz) / local_wlorentz + udotB *vlowz

          t00 = t00 + ( b2*( w2 - 0.5d0 ) &
                   - udotB**2 ) / localalp**2
          t0x = t0x + (b2*( w2*velxshift + 0.5d0*shiftx/localalp ) &
                            - udotB*bx)/localalp
          t0y = t0y + (b2*( w2*velyshift + 0.5d0*shifty/localalp ) &
                            - udotB*by)/localalp
          t0z = t0z + (b2*( w2*velzshift + 0.5d0*shiftz/localalp ) &
                            - udotB*bz)/localalp
          txx = txx + b2*( w2*velxshift**2 + 0.5d0*(uxx - shiftx*shiftx/localalp**2) ) - bx**2
          txy = txy + b2*( w2*velxshift*velyshift + 0.5d0*(uxy - shiftx*shifty/localalp**2) ) - bx*by
          txz = txz + b2*( w2*velxshift*velzshift + 0.5d0*(uxz - shiftx*shiftz/localalp**2) ) - bx*bz
          tyy = tyy + b2*( w2*velyshift**2 + 0.5d0*(uyy - shifty**2/localalp**2) ) - by**2
          tyz = tyz + b2*( w2*velyshift*velzshift + 0.5d0*(uyz - shifty*shiftz/localalp**2) ) - by*bz
          tzz = tzz + b2*( w2*velzshift**2 + 0.5d0*(uzz - shiftz**2/localalp**2) ) - bz**2

        end if

!!$        Derivatives of the lapse, metric and shift

        if (whisky_local_spatial_order .eq. 2) then

          dx_alp = DIFF_X_2(alp)
          dy_alp = DIFF_Y_2(alp)
          dz_alp = DIFF_Z_2(alp)

          dx_gxx = DIFF_X_2(gxx)
          dx_gxy = DIFF_X_2(gxy)
          dx_gxz = DIFF_X_2(gxz)
          dx_gyy = DIFF_X_2(gyy)
          dx_gyz = DIFF_X_2(gyz)
          dx_gzz = DIFF_X_2(gzz)
          dy_gxx = DIFF_Y_2(gxx)
          dy_gxy = DIFF_Y_2(gxy)
          dy_gxz = DIFF_Y_2(gxz)
          dy_gyy = DIFF_Y_2(gyy)
          dy_gyz = DIFF_Y_2(gyz)
          dy_gzz = DIFF_Y_2(gzz)
          dz_gxx = DIFF_Z_2(gxx)
          dz_gxy = DIFF_Z_2(gxy)
          dz_gxz = DIFF_Z_2(gxz)
          dz_gyy = DIFF_Z_2(gyy)
          dz_gyz = DIFF_Z_2(gyz)
          dz_gzz = DIFF_Z_2(gzz)

        else

          dx_alp = DIFF_X_4(alp)
          dy_alp = DIFF_Y_4(alp)
          dz_alp = DIFF_Z_4(alp)

          dx_gxx = DIFF_X_4(gxx)
          dx_gxy = DIFF_X_4(gxy)
          dx_gxz = DIFF_X_4(gxz)
          dx_gyy = DIFF_X_4(gyy)
          dx_gyz = DIFF_X_4(gyz)
          dx_gzz = DIFF_X_4(gzz)
          dy_gxx = DIFF_Y_4(gxx)
          dy_gxy = DIFF_Y_4(gxy)
          dy_gxz = DIFF_Y_4(gxz)
          dy_gyy = DIFF_Y_4(gyy)
          dy_gyz = DIFF_Y_4(gyz)
          dy_gzz = DIFF_Y_4(gzz)
          dz_gxx = DIFF_Z_4(gxx)
          dz_gxy = DIFF_Z_4(gxy)
          dz_gxz = DIFF_Z_4(gxz)
          dz_gyy = DIFF_Z_4(gyy)
          dz_gyz = DIFF_Z_4(gyz)
          dz_gzz = DIFF_Z_4(gzz)

        end if
          
!!$        Contract the shift with the extrinsic curvature

        shiftshiftk = shiftx*shiftx*kxx(i,j,k) + &
                      shifty*shifty*kyy(i,j,k) + &
                      shiftz*shiftz*kzz(i,j,k) + &
             two*(shiftx*shifty*kxy(i,j,k) + &
                  shiftx*shiftz*kxz(i,j,k) + &
                  shifty*shiftz*kyz(i,j,k))

        shiftkx = shiftx*kxx(i,j,k) + shifty*kxy(i,j,k) + shiftz*kxz(i,j,k)
        shiftky = shiftx*kxy(i,j,k) + shifty*kyy(i,j,k) + shiftz*kyz(i,j,k)
        shiftkz = shiftx*kxz(i,j,k) + shifty*kyz(i,j,k) + shiftz*kzz(i,j,k)

!!$        Contract the matter terms with the extrinsic curvature

        sumTK = txx*kxx(i,j,k) + tyy*kyy(i,j,k) + tzz*kzz(i,j,k) &
             + two*(txy*kxy(i,j,k) + txz*kxz(i,j,k) + tyz*kyz(i,j,k))

!!$        Update term for tau
        
        tau_source = t00* &
             (shiftshiftk - (shiftx*dx_alp + shifty*dy_alp + shiftz*dz_alp) )&
             + t0x*(-dx_alp + two*shiftkx) &
             + t0y*(-dy_alp + two*shiftky) &
             + t0z*(-dz_alp + two*shiftkz) &
             + sumTK

!!$        The following looks very little like the terms in the
!!$        standard papers. Take a look in the ThornGuide to see why
!!$        it is really the same thing.

!!$        Contract the shift with derivatives of the metric

        halfshiftdgx = half*(shiftx*shiftx*dx_gxx + &
             shifty*shifty*dx_gyy + shiftz*shiftz*dx_gzz) + &
             shiftx*shifty*dx_gxy + shiftx*shiftz*dx_gxz + &
             shifty*shiftz*dx_gyz
        halfshiftdgy = half*(shiftx*shiftx*dy_gxx + &
             shifty*shifty*dy_gyy + shiftz*shiftz*dy_gzz) + &
             shiftx*shifty*dy_gxy + shiftx*shiftz*dy_gxz + &
             shifty*shiftz*dy_gyz
        halfshiftdgz = half*(shiftx*shiftx*dz_gxx + &
             shifty*shifty*dz_gyy + shiftz*shiftz*dz_gzz) + &
             shiftx*shifty*dz_gxy + shiftx*shiftz*dz_gxz + &
             shifty*shiftz*dz_gyz

!!$        Contract the matter with derivatives of the metric

        halfTdgx = half*(txx*dx_gxx + tyy*dx_gyy + tzz*dx_gzz) +&
             txy*dx_gxy + txz*dx_gxz + tyz*dx_gyz
        halfTdgy = half*(txx*dy_gxx + tyy*dy_gyy + tzz*dy_gzz) +&
             txy*dy_gxy + txz*dy_gxz + tyz*dy_gyz
        halfTdgz = half*(txx*dz_gxx + tyy*dz_gyy + tzz*dz_gzz) +&
             txy*dz_gxy + txz*dz_gxz + tyz*dz_gyz

        sx_source = t00*&
             (halfshiftdgx - localalp*dx_alp) +&
             t0x*(shiftx*dx_gxx + shifty*dx_gxy + shiftz*dx_gxz) +&
             t0y*(shiftx*dx_gxy + shifty*dx_gyy + shiftz*dx_gyz) +&
             t0z*(shiftx*dx_gxz + shifty*dx_gyz + shiftz*dx_gzz) +&
             halfTdgx + rhoenthalpyW2*&
             (vlowx*dx_betax + vlowy*dx_betay + vlowz*dx_betaz)/&
             localalp
        sy_source = t00*&
             (halfshiftdgy - localalp*dy_alp) +&
             t0x*(shiftx*dy_gxx + shifty*dy_gxy + shiftz*dy_gxz) +&
             t0y*(shiftx*dy_gxy + shifty*dy_gyy + shiftz*dy_gyz) +&
             t0z*(shiftx*dy_gxz + shifty*dy_gyz + shiftz*dy_gzz) +&
             halfTdgy + rhoenthalpyW2*&
             (vlowx*dy_betax + vlowy*dy_betay + vlowz*dy_betaz)/&
             localalp
        sz_source = t00*&
             (halfshiftdgz - localalp*dz_alp) +&
             t0x*(shiftx*dz_gxx + shifty*dz_gxy + shiftz*dz_gxz) +&
             t0y*(shiftx*dz_gxy + shifty*dz_gyy + shiftz*dz_gyz) +&
             t0z*(shiftx*dz_gxz + shifty*dz_gyz + shiftz*dz_gzz) +&
             halfTdgz + rhoenthalpyW2*&
             (vlowx*dz_betax + vlowy*dz_betay + vlowz*dz_betaz)/&
             localalp

        if ( whisky_mhd_handle.gt.0 ) then !! T^0_i d_j beta^i for B-field
             ! can be simplified a bit by using Bvxlow and B^2 directly 
             sx_source = sx_source + b2*w2*(vlowx*dx_betax + vlowy*dx_betay &
                + vlowz*dx_betaz)/localalp - udotB*(blowx*dx_betax &
                + blowy*dx_betay + blowz*dx_betaz)/localalp
             sy_source = sy_source + b2*w2*(vlowx*dy_betax + vlowy*dy_betay &
                + vlowz*dy_betaz)/localalp - udotB*(blowx*dy_betax &
                + blowy*dy_betay + blowz*dy_betaz)/localalp
             sz_source = sz_source + b2*w2*(vlowx*dz_betax + vlowy*dz_betay &
                + vlowz*dz_betaz)/localalp - udotB*(blowx*dz_betax &
                + blowy*dz_betay + blowz*dz_betaz)/localalp
        end if

        densrhs(i,j,k) = 0.d0
        sxrhs(i,j,k)  = localalp*sqrtdet*sx_source
        syrhs(i,j,k)  = localalp*sqrtdet*sy_source
        szrhs(i,j,k)  = localalp*sqrtdet*sz_source
        taurhs(i,j,k) = localalp*sqrtdet*tau_source

        if ( whisky_mhd_handle.eq.2) then

            !! Extra terms, Roland put in. Not sure where this came from.
            if (mhd_zeta.ne.0.d0) then

               if ( whisky_local_spatial_order.eq.2 ) then
                  djBnj = DIFF_X_2(Bnx) + DIFF_Y_2(Bny) + DIFF_Z_2(Bnz)
               else
                  djBnj = DIFF_X_4(Bnx) + DIFF_Y_4(Bny) + DIFF_Z_4(Bnz)
               end if

               Bnxrhs(i,j,k) = - mhd_zeta*localalp*localvx*djBnj
               Bnyrhs(i,j,k) = - mhd_zeta*localalp*localvy*djBnj
               Bnzrhs(i,j,k) = - mhd_zeta*localalp*localvz*djBnj
               taurhs(i,j,k) = taurhs(i,j,k) - mhd_zeta*localalp*djBnj* &
                  (vlowx*localBnx + vlowy*localBny + vlowz*localBnz)
               sxrhs(i,j,k) = sxrhs(i,j,k) - mhd_zeta*localalp*djBnj* &
                  ( (localgxx*localBnx+localgxy*localBny+localgxz*localBnz)/w2 &
                    + vlowx*(vlowx*localBnx+vlowy*localBny+vlowz*localBnz ) )
               syrhs(i,j,k) = syrhs(i,j,k) - mhd_zeta*localalp*djBnj* &
                  ( (localgxy*localBnx+localgyy*localBny+localgyz*localBnz)/w2 &
                    + vlowy*(vlowx*localBnx+vlowy*localBny+vlowz*localBnz ) )
               szrhs(i,j,k) = szrhs(i,j,k) - mhd_zeta*localalp*djBnj* &
                  ( (localgxz*localBnx+localgyz*localBny+localgzz*localBnz)/w2 &
                    + vlowz*(vlowx*localBnx+vlowy*localBny+vlowz*localBnz ) )
            end if

            if ( clean_divergence.ne.0 ) then

                !! g^{jk} d_i g_{kj} = d_i (g) / g
                gdxg = uxx*dx_gxx + uyy*dx_gyy + uzz*dx_gzz + 2.d0*uxy*dx_gxy + 2.d0*uxz*dx_gxz + 2.d0*uyz*dx_gyz
                gdyg = uxx*dy_gxx + uyy*dy_gyy + uzz*dy_gzz + 2.d0*uxy*dy_gxy + 2.d0*uxz*dy_gxz + 2.d0*uyz*dy_gyz
                gdzg = uxx*dz_gxx + uyy*dz_gyy + uzz*dz_gzz + 2.d0*uxy*dz_gxy + 2.d0*uxz*dz_gxz + 2.d0*uyz*dz_gyz

                !! g^{ik} d_k g_{li} --> raise index to get d_i(g^{ij})
                div_ug_lx = uxx*dx_gxx + uxy*dy_gxx + uxz*dz_gxx + &
                           uxy*dx_gxy + uyy*dy_gxy + uyz*dz_gxy + &
                           uxz*dx_gxz + uyz*dy_gxz + uzz*dz_gxz

                div_ug_ly = uxx*dx_gxy + uxy*dy_gxy + uxz*dz_gxy + &
                           uxy*dx_gyy + uyy*dy_gyy + uyz*dz_gyy + &
                           uxz*dx_gyz + uyz*dy_gyz + uzz*dz_gyz

                div_ug_lz = uxx*dx_gxz + uxy*dy_gxz + uxz*dz_gxz + &
                           uxy*dx_gyz + uyy*dy_gyz + uyz*dz_gyz + &
                           uxz*dx_gzz + uyz*dy_gzz + uzz*dz_gzz

                local_dcpsi = divclean_psi(i,j,k)

                !! We use the divergence cleaning derived via Penner's PhD thesis. This is in GRHydro's MHD
                !! Note: To get the proper Bfield evolution via Penner, mhd_zeta=>0.
                !! We put the d_i(psi) terms into the flux calculation for the B-field.
                Bnxrhs(i,j,k) = Bnxrhs(i,j,k) &
                   - sqrtdet*( localBnx*dx_betax + localBny*dy_betax + localBnz*dz_betax ) &
                   + local_dcpsi*sqrtdet*( uxx*dx_alp + uxy*dy_alp + uxz*dz_alp ) &
                   + half*local_dcpsi*localalp*( uxx*gdxg + uxy*gdyg + uxz*gdzg ) &
                   + local_dcpsi*localalp*sqrtdet*( uxx*div_ug_lx + &
                       uxy*div_ug_ly + uxz*div_ug_lz )

                Bnyrhs(i,j,k) = Bnyrhs(i,j,k) &
                   - sqrtdet*( localBnx*dx_betay + localBny*dy_betay + localBnz*dz_betay ) &
                   + local_dcpsi*sqrtdet*( uxy*dx_alp + uyy*dy_alp + uyz*dz_alp ) &
                   + half*local_dcpsi*localalp*( uxy*gdxg + uyy*gdyg + uyz*gdzg ) &
                   + local_dcpsi*localalp*sqrtdet*( uxy*div_ug_lx + &
                       uyy*div_ug_ly + uyz*div_ug_lz )

                Bnzrhs(i,j,k) = Bnzrhs(i,j,k) &
                   - sqrtdet*( localBnx*dx_betaz + localBny*dy_betaz + localBnz*dz_betaz ) &
                   + local_dcpsi*sqrtdet*( uxz*dx_alp + uyz*dy_alp + uzz*dz_alp ) &
                   + half*local_dcpsi*localalp*( uxz*gdxg + uyz*gdyg + uzz*gdzg ) &
                   + local_dcpsi*localalp*sqrtdet*( uxz*div_ug_lx + &
                       uyz*div_ug_ly + uzz*div_ug_lz )

                divclean_psirhs(i,j,k) = - 1.d0 * local_dcpsi * ( divclean_cr*localalp &
                   + dx_betax + dy_betay + dz_betaz ) &
                   + localBnx*( dx_alp - half*localalp * &
                     ( uxx*dx_gxx + uyy*dx_gyy + uzz*dx_gzz + 2.d0*uxy*dx_gxy + &
                       2.d0*uxz*dx_gxz + 2.d0*uyz*dx_gyz ) )/ sqrtdet &
                   + localBny*( dy_alp - half*localalp * &
                     ( uxx*dy_gxx + uyy*dy_gyy + uzz*dy_gzz + 2.d0*uxy*dy_gxy + &
                       2.d0*uxz*dy_gxz + 2.d0*uyz*dy_gyz ) )/ sqrtdet &
                   + localBnz*( dz_alp - half*localalp * &
                     ( uxx*dz_gxx + uyy*dy_gyy + uzz*dy_gzz + 2.d0*uxy*dy_gxy + &
                       2.d0*uxz*dz_gxz + 2.d0*uyz*dy_gyz ) )/ sqrtdet
                   
               if ( (Bnxrhs(i,j,k).gt.1.d100) .or. (Bnyrhs(i,j,k).gt.1.d100) .or. (Bnzrhs(i,j,k).gt.1.d100) ) then
                   write(warnline,'("Bad Bfield sources. Bn=(",3(g15.6,","),") and divclean_psi=",g15.6)') &
                      localBnx,localBny,localBnz,local_dcpsi
                   call CCTK_WARN(2,warnline)
               end if 
               if ( divclean_psirhs(i,j,k).gt.1.d100) then
                   write(warnline,'("Bad divclean_psi source. Bn=(",3(g15.6,","),") and divclean_psi=",g15.6)') &
                      localBnx,localBny,localBnz,local_dcpsi
                   call CCTK_WARN(2,warnline)
               end if 

                !! Below are the sources from Lehner arXiv:1001.0575
                !! Bnxrhs(i,j,k) = Bnxrhs(i,j,k) - localalp*sqrtdet*( uxx*dx_psi + uxy*dy_psi + uxz*dz_psi)
                !! Bnyrhs(i,j,k) = Bnyrhs(i,j,k) - localalp*sqrtdet*( uxy*dx_psi + uyy*dy_psi + uyz*dz_psi)
                !! Bnzrhs(i,j,k) = Bnzrhs(i,j,k) - localalp*sqrtdet*( uxz*dx_psi + uyz*dy_psi + uzz*dz_psi)
                !! divclean_psirhs(i,j,k) = - 2.d0*localalp*divclean_cr*local_dcpsi &
                !!   + (divclean_ch/sqrtdet)*(localBnx*dx_alp + localBny*dy_alp + localBnz*dz_alp  &
                !!   - 0.5d0*localalp*(localBnx*gdxg + localBny*gdyg + localBnz*gdzg) ) &
                !!   + localalp*(velxshift*dx_psi+velyshift*dy_psi+velzshift*dz_psi)

            end if !! End divergence cleaning sources

        else if (whisky_mhd_handle.eq.3) then !! Avec, algebraic gauge

            !! Just use conservative Bfield?
            !! Check factor of alpha in vel
            Axrhs(i,j,k) = localalp*sqrtdet*( velyshift*localBz - velzshift*localBy )
            Ayrhs(i,j,k) = localalp*sqrtdet*( velzshift*localBx - velxshift*localBz )
            Azrhs(i,j,k) = localalp*sqrtdet*( velxshift*localBy - velyshift*localBx )

        else if (whisky_mhd_handle.eq.4) then !! Avec, Lorenz gauge


            localAx = Ax(i,j,k)
            localAy = Ay(i,j,k)
            localAz = Az(i,j,k)
            localAphi = Aphi(i,j,k)

            upAx = uxx*localAx + uxy*localAy + uxz*localAz
            upAy = uxy*localAx + uyy*localAy + uyz*localAz
            upAz = uxz*localAx + uyz*localAy + uzz*localAz

            !! g^{jk} d_i g_{kj} = d_i (g) / g
            gdxg = uxx*dx_gxx + uyy*dx_gyy + uzz*dx_gzz + 2.d0*uxy*dx_gxy + 2.d0*uxz*dx_gxz + 2.d0*uyz*dx_gyz
            gdyg = uxx*dy_gxx + uyy*dy_gyy + uzz*dy_gzz + 2.d0*uxy*dy_gxy + 2.d0*uxz*dy_gxz + 2.d0*uyz*dy_gyz
            gdzg = uxx*dz_gxx + uyy*dz_gyy + uzz*dz_gzz + 2.d0*uxy*dz_gxy + 2.d0*uxz*dz_gxz + 2.d0*uyz*dz_gyz

            !! g^{ik} d_k g_{li} --> raise index to get d_i(g^{ij})
            div_ug_lx = uxx*dx_gxx + uxy*dy_gxx + uxz*dz_gxx + &
                       uxy*dx_gxy + uyy*dy_gxy + uyz*dz_gxy + &
                       uxz*dx_gxz + uyz*dy_gxz + uzz*dz_gxz

            div_ug_ly = uxx*dx_gxy + uxy*dy_gxy + uxz*dz_gxy + &
                       uxy*dx_gyy + uyy*dy_gyy + uyz*dz_gyy + &
                       uxz*dx_gyz + uyz*dy_gyz + uzz*dz_gyz

            div_ug_lz = uxx*dx_gxz + uxy*dy_gxz + uxz*dz_gxz + &
                       uxy*dx_gyz + uyy*dy_gyz + uyz*dz_gyz + &
                       uxz*dx_gzz + uyz*dy_gzz + uzz*dz_gzz

            if (whisky_local_spatial_order .eq. 2) then

                dx_Ax = DIFF_X_2(Ax) 
                dy_Ax = DIFF_Y_2(Ax) 
                dz_Ax = DIFF_Z_2(Ax) 
                dx_Ay = DIFF_X_2(Ay) 
                dy_Ay = DIFF_Y_2(Ay) 
                dz_Ay = DIFF_Z_2(Ay) 
                dx_Az = DIFF_X_2(Az) 
                dy_Az = DIFF_Y_2(Az) 
                dz_Az = DIFF_Z_2(Az) 

                dx_Aphi = DIFF_X_2(Aphi)
                dy_Aphi = DIFF_Y_2(Aphi)
                dz_Aphi = DIFF_Z_2(Aphi)

             else

                dx_Ax = DIFF_X_4(Ax) 
                dy_Ax = DIFF_Y_4(Ax) 
                dz_Ax = DIFF_Z_4(Ax) 
                dx_Ay = DIFF_X_4(Ay) 
                dy_Ay = DIFF_Y_4(Ay) 
                dz_Ay = DIFF_Z_4(Ay) 
                dx_Az = DIFF_X_4(Az) 
                dy_Az = DIFF_Y_4(Az) 
                dz_Az = DIFF_Z_4(Az) 

                dx_Aphi = DIFF_X_4(Aphi)
                dy_Aphi = DIFF_Y_4(Aphi)
                dz_Aphi = DIFF_Z_4(Aphi)

             end if

            !!betaDotA = shiftx*Ax + shifty*Ay + shiftz*Az 

            !! Just use conservative Bfield?
            !! Check factor of alpha in vel
            !! d_t A_i = epsilon_ijk v^j B^k - d_i(alp*Aphi - beta.A)
            !! d_t (sqrtdet*Aphi) = - d_j ( alp*sqrtdet*A^j - sqrtdet*Aphi*beta^j )

            Axrhs(i,j,k) = localalp*sqrtdet*( velyshift*localBz - velzshift*localBy ) &
                           - (1.d0/sqrtdet)*( localalp*dx_Aphi + localAphi*dx_alp ) &
                           - sqrtdet*localalp*localAphi*gdxg &
                           + shiftx*dx_Ax + shifty*dx_Ay + shiftz*dx_Az &
                           + localAx*dx_betax + localAy*dx_betay + localAz*dx_betaz
            Ayrhs(i,j,k) = localalp*sqrtdet*( velzshift*localBx - velxshift*localBz ) &
                           - (1.d0/sqrtdet)*( localalp*dy_Aphi + localAphi*dy_alp ) &
                           - sqrtdet*localalp*localAphi*gdyg &
                           + shiftx*dy_Ax + shifty*dy_Ay + shiftz*dy_Az &
                           + localAx*dy_betax + localAy*dy_betay + localAz*dy_betaz
            Azrhs(i,j,k) = localalp*sqrtdet*( velxshift*localBy - velyshift*localBx ) &
                           - (1.d0/sqrtdet)*( localalp*dz_Aphi + localAphi*dz_alp ) &
                           - sqrtdet*localalp*localAphi*gdzg &
                           + shiftx*dz_Ax + shifty*dz_Ay + shiftz*dz_Az &
                           + localAx*dz_betax + localAy*dz_betay + localAz*dz_betaz

            Aphirhs(i,j,k) = - sqrtdet*( upAx*dx_alp + upAy*dy_alp + upAz*dz_alp ) &
                 - 0.5d0*localalp*sqrtdet*( upAx*gdxg + upAy*gdyg + upAz*gdzg ) &
                 - localalp*sqrtdet*( localAx*div_ug_lx + localAy*div_ug_ly &
                                        + localAz*div_ug_lz ) &
                 - localalp*sqrtdet*( uxx*dx_Ax + uxy*dx_Ay + uxz*dx_Az &
                                     + uxy*dy_Ax + uyy*dy_Ay + uyz*dy_Az &
                                     + uxz*dz_Ax + uyz*dz_Ay + uzz*dz_Az ) &
                 + shiftx*dx_Aphi + shifty*dy_Aphi + shiftz*dz_Aphi &
                 + localAphi*( dx_betax + dy_betay + dz_betaz ) 

 
        end if !! End MHD rhs

        
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO

  deallocate(force_spatial_second_order)

end subroutine SourceTerms



