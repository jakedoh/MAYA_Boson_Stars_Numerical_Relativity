#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

c     ========================================================================

      SUBROUTINE ADMmass_integrand3D(origin,Dx,Dy,Dz,x,y,z,gxx,gxy,gxz,
     &     gyy,gyz,gzz,ADMmass_int,Psi,Psi_power,conformal_state)

c     ------------------------------------------------------------------------
c
c     Estimates the ADM mass at a given radius using Equation (7) from
c     O Murchadha and York, Phys Rev D, 10, 1974 page 2345
c
c     ------------------------------------------------------------------------
      
      IMPLICIT NONE

c     Input variables
      INTEGER,INTENT(IN) ::
     &     Psi_power,conformal_state
      CCTK_REAL,INTENT(IN) :: 
     &     Dx,Dy,Dz,origin(3)
      CCTK_REAL,DIMENSION(:),INTENT(IN) :: 
     &     x,y,z
      CCTK_REAL,DIMENSION(:,:,:),INTENT(IN) :: 
     &     gxx,gxy,gxz,gyy,gyz,gzz,Psi
      
c     Output variables
      CCTK_REAL,DIMENSION(:,:,:),INTENT(OUT) :: 
     &     ADMmass_int
      
c     Local variables, here only
      INTEGER :: 
     &     i,j,k,ip
      CCTK_REAL,PARAMETER ::
     &     half = 0.5D0
      CCTK_REAL :: 
     &     rad,ux,uy,uz,det,dxx,dxy,dxz,dyy,dyz,dzz,uxx,uxy,uxz,uyy,
     &     uyz,uzz,term1,term2,term3,dxx_y,dxx_z,dxy_x,dxy_y,dxy_z,
     &     dyy_x,dxz_x,dxz_y,dxz_z,dyz_x,dzz_x,dyy_z,dyz_y,dyz_z,dzz_y,
     &     Pi,idet,p,pip,pim,pjp,pjm,pkp,pkm

      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

c     ------------------------------------------------------------------------

      Pi = ACOS(-1D0)

c     Because other codes evolve Psi**4
      SELECT CASE (Psi_power)

      CASE (1)
         ip = 4

      CASE (4)
         ip = 1

      CASE DEFAULT
         WRITE(*,*) "This value of Psi_power is not supported"

      END SELECT   
      

      DO k = 2, SIZE(z)-1
         DO j = 2, SIZE(y)-1
            DO i = 2, SIZE(x)-1
               
               rad = SQRT((x(i)-origin(1))**2 
     &              +(y(j)-origin(2))**2
     &              +(z(k)-origin(3))**2)
               
               IF (rad.NE.0) THEN

                  ux = (x(i)-origin(1))/rad
                  uy = (y(j)-origin(2))/rad
                  uz = (z(k)-origin(3))/rad

c                 Abbreviations for metric coefficients
c                 -------------------------------------

                  if (conformal_state > 0) then
  
                    p = psi(i,j,k)**ip

                    dxx = p*gxx(i,j,k); dxy = p*gxy(i,j,k)
                    dxz = p*gxz(i,j,k); dyy = p*gyy(i,j,k)
                    dyz = p*gyz(i,j,k); dzz = p*gzz(i,j,k)

                 else 

                    p = 1.0d0

                    dxx = gxx(i,j,k); dxy = gxy(i,j,k)
                    dxz = gxz(i,j,k); dyy = gyy(i,j,k)
                    dyz = gyz(i,j,k); dzz = gzz(i,j,k)

                 end if

c                 Determinant of 3-metric
c                 -----------------------
                  det = (dxx*dyy*dzz + 2.0D0*dxy*dxz*dyz
     &                 - (dxx*dyz**2 + dyy*dxz**2 + dzz*dxy**2))

                  idet = 1.0/det

c                 Inverse 3-metric
c                 ----------------
                  uxx = idet*(dyy*dzz - dyz**2)
                  uyy = idet*(dxx*dzz - dxz**2)
                  uzz = idet*(dxx*dyy - dxy**2)
                  uxy = idet*(dxz*dyz - dzz*dxy)
                  uxz = idet*(dxy*dyz - dyy*dxz)
                  uyz = idet*(dxy*dxz - dxx*dyz)

c                 Derivatives of the 3-metric
c                 ---------------------------

                  IF (conformal_state > 0) THEN

                    pip = psi(i+1,j,k)**ip
                    pim = psi(i-1,j,k)**ip
                    pjp = psi(i,j+1,k)**ip
                    pjm = psi(i,j-1,k)**ip
                    pkp = psi(i,j,k+1)**ip
                    pkm = psi(i,j,k-1)**ip

                 else

                    pip = 1.0d0
                    pim = 1.0d0
                    pjp = 1.0d0
                    pjm = 1.0d0
                    pkp = 1.0d0
                    pkm = 1.0d0

                 end if

                  dxx_y = half/Dy*(pjp*gxx(i,j+1,k)-pjm*gxx(i,j-1,k))
                  dxx_z = half/Dz*(pkp*gxx(i,j,k+1)-pkm*gxx(i,j,k-1))
                  dxy_x = half/Dx*(pip*gxy(i+1,j,k)-pim*gxy(i-1,j,k))
                  dxy_y = half/Dy*(pjp*gxy(i,j+1,k)-pjm*gxy(i,j-1,k))
                  dxy_z = half/Dz*(pkp*gxy(i,j,k+1)-pkm*gxy(i,j,k-1))
                  dyy_x = half/Dx*(pip*gyy(i+1,j,k)-pim*gyy(i-1,j,k))
                  dyy_z = half/Dz*(pkp*gyy(i,j,k+1)-pkm*gyy(i,j,k-1))
                  dxz_x = half/Dx*(pip*gxz(i+1,j,k)-pim*gxz(i-1,j,k))
                  dxz_y = half/Dy*(pjp*gxz(i,j+1,k)-pjm*gxz(i,j-1,k))
                  dxz_z = half/Dz*(pkp*gxz(i,j,k+1)-pkm*gxz(i,j,k-1))
                  dyz_x = half/Dx*(pip*gyz(i+1,j,k)-pim*gyz(i-1,j,k))
                  dyz_y = half/Dy*(pjp*gyz(i,j+1,k)-pjm*gyz(i,j-1,k))
                  dyz_z = half/Dz*(pkp*gyz(i,j,k+1)-pkm*gyz(i,j,k-1))
                  dzz_x = half/Dx*(pip*gzz(i+1,j,k)-pim*gzz(i-1,j,k))
                  dzz_y = half/Dy*(pjp*gzz(i,j+1,k)-pjm*gzz(i,j-1,k))

                  term1 = uxy*(dxx_y-dxy_x)+uxz*(dxx_z-dxz_x)
     &                 +uyy*(dxy_y-dyy_x)+uyz*(dxy_z-dyz_x)
     &                 +uyz*(dxz_y-dyz_x)+uzz*(dxz_z-dzz_x)

                  term2 = uyz*(dyy_z-dyz_y)+uxy*(dyy_x-dxy_y)
     &                 +uzz*(dyz_z-dzz_y)+uxz*(dyz_x-dxz_y)
     &                 +uxz*(dxy_z-dxz_y)+uxx*(dxy_x-dxx_y)

                  term3 = uxz*(dzz_x-dxz_z)+uyz*(dzz_y-dyz_z)
     &                 +uxx*(dxz_x-dxx_z)+uxy*(dxz_y-dxy_z)
     &                 +uxy*(dyz_x-dxy_z)+uyy*(dyz_y-dyy_z)

                  ADMmass_int(i,j,k) = 1.0D0/16.0D0/Pi*(ux*term1+
     &                 uy*term2+uz*term3)*SQRT(det)*rad**2

               ELSE
                  
                  ADMmass_int(i,j,k) = 0.0D0
                  
               ENDIF

          ENDDO
        ENDDO
      ENDDO


c     This is needed when the grid is an octant, but it does not hurt
c     if it is not

      DO k = 2, size(z)-1
         DO j = 2, size(y)-1
            ADMmass_int(1,j,k) = ADMmass_int(2,j,k)
         ENDDO
      ENDDO
      DO k = 2, size(z)-1
         DO i = 1, size(x)-1
            ADMmass_int(i,1,k) = ADMmass_int(i,2,k)
         ENDDO
      ENDDO
      DO j = 1, size(y)-1
         DO i = 1, size(x)-1
            ADMmass_int(i,j,1) = ADMmass_int(i,j,2)
         ENDDO
      ENDDO


      END SUBROUTINE ADMmass_integrand3D
