
#include "cctk.h"

c     ========================================================================

      SUBROUTINE spin_integrand3D(origin,x,y,z,gxx,gxy,gxz,
     &     gyy,gyz,gzz,hxx,hxy,hxz,hyy,hyz,hzz,
     &     spin_int,Psi,Psi_power,conformal_state)

c     ------------------------------------------------------------------------
c
c     Estimates the spin at a given radius using Equation (12) from
c     <Time-asymmetric initial data for black holes and black-hole 
c     collisions>, Bowen and York, Phys. Rev. D, 21 2047(1980)
c
c     -----------------------------------------------------------------------
      
      IMPLICIT NONE

c     Input variables
      INTEGER,INTENT(IN) ::
     &     Psi_power,conformal_state
      CCTK_REAL,INTENT(IN) :: 
     &     origin(3)
      CCTK_REAL,DIMENSION(:),INTENT(IN) :: 
     &     x,y,z
      CCTK_REAL,DIMENSION(:,:,:),INTENT(IN) :: 
     &     gxx,gxy,gxz,gyy,gyz,gzz,Psi,
     &     hxx,hxy,hxz,hyy,hyz,hzz

      
c     Output variables
      CCTK_REAL,DIMENSION(:,:,:),INTENT(OUT) :: 
     &     spin_int
      
c     Local variables, here only
      INTEGER :: 
     &     i,j,k,ip,count

      CCTK_REAL,PARAMETER ::
     &     half = 0.5D0
      CCTK_REAL :: 
     &     rad,ux,uy,uz,det,dxx,dxy,dxz,dyy,dyz,dzz,uxx,uxy,uxz,uyy,
     &     uyz,uzz,
     &     Pi,idet,p,
     &     term1,term2

      data count / 1 /
      save count
c     ------------------------------------------------------------------------

      Pi = ACOS(-1D0)

      SELECT CASE (Psi_power)

      CASE (1)
         ip = 4

      CASE (4)
         ip = 1

      CASE DEFAULT
         CALL CCTK_ERROR("This value of Psi_power is not supported")

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


c                 Integrands
c                 ----------

                  if (count.eq.1) then
              term1 = uxz*(ux*hxx(i,j,k)+uy*hxy(i,j,k)+uz*hxz(i,j,k))
     &                +uyz*(ux*hxy(i,j,k)+uy*hyy(i,j,k))
     &                +uzz*(ux*hxz(i,j,k)+uy*hyz(i,j,k)+uz*hzz(i,j,k))

              term2 = uxy*(ux*hxx(i,j,k)+uy*hxy(i,j,k)+uz*hxz(i,j,k))
     &                +uyy*(ux*hxy(i,j,k)+uy*hyy(i,j,k)+uz*hyz(i,j,k))
     &                +uyz*(ux*hxz(i,j,k)+uz*hzz(i,j,k))

                     spin_int(i,j,k) = 1.0D0/8.0D0/Pi*
     &                    (uy*term1-uz*term2)*rad**3

                  else if (count.eq.2) then
              term1 = uxx*(ux*hxx(i,j,k)+uy*hxy(i,j,k)+uz*hxz(i,j,k))
     &                +uxy*(ux*hxy(i,j,k)+uy*hyy(i,j,k)+uz*hyz(i,j,k))
     &                +uxz*(uy*hyz(i,j,k)+uz*hzz(i,j,k))

              term2 = uxz*(ux*hxx(i,j,k)+uy*hxy(i,j,k))
     &                +uyz*(ux*hxy(i,j,k)+uy*hyy(i,j,k)+uz*hyz(i,j,k))
     &                +uzz*(ux*hxz(i,j,k)+uy*hyz(i,j,k)+uz*hzz(i,j,k))

                     spin_int(i,j,k) = 1.0D0/8.0D0/Pi*
     &                    (uz*term1-ux*term2)*rad**3

                  else if (count.eq.3) then
              term1 = uxy*(ux*hxx(i,j,k)+uz*hxz(i,j,k))
     &                +uyy*(ux*hxy(i,j,k)+uy*hyy(i,j,k)+uz*hyz(i,j,k))
     &                +uyz*(ux*hxz(i,j,k)+uy*hyz(i,j,k)+uz*hzz(i,j,k))

              term2 = uxx*(ux*hxx(i,j,k)+uy*hxy(i,j,k)+uz*hxz(i,j,k))
     &                +uxy*(uy*hyy(i,j,k)+uz*hyz(i,j,k))
     &                +uxz*(ux*hxz(i,j,k)+uy*hyz(i,j,k)+uz*hzz(i,j,k))

                     spin_int(i,j,k) = 1.0D0/8.0D0/Pi*
     &                    (ux*term1-uy*term2)*rad**3

                   end if


               ELSE
                  
                  spin_int(i,j,k) = 0.0D0
                  
               ENDIF
               
          ENDDO
        ENDDO
      ENDDO


c     This is needed when the grid is an octant, but it does not hurt
c     if it is not

      DO k = 2, size(z)-1
         DO j = 2, size(y)-1
            spin_int(1,j,k) = spin_int(2,j,k)
         ENDDO
      ENDDO
      DO k = 2, size(z)-1
         DO i = 1, size(x)-1
            spin_int(i,1,k) = spin_int(i,2,k)
         ENDDO
      ENDDO
      DO j = 1, size(y)-1
         DO i = 1, size(x)-1
            spin_int(i,j,1) = spin_int(i,j,2)
         ENDDO
      ENDDO


c     setting counter
c     ---------------

      if (count.eq.3) then
         count = 1
      else
         count = count + 1
      end if


c     end
c     ---

      END SUBROUTINE spin_integrand3D
