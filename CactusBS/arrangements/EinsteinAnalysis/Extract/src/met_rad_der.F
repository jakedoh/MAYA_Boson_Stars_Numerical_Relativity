
#include "cctk.h"

c     ========================================================================

      SUBROUTINE met_rad_der(origin,Dx,Dy,Dz,x,y,z,g,dg)

c     ------------------------------------------------------------------------
c
c     Calculate isotropic radial derivative of unphysical metric,
c     component or conformal factor at each point on the 3D grid, 
c     using derivatives in the Cartesian directions.
c
c     ------------------------------------------------------------------------
      
      IMPLICIT NONE

c     Input variables
      CCTK_REAL,INTENT(IN) :: 
     &     Dx,Dy,Dz,origin(3)
      CCTK_REAL,DIMENSION(:),INTENT(IN) :: 
     &     x,y,z
      CCTK_REAL,DIMENSION(:,:,:),INTENT(IN) :: 
     &     g
      
c     Output variables
      CCTK_REAL,DIMENSION(:,:,:),INTENT(OUT) :: 
     &     dg
      
c     Local variables, here only
      INTEGER :: 
     &     i,j,k
      CCTK_REAL,PARAMETER ::
     &     zero = 0.0D0,
     &     half = 0.5D0
      CCTK_REAL :: 
     &     rad

c     ------------------------------------------------------------------------

      DO k = 2, SIZE(z)-1
         DO j = 2, SIZE(y)-1
            DO i = 2, SIZE(x)-1
               
               rad = SQRT((x(i)-origin(1))**2 
     &              +(y(j)-origin(2))**2
     &              +(z(k)-origin(3))**2)
               
               IF (rad.NE.0) THEN
                  
                  dg(i,j,k) = half/rad*
     &                 ((x(i)-origin(1))/Dx*(g(i+1,j,k)-g(i-1,j,k))
     &                 +(y(j)-origin(2))/Dy*(g(i,j+1,k)-g(i,j-1,k))
     &                 +(z(k)-origin(3))/Dz*(g(i,j,k+1)-g(i,j,k-1)))
                  
               ELSE
                  
                  dg(i,j,k) = zero
                  
               ENDIF

          ENDDO
        ENDDO
      ENDDO

c     This is needed when the grid is an octant, but it does not hurt
c     if it is not

      DO k = 2, size(z)-1
         DO j = 2, size(y)-1
            dg(1,j,k) = dg(2,j,k)
         ENDDO
      ENDDO
      DO k = 2, size(z)-1
         DO i = 1, size(x)-1
            dg(i,1,k) = dg(i,2,k)
         ENDDO
      ENDDO
      DO j = 1, size(y)-1
         DO i = 1, size(x)-1
            dg(i,j,1) = dg(i,j,2)
         ENDDO
      ENDDO

      END SUBROUTINE met_rad_der






