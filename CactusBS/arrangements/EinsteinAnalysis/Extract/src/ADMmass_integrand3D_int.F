
#include "cctk.h"
 
      MODULE ADMmass_integrand3D_int

c     ------------------------------------------------------------------

      INTERFACE

          SUBROUTINE  ADMmass_integrand3D(origin,Dx,Dy,Dz,x,y,z,gxx,gxy,
     &        gxz,gyy,gyz,gzz,ADMmass_int,Psi,Psi_power,conformal_state)
          IMPLICIT NONE
          CCTK_INT,INTENT(IN) ::
     &         conformal_state
          INTEGER,INTENT(IN) ::
     &         Psi_power
          CCTK_REAL,INTENT(IN) :: 
     &        Dx,Dy,Dz,origin(3)
          CCTK_REAL,DIMENSION(:),INTENT(IN) :: 
     &        x,y,z
          CCTK_REAL,DIMENSION(:,:,:),INTENT(IN) :: 
     &        gxx,gxy,gxz,gyy,gyz,gzz,Psi
          CCTK_REAL,DIMENSION(:,:,:),INTENT(OUT) :: 
     &        ADMmass_int
          END SUBROUTINE  ADMmass_integrand3D

      END INTERFACE

      END MODULE ADMmass_integrand3D_int
