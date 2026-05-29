
#include "cctk.h"
 
      MODULE spin_integrand3D_int

c     ------------------------------------------------------------------

      INTERFACE

          SUBROUTINE  spin_integrand3D(origin,x,y,z,gxx,
     &        gxy,gxz,gyy,gyz,gzz,hxx,hxy,hxz,hyy,hyz,hzz,spin_int,
     &        Psi,Psi_power,conformal_state)

          IMPLICIT NONE
          CCTK_INT,INTENT(IN) ::
     &         conformal_state
          INTEGER,INTENT(IN) ::
     &         Psi_power
          CCTK_REAL,INTENT(IN) :: 
     &        origin(3)
          CCTK_REAL,DIMENSION(:),INTENT(IN) :: 
     &        x,y,z
          CCTK_REAL,DIMENSION(:,:,:),INTENT(IN) :: 
     &         gxx,gxy,gxz,gyy,gyz,gzz,Psi,
     &         hxx,hxy,hxz,hyy,hyz,hzz
          CCTK_REAL,DIMENSION(:,:,:),INTENT(OUT) :: 
     &        spin_int
          END SUBROUTINE  spin_integrand3D

      END INTERFACE

      END MODULE spin_integrand3D_int
