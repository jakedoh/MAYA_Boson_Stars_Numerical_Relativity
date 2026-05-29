
#include "cctk.h"
 
      MODULE met_rad_der_int

c     ------------------------------------------------------------------

      INTERFACE

          SUBROUTINE met_rad_der(origin,Dx,Dy,Dz,x,y,z,g,dg)
          IMPLICIT NONE
          CCTK_REAL,INTENT(IN) :: 
     &        Dx,Dy,Dz,origin(3)
          CCTK_REAL,DIMENSION(:),INTENT(IN) :: 
     &        x,y,z
          CCTK_REAL,DIMENSION(:,:,:),INTENT(IN) :: 
     &        g
          CCTK_REAL,DIMENSION(:,:,:),INTENT(OUT) :: 
     &        dg
          END SUBROUTINE met_rad_der

      END INTERFACE

      END MODULE met_rad_der_int
