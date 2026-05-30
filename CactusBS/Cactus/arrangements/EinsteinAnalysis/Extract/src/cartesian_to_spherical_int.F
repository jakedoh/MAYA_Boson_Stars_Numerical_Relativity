
#include "cctk.h"

      MODULE cartesian_to_spherical_int

c     ------------------------------------------------------------------

      INTERFACE

      SUBROUTINE cartesian_to_spherical(theta,phi,rad,gxxs,gxys,gxzs,
     &     gyys,gyzs,gzzs,grrs,grts,grps,gtts,gtps,gpps,dgxxs,dgxys,
     &     dgxzs,dgyys,dgyzs,dgzzs,dgtts,dgtps,dgpps)
      IMPLICIT NONE
      CCTK_REAL  :: 
     &     rad,theta(:),phi(:)
      CCTK_REAL, DIMENSION (:,:) :: 
     &     gxxs,gxys,gxzs,gyys,gyzs,gzzs,dgxxs,dgxys,dgxzs,dgyys,
     &     dgyzs,dgzzs,grrs,grts,grps,gtts,gtps,gpps,dgtts,dgtps,
     &     dgpps
      END SUBROUTINE

      END INTERFACE

      END MODULE cartesian_to_spherical_int
