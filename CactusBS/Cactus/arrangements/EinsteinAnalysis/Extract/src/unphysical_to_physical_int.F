
#include "cctk.h"

      MODULE unphysical_to_physical_int

c     ------------------------------------------------------------------

      INTERFACE

      SUBROUTINE unphysical_to_physical(grr,grt,grp,gtt,gtp,gpp,dgtt,
     &     dgtp,dgpp,Psis,dPsis,Psi_power)
      IMPLICIT NONE
      INTEGER :: 
     &     Psi_power
      CCTK_REAL, DIMENSION (:,:) :: 
     &     grr,grt,grp,gtt,gtp,gpp,dgtt,dgtp,dgpp,Psis,dPsis
      END SUBROUTINE
      
      END INTERFACE
      
      END MODULE unphysical_to_physical_int
