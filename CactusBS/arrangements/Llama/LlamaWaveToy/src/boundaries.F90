#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"



subroutine LWT_boundaries (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS  
  
  character :: fbound*1000
  integer   :: fboundlen
  integer   :: ierr
  
   call CCTK_FortranString (fboundlen, bound, fbound)
   if (fboundlen > len(fbound)) call CCTK_WARN (0, "internal error")
   
    ierr = Boundary_SelectGroupForBC &
         (cctkGH, CCTK_ALL_FACES, +1, -1, "LlamaWaveToy::scalar", fbound)
    if (ierr/=0) call CCTK_WARN (0, "internal error")
 
    ierr = Boundary_SelectGroupForBC &
        (cctkGH, CCTK_ALL_FACES, +1, -1, "LlamaWaveToy::density", fbound)
   if (ierr/=0) call CCTK_WARN (0, "internal error")
   
end subroutine LWT_boundaries
