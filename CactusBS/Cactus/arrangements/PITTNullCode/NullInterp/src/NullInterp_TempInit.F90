! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

subroutine NullInterp_TempInit(CCTK_ARGUMENTS)
   use cctk
   implicit none
   DECLARE_CCTK_PARAMETERS
   DECLARE_CCTK_ARGUMENTS
   DECLARE_CCTK_FUNCTIONS

   integer levels_tmp_cgf, levels_tmp_rgf, levels_tmp_cgf3

! write (*,*) "SETTING TEMPORARIES TO ZERO ...."

! tmp_cgfn = 0; tmp_cgfs = 0; tmp_rgfn = 0; tmp_rgfs = 0; tmp_cgfn1 = 0; tmp_cgfs1 = 0; tmp_cgfn2 = 0; tmp_cgfs2 = 0; tmp_cgfn3 = 0; tmp_cgfs3 = 0
!return

  call CCTK_ActiveTimeLevels(levels_tmp_cgf, cctkGH, "NullInterp::tmp_cgf")
  call CCTK_ActiveTimeLevels(levels_tmp_rgf, cctkGH, "NullInterp::tmp_rgf")
  call CCTK_ActiveTimeLevels(levels_tmp_cgf3, cctkGH, "NullInterp::tmp_cgf3")
  
  if (levels_tmp_cgf.ge.1) then
  tmp_cgfn =1.e+0
  tmp_cgfs =1.e+0
  end if
 
  if (levels_tmp_rgf.ge.1) then
  tmp_rgfn=1.e+0
  tmp_rgfs=1.e+0
  end if
  
  if (levels_tmp_cgf3.ge.1) then
    tmp_cgfn1=1.e+0
    tmp_cgfs1=1.e+0
    tmp_cgfn2=1.e+0
    tmp_cgfs2=1.e+0
    tmp_cgfn3=1.e+0
    tmp_cgfs3=1.e+0
  end if


end subroutine NullInterp_TempInit
