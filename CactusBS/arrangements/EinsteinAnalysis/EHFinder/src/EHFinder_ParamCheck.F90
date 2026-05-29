!  @@
!  @file      EHFinder_ParamCheck.F90
!  @date      Tue May 21 21:26:45 2002
!  @author    Peter Diener
!  @desc 
!  Parameter checking stuff for EHFinder
!  @enddesc
!  @version $Header$
! @@

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine EHFinder_ParamCheck(CCTK_ARGUMENTS)
 
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  character(len=128) :: param_warn

  ! Check if the metric is one of the known types: physical or
  ! static conformal
  if ( ( .not. CCTK_EQUALS(metric_type,'physical') ) .and. &
       ( .not. CCTK_EQUALS(metric_type,'static conformal') ) ) then
    param_warn = 'Unknown ADMBase::metric_type - known types are "physical" '
    param_warn = trim(param_warn)//'and "static conformal"'
    call CCTK_PARAMWARN ( trim(param_warn) )
  end if

end subroutine EHFinder_ParamCheck
