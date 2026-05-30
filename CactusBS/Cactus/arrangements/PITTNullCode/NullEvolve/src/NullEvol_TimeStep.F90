! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "gauge.h"

subroutine NullEvol_TimeStep(CCTK_ARGUMENTS)
 use NullEvol_cfl_test
  implicit none
  
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

#define SKIP_THIS_ROUTINE
#ifndef SKIP_THIS_ROUTINE

!  call CCTK_INFO("NewsTest update time")
   if ( FirstTime ) then
     FirstTime = .false.
     call CCTK_ReductionArrayHandle(reduce_handle, "minimum");
     if (reduce_handle .lt. 0 ) then
       call CCTK_WARN(0,"Could not get reduction handle")
     endif
  endif
 
  !! REPLACE ME WITH A REAL CFL ROUTINE
#ifdef HORIZON_GAUGE
  null_dt  =  -cctk_time * dt_fact
write(*,*)'TimeStep HORIZON_GAUGE', null_dt
#endif

  if (mod(it, it_cfl) == 0 ) then
    if (it .eq. it_start) then 
      cfl_dt(1) = null_dt /cfl
      cfl_dt(2) = null_dt /cfl
    endif
    call null_cfl0(null_delta(1), null_delta(2),&
       null_lsh(1), null_lsh(2), N_radial_pts, jcn, bcn, ucn, wcn, cfl_dt(1),&
       boundary_maskn, stereo_q, stereo_p, null_xb, null_rb, null_dx, cfl, null_rwt)

    call null_cfl0(null_delta(1), null_delta(2),&
       null_lsh(1), null_lsh(2), N_radial_pts, jcs, bcs, ucs, wcs, cfl_dt(2),&
       boundary_masks, stereo_q, stereo_p, null_xb, null_rb, null_dx, cfl, null_rwt)

!Communication replacement
     cfl_dt(3) = min(cfl_dt(1),cfl_dt(2))
     
     call CCTK_ReduceLocScalar(reval, cctkGH, -1, reduce_handle,&
        cfl_dt(3), gl_dt, CCTK_VARIABLE_REAL)
     cfl_dt(3) = cfl * gl_dt

!     call CCTK_ReduceLocScalar(reval, cctkGH, -1, reduce_handle,&
!        cfl_dt(1), gl_dt, CCTK_VARIABLE_REAL)
!    cfl_dt(1) = gl_dt

!     call CCTK_ReduceLocScalar(reval, cctkGH, -1, reduce_handle,&
!        cfl_dt(2), gl_dt, CCTK_VARIABLE_REAL)
!    cfl_dt(2) = gl_dt

!    cfl_dt(3) = cfl * min(cfl_dt(1), cfl_dt(2))
    tmp_dt = null_dt
    null_dt = min(null_dt, cfl_dt(3))
 
    if (null_dt .ne. tmp_dt) then
      write(message,*) "CFL ADJUST !!! ", it, cctk_time, tmp_dt, null_dt
      call CCTK_INFO(trim(message))
    endif
  endif

  cctk_delta_time = null_dt

  cctk_time = cctk_time + null_dt
  it = it + 1   ! time counter

#endif

end subroutine NullEvol_TimeStep


