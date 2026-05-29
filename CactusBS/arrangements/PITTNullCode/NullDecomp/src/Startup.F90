/*
  vim: syntax=fortran
*/
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine NullDecomp_Startup(CCTK_ARGUMENTS)
    use NullDecomp_Vars
    use NullDecomp_SpinDecomp

    implicit none

    TARGET:: area, kern

    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_ARGUMENTS


    DECLARE_CCTK_PARAMETERS

    MyProc = CCTK_MyProc( cctkGH )
    NumProcs = CCTK_nProcs( cctkGH )

    lmax = l_max

    call CCTK_ReductionArrayHandle(min_handle, "minimum")
    call CCTK_ReductionArrayHandle(max_handle, "maximum")
    call CCTK_ReductionArrayHandle(sum_handle, "sum")

    if (min_handle .lt. 0 .OR. max_handle .lt. 0.OR.sum_handle.lt.0) then
      call CCTK_WARN(0, "Could not get reduction hadles")
    endif

  
    kern(:,:,1) = 4.0d0 / (1.0d0 + zeta*conjg(zeta))**2
    kern(:,:,2) = 4.0d0 / (1.0d0 + zeta*conjg(zeta))**2

    if (null_lbnd(1) .eq. 0) then
      lq = 1
    else
      lq = 1+N_ang_ghost_pts
    endif

    if (null_lbnd(2) .eq. 0) then
      lp = 1
    else
      lp = 1+N_ang_ghost_pts
    endif

    if (null_ubnd(1) .eq. null_gsh(1) - 1) then
      uq = null_lsh(1)
    else
      uq = null_lsh(1) - N_ang_ghost_pts + 1
    endif

    if (null_ubnd(2) .eq. null_gsh(2) - 1) then
      up = null_lsh(2)
    else
      up = null_lsh(2) - N_ang_ghost_pts + 1
    endif

    call wt_intarea(null_lsh(1), null_lsh(2), dble(zeta(:,1)), dimag(zeta(1,:)), area)

    area_p => area
    kern_p => kern

contains

subroutine wt_intarea(nq, np, q, p, area)

  implicit none
  CCTK_INT,          intent (in)  :: nq, np
  CCTK_REAL, intent (out) :: area(nq,np)
  CCTK_REAL, dimension(nq), intent (in) :: q
  CCTK_REAL, dimension(np), intent (in) :: p

  CCTK_INT ::  i, j
  CCTK_REAL :: dq, dp

  CCTK_REAL :: qp, qm, pp, pm
  CCTK_REAL :: rpp, rpm, rmp, rmm
  CCTK_REAL :: topd, bottomd, leftd, rightd



  dq = q(2) - q(1)   !!! Assume uniform grid
  dp = p(2) - p(1)   !!! Assume uniform grid

  area = 0.

  do j = 1, np-1
    do i = 1, nq-1
      qp = q(i+1)
      qm = q(i)
      pp = p(j+1)
      pm = p(j)

      !r** < 0  => corner point within unit circle

    
      rpp = (qp*qp + pp*pp) - 1.0d0
      rpm = (qp*qp + pm*pm) - 1.0d0
      rmp = (qm*qm + pp*pp) - 1.0d0
      rmm = (qm*qm + pm*pm) - 1.0d0

      ! we'll fudge if the corners lie on the unit circle
      if (abs(rpp) .lt. 1.0d-16) rpp = -1.0d-16
      if (abs(rpm) .lt. 1.0d-16) rpm = -1.0d-16
      if (abs(rmp) .lt. 1.0d-16) rmp = -1.0d-16
      if (abs(rmm) .lt. 1.0d-16) rmm = -1.0d-16

      if (rpp < 0.0 .AND. rpm < 0.0 .AND. rmp < 0.0 .AND. rmm < 0.0) then
        ! if the box is entirely within the unit circle then it's easy
        area(i,j) = dq*dp
      else if(rpp > 0.0 .AND. rpm > 0.0 .AND. rmp > 0.0 .AND. rmm > 0.0) then
         ! we don't integrate outside of the unit circle
         ! (that is covered by the other patch)
        area(i,j) = 0
      else
        ! now for the nasty case of the box only partially within 
        ! the unit circle

        topd = 0
        bottomd=0
        leftd = 0
        rightd = 0
 
        if (rmp*rpp < 0.0) then ! circle crossed top line
          !topd = -rmp/(rpp - rmp)
          ! the 1.0d-30 term ensures that topd is not zero. otherwise
          ! the logic below would get the wrong answer
          topd = abs(abs(qm) - dsqrt(1.0-pp*pp))/dq + 1.0d-30
        endif
        if (rmm*rpm < 0.0) then ! circle crossed bottom line
          !bottomd = -rmm/(rpm-rmm)
          bottomd = abs(abs(qm) - dsqrt(1.0-pm*pm))/dq + 1.0d-30
        endif
        if (rmm*rmp < 0.0) then ! circle crossed left line
          !leftd = -rmm/(rmp - rmm)
          leftd = abs(abs(pm) - dsqrt(1.0-qm*qm))/dp + 1.0d-30
        endif
        if (rpm*rpp < 0.0) then ! circle crossed right line
          !rightd = -rpm/ (rpp - rpm)
          rightd = abs(abs(pm) - dsqrt(1.0-qp*qp))/dp + 1.0d-30
        endif

        if (bottomd .ne. 0 .AND. topd .ne. 0) then
          if (rmp < 0.0) then
            area(i,j) = dp*dq*topd -.5 * dp *dq* (topd - bottomd)
          else 
            area(i,j) = dp*dq*(1.-topd) +.5 * dp * dq* (topd - bottomd)
          endif
        else if(leftd .ne. 0 .AND. rightd .ne. 0) then
          if (rmm < 0.0) then
            area(i,j) = dq*dp*rightd -.5 * dq *dp* (rightd - leftd)
          else 
            area(i,j) = dq*dp*(1.0-rightd) +.5 * dq * dp* (rightd - leftd)
          endif
        else if(topd .ne. 0 .AND. leftd .ne. 0) then
          if (rmp < 0.0) then
            area(i,j) = .5*topd*dq*dp*(1.0 - leftd)
          else
            area(i,j) = dq*dp - .5*topd*dq*dp*(1.0 - leftd)
          endif
        else if(topd .ne. 0 .AND. rightd .ne. 0) then
          if (rpp < 0.0) then
            area(i,j) = .5*dq*(1.0-topd)*dp*(1.0 - rightd)
          else
            area(i,j) = dq*dp - .5*dq*(1.0-topd)*dp*(1.0 - rightd)
          endif
        else if(bottomd .ne. 0 .AND. leftd .ne. 0) then
          if (rmm < 0.0 ) then
            area(i,j) = .5*bottomd*leftd *dq*dp
          else
            area(i,j) = dq*dp - .5*bottomd*leftd*dq*dp
          endif
        else if (bottomd .ne. 0 .AND. rightd .ne. 0) then! bottom and right
          if (rpm < 0.0) then
            area(i,j) = .5*dq*(1.0 - bottomd)*dp*rightd
          else
            area(i,j) = dq*dp - .5*dq*(1.0 - bottomd)*dp*rightd
          endif
         else
            call CCTK_WARN(0,"something went horribly wrong")
        endif
      endif
    end do 
  end do 

end subroutine wt_intarea


end subroutine NullDecomp_Startup


