! vim: syntax=fortran
!ifndef INTERFACE__
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#define LY 76
#define EPS 1.0d-8
#undef DEBUG__
module NullEvol_cfl_test
 contains 
! Revised CFL. To remove all possible complications, let s
! think about it this way. Being at the point p = (n+1,j,k,i)
! we want to know if dt was too large. So we  find the null rays 
! connecting p to those in the previous level that were
! involved in the calculation of the fields at p. The ones we care
! about are (n,j,k,i+1) and (n,j+A,k+B,i-1) where A,B=-1,1
! Certainly we can t evaluate everything cause we do not have the
! metric complete at n+1. We ll proceed for now doing a  first order
! extrapolation, ie setting metric at n+1 = metric at n, and then
! we ll take a smaller value for dt to be safe....
! also note that du is  negative  in the calculation, we ll absorb
! the sign in the metric quantities to deal with a positive dt
!endif
  subroutine null_cfl0(ddq, ddp, nq, np, nx,&
        jos, bos, uos, wos, dt_old, maskb, qs, ps, x, rb,&
        dx, cfl, rwt)


    implicit none

    CCTK_REAL, intent(in) :: ddq, ddp
    CCTK_INT, intent(in) :: nq, np, nx
    CCTK_COMPLEX,   dimension (nq,np,nx), intent(in)  :: jos, uos
    CCTK_REAL, dimension (nq,np,nx), intent(in)  :: bos, wos
    CCTK_REAL,                       intent(inout) :: dt_old
    CCTK_INT, dimension(nq,np), intent(in) :: maskb
    CCTK_REAL, dimension(nq, np), intent(in) :: qs, ps
    CCTK_REAL, dimension(nx), intent(in) :: x, rb
    CCTK_REAL, intent(in) :: dx, cfl, rwt

    CCTK_REAL, dimension (9)     :: dtry1, dtry2, dis, dq, dp

    CCTK_REAL :: dfor, dforU, dfor2
    CCTK_REAL :: guu, gup, guq, gur, gpp, gqq, gpq
    CCTK_COMPLEX   :: JJ, Jb, U, Ub

    CCTK_REAL :: r, P0, dr1, dr, KK, B, W, q, p, dmin1, dmin2, &
         dutime
    CCTK_INT:: i,j,k,c

#ifdef DEBUG__
character(len=500) :: message
#endif

    dmin2 = 100.

    !this condition is to have a  floor  for how small
    !we ll be willing to go. I m setting it to 1e-10 but
    ! if we get to even 1e-5 we re in trouble if we want
    ! to evolve for times of order 1....

    if(dt_old.lt.1e-10) then
       call CCTK_WARN(0, "WARNING!!!! TOO SMALL DT!!!!")
    end if

    !now set dutime to something big, it will be redefined
    dutime = 100. * dt_old

    !loop over all angular points involved in the algorithm

    dq = (/ ddq, ddq,  ddq, 0.0d0, 0.0d0,  0.0d0, -ddq, -ddq, -ddq/)
    dp = (/ ddp, 0.0d0, -ddp, ddp, 0.0d0, -ddp,  ddp,  0.0d0, -ddp/)

    !set all these equal to dt, ie the previous one obtained in this routine.
    !however we need to divide by cfl as when we leave this routine the main
    !code multiplies the dt obtained here by a  safe factor  cfl. Since we ll
    !require never to go  larger  than the dt previously used, if we don t 
    ! recover  the old dt obtained here by cfl, the net effect will be to
    ! keep multimpliying by "cfl factors" everytime we leave and we ll have
    ! dt converge to zero.

    dforU = dt_old/cfl
    dfor2 = dt_old/cfl
    dfor  = dt_old/cfl
    dmin1 = dt_old/cfl

    !evaluate the null  distances  over all points
    !but masking the first 4 and the last one. just
    !in case, soon we ll use the  real mask  from the
    !code to go over valid points only.

    do k = 2, nx - 1

       r = rb(k)
       dr1 = dx * rwt * 1. / (1.-x(k) )**2
       dr  = -dx * rwt * 1. / (1.-x(k) )**2 !note sign

       do j = 2, np-1
          do i = 2, nq-1

             !use the mask here to evaluate only valid points
             if ( k .gt. maskb(i,j) ) then
               JJ = jos(i,j,k)
                Jb = conjg(JJ)
                B = bos(i,j,k)
                U = 0.5 * ( uos(i,j,k) + uos(i,j,k-1) )
                Ub = conjg(U)
                KK = sqrt( 1. + JJ*Jb )
                W = wos(i,j,k)

                q = qs(i,j)
                p = ps(i,j)
                P0 = 1. / (1. + q**2 + p**2)

                !define metric components
                guu = -exp(2.*B) * (1. + r * W )+r**2 &
                     * ( KK * U * Ub + dble(JJ*Ub**2) )
                gur = exp(2.*B) !note sign changed!
                guq = r**2 * 2.* P0 * dble( JJ*Ub + KK * U ) !sign changed
                gup = r**2 * 2.* P0 * dimag((JJ*Ub+KK*U) ) !sign changed
                gqq = (r*P0)**2 * 4. * ( KK + dble( JJ) )
                gpp = (r*P0)**2 * 4. * ( KK - dble( JJ) )
                gpq = (r*P0)**2 * 4. * dimag(JJ)

                !evaluate  distance  of point in the radial direction. ignoring
                ! U for the moment
                dfor2 = 2. * dr1 * gur  / abs( exp(2.*B) * (1. + r * W) )

                !set to 10^-8 to avoid taking insanely small steps1
                if (dfor2 .lt. EPS) then
                   dfor2 = 100.
                end if
#ifdef DEBUG__
if (i.eq.76 .and. j.eq. LY .and.k.eq.3)then
  write(message,*) "dfor2 ", dfor2
  call CCTK_WARN(1, trim(message))
endif
#endif

                !calculate  distance  i with the full guu

                if (abs(guu) .ge. EPS) then
                   !radial driection only 
                   dfor = 2. * dr1 * gur  / (-guu)

                   if(dfor .le. 0.) then
                      dfor = 100.
                   end if
#ifdef DEBUG__
if (i.eq.76 .and. j.eq. LY .and.k.eq.3)then
  write(message,*) "dfor ", dfor
  call CCTK_WARN(1, trim(message))
endif
#endif

                   !with all directions.
                 !  dr = -dr 

                   dis(:) =  ( gur * dr + guq *dq(:) + gup*dp(:) )**2 &
                        - guu * ( 2.* gpq * dp(:) * dq(:) + gpp * dp(:)**2 &
                        + gqq * dq(:)**2 )

#ifdef DEBUG__
if (i.eq.76 .and. j.eq. LY .and.k.eq.3)then
  write(message,*) "dis ", dis
  call CCTK_WARN(1, trim(message))
  write(message,*) "dq ", dq
  call CCTK_WARN(1, trim(message))
  write(message,*) "dp ", dp
  call CCTK_WARN(1, trim(message))
  write(message,*) "dr ", dr
  call CCTK_WARN(1, trim(message))
endif
#endif
                 !check that we take square root of a positive quantity
                  do c = 1,9
                      if (dis(c) .ge. 0.0d0) then
                         dtry1(c) = (- ( gur * dr + guq * dq(c) + gup * dp(c) ) &
                              - dsqrt(dis(c) )) / guu

                         dtry2(c) = (- ( gur * dr + guq * dq(c) + gup * dp(c) ) &
                              + dsqrt(dis(c) ))/ guu
                      else
                         dtry1(c) = 100.
                         dtry2(c) = 100.
                      end if

                      if (dtry1(c).lt.-EPS) then
                         dtry1(c) = -dtry1(c)
                      else
                         dtry1(c) = 100.
                      end if

                      if (dtry2(c).lt.-EPS) then
                         dtry2(c) = -dtry2(c)
                      else
                         dtry2(c) = 100.
                      end if
                  end do
#ifdef DEBUG__
if (i.eq.76 .and. j.eq. LY .and.k.eq.3)then
  write(message,*) "dtry1 ", dtry1
  call CCTK_WARN(1, trim(message))
  write(message,*) "dtry2 ", dtry2
  call CCTK_WARN(1, trim(message))
endif
#endif
                   dmin1 = min( minval(dtry1), minval(dtry2))

                else  !(i.e. guu = 0.0)
                  ! dr = -dr 
                  do c = 1, 9
                      if (abs(gur * dr + guq *dq(c) + gup*dp(c)).gt.EPS) then
                      dtry1(c) = (2.* gpq * dp(c) * dq(c) + gpp * dp(c)**2  &
                              + gqq * dq(c)**2 ) / &
                              ( gur * dr + guq *dq(c) + gup*dp(c) ) * 0.5
                      else
                      dtry1(c) = EPS
                      end if
                  end do

                   dmin1 = max(EPS,minval(dtry1))

                end if

#ifdef DEBUG__
if (min(dmin1, dfor, dfor2) < 6.2d-3) then
   write(message,*) "MIN AT ", i,j, k
   call CCTK_WARN(1, trim(message))
   write(message,*) dmin1, dfor, dfor2
   call CCTK_WARN(1, trim(message))
endif
!if (i.eq.76 .and. j.eq.39.and.k.eq.3)then
if (i.eq.76 .and. j.eq. LY .and.k.eq.3)then
  ! write(message,*) JJ, U, B, W, q, p
   write(message,*) guu, gur, guq, gup, gqq, gpq, gpp
   call CCTK_WARN(1, trim(message))
   write(message,*) dmin1, dfor, dfor2
   call CCTK_WARN(1, trim(message))
endif
#endif
                dutime = min(dmin1, dfor, dfor2, dutime)
             endif
          end do
       end do

    end do

    if(dutime.lt.dt_old) dt_old =  dutime

    !FINAL WARNING IN CASE IS TOO SMALL

    if(dt_old.le.EPS) then
       call CCTK_INFO("WAY TOO SMALL DT FROM CFL ROUTINE!")
    endif

  end subroutine null_cfl0

end module NullEvol_cfl_test
