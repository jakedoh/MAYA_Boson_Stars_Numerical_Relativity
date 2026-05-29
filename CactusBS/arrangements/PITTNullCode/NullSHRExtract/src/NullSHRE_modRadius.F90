! vim: syntax=fortran
#include "cctk.h"

module NullSHRE_modRadius

  use cctk
  implicit none

contains

  subroutine wt_r(nn1, nn2, pp, eta0, r0)
! the bondi radius eq. (53) 
    use NullSHRE_modGFdef
    implicit none

    type (gf2d), dimension (4,4),  intent (in)    :: eta0
    type (gf2d),                   intent (inout) :: r0
    CCTK_INT,                      intent (in)    :: nn1,nn2
    CCTK_REAL, dimension(nn1,nn2), intent (in)    :: pp
    
    r0%d = 0.5D0 * (eta0(2,2)%d * eta0(3,3)%d - eta0(2,3)%d ** 2) ** (1.D0/4.D0) * pp

  end subroutine wt_r


  subroutine wt_rl(etaup0, eta1, r0, temp, dr0, halt_on_negative_rl, rl_min)
!lambda derivative of bondi radius eq. (54)
    use NullSHRE_modGFdef
    implicit none

    type (gf2d), dimension (4,4), intent (in)    :: etaup0, eta1
    type (gf2d),                  intent (in)    :: r0
    type (gf2d),                  intent (inout) :: temp
    type (gf2d), dimension (4),   intent (inout) :: dr0
    CCTK_INT,                     intent (in)    :: halt_on_negative_rl
    CCTK_REAL,                    intent (in)    :: rl_min ! minimum allowed value
    CCTK_INT :: a, b

    character(200) message

    ! lambda derivative of bondi r

    temp%d = 0.d0
    do a = 2, 3
       do b = 2, 3
         temp%d = temp%d + etaup0(a,b)%d * eta1(a,b)%d
       end do
    end do
    dr0(1)%d = (r0%d / 4.d0) * temp%d

    if (halt_on_negative_rl.ne.0) then
       if (any(dr0(1)%d.lt.0)) call CCTK_WARN(0, "negative value for r1_wt")
    else
       if (any(dr0(1)%d.lt.rl_min)) then
          write (message,'("r_{,lambda} is less than minimum accepted value: ",g15.5," < ",g15.5)') minval(dr0(1)%d), rl_min
          call CCTK_INFO(message)
          call CCTK_INFO("adding artificial correction term to r1_wt")
          dr0(1)%d = sqrt( dr0(1)%d**2 + rl_min**2 )
       end if
    end if

  end subroutine wt_rl


  subroutine wt_ra_long (nn1, nn2, qs, ps, pp,&
             &g, dg, j0, dj0, etaup0, r0, temp, deta0, dr0)
!algebraic computation of bondi radius angular derivatives eq. (55)
    
    use NullSHRE_modGFdef
    implicit none

    CCTK_INT,                           intent (in)       :: nn1, nn2
    CCTK_REAL,   dimension (nn1, nn2),  intent (in)       :: qs, ps, pp
    type (gf2d), dimension (4,4),       intent (in)       :: g
    type (gf2d), dimension (4,4,4),     intent (in)       :: dg
    type (gf2d), dimension (4,4),       intent (in)       :: j0
    type (gf2d), dimension (3,2:3,2:3), intent (in)       :: dj0
    type (gf2d), dimension (4,4),       intent (in)       :: etaup0
    type (gf2d),                        intent (in)       :: r0
    type (gf2d),                        intent (inout)    :: temp
    type (gf2d), dimension (2:3,2:3,2:3), intent (inout)  :: deta0
    type (gf2d), dimension (4),         intent (inout) :: dr0

    CCTK_INT :: i, j, a, b, c

    ! angular derivatives of bondi r

    do c = 2, 3

       dr0(c)%d = 0.d0

       do a = 2, 3
          do b = 2, 3

! angular derivatives of eta_{ab} -- eq. (56)
             ! first we add the g contribution

             temp%d = 0.d0

             do i = 1, 3
                do j = 1, 3
                   temp%d = temp%d + (dj0(i,a,c)%d * j0(j,b)%d &
                          + j0(i,a)%d * dj0(j,b,c)%d) * g(i,j)%d
                end do
             end do

             ! and then we add the dg part

             do i = 1, 3
                do j = 1, 3
                      temp%d = temp%d + j0(i,a)%d * j0(j,b)%d * dg(i,j,c)%d
                end do
             end do

             deta0(a,b,c)%d = temp%d

! then we multiply by eta^{ab}, and accumulate

             dr0(c)%d = dr0(c)%d + etaup0(a,b)%d * temp%d

          end do
       end do
    end do

    ! then we substract the angular derivative of |q_{ab}| * r^4,
    ! and divide by |q_{ab}|
! eq. (55)
    dr0(2)%d = 0.25 * r0%d * (dr0(2)%d + 8.d0 / pp * qs)
    dr0(3)%d = 0.25 * r0%d * (dr0(3)%d + 8.d0 / pp * ps)

  end subroutine wt_ra_long


  subroutine wt_ru(dg, j0, etaup0, r0, temp, dr0)
! time derivative of bondi r eq. (57)

    use NullSHRE_modGFdef
    implicit none

    type (gf2d), dimension (4,4,4), intent (in)    :: dg
    type (gf2d), dimension (4,4),   intent (in)    :: j0, etaup0
    type (gf2d),                    intent (in)    :: r0
    type (gf2d),                    intent (inout) :: temp
    type (gf2d), dimension (4),     intent (inout) :: dr0

    CCTK_INT :: i, j, a, b


    dr0(4)%d = 0.d0
    do a = 2, 3
       do b = 2, 3
          temp%d = 0.d0
          do i = 1, 3
             do j = 1, 3
!eq. (58)
                temp%d = temp%d + j0(i,a)%d * j0(j,b)%d * dg(i,j,4)%d
             end do
          end do
          dr0(4)%d = dr0(4)%d + etaup0(a,b)%d * temp%d
       end do
    end do
    dr0(4)%d = (r0%d / 4.d0) * dr0(4)%d

  end subroutine wt_ru


  subroutine wt_rll (beta_l, dr0, dr1)

    use NullSHRE_modGFdef
    implicit none

    type (gf2d),                intent (in)    :: beta_l
    type (gf2d), dimension (4), intent (in)    :: dr0
    type (gf2d), dimension (4), intent (inout) :: dr1

      dr1(1)%d = -2.d0 * beta_l%d * dr0(1)%d

   end subroutine wt_rll


   subroutine wt_rlu (dt, dr0, rl_old, dr1)

    use NullSHRE_modGFdef
    implicit none

      CCTK_REAL,                  intent (in)    :: dt
      type (gf2d), dimension (4), intent (in)    :: dr0
      type (gf2d), dimension (2), intent (inout) :: rl_old
      type (gf2d), dimension (4), intent (inout) :: dr1


         dr1(4)%d = (3.d0 * dr0(1)%d - 4.d0 * rl_old(1)%d + rl_old(2)%d) / (2.d0 * dt)

   end subroutine wt_rlu

end module NullSHRE_modRadius
