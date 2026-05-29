! vim: syntax=fortran
#include "cctk.h"

module NullSHRE_modAnaCoord

  use cctk
  implicit none

   contains

   subroutine wt_coord (nn1, nn2, qs, ps, pp, ip, cr, wt) ! eq. (17)

    use NullSHRE_modGFdef
    implicit none
      
      CCTK_INT,                        intent (in) :: nn1, nn2, ip
      CCTK_REAL, dimension (nn1, nn2), intent (in) :: qs, ps, pp
      CCTK_REAL,                       intent (in) :: cr
      type(gf2d), dimension (3),       intent (inout) :: wt

      ! world-tube coordinates

      wt(1)%d = 2.d0 * cr * qs / pp
      wt(2)%d = (3 - 2 * ip) * 2.d0 * cr * ps / pp
      wt(3)%d = (3 - 2 * ip) * cr * (1.d0 - qs ** 2 - ps ** 2) / pp

   end subroutine wt_coord


   subroutine wt_sigma (nn1, nn2, qs, ps, pp, ip, cr, sigma) ! eq. (23)

     use NullSHRE_modGFdef
     implicit none

      CCTK_INT,                        intent (in) :: nn1, nn2, ip
      CCTK_REAL, dimension (nn1, nn2), intent (in) :: qs, ps, pp
      CCTK_REAL,                    intent (in)    :: cr
      type(gf2d), dimension (3),     intent (inout) :: sigma

      ! 1-form normal to the world-tube

      sigma(1)%d = 8 * cr ** 2 * qs / pp ** 3
      sigma(2)%d = (3 - 2 * ip) * 8 * cr ** 2 * ps / pp ** 3
      sigma(3)%d = (3 - 2 * ip) * 4 * cr ** 2 * (1 - qs ** 2 - ps ** 2) / pp ** 3

   end subroutine wt_sigma


   subroutine wt_dsigma (nn1, nn2, qs, ps, pp, ip, cr, dsigma)

     use NullSHRE_modGFdef
     implicit none

      CCTK_INT,                        intent (in) :: nn1, nn2, ip
      CCTK_REAL, dimension (nn1, nn2), intent (in) :: qs, ps, pp
      CCTK_REAL,                    intent (in)    :: cr
      type(gf2d), dimension (3,2:3), intent (inout) :: dsigma

      ! angular derivatives of sigma

      dsigma(1,2)%d = - 8 * cr ** 2 * (- 1 + 5 * qs ** 2 - ps ** 2) /&
         & pp ** 4
      dsigma(1,3)%d =  - 48 * cr ** 2 * qs * ps / pp ** 4
      dsigma(2,2)%d = 48 * ( - 3 + 2 * ip) * cr ** 2 * ps * qs / pp ** 4
      dsigma(2,3)%d =  - 8 * ( - 3 + 2 * ip) * cr ** 2 * (1 + qs ** 2 -&
         & 5 * ps ** 2) / pp ** 4
      dsigma(3,2)%d =  - 16 * ( - 3 + 2 * ip) * cr ** 2 * qs * ( - 2 +&
         & qs ** 2 + ps ** 2) / pp ** 4
      dsigma(3,3)%d =  - 16 * ( - 3 + 2 * ip) * cr ** 2 * ps * ( - 2 +&
         & qs ** 2 + ps ** 2) / pp ** 4

   end subroutine wt_dsigma


   subroutine wt_j0 (nn1, nn2, qs, ps, pp, ip, cr, ell, j0)  ! eq. (33)

     use NullSHRE_modGFdef
     implicit none

      CCTK_INT,                        intent (in) :: nn1, nn2, ip
      CCTK_REAL, dimension (nn1, nn2), intent (in) :: qs, ps, pp
      CCTK_REAL,                       intent (in)    :: cr
      type(gf2d), dimension (4),       intent (in)    :: ell
      type(gf2d), dimension (4,4),     intent (inout) :: j0

      CCTK_INT :: i

      ! lambda derivatives

      do i = 1, 4
         j0(i,1)%d = ell(i)%d
      end do

      ! angular derivatives

      j0(1,2)%d = - 2 * cr * (-1.d0 + qs ** 2 - ps ** 2) / pp ** 2
      j0(2,2)%d = - 4 * (3 - 2 * ip) * cr * qs * ps / pp ** 2
      j0(3,2)%d = - 4 * (3 - 2 * ip) * cr * qs / pp ** 2
      j0(4,2)%d = 0.d0

      j0(1,3)%d = - 4 * cr * qs * ps / pp ** 2
      j0(2,3)%d = 2 * (3 - 2 * ip) * cr * (1.d0 + qs ** 2 - ps ** 2) / pp ** 2
      j0(3,3)%d = - 4 * (3 - 2 * ip) * cr * ps / pp ** 2
      j0(4,3)%d = 0.d0

      ! time derivatives

      j0(1,4)%d = 0.d0
      j0(2,4)%d = 0.d0
      j0(3,4)%d = 0.d0
      j0(4,4)%d = 1.d0

   end subroutine wt_j0


   subroutine wt_dj0 (nn1, nn2, qs, ps, pp, ip, cr, dj0)

      use NullSHRE_modGFdef
      implicit none

      CCTK_INT,                        intent (in) :: nn1, nn2, ip
      CCTK_REAL, dimension (nn1, nn2), intent (in) :: qs, ps, pp
      CCTK_REAL,                        intent (in)    :: cr
      type(gf2d), dimension (3,2:3,2:3), intent (inout) :: dj0

  ! angular derivatives of the o(1) jacobian

       dj0(1,2,2)%d = qs * ( - 3 + qs ** 2 - 3 * ps ** 2) * (4.d0 * cr) / pp ** 3
       dj0(1,2,3)%d = ps * ( - 1 + 3 * qs ** 2 - ps ** 2) * (4.d0 * cr) / pp ** 3
       dj0(1,3,3)%d =  - qs * (1 + qs ** 2 - 3 * ps ** 2) * (4.d0 * cr) / pp ** 3

       dj0(2,2,2)%d =  -( - 3 + 2 * ip) * ps * ( - 1 + 3 * qs ** 2 - ps ** 2) * (4.d0 * cr) / pp ** 3
       dj0(2,2,3)%d = ( - 3 + 2 * ip) * qs * (1 + qs ** 2 - 3 * ps ** 2) * (4.d0 * cr) / pp ** 3
       dj0(2,3,3)%d = ( - 3 + 2 * ip) * ps * (3 + 3 * qs ** 2 - ps ** 2) * (4.d0 * cr) / pp ** 3

       dj0(3,2,2)%d =  - ( - 3 + 2 * ip) * ( - 1 + 3 * qs ** 2 - ps ** 2) * (4.d0 * cr) / pp ** 3
       dj0(3,2,3)%d =  - 4 * ( - 3 + 2 * ip) * qs * ps * (4.d0 * cr) / pp ** 3
       dj0(3,3,3)%d = ( - 3 + 2 * ip) * (1 + qs ** 2 - 3 * ps ** 2) * (4.d0 * cr) / pp ** 3
      

   end subroutine wt_dj0

   subroutine wt_j0inv (nn1, nn2, qs, ps, pp, ip, cr, j0inv)

     use NullSHRE_modGFdef
     implicit none

      CCTK_INT,                        intent (in) :: nn1, nn2, ip
      CCTK_REAL, dimension (nn1, nn2), intent (in) :: qs, ps, pp
      CCTK_REAL,                       intent (in)    :: cr
      type(gf2d), dimension (4,4),     intent (inout) :: j0inv

      ! x,y,z derivatives of (r, q, p) coordinates

      j0inv(1,1)%d = 2*qs/pp
      j0inv(2,1)%d = -(-1+qs**2-ps**2)/cr/2
      j0inv(3,1)%d = -qs*ps/cr
      j0inv(4,1)%d = 0.d0
 
      j0inv(1,2)%d = 2*(3-2*ip)*ps/pp
      j0inv(2,2)%d = -ps*(3- 2*ip)*qs/cr
      j0inv(3,2)%d = (3 - 2*ip)/cr*(1+qs**2-ps**2)/2
      j0inv(4,2)%d = 0.d0

      j0inv(1,3)%d = -(3-2*ip)*(-1+qs**2+ps**2)/pp
      j0inv(2,3)%d = -(3-2*ip)*qs/cr
      j0inv(3,3)%d = -(3-2*ip)*ps/cr
      j0inv(4,3)%d = 0.d0

      ! time derivatives

      j0inv(1,4)%d = 0.d0
      j0inv(2,4)%d = 0.d0
      j0inv(3,4)%d = 0.d0
      j0inv(4,4)%d = 1.d0

   end subroutine wt_j0inv


end module NullSHRE_modAnaCoord
