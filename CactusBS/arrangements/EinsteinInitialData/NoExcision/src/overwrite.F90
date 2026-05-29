! $Header$

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine NoExcision_Overwrite (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL, parameter :: zero=0
  CCTK_REAL :: dist2 (cctk_lsh(1), cctk_lsh(2), cctk_lsh(3))
  CCTK_REAL :: cx, cy, cz, radx, rady, radz, width
  integer   :: n

  integer, parameter :: sm_linear = 1
  integer, parameter :: sm_spline = 2
  integer, parameter :: sm_cosine = 3
  integer :: sm_type

  do n = 1, num_regions

     cx = centre_x(n)
     cy = centre_y(n)
     cz = centre_z(n)
     if (CCTK_EQUALS(region_shape(n), "sphere")) then
        radx = radius(n)
        rady = radius(n)
        radz = radius(n)
     else if (CCTK_EQUALS(region_shape(n), "ellipsoid")) then
        radx = radius_x(n)
        rady = radius_y(n)
        radz = radius_z(n)
     else
        call CCTK_WARN (0, "internal error")
     end if

     if (CCTK_EQUALS (smoothing_function(n), "linear")) then
        sm_type = sm_linear
     else if (CCTK_EQUALS (smoothing_function(n), "spline")) then
        sm_type = sm_spline
     else if (CCTK_EQUALS (smoothing_function(n), "cosine")) then
        sm_type = sm_cosine
     else
        call CCTK_WARN (0, "internal error")
     end if
     width = smoothing_zone_width(n)

     dist2 =   ((x - cx) / radx)**2 + ((y - cy) / rady)**2 &
          &  + ((z - cz) / radz)**2

     if (overwrite_geometry(n) /= 0) then

        if (conformal_state >= 1) then
           where (dist2 <= 1)
              psi = 1
           end where
        end if
        if (conformal_state >= 2) then
           where (dist2 <= 1)
              psix = 0
              psiy = 0
              psiz = 0
           end where
        end if
        if (conformal_state >= 3) then
           where (dist2 <= 1)
              psixx = 0
              psixy = 0
              psixz = 0
              psiyy = 0
              psiyz = 0
              psizz = 0
           end where
        end if

        where (dist2 <= 1)
           gxx = smooth (gxx,  Minkowski_scale(n), dist2)
           gxy = smooth (gxy,  zero              , dist2)
           gxz = smooth (gxz,  zero              , dist2)
           gyy = smooth (gyy,  Minkowski_scale(n), dist2)
           gyz = smooth (gyz,  zero              , dist2)
           gzz = smooth (gzz,  Minkowski_scale(n), dist2)
           kxx = smooth (kxx,  zero              , dist2)
           kxy = smooth (kxy,  zero              , dist2)
           kxz = smooth (kxz,  zero              , dist2)
           kyy = smooth (kyy,  zero              , dist2)
           kyz = smooth (kyz,  zero              , dist2)
           kzz = smooth (kzz,  zero              , dist2)
        end where

     end if

     if (overwrite_lapse(n) /= 0) then

        where (dist2 <= 1)
           alp = smooth (alp, lapse_scale(n), dist2)
        end where

     end if

     if (overwrite_shift(n) /= 0) then

        if (shift_state /= 0) then
           where (dist2 <= 1)
              betax = smooth (betax, zero, dist2)
              betay = smooth (betay, zero, dist2)
              betaz = smooth (betaz, zero, dist2)
           end where
        end if

     end if

  end do
  
contains
  
  elemental function smooth (oldval, goodval, dist2) result (newval)
    CCTK_REAL             :: newval
    CCTK_REAL, intent(in) :: oldval
    CCTK_REAL, intent(in) :: goodval
    CCTK_REAL, intent(in) :: dist2
    
    CCTK_REAL, parameter :: zero=0, two=2, half=1/two
    integer,   parameter :: rk = kind(zero)
    CCTK_REAL, parameter :: pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068_rk
    CCTK_REAL :: dist
    CCTK_REAL :: x, f
    
    if (dist2 <= (1 - width)**2) then
       newval = goodval
    else if (width == 0 .or. dist2 >= 1) then
       newval = oldval
    else
       dist = sqrt(dist2)
       x = 1 - (1 - dist) / width
       select case (sm_type)
       case (sm_linear)
          f = x
       case (sm_spline)
          if (x < half) then
             f = 2 * x**2
          else
             f = 1 - 2 * (1-x)**2
          end if
       case (sm_cosine)
          f = (1 - cos (pi * x)) / 2
       end select
       newval = oldval * f + goodval * (1 - f)
    end if
    
  end function smooth
  
end subroutine NoExcision_Overwrite
