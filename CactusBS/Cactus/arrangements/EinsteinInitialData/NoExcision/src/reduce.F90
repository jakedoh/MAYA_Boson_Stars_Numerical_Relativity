! $Header$

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine NoExcision_Reduce (cctk_iteration, cctk_lsh, rhs, x, y, z)
  implicit none
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  integer :: cctk_iteration
  integer :: cctk_lsh(3)
  CCTK_REAL, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)) :: rhs, x, y, z
  
  CCTK_REAL, parameter :: zero=0
  CCTK_REAL, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)) :: dist2
  CCTK_REAL :: cx, cy, cz, radx, rady, radz, width
  integer   :: n
  
  integer, parameter :: sm_linear = 1
  integer, parameter :: sm_spline = 2
  integer :: sm_type
  
  do n = 1, num_regions
     if (reduce_rhs(n) /= 0) then

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
        else
           call CCTK_WARN (0, "internal error")
        end if
        width = smoothing_zone_width(n)

        dist2 = (x - cx)**2 + (y - cy)**2 + (z - cz)**2

        ! TODO: reduce not only the RHS, but also drive the variables
        ! towards Minkowski
        where (dist2 <= 1)
           rhs = smooth (rhs, rhs * reduction_factor(n), dist2)
        end where
        
     end if
  end do
  
contains
  
  elemental function smooth (oldval, goodval, dist2) result (newval)
    CCTK_REAL             :: newval
    CCTK_REAL, intent(in) :: oldval
    CCTK_REAL, intent(in) :: goodval
    CCTK_REAL, intent(in) :: dist2
    
    CCTK_REAL, parameter :: two=2, half=1/two
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
       end select
       newval = oldval * f + goodval * (1 - f)
    end if
    
  end function smooth
  
end subroutine NoExcision_Reduce
