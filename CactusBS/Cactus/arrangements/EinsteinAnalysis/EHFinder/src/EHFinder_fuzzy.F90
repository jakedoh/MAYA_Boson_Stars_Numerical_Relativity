! Module with routines for fuzzy comparisons.
! $Header$
#include "cctk.h"

module EHFinder_fuzzy

  use EHFinder_Constants

  interface fuzzy_eq
    module procedure fuzzy_eq_0d, fuzzy_eq_1d, fuzzy_eq_2d, fuzzy_eq_3d
  end interface

  interface fuzzy_ne
    module procedure fuzzy_ne_0d, fuzzy_ne_1d, fuzzy_ne_2d, fuzzy_ne_3d
  end interface

contains

  function fuzzy_eq_0d ( x, y, eps )
    CCTK_REAL, intent(IN) :: x
    CCTK_REAL, intent(IN) :: y, eps
    logical :: fuzzy_eq_0d

    fuzzy_eq_0d = ( ( x .gt. y-eps ) .and. ( x .lt. y+eps ) )
  end function fuzzy_eq_0d 

  function fuzzy_eq_1d ( x, y, eps )
    CCTK_REAL, intent(IN), dimension(:) :: x
    CCTK_REAL, intent(IN) :: y, eps
    logical, dimension(size(x,1)) :: fuzzy_eq_1d

    fuzzy_eq_1d = ( ( x .gt. y-eps ) .and. ( x .lt. y+eps ) )
  end function fuzzy_eq_1d 

  function fuzzy_eq_2d ( x, y, eps )
    CCTK_REAL, intent(IN), dimension(:,:) :: x
    CCTK_REAL, intent(IN) :: y, eps
    logical, dimension(size(x,1),size(x,2)) :: fuzzy_eq_2d

    fuzzy_eq_2d = ( ( x .gt. y-eps ) .and. ( x .lt. y+eps ) )
  end function fuzzy_eq_2d 

  function fuzzy_eq_3d ( x, y, eps )
    CCTK_REAL, intent(IN), dimension(:,:,:) :: x
    CCTK_REAL, intent(IN) :: y, eps
    logical, dimension(size(x,1),size(x,2),size(x,3)) :: fuzzy_eq_3d

    fuzzy_eq_3d = ( ( x .gt. y-eps ) .and. ( x .lt. y+eps ) )
  end function fuzzy_eq_3d 

  function fuzzy_ne_0d ( x, y, eps )
    CCTK_REAL, intent(IN) :: x
    CCTK_REAL, intent(IN) :: y, eps
    logical :: fuzzy_ne_0d

    fuzzy_ne_0d = ( ( x .lt. y-eps ) .or. ( x .gt. y+eps ) )
  end function fuzzy_ne_0d 

  function fuzzy_ne_1d ( x, y, eps )
    CCTK_REAL, intent(IN), dimension(:) :: x
    CCTK_REAL, intent(IN) :: y, eps
    logical, dimension(size(x,1)) :: fuzzy_ne_1d

    fuzzy_ne_1d = ( ( x .lt. y-eps ) .or. ( x .gt. y+eps ) )
  end function fuzzy_ne_1d 

  function fuzzy_ne_2d ( x, y, eps )
    CCTK_REAL, intent(IN), dimension(:,:) :: x
    CCTK_REAL, intent(IN) :: y, eps
    logical, dimension(size(x,1),size(x,2)) :: fuzzy_ne_2d

    fuzzy_ne_2d = ( ( x .lt. y-eps ) .or. ( x .gt. y+eps ) )
  end function fuzzy_ne_2d 

  function fuzzy_ne_3d ( x, y, eps )
    CCTK_REAL, intent(IN), dimension(:,:,:) :: x
    CCTK_REAL, intent(IN) :: y, eps
    logical, dimension(size(x,1),size(x,2),size(x,3)) :: fuzzy_ne_3d

    fuzzy_ne_3d = ( ( x .lt. y-eps ) .or. ( x .gt. y+eps ) )
  end function fuzzy_ne_3d 

end module EHFinder_fuzzy
