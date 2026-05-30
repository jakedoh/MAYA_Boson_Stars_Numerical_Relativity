! Module with various common variable.
! $Header$
#include "cctk.h"

module IsoSurface_mod

  use EHFinder_Constants

  CCTK_REAL,dimension(:),allocatable:: coord_x,coord_y,coord_z
  CCTK_INT,dimension (:),allocatable::triangle_index
  CCTK_INT,dimension(:),allocatable::resize_triangle_index
  CCTK_REAL,dimension(:),allocatable::resize_coord_x,resize_coord_y,resize_coord_z
end module IsoSurface_mod
