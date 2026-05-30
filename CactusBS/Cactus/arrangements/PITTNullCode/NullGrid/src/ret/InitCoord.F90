#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

subroutine InitCoord(CCTK_ARGUMENTS)
  implicit none

   CCTK_INT :: i, j, nx, ny
   CCTK_REAL :: dx, dy

   DECLARE_CCTK_ARGUMENTS
   DECLARE_CCTK_PARAMETER


   nx = -1
   ny = -1
   call getset(nx, ny)
   if (nx == -1 .OR. ny == -1) then
     call CCTK_WARN("nx, ny not setup")
   endif
   dx = 2 * mygs / (mynx -1)
   dy = 2 * mygs / (myny -1)

   do i = 1, nx
     do j = 1, ny
       
