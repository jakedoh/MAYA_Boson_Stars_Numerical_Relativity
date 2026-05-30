#include "cctk.h"

subroutine TestFortranCrayPointers_sub(pointers, n)
  implicit none
  
  ! Find the integer type used for pointers
  ! (This should be integer*SIZEOF_CHAR_P)
  integer dummy
  pointer (pdummy, dummy)
  integer, parameter :: pk = kind(pdummy)
  
  CCTK_POINTER p
  
  ! Ensure that the pointer kinds (pointer sizes) match
  integer, parameter :: compare_sizes(1+abs(pk-kind(p))) = [0]
  
  ! An array of pointers that is passed in
  integer(pk) pointers(3)
  
  ! The array size
  integer n
  
  ! Declare local variables which are pointers, using the Fortran Cray
  ! pointer extension

  ! Explanation: The variables "a" and "pa" belong together. Only "pa"
  ! is a variable. Whenever "a" is used, the pointer which is stored
  ! in "pa" is dereferenced. In C, one would write "*pa" instead of
  ! "a".
  CCTK_REAL a(n,n), b(n,n), c(n,n)
  pointer (pa, a)
  pointer (pb, b)
  pointer (pc, c)
  
  ! Local loop indices
  integer i, j, k
  
  ! Set the pointers from the array which is passed in
  pa = pointers(1)
  pb = pointers(2)
  pc = pointers(3)
  
  ! Do some work on the arrays, as if they were not pointers
  do i = 1, n
     do j = 1, n
        a(i,j) = 0
        do k = 1, n
           a(i,j) = a(i,j) + b(i,k) * c(k,j)
        end do
     end do
  end do
  
end subroutine TestFortranCrayPointers_sub
