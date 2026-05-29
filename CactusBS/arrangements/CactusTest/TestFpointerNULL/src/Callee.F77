c @@
c @file       Callee.c
c  @date      Thu 13 June 2002
c  @author    Thomas Radke
c  @desc
c             Fortran interface of the TestFpointerNULL thorn.
c  @enddesc
c  @version   $Id$
c@@

#include "cctk.h"


c  @@
c  @routine    TestFpointerNULL_Callee
c  @author     Thomas Radke
c  @date       Thu 13 June 2002
c  @desc
c              Fortran routine which gets passed a scalar and a 2D array
c              which will be set to a value != 0.
c              In addition, a dummy scalar and array are passed in as NULL
c              pointers by the C caller.
c              These variables should not be touched by this routine.
c  @enddesc
c
c  @var        null_scalar
c  @vdesc      invalid reference to a scalar variable
c  @vtype      CCTK_INT
c  @vio        in
c  @endvar
c  @var        scalar
c  @vdesc      valid reference to a scalar variable
c  @vtype      CCTK_INT
c  @vio        inout
c  @endvar
c  @var        null_array
c  @vdesc      invalid reference to an array variable
c  @vtype      CCTK_REAL
c  @vio        in
c  @endvar
c  @var        array
c  @vdesc      valid reference to an array variable
c  @vtype      CCTK_REAL
c  @vio        inout
c  @endvar
c  @var        lsh
c  @vdesc      dimensions of the array variables
c  @vtype      int *
c  @vio        in
c  @endvar
c@@
      subroutine TestFpointerNULL_Callee (null_scalar, scalar,
     &                                    null_array, array, lsh)

      implicit none

      integer lsh(2)
      CCTK_INT null_scalar, scalar
      CCTK_REAL null_array(lsh(1), lsh(2)), array(lsh(1), lsh(2))

      integer i, j


c     set the variables to some value != 0

      scalar = 42

      do j = 1, lsh(2)
        do i = 1, lsh(1)
   
          array(i, j) = 42

        end do
      end do

c     Uncomment the following lines to make the program crash.
c     null_scalar = 42
c     null_array  = 42

      end
