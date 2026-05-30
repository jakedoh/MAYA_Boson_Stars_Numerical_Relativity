#include "cctk.h"
#include "cctk_Arguments.h"

      subroutine PrintThoseFortranStrings(CCTK_ARGUMENTS)

      DECLARE_CCTK_ARGUMENTS
      
      integer num1, num2

      num1 = 1
      num2 = 2
      
      call PrintOneString(num1, num2, "Just One")

      call PrintTwoString(num1, num2, "Just One", "Now Two")

      call PrintThreeString(num1, num2, "Just One", "Now Two", "Full Three")

      return
      end
  
