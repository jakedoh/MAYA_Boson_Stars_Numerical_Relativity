/*@@
  @file      tfd_use_module_far.F90
  @author    Erik Schnetter
  @date      2004/04/01
  @desc
             Some code using the sample Fortran module from a different thorn
  @version   $Header$
  @enddesc
@@*/


subroutine tfd_use_module_far
  use tfd_module
  implicit none
  some_variable = 2
end subroutine tfd_use_module_far
