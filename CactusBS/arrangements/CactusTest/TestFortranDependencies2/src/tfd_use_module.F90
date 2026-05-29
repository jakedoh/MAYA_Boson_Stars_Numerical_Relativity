/*@@
  @file      tfd_use_module.F90
  @author    Erik Schnetter
  @date      2004/04/01
  @desc
             Some code using the sample Fortran module
  @version   $Header$
  @enddesc
@@*/


subroutine tfd_use_module
  use tfd_module
  implicit none
  some_variable = 1
end subroutine tfd_use_module
