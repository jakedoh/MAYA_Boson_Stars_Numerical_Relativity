! vim: syntax=fortran
#include "cctk.h"


module NullSHRE_modGFdef

  ! define derived types

  type gf2d
     CCTK_REAL, dimension (:,:), pointer :: d
  end type gf2d

   type gf2dc
      CCTK_COMPLEX, dimension (:,:), pointer :: d
   end type gf2dc

end module NullSHRE_modGFdef
