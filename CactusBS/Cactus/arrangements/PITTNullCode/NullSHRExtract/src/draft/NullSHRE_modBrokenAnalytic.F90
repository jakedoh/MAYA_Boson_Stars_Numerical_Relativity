! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Parameters.h"

 module NullSHRE_modAnalytic

  use cctk
  implicit none

 contains

 !-----------------------------------------------------------------
 ! the analytic expressions for the Cartesian 4-metric in (p,q) coords
 ! the analytic expressions for the lambda derivative of covariant metric
 !-----------------------------------------------------------------

  subroutine modAna_g(n1, n2, qs, ps, pp, ip, cr, mass, ana_g)

    use NullSHRE_modGFdef
    implicit none

    CCTK_INT,                      intent (in) :: n1, n2, ip
    CCTK_REAL, dimension (n1, n2), intent (in) :: qs, ps, pp
    CCTK_REAL,                     intent (in) :: cr, mass
    type (gf2d), dimension (4,4),  intent(out) :: ana_g

    CCTK_INT :: sign
    sign = (3 - 2 * ip) 

    ana_g(1,1)%d = 1+8*mass*qs**2/pp**2/cr
    ana_g(2,1)%d = 8*mass/pp**2*qs*sign*ps/cr
    ana_g(3,1)%d = -4*mass/pp**2*qs*sign*(-1+qs**2+ps**2)/cr
    ana_g(4,1)%d = -4*qs*mass/cr/pp 
         
    ana_g(1,2)%d = 8*mass/pp**2*qs*sign*ps/cr
    ana_g(2,2)%d = 1+8*mass*ps**2/pp**2/cr
    ana_g(3,2)%d = -4*mass/pp**2*ps*(-1+qs**2+ps**2)/cr
    ana_g(4,2)%d = -4*mass*sign*ps/cr/pp
         
    ana_g(1,3)%d = -4*mass/pp**2*qs*sign*(-1+qs**2+ps**2)/cr
    ana_g(2,3)%d = -4*mass/pp**2*ps*(-1+qs**2+ps**2)/cr
    ana_g(3,3)%d = 1+2*mass*(-1+qs**2+ps**2)**2/pp**2/cr
    ana_g(4,3)%d = 2*mass*sign*(-1+qs**2+ps**2)/cr/pp
         
    ana_g(1,4)%d = -4*qs*mass/cr/pp
    ana_g(2,4)%d = -4*mass*sign*ps/cr/pp
    ana_g(3,4)%d = 2*mass*sign*(-1+qs**2+ps**2)/cr/pp
    ana_g(4,4)%d = (2*mass-cr)/cr
         
  end subroutine modAna_g


  subroutine modAna_dlg(n1, n2, qs, ps, pp, ip, cr, mass, ana_dlg)

    use NullSHRE_modGFdef
    implicit none

    CCTK_INT,                      intent (in) :: n1, n2, ip
    CCTK_REAL, dimension (n1, n2), intent (in) :: qs, ps, pp
    CCTK_REAL,                     intent (in) :: cr, mass
    type (gf2d), dimension (4,4),  intent(out) :: ana_dlg

    CCTK_INT :: sign
    sign = (3 - 2 * ip) 

    ana_dlg(1,1)%d = -8*qs**2*mass*(2*mass-cr)/cr**3/pp**2 
    ana_dlg(2,1)%d = -8*sign*ps*mass*qs*(2*mass-cr)/cr**3/pp**2
    ana_dlg(3,1)%d = 4*sign*(-1+qs**2+ps**2)*mass*qs*(2*mass-cr)/cr**3/pp**2
    ana_dlg(4,1)%d = 4*mass*qs*(2*mass-cr)/cr**3/pp
         
    ana_dlg(1,2)%d = -8*sign*ps*mass*qs*(2*mass-cr)/cr**3/pp**2 
    ana_dlg(2,2)%d = -8*ps**2*mass*(2*mass-cr)/cr**3/pp**2
    ana_dlg(3,2)%d = 4*mass*ps*(-1+qs**2+ps**2)*(2*mass-cr)/cr**3/pp**2
    ana_dlg(4,2)%d = 4*mass*sign*ps*(2*mass-cr)/cr**3/pp
         
    ana_dlg(1,3)%d = 4*sign*(-1+qs**2+ps**2)*mass*qs*(2*mass-cr)/cr**3/pp**2
    ana_dlg(2,3)%d = 4*mass*ps*(-1+qs**2+ps**2)*(2*mass-cr)/cr**3/pp**2
    ana_dlg(3,3)%d = -2*(-1+qs**2+ps**2)**2*mass*(2*mass-cr)/cr**3/pp**2
    ana_dlg(4,3)%d = -2*mass*sign*(-1+qs**2+ps**2)*(2*mass-cr)/cr**3/pp
         
    ana_dlg(1,4)%d = 4*mass*qs*(2*mass-cr)/cr**3/pp
    ana_dlg(2,4)%d = 4*mass*sign*ps*(2*mass-cr)/cr**3/pp
    ana_dlg(3,4)%d = -2*mass*sign*(-1+qs**2+ps**2)*(2*mass-cr)/cr**3/pp
    ana_dlg(4,4)%d = -2*mass*(2*mass-cr)/cr**3
         
  end subroutine modAna_dlg


 end module NullSHRE_modAnalytic
