! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

module NullEvol_Pulse
  implicit none

contains

  double complex function WKB_pulse(q, p, x, rwt, qmin, qmax, pmin, pmax,&
                      xmin, xmax, wrot, time, ID_AMP)
     implicit none
     double precision, intent(in) :: q, p, x, rwt, wrot, time
     double precision, intent(in) :: qmin, qmax, pmin, pmax, xmin, xmax


     double precision :: phase, r, pp, ID_AMP
     double complex :: ii = (0, 1.0)

     if ( x .gt. xmin .AND. x .lt. xmax .AND. q .gt. qmin&
           .AND. q .lt. qmax .AND. p .gt. pmin .AND. p .lt. pmax) then

            if ( x .eq. 1.0) then
               r = 1.0d10
            else
               r = rwt * x / ( 1.0d0 - x)
            endif

            pp = 1.0 + q**2 + p**2

            phase = wrot * (-time * (1 - 2./r)  - r &
                           +dsqrt(1 - 2./r)*q*2*r/pp)

            WKB_pulse = ID_amp*Polyr(q,qmin,qmax,5) *&
                           Polyr(p, pmin, pmax, 5) * &
                           Polyr(x, xmin, xmax,3) * (cos(phase) + &
                           ii * sin(phase))
      else
            WKB_pulse = 0.0d0
      endif
      return
   end function WKB_pulse
   double precision function Polyr(x, x1, x2,p)
    implicit none
    double precision, intent(in) :: x, x1, x2
    integer, intent(in) :: p

    if ( x < x1 .OR. x > x2 ) then
       Polyr = 0
    else
    Polyr = (x - x1)**p * (x2 - x)**p * 2**(2*p)/(x2 - x1)**(2*p)
    endif
    return
  end function Polyr
  double precision function Poly(x, x1, x2,p)
    implicit none
    double precision, intent(in) :: x, x1, x2
    integer, intent(in) :: p

    !Poly = (x - x1)**p * (x2 - x)**p * 2**(2*p)/(x2 - x1)**(2*p)

     if ( x < x2 ) then
       Poly = 1
     else if (x > x1) then
       Poly = 0
     else
     Poly =  -(((x - x1)**5*(70*x**4 + 35*x**3*x1 + 15*x**2*x1**2 +&
       5*x*x1**3 +&
     x1**4 - 9*(35*x**3 + 15*x**2*x1 + 5*x*x1**2 + x1**3)*x2 +&
     36*(15*x**2 + 5*x*x1 + x1**2)*x2**2 - 84*(5*x + x1)*x2**3 +&
      126*x2**4))/(x1 - x2)**9)
     endif
    return
  end function Poly

end module NullEvol_Pulse

