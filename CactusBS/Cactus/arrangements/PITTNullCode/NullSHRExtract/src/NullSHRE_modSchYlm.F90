! vim: syntax=fortran
#include "cctk.h"

module NullSHRE_modSchYlm

  contains

  subroutine  SchYlm (l, m, nq, np, qs, ps, pp, Ylm)

     implicit none

     CCTK_INT, intent(in) :: l, m, nq, np
     CCTK_REAL, dimension(nq, np), intent(in) :: qs, ps, pp
     CCTK_COMPLEX, dimension(nq, np, 2) :: Ylm
     CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
     CCTK_COMPLEX :: ii
   
     ii = dcmplx(0.,1.)

     if (l.lt.0 ) then
       call CCTK_WARN(0, "Ylm ERROR: l < 0 ")
     endif
   
     if (l.gt.3 ) then
       call CCTK_WARN(0, "Ylm not implemented: l > 3 ")
     endif

     if (abs(m).gt.l)  then
       call CCTK_WARN(0, "Ylm ERROR: abs(m) > l")
     endif

!l=0, m=0
     if (l.eq.0) then      
        Ylm(:,:,1) = 1/sqrt(pi)/2
        Ylm(:,:,2) = 1/sqrt(pi)/2
     endif

!l=1
     if (l.eq.1) then
!m=-1
       if (m.eq.-1) then       
          Ylm(:,:,1) = -sqrt(6.D0)*(-qs+ii*ps)/pp/sqrt(pi)/2
          Ylm(:,:,2) = sqrt(6.D0)/sqrt(pi)/pp*(qs+ii*ps)/2
!m=0
       else if (m.eq.0) then
          Ylm(:,:,1) = -sqrt(3.D0)*(-1+qs**2+ps**2)/sqrt(pi)/pp/2
          Ylm(:,:,2) = sqrt(3.D0)*(-1+qs**2+ps**2)/sqrt(pi)/pp/2
!m=1
       else if (m.eq.1) then
          Ylm(:,:,1) = -sqrt(6.D0)/sqrt(pi)/pp*(qs+ii*ps)/2
          Ylm(:,:,2) = sqrt(6.D0)*(-qs+ii*ps)/sqrt(pi)/pp/2
       endif
     endif

!l=2
     if (l.eq.2) then
!m=-2
       if (m.eq.-2) then
          Ylm(:,:,1) = sqrt(30.D0)*(-qs+ii*ps)**2/pp**2/sqrt(pi)/2
          Ylm(:,:,2) = sqrt(30.D0)/sqrt(pi)/pp**2*(qs+ii*ps)**2/2
!m=-1
       else if (m.eq.-1) then
          Ylm(:,:,1) = sqrt(30.D0)*(-1+qs**2+ps**2)&
                        *(-qs+ii*ps)/pp**2/sqrt(pi)/2 
          Ylm(:,:,2) = sqrt(30.D0)/sqrt(pi)/pp**2&
                        *(-1+qs**2+ps**2)*(qs+ii*ps)/2
!m=0
       else if (m.eq.0) then
          Ylm(:,:,1) = sqrt(5.D0)/sqrt(pi)&
                *(1-4*qs**2-4*ps**2+qs**4+2*qs**2*ps**2+ps**4)/pp**2/2
          Ylm(:,:,2) = sqrt(5.D0)/sqrt(pi)&
                *(1-4*qs**2-4*ps**2+qs**4+2*qs**2*ps**2+ps**4)/pp**2/2
!m=1
       else if (m.eq.1) then 
          Ylm(:,:,1) = sqrt(30.D0)*(-1+qs**2+ps**2)&
                           *(qs+ii*ps)/sqrt(pi)/pp**2/2
          Ylm(:,:,2) = sqrt(30.D0)*(-1+qs**2+ps**2)&
                           *(-qs+ii*ps)/sqrt(pi)/pp**2/2
!m=2
       else if (m.eq.2) then
          Ylm(:,:,1) = sqrt(30.D0)/sqrt(pi)/pp**2*(qs+ii*ps)**2/2
          Ylm(:,:,2) = sqrt(30.D0)*(-qs+ii*ps)**2/sqrt(pi)/pp**2/2
       end if
     end if
	 
  end subroutine SchYlm	 


  subroutine  Sch1Ylm (l, m, nq, np, qs, ps, pp, Y1lm)

     implicit none

     CCTK_INT, intent(in) :: l, m, nq, np
     CCTK_REAL, dimension(nq, np), intent(in) :: qs, ps, pp
     CCTK_COMPLEX, dimension(nq, np, 2) :: Y1lm
     CCTK_REAL, parameter ::&
        pi = 3.1415926535897932384626433832795d0
     CCTK_COMPLEX :: ii
   
     ii = dcmplx(0.,1.)

     if (l.lt.0 ) then
       call CCTK_WARN(0, "Y1lm ERROR: l < 0 ")
     endif
   
     if (l.gt.3 ) then
       call CCTK_WARN(0, "Y1lm not implemented: l > 3 ")
     endif

     if (abs(m).gt.l)  then
       call CCTK_WARN(0, "Y1lm ERROR: abs(m) > l")
     endif

!l=0, m=0
     if (l.eq.0) then      
        Y1lm(:,:,1) = 0.d0 
        Y1lm(:,:,2) = 0.d0
     endif

!l=1
     if (l.eq.1) then
!m=-1
       if (m.eq.-1) then       
          Y1lm(:,:,1) = sqrt(3.D0)/pp/sqrt(pi)/2 
          Y1lm(:,:,2) = -sqrt(3.D0)*(qs+dcmplx(0.D0,1.D0)*ps)**2/pp/sqrt(pi)/2
!m=0
       else if (m.eq.0) then
          Y1lm(:,:,1) = -sqrt(3.D0)*sqrt(2.D0)*(qs+dcmplx(0.D0,1.D0)*ps)&
                       /pp/sqrt(pi)/2 
          Y1lm(:,:,2) = sqrt(3.D0)*sqrt(2.D0)*(qs+dcmplx(0.D0,1.D0)*ps)&
                       /pp/sqrt(pi)/2
!m=1
       else if (m.eq.1) then
          Y1lm(:,:,1) = sqrt(3.D0)*(qs+dcmplx(0.D0,1.D0)*ps)**2/pp/sqrt(pi)/2 
          Y1lm(:,:,2) = -sqrt(3.D0)/pp/sqrt(pi)/2
       endif
     endif

!l=2
     if (l.eq.2) then
!m=-2
       if (m.eq.-2) then
          Y1lm(:,:,1) = sqrt(5.D0)*(qs+dcmplx(0.D0,-1.D0)*ps)/pp**2/sqrt(pi) 
          Y1lm(:,:,2) = -sqrt(5.D0)*(qs+dcmplx(0.D0,1.D0)*ps)**3/pp**2/sqrt(pi)
!m=-1
       else if (m.eq.-1) then
          Y1lm(:,:,1) = -sqrt(5.D0)*(-1+3*qs**2+3*ps**2)/pp**2/sqrt(pi)/2
          Y1lm(:,:,2) = -sqrt(5.D0)*(qs+dcmplx(0.D0,1.D0)*ps)**2&
                        *(qs**2+ps**2-3)/pp**2/sqrt(pi)/2
!m=0
       else if (m.eq.0) then
          Y1lm(:,:,1) = (-1+qs**2+ps**2)*(qs+dcmplx(0.D0,1.D0)*ps)&
                       *sqrt(5.D0)*sqrt(6.D0)/pp**2/sqrt(pi)/2
          Y1lm(:,:,2) = (-1+qs**2+ps**2)*(qs+dcmplx(0.D0,1.D0)*ps)&
                       *sqrt(5.D0)*sqrt(6.D0)/pp**2/sqrt(pi)/2
!m=1
       else if (m.eq.1) then 
          Y1lm(:,:,1) = -sqrt(5.D0)*(qs+dcmplx(0.D0,1.D0)*ps)**2&
                        *(qs**2+ps**2-3)/pp**2/sqrt(pi)/2
          Y1lm(:,:,2) = -sqrt(5.D0)*(-1+3*qs**2+3*ps**2)/pp**2/sqrt(pi)/2 
!m=2
       else if (m.eq.2) then
          Y1lm(:,:,1) = -sqrt(5.D0)*(qs+dcmplx(0.D0,1.D0)*ps)**3/pp**2/sqrt(pi)
          Y1lm(:,:,2) = sqrt(5.D0)*(qs+dcmplx(0.D0,-1.D0)*ps)/pp**2/sqrt(pi)
       end if
     end if
	 
   end subroutine Sch1Ylm	 

end  module NullSHRE_modSchYlm
