/*@@
  @file      brilldata.F
  @date
  @author    Carsten Gundlach, Miguel Alcubierre.
  @desc
             Construct Brill wave q function.
  @enddesc
  @version   $Header$
@@*/

#include "cctk.h" 
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

      function brillq(rho,z,phi)

c     Calculates the function q that appear in the conformal
c     metric for Brill waves:
c
c     ds^2  =  psi^4 ( e^(2q) (drho^2 + dz^2) + rho^2 dphi^2 )
c
c     There are three different choices for the form of q depending
c     on the value of the parameter "brill_q":
c
c     brill_q = 0:
c                                           2
c                           2+b  - (rho-rho0)      2     2
c                  q = a rho    e            (z/sz)  /  r
c
c
c     brill_q = 1:  (includes Eppleys form as special case)
c
c
c                                  b               2    2      2  c/2
c                  q = a (rho/srho)  /  { 1 + [ ( r - r0 ) / sr  ]    }
c
c
c     brill_q = 2:  (includes Holz et al form as special case)
c
c                                            2    2      2  c/2 
c                                  b  - [ ( r - r0 ) / sr  ]
c                  q = a (rho/srho)  e
c
c
c     If we want 3D initial data, the function q is multiplied by an
c     additional factor:
c
c                           m    2                              m
c         q -> q [ 1 + d rho  cos  (n (phi)+phi0 ) / ( 1 + e rho ) ]

      implicit none

      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

      logical firstcall

      integer qtype

      CCTK_REAL brillq,rho,z,phi

      data firstcall /.true./

      save firstcall,qtype

c     Get parameters at first call.

      if (firstcall) then

         if (CCTK_EQUALS(q_function,"exp")) then
            qtype = 0
         else if (CCTK_EQUALS(q_function,"eppley")) then
            qtype = 1
         else if (CCTK_EQUALS(q_function,"gundlach")) then
            qtype = 2
         else
            call CCTK_WARN(0,"Brill wave data type not recognised")
         end if

         firstcall = .false.

      end if

      if (rho < 0) then
         rho = -rho
      end if

c     Calculate q(rho,z) from a choice of forms.

c     "exp" brill data

      if (qtype.eq.0) then

         brillq = exp_a*abs(rho)**(2 + exp_b)
     .          / exp((rho - exp_rho0)**2)/(rho**2 + z**2)

         if (exp_sigmaz.ne.0.0D0) then
            brillq = brillq/exp((z/exp_sigmaz)**2)
         end if

      else if (qtype.eq.1) then

c     "Eppley" brill data.  This includes Eppleys choice of q. 
c     But note that q(Eppley) = 2q(Cactus).

         brillq = eppley_a*(rho/eppley_sigmarho)**eppley_b 
     &          / ( 1.0D0 + ((rho**2 + z**2 - eppley_r0**2)/
     &          eppley_sigmar**2)**(eppley_c/2) )

      else if (qtype.eq.2) then

c     "Gundlach" brill data   This includes my (Carstens) notion of what a
c     smooth "pure quadrupole" q should be, plus the choice of
c     Holz et al.

         brillq = gundlach_a*(rho/gundlach_sigmarho)**gundlach_b
     &          / exp(((rho**2 + z**2 - gundlach_r0**2)/
     &            gundlach_sigmar**2)**(gundlach_c/2))

      else

c     Unknown type for q function.

         call CCTK_WARN(0,"brillq: Unknown type of Brill function q")

      end if

c     If desired, multiply with a phi dependent factor.

      if (CCTK_EQUALS(initial_data,"brilldata")) then
         brillq = brillq*(1.0D0 + brill3d_d*rho**brill3d_m*
     &        cos(brill3D_n*(phi-brill3d_phi0))**2
     &        / (1.0D0 + brill3d_e*rho**brill3d_m))
      end if

      return
      end
