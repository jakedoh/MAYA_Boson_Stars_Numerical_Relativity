 /*@@
   @file      HLLE_Eigenproblem.F90
   @date      Sat Jan 26 01:25:44 2002
   @author    Ian Hawke, Pedro Montero, Joachim Frieben
   @desc
   Computes the spectral decomposition of a given state.
   Implements the analytical scheme devised by J. M. Ibanez
   et al., "Godunov Methods: Theory and Applications", New
   York, 2001, 485-503. The optimized method for computing
   the Roe flux in the special relativistic case is due to
   M. A. Aloy et al., Comput. Phys. Commun. 120 (1999)
   115-121, and has been extended to the general relativistic
   case as employed in this subroutine by J. Frieben, J. M.
   Ibanez, and J. Pons (in preparation).
   @enddesc
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

 /*@@
   @routine    eigenvalues_hlle
   @date       Sat Jan 26 01:26:20 2002
   @author     Ian Hawke
   @desc
   Computes the eigenvalues of the Jacobian matrix evaluated
   at the given state.
   @enddesc
   @calls
   @calledby
   @history
   Culled from the routines in GR3D, author Mark Miller.
   @endhistory

@@*/

subroutine eigenvalues_hlle(handle,rho,velx,vely,velz,eps,w_lorentz,&
     lam,gxx,gxy,gxz,gyy,gyz,gzz,u,det,alp,beta)

  implicit none

  DECLARE_CCTK_PARAMETERS

  CCTK_REAL rho,velx,vely,velz,eps,w_lorentz
  CCTK_REAL lam(5)
  CCTK_REAL gxx,gxy,gxz,gyy,gyz,gzz,det
  CCTK_REAL alp,beta,u

  CCTK_REAL cs2,one,two
  CCTK_REAL vlowx,vlowy,vlowz,v2,w
  CCTK_REAL lam1,lam2,lam3,lamm,lamp,lamm_nobeta,lamp_nobeta
  CCTK_INT handle
  CCTK_REAL dpdrho,dpdeps,press

#ifdef _EOS_BASE_INC_
#undef _EOS_BASE_INC_
#endif
#include "EOS_Base.inc"

  one = 1.0d0
  two = 2.0d0

!!$  Set required fluid quantities

  press = EOS_Pressure(handle,rho,eps)
  dpdrho = EOS_DPressByDRho(handle,rho,eps)
  dpdeps = EOS_DPressByDEps(handle,rho,eps)
  cs2 = (dpdrho + press * dpdeps / (rho**2))/ &
       (1.0d0 + eps + press/rho)

  vlowx = gxx*velx + gxy*vely + gxz*velz
  vlowy = gxy*velx + gyy*vely + gyz*velz
  vlowz = gxz*velx + gyz*vely + gzz*velz
  v2 = vlowx*velx + vlowy*vely + vlowz*velz

  w = w_lorentz

!!$  Calculate eigenvalues

  lam1 = velx - beta/alp
  lam2 = velx - beta/alp
  lam3 = velx - beta/alp
  lamp_nobeta = (velx*(one-cs2) + sqrt(cs2*(one-v2)*&
       (u*(one-v2*cs2) - velx**2*(one-cs2))))/(one-v2*cs2)
  lamm_nobeta = (velx*(one-cs2) - sqrt(cs2*(one-v2)*&
       (u*(one-v2*cs2) - velx**2*(one-cs2))))/(one-v2*cs2)

  lamp = lamp_nobeta - beta/alp
  lamm = lamm_nobeta - beta/alp

  lam(1) = lamm
  lam(2) = lam1
  lam(3) = lam2
  lam(4) = lam3
  lam(5) = lamp

end subroutine eigenvalues_hlle

 /*@@
   @routine    eigenvalues_hlle
   @date       Sat Jan 26 01:26:20 2002
   @author     Ian Hawke
   @desc
   Computes the eigenvalues of the Jacobian matrix evaluated
   at the given state.
   @enddesc
   @calls
   @calledby
   @history
   Culled from the routines in GR3D, author Mark Miller.
   @endhistory

@@*/

subroutine eigenvalues_general_hlle(&
     rho,velx,vely,velz,eps,press,cs2,&
     lam,&
     gxx,gxy,gxz,gyy,gyz,gzz,&
     u,det,alp,beta)

  implicit none

  DECLARE_CCTK_PARAMETERS

  CCTK_REAL rho,velx,vely,velz,eps
  CCTK_REAL lam(5)
  CCTK_REAL gxx,gxy,gxz,gyy,gyz,gzz,det
  CCTK_REAL alp,beta,u

  CCTK_REAL cs2,one,two
  CCTK_REAL vlowx,vlowy,vlowz,v2,w
  CCTK_REAL lam1,lam2,lam3,lamm,lamp,lamm_nobeta,lamp_nobeta
  CCTK_REAL press

  one = 1.0d0
  two = 2.0d0

!!$  Set required fluid quantities

  vlowx = gxx*velx + gxy*vely + gxz*velz
  vlowy = gxy*velx + gyy*vely + gyz*velz
  vlowz = gxz*velx + gyz*vely + gzz*velz
  v2 = vlowx*velx + vlowy*vely + vlowz*velz

  w = one / sqrt(one - v2)

!!$  Calculate eigenvalues

  lam1 = velx - beta/alp
  lam2 = velx - beta/alp
  lam3 = velx - beta/alp
  lamp_nobeta = (velx*(one-cs2) + sqrt(cs2*(one-v2)*&
       (u*(one-v2*cs2) - velx**2*(one-cs2))))/(one-v2*cs2)
  lamm_nobeta = (velx*(one-cs2) - sqrt(cs2*(one-v2)*&
       (u*(one-v2*cs2) - velx**2*(one-cs2))))/(one-v2*cs2)

  lamp = lamp_nobeta - beta/alp
  lamm = lamm_nobeta - beta/alp

  lam(1) = lamm
  lam(2) = lam1
  lam(3) = lam2
  lam(4) = lam3
  lam(5) = lamp

end subroutine eigenvalues_general_hlle
