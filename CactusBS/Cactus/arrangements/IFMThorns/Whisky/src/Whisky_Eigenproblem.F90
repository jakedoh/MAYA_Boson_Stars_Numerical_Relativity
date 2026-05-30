 /*@@
   @file      Whisky_Eigenproblem.F90
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
   @routine    eigenvalues
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

subroutine eigenvalues(handle,rho,velx,vely,velz,eps,w_lorentz,&
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

end subroutine eigenvalues

 /*@@
   @routine    eigenproblem_leftright
   @date       Sat Jan 26 01:27:59 2002
   @author     Ian Hawke, Pedro Montero, Joachim Frieben
   @desc
   Returns the left and right eigenvectors.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory

@@*/

subroutine eigenproblem_leftright(handle,rho,velx,vely,velz,eps,&
     w_lorentz,gxx,gxy,gxz,gyy,gyz,gzz,u,det,&
     alp,beta,lambda,levec,revec)

  USE Whisky_Scalars

  implicit none

  DECLARE_CCTK_PARAMETERS

  CCTK_REAL rho,velx,vely,velz,eps,w_lorentz
  CCTK_REAL lambda(5),levec(5,5),revec(5,5)
  CCTK_REAL gxx,gxy,gxz,gyy,gyz,gzz,u,det
  CCTK_REAL alp,beta

  integer i,j,k,l
  CCTK_REAL paug(5,6),tmp1,tmp2,sump,summ,f_du(5)
  integer ii,jj,kk
  CCTK_REAL leivec1(5),leivec2(5),leivec3(5),leivecp(5),leivecm(5)
  CCTK_REAL reivec1(5),reivec2(5),reivec3(5),reivecp(5),reivecm(5)
  CCTK_REAL lam1,lam2,lam3,lamm,lamp,lamm_nobeta,lamp_nobeta
  CCTK_REAL cs2,one,two
  CCTK_REAL vlowx,vlowy,vlowz,v2,w
  CCTK_REAL press,dpdrho,dpdeps,enthalpy,kappa
  CCTK_REAL axp,axm,vxp,vxm,cxx,cxy,cxz,gam,xsi,dlt,vxa,vxb
  CCTK_INT handle

!!$  Warning, warning. Nasty hack follows

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
!  if (cs2<0) cs2=0 ! this does not modify the roe crashing problem with shocktube
  enthalpy = one + eps + press / rho

  vlowx = gxx*velx + gxy*vely + gxz*velz
  vlowy = gxy*velx + gyy*vely + gyz*velz
  vlowz = gxz*velx + gyz*vely + gzz*velz
  v2 = vlowx*velx + vlowy*vely + vlowz*velz

  w = w_lorentz

!!$Calculate eigenvalues and put them in conventional order

  lam1 = velx - beta/alp
  lam2 = velx - beta/alp
  lam3 = velx - beta/alp

  lamp_nobeta = (velx*(one-cs2) + sqrt(cs2*(one-v2)*&
       (u*(one-v2*cs2) - velx**2*(one-cs2))))/(one-v2*cs2)
  lamm_nobeta = (velx*(one-cs2) - sqrt(cs2*(one-v2)*&
       (u*(one-v2*cs2) - velx**2*(one-cs2))))/(one-v2*cs2)

  lamp = lamp_nobeta - beta/alp
  lamm = lamm_nobeta - beta/alp

!!$  lam(1) = lamm
!!$  lam(2) = lam1
!!$  lam(3) = lam2
!!$  lam(4) = lam3
!!$  lam(5) = lamp

  lambda(1) = lamm
  lambda(2) = lam1
  lambda(3) = lam2
  lambda(4) = lam3
  lambda(5) = lamp

!!$Compute some auxiliary quantities

  axp = (u - velx*velx)/(u - velx*lamp_nobeta)
  axm = (u - velx*velx)/(u - velx*lamm_nobeta)
  vxp = (velx - lamp_nobeta)/(u - velx * lamp_nobeta)
  vxm = (velx - lamm_nobeta)/(u - velx * lamm_nobeta)

!!$Calculate associated right eigenvectors

  kappa = dpdeps / (dpdeps - rho * cs2)

  reivec1(1) = kappa / (enthalpy * w)
  reivec1(2) = vlowx
  reivec1(3) = vlowy
  reivec1(4) = vlowz
  reivec1(5) = one - reivec1(1)

  reivec2(1) = w * vlowy
  reivec2(2) = enthalpy * (gxy + two * w * w * vlowx * vlowy)
  reivec2(3) = enthalpy * (gyy + two * w * w * vlowy * vlowy)
  reivec2(4) = enthalpy * (gyz + two * w * w * vlowy * vlowz)
  reivec2(5) = vlowy * w * (two * w * enthalpy - one)

  reivec3(1) = w * vlowz
  reivec3(2) = enthalpy * (gxz + two * w * w * vlowx * vlowz)
  reivec3(3) = enthalpy * (gyz + two * w * w * vlowy * vlowz)
  reivec3(4) = enthalpy * (gzz + two * w * w * vlowz * vlowz)
  reivec3(5) = vlowz * w * (two * w * enthalpy - one)

  reivecp(1) = one
  reivecp(2) = enthalpy * w * (vlowx - vxp)
  reivecp(3) = enthalpy * w * vlowy
  reivecp(4) = enthalpy * w * vlowz
  reivecp(5) = enthalpy * w * axp - one

  reivecm(1) = one
  reivecm(2) = enthalpy * w * (vlowx - vxm)
  reivecm(3) = enthalpy * w * vlowy
  reivecm(4) = enthalpy * w * vlowz
  reivecm(5) = enthalpy * w * axm - one

  revec(1,:) = reivecm
  revec(2,:) = reivec1
  revec(3,:) = reivec2
  revec(4,:) = reivec3
  revec(5,:) = reivecp

!!$Calculate associated left eigenvectors if requested

  cxx = gyy * gzz - gyz * gyz
  cxy = gxz * gyz - gxy * gzz
  cxz = gxy * gyz - gxz * gyy
  gam = gxx * cxx + gxy * cxy + gxz * cxz
  xsi = cxx - gam * velx * velx
  dlt = enthalpy**3 * w * (kappa - one) * (vxm - vxp) * xsi
  
  tmp1 = w / (kappa - one)
  
  leivec1(1) = tmp1 * (enthalpy - w)
  leivec1(2) = tmp1 * w * velx
  leivec1(3) = tmp1 * w * vely
  leivec1(4) = tmp1 * w * velz
  leivec1(5) =-tmp1 * w
  
  tmp1 = one / (xsi * enthalpy)
  
  leivec2(1) = (gyz * vlowz - gzz * vlowy) * tmp1
  leivec2(2) = (gzz * vlowy - gyz * vlowz) * tmp1 * velx
  leivec2(3) = (gzz * (one - velx * vlowx) + gxz * vlowz * velx) * tmp1
  leivec2(4) = (gyz * (velx * vlowx - one) - gxz * velx * vlowy) * tmp1
  leivec2(5) = (gyz * vlowz - gzz * vlowy) * tmp1
  
  leivec3(1) = (gyz * vlowy - gyy * vlowz) * tmp1
  leivec3(2) = (gyy * vlowz - gyz * vlowy) * tmp1 * velx
  leivec3(3) = (gyz * (velx * vlowx - one) - gxy * velx * vlowz) * tmp1
  leivec3(4) = (gyy * (one - velx * vlowx) + gxy * velx * vlowy) * tmp1
  leivec3(5) = (gyz * vlowy - gyy * vlowz) * tmp1
  
  tmp1 = enthalpy * enthalpy / dlt
  tmp2 = w * w * xsi
  
  leivecp(1) = - (enthalpy * w * vxm * xsi + (one - kappa) * (vxm * &
       (tmp2 - cxx) - gam * velx) - kappa * tmp2 * vxm) * tmp1
  leivecp(2) = - (cxx * (one - kappa * axm) + (two * kappa - one) * vxm * &
       (tmp2 * velx - cxx * velx)) * tmp1
  leivecp(3) = - (cxy * (one - kappa * axm) + (two * kappa - one) * vxm * &
       (tmp2 * vely - cxy * velx)) * tmp1
  leivecp(4) = - (cxz * (one - kappa * axm) + (two * kappa - one) * vxm * &
       (tmp2 * velz - cxz * velx)) * tmp1
  leivecp(5) = - ((one - kappa) * (vxm * (tmp2 - cxx) - gam * velx) - &
       kappa * tmp2 * vxm) * tmp1
  
  leivecm(1) = (enthalpy * w * vxp * xsi + (one - kappa) * (vxp * &
       (tmp2 - cxx) - gam * velx) - kappa * tmp2 * vxp) * tmp1
  leivecm(2) = (cxx * (one - kappa * axp) + (two * kappa - one) * vxp * &
       (tmp2 * velx - cxx * velx)) * tmp1
  leivecm(3) = (cxy * (one - kappa * axp) + (two * kappa - one) * vxp * &
       (tmp2 * vely - cxy * velx)) * tmp1
  leivecm(4) = (cxz * (one - kappa * axp) + (two * kappa - one) * vxp * &
       (tmp2 * velz - cxz * velx)) * tmp1
  leivecm(5) = ((one - kappa) * (vxp * (tmp2 - cxx) - gam * velx) - &
       kappa * tmp2 * vxp) * tmp1

  levec(1,:) = leivecm
  levec(2,:) = leivec1
  levec(3,:) = leivec2
  levec(4,:) = leivec3
  levec(5,:) = leivecp

end subroutine eigenproblem_leftright

 /*@@
   @routine    eigenvalues
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

subroutine eigenvalues_general(&
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

end subroutine eigenvalues_general
