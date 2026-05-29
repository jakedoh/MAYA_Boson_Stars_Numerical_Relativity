/*@@
  @file    Whisky_Eigenproblem_Marquina.F90
  @date    Wed Feb 13 12:30 2002
  @author  Pedro Montero, Toni Font, Joachim Frieben
  @desc
  Computes the eigenvectors (p) and eigenvalues (lam) for Marquina flux
  formula for the input primitive state. Note that this is the MODIFIED
  Marquina formula as given by M. A. Aloy et al., Astrophys. J. Suppl.
  122 (1999) 151, and not the full Marquina flux of Donat and Marquina.
  The optimized method for computing the Marquina flux in the special
  relativistic case is due to M. A. Aloy et al., Comput. Phys. Commun.
  120 (1999) 115-121, and has been extended to the general relativistic
  case as employed in this subroutine by J. Frieben, J. M. Ibanez, and
  J. Pons (in preparation).
  @enddesc
@@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

 /*@@
   @routine    eigenproblem_marquina
   @date       Wed Feb 13 12:27:59 2002
   @author     Pedro Montero, Toni Font, Joachim Frieben
   @desc
   Computes the eigenvectors (p) and eigenvalues (lam) for
   Marquina flux formula
   for the input primitive state.
   @enddesc
   @calls
   @calledby
   @history
   Follows the routines by Toni Font and GR3D of Mark Miller.
   @endhistory

@@*/

subroutine eigenproblem_marquina(handle,rhor,velxr,velyr,&
     velzr,epsr,rhol,velxl,velyl,velzl,epsl,&
     gxxr,gxyr,gxzr,gyyr,gyzr,gzzr,ur,detr,alpr,betar,&
     gxxl,gxyl,gxzl,gyyl,gyzl,gzzl,ul,detl,alpl,betal,&
     densl,sxl,syl,szl,taul,&
     densr,sxr,syr,szr,taur,flux1,flux2,flux3,flux4,flux5,&
     weight)

  USE Marquina_Scalars

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL rhor,velxr,velyr,velzr,epsr,w_lorentzr
  CCTK_REAL rhol,velxl,velyl,velzl,epsl,w_lorentzl
  CCTK_REAL densl,sxl,syl,szl,taul
  CCTK_REAL densr,sxr,syr,szr,taur

  CCTK_REAL lam(5),p(5,5),q(5,5),dw(5)
  CCTK_REAL rfluxr(5),rfluxl(5)
  CCTK_REAL gxxr,gxyr,gxzr,gyyr,gyzr,gzzr,ur,detr,alpr,betar
  CCTK_REAL gxxl,gxyl,gxzl,gyyl,gyzl,gzzl,ul,detl,alpl,betal
  CCTK_REAL flux1r,flux2r,flux3r,flux4r,flux5r
  CCTK_REAL flux1l,flux2l,flux3l,flux4l,flux5l
  CCTK_REAL flux1,flux2,flux3,flux4,flux5
  CCTK_REAL weight

!!$  LOCAL VARS

  CCTK_REAL du(5),aa(5,5),qdiff1,qdiff2,qdiff3,qdiff4,qdiff5
  integer i,j,k,l
  CCTK_REAL paug(5,6),tmp1,tmp2,sump,summ,f_du(5)
  CCTK_REAL leivec1l(5),leivec2l(5),leivec3l(5),leivecpl(5),leivecml(5)
  CCTK_REAL leivec1r(5),leivec2r(5),leivec3r(5),leivecpr(5),leivecmr(5)
  CCTK_REAL reivec1l(5),reivec2l(5),reivec3l(5),reivecpl(5),reivecml(5)
  CCTK_REAL reivec1r(5),reivec2r(5),reivec3r(5),reivecpr(5),reivecmr(5)
  CCTK_REAL lam1l,lam2l,lam3l,lamml,lampl,lamm_nobetal,lamp_nobetal
  CCTK_REAL lam1r,lam2r,lam3r,lammr,lampr,lamm_nobetar,lamp_nobetar
  CCTK_REAL lam1,lam2,lam3,lamm,lamp

  CCTK_REAL cs2l,cs2r,one,two
  CCTK_REAL vlowxr,vlowyr,vlowzr,v2r,wr
  CCTK_REAL vlowxl,vlowyl,vlowzl,v2l,wl

  CCTK_REAL lamp_nobeta,lamm_nobeta

  CCTK_REAL pressl,dpdrhol,dpdepsl,enthalpyl,kappal
  CCTK_REAL pressr,dpdrhor,dpdepsr,enthalpyr,kappar
  CCTK_REAL axpl,axml,vxpl,vxml,xsil,dltl
  CCTK_REAL axpr,axmr,vxpr,vxmr,xsir,dltr
  CCTK_REAL cxx,cxy,cxz,gam,vxa,vxb
  CCTK_INT handle
  character(len=400) warnline

!!$  Warning, warning. Nasty hack follows

#ifdef _EOS_BASE_INC_
#undef _EOS_BASE_INC_
#endif
#include "EOS_Base.inc"

  one = 1.0d0
  two = 2.0d0

!!$  LEFT

!!$  Set required fluid quantities

  pressl = EOS_Pressure(handle,rhol,epsl)
  dpdrhol = EOS_DPressByDRho(handle,rhol,epsl)
  dpdepsl = EOS_DPressByDEps(handle,rhol,epsl)
  cs2l = (dpdrhol + pressl * dpdepsl / (rhol**2))/ &
       (1.0d0 + epsl + pressl/rhol)
  if((cs2l .lt. 0.d0 .or. cs2l .gt. one .or. cs2l.ne.cs2l) &
     .and. weight.eq.1.d0) then
    write(warnline,*) cs2l,rhol,epsl
    call CCTK_WARN(1, warnline)
    call CCTK_WARN(1, "cs2l unphysical marquina_eigenproblem")
  end if

  enthalpyl = one + epsl + pressl / rhol

  vlowxl = gxxl*velxl + gxyl*velyl + gxzl*velzl
  vlowyl = gxyl*velxl + gyyl*velyl + gyzl*velzl
  vlowzl = gxzl*velxl + gyzl*velyl + gzzl*velzl
  v2l = vlowxl*velxl + vlowyl*velyl + vlowzl*velzl

!!$  Assume consistent primitive data

  wl = one / sqrt(one - v2l)
  if (wl.ne.wl .and. weight.eq.1.d0) then
    call CCTK_WARN(1, "Left Lorentz factor is NaN in marquina_eigenproblem")
  end if 

!!$  EIGENVALUES

  lam1l = velxl - betal/alpl
  lam2l = lam1l
  lam3l = lam1l
  lamp_nobetal = (velxl*(one-cs2l) + sqrt(cs2l*(one-v2l)* &
       (ul*(one-v2l*cs2l) - velxl**2*(one-cs2l))))/(one-v2l*cs2l)
  lamm_nobetal = (velxl*(one-cs2l) - sqrt(cs2l*(one-v2l)* &
       (ul*(one-v2l*cs2l) - velxl**2*(one-cs2l))))/(one-v2l*cs2l)
  if (lamp_nobetal.ne.lamp_nobetal .and. weight.eq.1.d0) then
    write(warnline,*) velxl,cs2l,v2l,ul
    call CCTK_WARN(1, warnline)
    call CCTK_WARN(1, "lamp_nobetal is NaN in marquina_eigenproblem")
  end if 
  if (lamm_nobetal.ne.lamm_nobetal .and. weight.eq.1.d0) then
    write(warnline,*) velxl,cs2l,v2l,ul
    call CCTK_WARN(1, warnline)
    call CCTK_WARN(1, "lamm_nobetal is NaN in marquina_eigenproblem")
  end if 

  lampl = lamp_nobetal - betal/alpl
  lamml = lamm_nobetal - betal/alpl

!!$  RIGHT

!!$  Set required fluid quantities

  pressr = EOS_Pressure(handle,rhor,epsr)
  dpdrhor = EOS_DPressByDRho(handle,rhor,epsr)
  dpdepsr = EOS_DPressByDEps(handle,rhor,epsr)
  cs2r = (dpdrhor + pressr * dpdepsr / (rhor**2))/ &
       (1.0d0 + epsr + pressr/rhor)
  if((cs2r .lt. 0.d0 .or. cs2r .gt. one .or. cs2r.ne.cs2r) &
     .and. weight.eq.1.d0) then
    write(warnline,*) cs2r,rhor,epsr
    call CCTK_WARN(1, warnline)
    call CCTK_WARN(1, "cs2l unphysical marquina_eigenproblem")
  end if

  enthalpyr = one + epsr + pressr / rhor

  vlowxr = gxxr*velxr + gxyr*velyr + gxzr*velzr
  vlowyr = gxyr*velxr + gyyr*velyr + gyzr*velzr
  vlowzr = gxzr*velxr + gyzr*velyr + gzzr*velzr
  v2r = vlowxr*velxr + vlowyr*velyr + vlowzr*velzr

!!$  Assume consistent primitive data

  wr = one / sqrt(one - v2r)
  if (wr.ne.wr .and. weight.eq.1.d0) then
    call CCTK_WARN(1, "Right Lorentz factor is NaN in marquina_eigenproblem")
  end if 

!!$  EIGENVALUES

  lam1r = velxr - betar/alpr
  lam2r = lam1r
  lam3r = lam1r
  lamp_nobetar = (velxr*(one-cs2r) + sqrt(cs2r*(one-v2r)* &
       (ur*(one-v2r*cs2r) - velxr**2*(one-cs2r))))/(one-v2r*cs2r)
  lamm_nobetar = (velxr*(one-cs2r) - sqrt(cs2r*(one-v2r)* &
       (ur*(one-v2r*cs2r) - velxr**2*(one-cs2r))))/(one-v2r*cs2r)
  if (lamp_nobetar.ne.lamp_nobetar .and. weight.eq.1.d0) then
    write(warnline,*) velxr,cs2r,v2r,ur
    call CCTK_WARN(1, warnline)
    call CCTK_WARN(1, "lamp_nobetar is NaN in marquina_eigenproblem")
  end if 
  if (lamm_nobetar.ne.lamm_nobetar .and. weight.eq.1.d0) then
    write(warnline,*) velxr,cs2r,v2r,ur
    call CCTK_WARN(1, warnline)
    call CCTK_WARN(1, "lamm_nobetar is NaN in marquina_eigenproblem")
  end if 

  lampr = lamp_nobetar - betar/alpr
  lammr = lamm_nobetar - betar/alpr

!!$  FINAL

  lam1 = dmax1(dabs(lam1l),dabs(lam1r))
  lam2 = lam1
  lam3 = lam1
  lamp = dmax1(dabs(lampl),dabs(lampr))
  lamm = dmax1(dabs(lamml),dabs(lammr))

!!$  lam(1) = lamm
!!$  lam(2) = lam1
!!$  lam(3) = lam2
!!$  lam(4) = lam3
!!$  lam(5) = lamp

  lam(1) = lamm
  lam(5) = lam1
  lam(3) = lam2
  lam(4) = lam3
  lam(2) = lamp
  
!!$  LEFT

!!$  Compute some auxiliary quantities

  axpl = (ul - velxl*velxl)/(ul - velxl*lamp_nobetal)
  axml = (ul - velxl*velxl)/(ul - velxl*lamm_nobetal)
  vxpl = (velxl - lamp_nobetal)/(ul - velxl * lamp_nobetal)
  vxml = (velxl - lamm_nobetal)/(ul - velxl * lamm_nobetal)

!!$  Calculate associated right eigenvectors

  kappal = dpdepsl / (dpdepsl - rhol * cs2l)

!!$  Right eigenvector # 1

  reivec1l(1) = kappal / (enthalpyl * wl)
  reivec1l(2) = vlowxl
  reivec1l(3) = vlowyl
  reivec1l(4) = vlowzl
  reivec1l(5) = one - reivec1l(1)

!!$  Right eigenvector # 2

  reivec2l(1) = wl * vlowyl
  reivec2l(2) = enthalpyl * (gxyl + two * wl * wl * vlowxl * vlowyl)
  reivec2l(3) = enthalpyl * (gyyl + two * wl * wl * vlowyl * vlowyl)
  reivec2l(4) = enthalpyl * (gyzl + two * wl * wl * vlowyl * vlowzl)
  reivec2l(5) = vlowyl * wl * (two * wl * enthalpyl - one)

!!$  Right eigenvector # 3

  reivec3l(1) = wl * vlowzl
  reivec3l(2) = enthalpyl * (gxzl + two * wl * wl * vlowxl * vlowzl)
  reivec3l(3) = enthalpyl * (gyzl + two * wl * wl * vlowyl * vlowzl)
  reivec3l(4) = enthalpyl * (gzzl + two * wl * wl * vlowzl * vlowzl)
  reivec3l(5) = vlowzl * wl * (two * wl * enthalpyl - one)

!!$  Right + eigenvector

  reivecpl(1) = one
  reivecpl(2) = enthalpyl * wl * (vlowxl - vxpl)
  reivecpl(3) = enthalpyl * wl * vlowyl
  reivecpl(4) = enthalpyl * wl * vlowzl
  reivecpl(5) = enthalpyl * wl * axpl - one

!!$  Right - eigenvector

  reivecml(1) = one
  reivecml(2) = enthalpyl * wl * (vlowxl - vxml)
  reivecml(3) = enthalpyl * wl * vlowyl
  reivecml(4) = enthalpyl * wl * vlowzl
  reivecml(5) = enthalpyl * wl * axml - one

!!$  RIGHT

!!$  Compute some auxiliary quantities

  axpr = (ur - velxr*velxr)/(ur - velxr*lamp_nobetar)
  axmr = (ur - velxr*velxr)/(ur - velxr*lamm_nobetar)
  vxpr = (velxr - lamp_nobetar)/(ur - velxr * lamp_nobetar)
  vxmr = (velxr - lamm_nobetar)/(ur - velxr * lamm_nobetar)

!!$  Calculate associated right eigenvectors

  kappar = dpdepsr / (dpdepsr - rhor * cs2r)

!!$  Right eigenvector # 1

  reivec1r(1) = kappar / (enthalpyr * wr)
  reivec1r(2) = vlowxr
  reivec1r(3) = vlowyr
  reivec1r(4) = vlowzr
  reivec1r(5) = one - reivec1r(1)

!!$  Right eigenvector # 2

  reivec2r(1) = wr * vlowyr
  reivec2r(2) = enthalpyr * (gxyr + two * wr * wr * vlowxr * vlowyr)
  reivec2r(3) = enthalpyr * (gyyr + two * wr * wr * vlowyr * vlowyr)
  reivec2r(4) = enthalpyr * (gyzr + two * wr * wr * vlowyr * vlowzr)
  reivec2r(5) = vlowyr * wr * (two * wr * enthalpyr - one)

!!$  Right eigenvector # 3

  reivec3r(1) = wr * vlowzr
  reivec3r(2) = enthalpyr * (gxzr + two * wr * wr * vlowxr * vlowzr)
  reivec3r(3) = enthalpyr * (gyzr + two * wr * wr * vlowyr * vlowzr)
  reivec3r(4) = enthalpyr * (gzzr + two * wr * wr * vlowzr * vlowzr)
  reivec3r(5) = vlowzr * wr * (two * wr * enthalpyr - one)

!!$  Right + eigenvector

  reivecpr(1) = one
  reivecpr(2) = enthalpyr * wr * (vlowxr - vxpr)
  reivecpr(3) = enthalpyr * wr * vlowyr
  reivecpr(4) = enthalpyr * wr * vlowzr
  reivecpr(5) = enthalpyr * wr * axpr - one

!!$  Right - eigenvector

  reivecmr(1) = one
  reivecmr(2) = enthalpyr * wr * (vlowxr - vxmr)
  reivecmr(3) = enthalpyr * wr * vlowyr
  reivecmr(4) = enthalpyr * wr * vlowzr
  reivecmr(5) = enthalpyr * wr * axmr - one

!!$  Calculate associated left eigenvectors if requested

  if (ANALYTICAL) then

!!$  LEFT

    cxx = gyyl * gzzl - gyzl * gyzl
    cxy = gxzl * gyzl - gxyl * gzzl
    cxz = gxyl * gyzl - gxzl * gyyl
    gam = gxxl * cxx + gxyl * cxy + gxzl * cxz

    xsil = cxx - gam * velxl * velxl
    dltl = enthalpyl**3 * wl * (kappal - one) * (vxml - vxpl) * xsil

!!$  Left eigenvector # 1

    tmp1 = wl / (kappal - one)

    leivec1l(1) = tmp1 * (enthalpyl - wl)
    leivec1l(2) = tmp1 * wl * velxl
    leivec1l(3) = tmp1 * wl * velyl
    leivec1l(4) = tmp1 * wl * velzl
    leivec1l(5) =-tmp1 * wl

!!$  Left eigenvector # 2

    tmp1 = one / (xsil * enthalpyl)

    leivec2l(1) = (gyzl * vlowzl - gzzl * vlowyl) * tmp1
    leivec2l(2) = (gzzl * vlowyl - gyzl * vlowzl) * tmp1 * velxl
    leivec2l(3) = (gzzl * (one - velxl * vlowxl) + gxzl * vlowzl * velxl ) * tmp1
    leivec2l(4) = (gyzl * (velxl * vlowxl - one) - gxzl * velxl * vlowyl) * tmp1
    leivec2l(5) = (gyzl * vlowzl - gzzl * vlowyl) * tmp1

!!$  Left eigenvector # 3

    leivec3l(1) = (gyzl * vlowyl - gyyl * vlowzl) * tmp1
    leivec3l(2) = (gyyl * vlowzl - gyzl * vlowyl) * tmp1 * velxl
    leivec3l(3) = (gyzl * (velxl * vlowxl - one) - gxyl * velxl * vlowzl) * tmp1
    leivec3l(4) = (gyyl * (one - velxl * vlowxl) + gxyl * velxl * vlowyl) * tmp1
    leivec3l(5) = (gyzl * vlowyl - gyyl * vlowzl) * tmp1

!!$  Left + eigenvector

    tmp1 = enthalpyl * enthalpyl / dltl
    tmp2 = wl * wl * xsil

    leivecpl(1) = - (enthalpyl * wl * vxml * xsil + (one - kappal) * &
      (vxml * (tmp2 - cxx) - gam * velxl) - kappal * tmp2 * vxml) * tmp1
    leivecpl(2) = - (cxx * (one - kappal * axml) + (two * kappal - one) * &
      vxml * (tmp2 * velxl - cxx * velxl)) * tmp1
    leivecpl(3) = - (cxy * (one - kappal * axml) + (two * kappal - one) * &
      vxml * (tmp2 * velyl - cxy * velxl)) * tmp1
    leivecpl(4) = - (cxz * (one - kappal * axml) + (two * kappal - one) * &
      vxml * (tmp2 * velzl - cxz * velxl)) * tmp1
    leivecpl(5) = - ((one - kappal) * (vxml * (tmp2 - cxx) - gam * velxl) - &
      kappal * tmp2 * vxml) * tmp1

!!$  Left - eigenvector

    leivecml(1) = (enthalpyl * wl * vxpl * xsil + (one - kappal) * &
      (vxpl * (tmp2 - cxx) - gam * velxl) - kappal * tmp2 * vxpl) * tmp1
    leivecml(2) = (cxx * (one - kappal * axpl) + (two * kappal - one) * &
      vxpl * (tmp2 * velxl - cxx * velxl)) * tmp1
    leivecml(3) = (cxy * (one - kappal * axpl) + (two * kappal - one) * &
      vxpl * (tmp2 * velyl - cxy * velxl)) * tmp1
    leivecml(4) = (cxz * (one - kappal * axpl) + (two * kappal - one) * &
      vxpl * (tmp2 * velzl - cxz * velxl)) * tmp1
    leivecml(5) = ((one - kappal) * (vxpl * (tmp2 - cxx) - gam * &
      velxl) - kappal * tmp2 * vxpl) * tmp1

!!$  RIGHT

    cxx = gyyr * gzzr - gyzr * gyzr
    cxy = gxzr * gyzr - gxyr * gzzr
    cxz = gxyr * gyzr - gxzr * gyyr
    gam = gxxr * cxx + gxyr * cxy + gxzr * cxz

    xsir = cxx - gam * velxr * velxr
    dltr = enthalpyr**3 * wr * (kappar - one) * (vxmr - vxpr) * xsir

!!$  Left eigenvector # 1

    tmp1 = wr / (kappar - one)

    leivec1r(1) = tmp1 * (enthalpyr - wr)
    leivec1r(2) = tmp1 * wr * velxr
    leivec1r(3) = tmp1 * wr * velyr
    leivec1r(4) = tmp1 * wr * velzr
    leivec1r(5) =-tmp1 * wr

!!$  Left eigenvector # 2

    tmp1 = one / (xsir * enthalpyr)

    leivec2r(1) = (gyzr * vlowzr - gzzr * vlowyr) * tmp1
    leivec2r(2) = (gzzr * vlowyr - gyzr * vlowzr) * tmp1 * velxr
    leivec2r(3) = (gzzr * (one - velxr * vlowxr) + gxzr * vlowzr * velxr) * tmp1
    leivec2r(4) = (gyzr * (velxr * vlowxr - one) - gxzr * velxr * vlowyr) * tmp1
    leivec2r(5) = (gyzr * vlowzr - gzzr * vlowyr) * tmp1

!!$  Left eigenvector # 3

    leivec3r(1) = (gyzr * vlowyr - gyyr * vlowzr) * tmp1
    leivec3r(2) = (gyyr * vlowzr - gyzr * vlowyr) * tmp1 * velxr
    leivec3r(3) = (gyzr * (velxr * vlowxr - one) - gxyr * velxr * vlowzr) * tmp1
    leivec3r(4) = (gyyr * (one - velxr * vlowxr) + gxyr * velxr * vlowyr) * tmp1
    leivec3r(5) = (gyzr * vlowyr - gyyr * vlowzr) * tmp1

!!$  Left + eigenvector

    tmp1 = enthalpyr * enthalpyr / dltr
    tmp2 = wr * wr * xsir

    leivecpr(1) = - (enthalpyr * wr * vxmr * xsir + (one - kappar) * &
      (vxmr * (tmp2 - cxx) - gam * velxr) - kappar * tmp2 * vxmr) * tmp1
    leivecpr(2) = - (cxx * (one - kappar * axmr) + (two * kappar - one) * &
      vxmr * (tmp2 * velxr - cxx * velxr)) * tmp1
    leivecpr(3) = - (cxy * (one - kappar * axmr) + (two * kappar - one) * &
      vxmr * (tmp2 * velyr - cxy * velxr)) * tmp1
    leivecpr(4) = - (cxz * (one - kappar * axmr) + (two * kappar - one) * &
      vxmr * (tmp2 * velzr - cxz * velxr)) * tmp1
    leivecpr(5) = - ((one - kappar) * (vxmr * (tmp2 - cxx) - gam * velxr) - &
      kappar * tmp2 * vxmr) * tmp1

!!$  Left - eigenvector

    leivecmr(1) = (enthalpyr * wr * vxpr * xsir + (one - kappar) * &
      (vxpr * (tmp2 - cxx) - gam * velxr) - kappar * tmp2 * vxpr) * tmp1
    leivecmr(2) = (cxx * (one - kappar * axpr) + (two * kappar - one) * &
      vxpr * (tmp2 * velxr - cxx * velxr)) * tmp1
    leivecmr(3) = (cxy * (one - kappar * axpr) + (two * kappar - one) * &
      vxpr * (tmp2 * velyr - cxy * velxr)) * tmp1
    leivecmr(4) = (cxz * (one - kappar * axpr) + (two * kappar - one) * &
      vxpr * (tmp2 * velzr - cxz * velxr)) * tmp1
    leivecmr(5) = ((one - kappar) * (vxpr * (tmp2 - cxx) - gam * &
      velxr) - kappar * tmp2 * vxpr) * tmp1

  endif

!!$  LEFT
!!$  PUT RIGHT EIGENVECTORS IN THE P MATRIX

!!$  p(:,1) = reivecml(:)
!!$  p(:,2) = reivec1l(:)
!!$  p(:,3) = reivec2l(:)
!!$  p(:,4) = reivec3l(:)
!!$  p(:,5) = reivecpl(:)

  p(:,1) = reivecml(:)
  p(:,2) = reivecpl(:)
  p(:,3) = reivec2l(:)
  p(:,4) = reivec3l(:)
  p(:,5) = reivec1l(:)

!!$  Calculate change in u:

  du(1) = densl
  du(2) = sxl
  du(3) = syl
  du(4) = szl
  du(5) = taul

  if (ANALYTICAL) then

    if (FAST) then

      sump = 0.0d0
      summ = 0.0d0

      do i=1,5
         sump = sump + (lamp - lam1) * leivecpl(i) * du(i)
         summ = summ + (lamm - lam1) * leivecml(i) * du(i)
      enddo

      vxa = sump + summ
      vxb =-(sump * vxpl + summ * vxml)

      rfluxl(1) = lam1 * du(1) + vxa
      rfluxl(2) = lam1 * du(2) + enthalpyl * wl * (vlowxl * vxa + vxb)
      rfluxl(3) = lam1 * du(3) + enthalpyl * wl * (vlowyl * vxa)
      rfluxl(4) = lam1 * du(4) + enthalpyl * wl * (vlowzl * vxa)
      rfluxl(5) = lam1 * du(5) + enthalpyl * wl * (velxl  * vxb + vxa) - vxa

    else

!!$  PUT LEFT EIGENVECTORS IN THE Q MATRIX

      q(1,:) = leivecml(:)
      q(2,:) = leivecpl(:)
      q(3,:) = leivec2l(:)
      q(4,:) = leivec3l(:)
      q(5,:) = leivec1l(:)

      do i=1,5
        dw(i) = 0.0d0
        do j=1,5
          dw(i) = dw(i) + q(i,j) * du(j)
        enddo
      enddo

      do i = 1, 5
        rfluxl(i) = 0.d0
        do j = 1, 5
          rfluxl(i) = rfluxl(i) + p(i,j) * lam(j) * dw(j)
        end do
      end do

    endif

  else

!!$  Solve for characteristic variable change, dw

    dw=du
    aa = p

    do i=1,5
      paug(:,i) = p(:,i)
    enddo

!same, but in old F77 style!
!!$    do i=1,5
!!$      dw(i) = du(i)
!!$      do j=1,5
!!$        aa(i,j) = p(i,j)
!!$      end do
!!$    enddo
!!$
!!$    do i=1,5
!!$      paug(i,1) = p(i,1)
!!$      paug(i,2) = p(i,2)
!!$      paug(i,3) = p(i,3)
!!$      paug(i,4) = p(i,4)
!!$      paug(i,5) = p(i,5)
!!$    enddo

    paug(1,6) = du(1)
    paug(2,6) = du(2)
    paug(3,6) = du(3)
    paug(4,6) = du(4)
    paug(5,6) = du(5)

!!$  Get lower left triangle to be all zeros
    do i=1,5
!!$  First, make diagonal element 1
      tmp1 = paug(i,i)
      do j=i,6
        paug(i,j) = paug(i,j)/tmp1
      enddo
!!$  Now, get rid of everything below that diagonal
      do j=i+1,5
        tmp1 = - (paug(j,i))
        do k=i,6
          paug(j,k) = paug(j,k) + tmp1*paug(i,k)
        enddo
      enddo
    enddo
!!$  Back substitute
    f_du(5) = paug(5,6)
    do i=4,1,-1
      f_du(i) = paug(i,6)
      do j=i+1,5
        f_du(i) = f_du(i) - paug(i,j)*f_du(j)
      enddo
    enddo

    dw(1) = f_du(1)
    dw(2) = f_du(2)
    dw(3) = f_du(3)
    dw(4) = f_du(4)
    dw(5) = f_du(5)

    do i = 1, 5
      rfluxl(i) = 0.d0
      do j = 1, 5
        rfluxl(i) = rfluxl(i) + p(i,j) * lam(j) * dw(j)
      end do
    end do

  endif

    flux1l = rfluxl(1)
    flux2l = rfluxl(2)
    flux3l = rfluxl(3)
    flux4l = rfluxl(4)
    flux5l = rfluxl(5)

!!$  RIGHT

!!$  PUT RIGHT EIGENVECTORS IN THE P MATRIX

!!$  p(:,1) = reivecmr(:)
!!$  p(:,2) = reivec1r(:)
!!$  p(:,3) = reivec2r(:)
!!$  p(:,4) = reivec3r(:)
!!$  p(:,5) = reivecpr(:)

  p(:,1) = reivecmr(:)
  p(:,2) = reivecpr(:)
  p(:,3) = reivec2r(:)
  p(:,4) = reivec3r(:)
  p(:,5) = reivec1r(:)

  du(1) = densr
  du(2) = sxr
  du(3) = syr
  du(4) = szr
  du(5) = taur

  if (ANALYTICAL) then

    if (FAST) then

      sump = 0.0d0
      summ = 0.0d0

      do i=1,5
         sump = sump + (lamp - lam1) * leivecpr(i) * du(i)
         summ = summ + (lamm - lam1) * leivecmr(i) * du(i)
      enddo

      vxa = sump + summ
      vxb =-(sump * vxpr + summ * vxmr)

      rfluxr(1) = lam1 * du(1) + vxa
      rfluxr(2) = lam1 * du(2) + enthalpyr * wr * (vlowxr * vxa + vxb)
      rfluxr(3) = lam1 * du(3) + enthalpyr * wr * (vlowyr * vxa)
      rfluxr(4) = lam1 * du(4) + enthalpyr * wr * (vlowzr * vxa)
      rfluxr(5) = lam1 * du(5) + enthalpyr * wr * (velxr  * vxb + vxa) - vxa

    else

!!$  PUT LEFT EIGENVECTORS IN THE Q MATRIX

      q(1,:) = leivecmr(:)
      q(2,:) = leivecpr(:)
      q(3,:) = leivec2r(:)
      q(4,:) = leivec3r(:)
      q(5,:) = leivec1r(:)

      do i=1,5
        dw(i) = 0.0d0
        do j=1,5
          dw(i) = dw(i) + q(i,j) * du(j)
        enddo
      enddo

      do i = 1, 5
        rfluxr(i) = 0.d0
        do j = 1, 5
          rfluxr(i) = rfluxr(i) + p(i,j) * lam(j) * dw(j)
        end do
      end do

    endif

  else

!!$  Solve for characteristic variable change, dw

    do i=1,5
      dw(i) = du(i)
      do j=1,5
        aa(i,j) = p(i,j)
      end do
    enddo

    do i=1,5
      paug(i,1) = p(i,1)
      paug(i,2) = p(i,2)
      paug(i,3) = p(i,3)
      paug(i,4) = p(i,4)
      paug(i,5) = p(i,5)
    enddo

    paug(1,6) = du(1)
    paug(2,6) = du(2)
    paug(3,6) = du(3)
    paug(4,6) = du(4)
    paug(5,6) = du(5)

!!$  Get lower left triangle to be all zeros

    do i=1,5
!!$  First, make diagonal element 1
      tmp1 = paug(i,i)
      do j=i,6
        paug(i,j) = paug(i,j)/tmp1
      enddo
!!$  Now, get rid of everything below that diagonal
      do j=i+1,5
        tmp1 = - (paug(j,i))
        do k=i,6
          paug(j,k) = paug(j,k) + tmp1*paug(i,k)
        enddo
      enddo
    enddo
!!$  Back substitute

    f_du(5) = paug(5,6)
    do i=4,1,-1
      f_du(i) = paug(i,6)
      do j=i+1,5
        f_du(i) = f_du(i) - paug(i,j)*f_du(j)
      enddo
    enddo

    dw(1) = f_du(1)
    dw(2) = f_du(2)
    dw(3) = f_du(3)
    dw(4) = f_du(4)
    dw(5) = f_du(5)

    do i = 1, 5
      rfluxr(i) = 0.d0
      do j = 1, 5
        rfluxr(i) = rfluxr(i) + p(i,j) * lam(j) * dw(j)
      end do
    end do

  endif

  flux1r = rfluxr(1)
  flux2r = rfluxr(2)
  flux3r = rfluxr(3)
  flux4r = rfluxr(4)
  flux5r = rfluxr(5)

  flux1 = flux1r - flux1l 
  flux2 = flux2r - flux2l 
  flux3 = flux3r - flux3l 
  flux4 = flux4r - flux4l 
  flux5 = flux5r - flux5l 

  return

end subroutine eigenproblem_marquina

