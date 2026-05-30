  /*@@
   @file      Whisky_Prim2Con_Ptwise.F90
   @date      Thu Jan  11 11:03:32 2002
   @author    Tanja Bode 
   @desc 
   Pointwise primitive to conservative functions
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

 /*@@
   @routine    prim2con
   @date       Sat Jan 26 01:52:18 2002
   @author     Pedro Montero, Ian Hawke
   @desc 
   Converts from primitive to conservative at a single point
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine prim2con(handle, gxx, gxy, gxz, gyy, gyz, gzz, det, ddens, &
     dsx, dsy, dsz, dtau , drho, dvelx, dvely, dvelz, deps, dpress, w) 
  
  implicit none
  
  CCTK_REAL :: gxx, gxy, gxz, gyy, gyz, gzz, det
  CCTK_REAL :: ddens, dsx, dsy, dsz, dtau, drho, dvelx, dvely, dvelz,&
       deps, dpress, w, vlowx, vlowy, vlowz   
  CCTK_INT :: handle

#include "EOS_Base.inc"
  
  w = 1.d0 / sqrt(1.d0 - (gxx*dvelx*dvelx + gyy*dvely*dvely + gzz &
       *dvelz*dvelz + 2*gxy*dvelx*dvely + 2*gxz*dvelx *dvelz + 2*gyz&
       *dvely*dvelz))  

  dpress = EOS_Pressure(handle, drho, deps)

  vlowx = gxx*dvelx + gxy*dvely + gxz*dvelz
  vlowy = gxy*dvelx + gyy*dvely + gyz*dvelz
  vlowz = gxz*dvelx + gyz*dvely + gzz*dvelz

  ddens = sqrt(det) * drho * w 
  dsx = sqrt(det) * (drho*(1+deps)+dpress)*w*w * vlowx
  dsy = sqrt(det) * (drho*(1+deps)+dpress)*w*w * vlowy
  dsz = sqrt(det) * (drho*(1+deps)+dpress)*w*w * vlowz
  !dtau = sqrt(det) * ((drho*(1+deps)+dpress)*w*w - dpress) - ddens 
  ! hand-crafted expression to involve only positive terms
  dtau = sqrt(det) * (drho*(w-1.d0)*w + drho*deps*w**2 + dpress*(w**2-1.d0))
  !if ( w .gt. 1.d0 ) then
  !   dtau = sqrt(det) * (drho*(w-1.d0)*w + drho*deps*w**2 + dpress*(w**2-1.d0))
  !else
  !   dtau = sqrt(det)*drho*deps
  !end if

end subroutine prim2con

 /*@@
   @routine    prim2conpolytype
   @date       Sat Jan 26 01:52:18 2002
   @author     Pedro Montero, Ian Hawke
   @desc 
   Converts from primitive to conservative at a single point
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine prim2conpolytype(handle, gxx, gxy, gxz, gyy, gyz, &
     gzz, det, ddens, &
     dsx, dsy, dsz, dtau , drho, dvelx, dvely, dvelz, deps, dpress, w) 
  
  implicit none
  
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL :: gxx, gxy, gxz, gyy, gyz, gzz, det
  CCTK_REAL :: ddens, dsx, dsy, dsz, dtau, drho, dvelx, dvely, dvelz,&
       deps, dpress, w, vlowx, vlowy, vlowz   
  CCTK_INT :: handle

#ifdef _EOS_BASE_INC_
#undef _EOS_BASE_INC_
#endif
#include "EOS_Base.inc"
  
  w = 1.d0 / sqrt(1.d0 - (gxx*dvelx*dvelx + gyy*dvely*dvely + gzz &
       *dvelz*dvelz + 2*gxy*dvelx*dvely + 2*gxz*dvelx *dvelz + 2*gyz&
       *dvely*dvelz))  
  
  dpress = EOS_Pressure(handle, drho, deps)
  deps = EOS_SpecificIntEnergy(handle, drho, dpress)

  vlowx = gxx*dvelx + gxy*dvely + gxz*dvelz
  vlowy = gxy*dvelx + gyy*dvely + gyz*dvelz
  vlowz = gxz*dvelx + gyz*dvely + gzz*dvelz

  ddens = sqrt(det) * drho * w 
  dsx = sqrt(det) * (drho*(1+deps)+dpress)*w*w * vlowx
  dsy = sqrt(det) * (drho*(1+deps)+dpress)*w*w * vlowy
  dsz = sqrt(det) * (drho*(1+deps)+dpress)*w*w * vlowz
  !dtau = sqrt(det) * ((drho*(1+deps)+dpress)*w*w - dpress) - ddens 
  ! hand-crafted expression to involve only positive terms
  dtau = sqrt(det) * (drho*(w-1.d0)*w + drho*deps*w**2 + dpress*(w**2-1.d0))

end subroutine prim2conpolytype

 /*@@
   @routine    prim2conmhd
   @date       Sat Jan 26 01:52:18 2002
   @author     Tanja Bode, (based on prim2con by Pedro Montero, Ian Hawke)
   @desc 
   Converts from primitive to conservative at a single point
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine prim2conmhd(handle, gxx, gxy, gxz, gyy, gyz, gzz, det, ddens, &
     dsx, dsy, dsz, dtau , drho, dvelx, dvely, dvelz, deps, dpress, w, &
     bnx, bny, bnz, bvecx, bvecy, bvecz) 
  
  implicit none
  
  CCTK_REAL :: gxx, gxy, gxz, gyy, gyz, gzz, det
  CCTK_REAL :: ddens, dsx, dsy, dsz, dtau, drho, dvelx, dvely, dvelz,&
       deps, dpress, w, vlowx, vlowy, vlowz
  CCTK_REAL :: bnx, bny, bnz, bvecx, bvecy, bvecz 
  CCTK_REAL :: udotB, bvec2, b2, blowx, blowy, blowz ! blow[xyz] = F^*_munu u^mu
  CCTK_INT :: handle

#include "EOS_Base.inc"
  
  bnx = sqrt(det) * bvecx
  bny = sqrt(det) * bvecy
  bnz = sqrt(det) * bvecz

  w = 1.d0 / sqrt(1.d0 - (gxx*dvelx*dvelx + gyy*dvely*dvely + gzz &
       *dvelz*dvelz + 2*gxy*dvelx*dvely + 2*gxz*dvelx *dvelz + 2*gyz&
       *dvely*dvelz))  

  dpress = EOS_Pressure(handle, drho, deps)

  vlowx = gxx*dvelx + gxy*dvely + gxz*dvelz
  vlowy = gxy*dvelx + gyy*dvely + gyz*dvelz
  vlowz = gxz*dvelx + gyz*dvely + gzz*dvelz

  udotB = w *( vlowx*bvecx + vlowy*bvecy + vlowz*bvecz ) 
  blowx = ( gxx*bvecx + gxy*bvecy + gxz*bvecz ) / w + udotB * vlowx
  blowy = ( gxy*bvecx + gyy*bvecy + gyz*bvecz ) / w + udotB * vlowy
  blowz = ( gxz*bvecx + gyz*bvecy + gzz*bvecz ) / w + udotB * vlowz
  bvec2 = gxx*bvecx*bvecx + gyy*bvecy*bvecy + gzz*bvecz*bvecz + &
        2*gxy*bvecx*bvecy + 2*gxz*bvecx*bvecz + 2*gyz*bvecy*bvecz
  b2 = ( bvec2 + udotB**2 ) / (w**2)

  ddens = sqrt(det) * drho * w 
  dsx = sqrt(det) * ((drho*(1+deps)+dpress + b2)*w**2 * vlowx - udotB*blowx)
  dsy = sqrt(det) * ((drho*(1+deps)+dpress + b2)*w**2 * vlowy - udotB*blowy)
  dsz = sqrt(det) * ((drho*(1+deps)+dpress + b2)*w**2 * vlowz - udotB*blowz)
  dtau = sqrt(det) * ((drho*(w-1.d0)*w + drho*deps*w**2 + dpress*(w**2-1.d0)) &
         + b2*( w**2-0.5d0 ) - udotB**2)


end subroutine prim2conmhd

 /*@@
   @routine    prim2conmhdpolytype
   @date       Sat Jan 26 01:52:18 2002
   @author     Pedro Montero, Ian Hawke
   @desc 
   Converts from primitive to conservative at a single point
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine prim2conmhdpolytype(handle, gxx, gxy, gxz, gyy, gyz, &
     gzz, det, ddens, dsx, dsy, dsz, dtau , &
     drho, dvelx, dvely, dvelz, deps, dpress, w, &
     bnx, bny, bnz, bvecx, bvecy, bvecz) 
  
  implicit none
  
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL :: gxx, gxy, gxz, gyy, gyz, gzz, det
  CCTK_REAL :: ddens, dsx, dsy, dsz, dtau, drho, dvelx, dvely, dvelz,&
       deps, dpress, w, vlowx, vlowy, vlowz   
  CCTK_REAL :: bnx, bny, bnz, bvecx, bvecy, bvecz 
  CCTK_REAL :: udotB, bvec2, b2, blowx, blowy, blowz
  CCTK_INT :: handle

#ifdef _EOS_BASE_INC_
#undef _EOS_BASE_INC_
#endif
#include "EOS_Base.inc"
  
  bnx = sqrt(det) * bvecx
  bny = sqrt(det) * bvecy
  bnz = sqrt(det) * bvecz

  w = 1.d0 / sqrt(1.d0 - (gxx*dvelx*dvelx + gyy*dvely*dvely + gzz &
       *dvelz*dvelz + 2*gxy*dvelx*dvely + 2*gxz*dvelx *dvelz + 2*gyz&
       *dvely*dvelz))  
  
  dpress = EOS_Pressure(handle, drho, deps)
  deps = EOS_SpecificIntEnergy(handle, drho, dpress)

  vlowx = gxx*dvelx + gxy*dvely + gxz*dvelz
  vlowy = gxy*dvelx + gyy*dvely + gyz*dvelz
  vlowz = gxz*dvelx + gyz*dvely + gzz*dvelz

  udotB = w *( vlowx*bvecx + vlowy*bvecy + vlowz*bvecz ) 
  blowx = ( gxx*bvecx + gxy*bvecy + gxz*bvecz ) / w + udotB * vlowx
  blowy = ( gxy*bvecx + gyy*bvecy + gyz*bvecz ) / w + udotB * vlowy
  blowz = ( gxz*bvecx + gyz*bvecy + gzz*bvecz ) / w + udotB * vlowz
  bvec2 = gxx*bvecx*bvecx + gyy*bvecy*bvecy + gzz*bvecz*bvecz + &
        2*gxy*bvecx*bvecy + 2*gxz*bvecx*bvecz + 2*gyz*bvecy*bvecz
  b2 = ( bvec2 + udotB**2 ) / (w**2)

  ddens = sqrt(det) * drho * w 
  dsx = sqrt(det) * ((drho*(1+deps)+dpress + b2)*w**2 * vlowx - udotB*blowx)
  dsy = sqrt(det) * ((drho*(1+deps)+dpress + b2)*w**2 * vlowy - udotB*blowy)
  dsz = sqrt(det) * ((drho*(1+deps)+dpress + b2)*w**2 * vlowz - udotB*blowz)
  dtau = sqrt(det) * ((drho*(w-1.d0)*w + drho*deps*w**2 + dpress*(w**2-1.d0)) &
         + b2*( w**2-0.5d0 ) - udotB**2)

end subroutine prim2conmhdpolytype
