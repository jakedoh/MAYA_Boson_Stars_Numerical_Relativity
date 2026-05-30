  /*@@
   @file      primitive2conservative
   @date      Thu Jan  11 11:03:32 2002
   @author    Pedro Montero, Ian Hawke
   @desc 
   Primitive to conservative routine
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "Whisky_Utils.h"

 /*@@
   @routine   primitive2conservative.f90 
   @date       Thu Jan 11 11:03:32 2002
   @author     Pedro Montero, Ian Hawke
   @desc 
   Converts primitive to conserved variables for the boundary extended data.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine primitive2conservative(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  integer :: i, j, k
  CCTK_REAL :: gxxl,gxyl,gxzl,gyyl,gyzl,gzzl,avg_detl,&
       gxxr,gxyr,gxzr,gyyr,gyzr,gzzr,avg_detr
  CCTK_REAL :: localbxplus, localbyplus, localbzplus
  CCTK_REAL :: localbxminus, localbyminus, localbzminus
  CHARACTER(len=600) :: infoline
  
  CCTK_LOOP3_DECLARE(whisky_prim2con_loop)

  ! Intel 11.1 causes SEGFAULTs when declaring array whose size is unknown
  ! at compile time as FIRSTPRIVATE. We use SHARED insteadd.
  !$OMP PARALLEL PRIVATE (i,j,k,gxxl,gxyl,gxzl,gyyl,gyzl,gzzl,avg_detl,&
  !$OMP gxxr,gxyr,gxzr,gyyr,gyzr,gzzr,avg_detr,localbxplus,localbyplus,localbzplus,&
  !$OMP localbxminus,localbyminus,localbzminus,infoline)

  CCTK_LOOP3( whisky_prim2con_loop, i,j,k,
              whisky_stencil, whisky_stencil, whisky_stencil,
              cctk_lsh(1)-whisky_stencil+1,
              cctk_lsh(2)-whisky_stencil+1, 
              cctk_lsh(3)-whisky_stencil+1, 
              cctk_ash(1), cctk_ash(2), cctk_ash(3))
        
        gxxl = 0.5d0 * (gxx(i,j,k) + gxx(i-xoffset,j-yoffset,k-zoffset))
        gxyl = 0.5d0 * (gxy(i,j,k) + gxy(i-xoffset,j-yoffset,k-zoffset))
        gxzl = 0.5d0 * (gxz(i,j,k) + gxz(i-xoffset,j-yoffset,k-zoffset))
        gyyl = 0.5d0 * (gyy(i,j,k) + gyy(i-xoffset,j-yoffset,k-zoffset))
        gyzl = 0.5d0 * (gyz(i,j,k) + gyz(i-xoffset,j-yoffset,k-zoffset))
        gzzl = 0.5d0 * (gzz(i,j,k) + gzz(i-xoffset,j-yoffset,k-zoffset))
        gxxr = 0.5d0 * (gxx(i,j,k) + gxx(i+xoffset,j+yoffset,k+zoffset))
        gxyr = 0.5d0 * (gxy(i,j,k) + gxy(i+xoffset,j+yoffset,k+zoffset))
        gxzr = 0.5d0 * (gxz(i,j,k) + gxz(i+xoffset,j+yoffset,k+zoffset))
        gyyr = 0.5d0 * (gyy(i,j,k) + gyy(i+xoffset,j+yoffset,k+zoffset))
        gyzr = 0.5d0 * (gyz(i,j,k) + gyz(i+xoffset,j+yoffset,k+zoffset))
        gzzr = 0.5d0 * (gzz(i,j,k) + gzz(i+xoffset,j+yoffset,k+zoffset))

        avg_detl = SPATIAL_DET(gxxl,gxyl,gxzl,gyyl,gyzl,gzzl)
        avg_detr = SPATIAL_DET(gxxr,gxyr,gxzr,gyyr,gyzr,gzzr)

        !!if ( ((avg_detl.lt.0.d0) .or. (avg_detr.lt.0.d0)) .and. (Whisky_CarpetWeights(i,j,k).eq.1.d0) ) then
        !!   call CCTK_WARN(1,"Prim2Con(Bdies) will fail because det<0 on one or both sides");
        !!end if

        if ( whisky_mhd_handle.gt.1 ) then
            call prim2conmhd(whisky_eos_handle, gxxl,gxyl,gxzl,gyyl,gyzl,gzzl, &
                 avg_detl,densminus(i,j,k),sxminus(i,j,k),&
                 syminus(i,j,k),szminus(i,j,k),tauminus(i,j,k),rhominus(i,j,k), &
                 velxminus(i,j,k),velyminus(i,j,k),velzminus(i,j,k),&
                 epsminus(i,j,k),pressminus(i,j,k),w_lorentzminus(i, j, k), &
                 bnxminus(i,j,k), bnyminus(i,j,k), bnzminus(i,j,k), bxminus(i,j,k), &
                 byminus(i,j,k), bzminus(i,j,k))
            call prim2conmhd(whisky_eos_handle, gxxr,gxyr,gxzr,gyyr,gyzr,gzzr, &
                 avg_detr,densplus(i,j,k),sxplus(i,j,k),&
                 syplus(i,j,k),szplus(i,j ,k),tauplus(i,j,k),&
                 rhoplus(i,j,k),velxplus(i,j,k),velyplus(i,j,k),&
                 velzplus(i,j,k),epsplus(i,j,k),pressplus(i,j,k),&
                 w_lorentzplus(i,j,k), bnxplus(i,j,k), bnyplus(i,j,k), &
                 bnzplus(i,j,k), bxplus(i,j,k), byplus(i,j,k), bzplus(i,j,k))
            localbxplus=bxplus(i,j,k)
            localbyplus=byplus(i,j,k)
            localbzplus=bzplus(i,j,k)
            localbxminus=bxminus(i,j,k)
            localbyminus=byminus(i,j,k)
            localbzminus=bzminus(i,j,k)
        else
            call prim2con(whisky_eos_handle, gxxl,gxyl,gxzl,gyyl,gyzl,gzzl, &
                 avg_detl,densminus(i,j,k),sxminus(i,j,k),&
                 syminus(i,j,k),szminus(i,j,k),tauminus(i,j,k),rhominus(i,j,k), &
                 velxminus(i,j,k),velyminus(i,j,k),velzminus(i,j,k),&
                 epsminus(i,j,k),pressminus(i,j,k),w_lorentzminus(i, j, k))
            call prim2con(whisky_eos_handle, gxxr,gxyr,gxzr,gyyr,gyzr,gzzr, &
                 avg_detr, densplus(i,j,k),sxplus(i,j,k),&
                 syplus(i,j,k),szplus(i,j ,k),tauplus(i,j,k),&
                 rhoplus(i,j,k),velxplus(i,j,k),velyplus(i,j,k),&
                 velzplus(i,j,k),epsplus(i,j,k),pressplus(i,j,k),&
                 w_lorentzplus(i,j,k)) 
            localbxplus=0.d0
            localbyplus=0.d0
            localbzplus=0.d0
            localbxminus=0.d0
            localbyminus=0.d0
            localbzminus=0.d0
        end if

        if ((w_lorentzminus(i, j, k).lt.1.d0 .or. w_lorentzminus(i,j,k).ne.w_lorentzminus(i,j,k)) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
          !$OMP CRITICAL
          write(infoline,'("Primitive2Conservative: Unphysical Lorentz factor (l)",g15.6,&
            &" for data g = [",5(g15.6,","),g15.6,"] vel = [",&
            &2(g15.6,","),g15.6,"] Bvec=[",2(g15.6,","),g15.6,"], det=",g15.6," occured in iteration ",i6," at location [",&
            &2(g15.6,","),g15.6,"] on reflevel ",i2," weight ",f4.2)') &
            w_lorentzminus(i,j,k), gxxl,gxyl,gxzl,gyyl,gyzl,gzzl,&
            velxminus(i,j,k),velyminus(i,j,k),velzminus(i,j,k),&
            localbxminus,localbyminus,localbzminus,avg_detl,&
            cctk_iteration,&
            x(i,j,k)-0.5d0*xoffset*CCTK_DELTA_SPACE(1),&
            y(i,j,k)-0.5d0*yoffset*CCTK_DELTA_SPACE(2),z(i,j,k)-0.5d0*zoffset*CCTK_DELTA_SPACE(3),&
            whisky_reflevel, Whisky_CarpetWeights(i,j,k)
          call CCTK_WARN(1,infoline)
          !$OMP END CRITICAL
        end if
        if ((tauminus(i,j,k).le.0 .or. tauminus(i,j,k).ne.tauminus(i,j,k)) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
          !$OMP CRITICAL
          write(infoline,'("Primitive2Conservative: Unphysical tau (l) = ",g15.6," found for data rho = ",&
            &g15.6,", eps = ",g15.6,", press = ",g15.6,", w = ",g15.6,&
            &" Bvec=[",2(g15.6,","),g15.6,"], det=",g15.6," occured in iteration ",i6," at location [",&
            &2(g15.6,","),g15.6,"] on reflevel ",i2," weight ",f4.2)') &
            tauminus(i,j,k), rhominus(i,j,k),epsminus(i,j,k),pressminus(i,j,k),&
            w_lorentzminus(i,j,k), localbxminus, localbyminus, localbzminus,&
            avg_detl, cctk_iteration, &
            x(i,j,k)-0.5d0*xoffset*CCTK_DELTA_SPACE(1),&
            y(i,j,k)-0.5d0*yoffset*CCTK_DELTA_SPACE(2),z(i,j,k)-0.5d0*zoffset*CCTK_DELTA_SPACE(3),&
            whisky_reflevel, &
            Whisky_CarpetWeights(i,j,k)
          call CCTK_WARN(1,infoline)
          !$OMP END CRITICAL
        end if

        if ((w_lorentzplus(i, j, k).lt.1.d0 .or. w_lorentzplus(i,j,k).ne.w_lorentzplus(i,j,k)) &
            .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
          !$OMP CRITICAL
          write(infoline,'("Primitive2Conservative: Unphysical Lorentz factor (r)",g15.6,&
            &" for data g = [",5(g15.6,","),g15.6,"] vel = [",&
            &2(g15.6,","),g15.6,"] Bvec=[",2(g15.6,","),g15.6,"], det=",g15.6," occured in iteration ",i6," at location [",&
            &2(g15.6,","),g15.6,"] on reflevel ",i2," weight ",f4.2)') &
            w_lorentzplus(i,j,k), gxxr,gxyr,gxzr,gyyr,gyzr,gzzr,&
            velxplus(i,j,k),velyplus(i,j,k),velzplus(i,j,k),&
            localbxplus,localbyplus,localbzplus,avg_detr,&
            cctk_iteration,&
            x(i,j,k)+0.5d0*xoffset*CCTK_DELTA_SPACE(1),&
            y(i,j,k)+0.5d0*yoffset*CCTK_DELTA_SPACE(2),z(i,j,k)+0.5d0*zoffset*CCTK_DELTA_SPACE(3),&
            whisky_reflevel, &
            Whisky_CarpetWeights(i,j,k)
          call CCTK_WARN(1,infoline)
          !$OMP END CRITICAL
        end if
        if ((tauplus(i,j,k).le.0 .or. tauplus(i,j,k).ne.tauplus(i,j,k)) &
            .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
          !$OMP CRITICAL
          write(infoline,'("Unphysical tau (r) = ",g15.6," found for data rho = ",&
            &g15.6,", eps = ",g15.6,", press = ",g15.6,", w = ",g15.6,&
            &" Bvec=[",2(g15.6,","),g15.6,"], det=",g15.6," occured in iteration ",i6," at location [",&
            &2(g15.6,","),g15.6,"] on reflevel ",i2," weight ",f4.2)') &
            tauplus(i,j,k), rhoplus(i,j,k),epsplus(i,j,k),pressplus(i,j,k),&
            w_lorentzplus(i,j,k), localbxplus, localbyplus, localbzplus, &
            avg_detr, cctk_iteration,&
            x(i,j,k)+0.5d0*xoffset*CCTK_DELTA_SPACE(1),&
            y(i,j,k)+0.5d0*yoffset*CCTK_DELTA_SPACE(2),z(i,j,k)+0.5d0*zoffset*CCTK_DELTA_SPACE(3),&
            whisky_reflevel, &
            Whisky_CarpetWeights(i,j,k)
          call CCTK_WARN(1,infoline)
          !$OMP END CRITICAL
        end if

  CCTK_ENDLOOP3(whisky_prim2con_loop) 

  !$OMP END PARALLEL

end subroutine primitive2conservative



 /*@@
   @routine    Primitive2ConservativeCells
   @date       Sun Mar 10 21:16:20 2002
   @author     
   @desc 
   Wrapper function that converts primitive to conservative at the 
     cell centres. 
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/


subroutine Primitive2ConservativeCells(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  CCTK_INT :: i, j, k
  CCTK_REAL :: det
  CCTK_REAL :: localbx, localby, localbz
  CHARACTER(len=600) :: infoline
  
  !$OMP PARALLEL DO PRIVATE (i,j,k,det,infoline,localbx,localby,localbz)
  do k = whisky_stencil,cctk_lsh(3)-whisky_stencil+1
    do j = whisky_stencil,cctk_lsh(2)-whisky_stencil+1
      do i = whisky_stencil,cctk_lsh(1)-whisky_stencil+1
        
        det = SPATIAL_DET(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))

        if (Whisky_CarpetWeights(i,j,k).eq.1.d0.and.(det .ne. det .or. det .lt. 0d0)) then
             !$OMP CRITICAL
            write(infoline,'("Metric determinant in iteration ",i6," at [",g15.6,",",g15.6,",",&
              &g15.6,"] on level ",i2," is : ",g15.6)') cctk_iteration,&
              x(i,j,k),y(i,j,k),z(i,j,k),whisky_reflevel,det
             call CCTK_WARN(1, infoline)
             !$OMP END CRITICAL
        end if

        if ( whisky_mhd_handle.gt.1 ) then
            call prim2conmhd(whisky_eos_handle,gxx(i,j,k),&
                 gxy(i,j,k),gxz(i,j,k),&
                 gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
                 det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
                 tau(i,j,k),rho(i,j,k),velx(i,j,k),vely(i,j,k),velz(i,j,k),&
                 eps(i,j,k),press(i,j,k),w_lorentz(i,j,k),bnx(i,j,k), &
                 bny(i,j,k),bnz(i,j,k),bvec(i,j,k,1),bvec(i,j,k,2),bvec(i,j,k,3))
            localbx = bvec(i,j,k,1)
            localby = bvec(i,j,k,2)
            localbz = bvec(i,j,k,3)
        else
            call prim2con(whisky_eos_handle,gxx(i,j,k),&
                 gxy(i,j,k),gxz(i,j,k),&
                 gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
                 det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
                 tau(i,j,k),rho(i,j,k),velx(i,j,k),vely(i,j,k),velz(i,j,k),&
                 eps(i,j,k),press(i,j,k),w_lorentz(i,j,k))
            localbx = 0
            localby = 0
            localbz = 0
        end if
        if ((w_lorentz(i, j, k).lt.1.d0 .or. w_lorentz(i,j,k).ne.w_lorentz(i,j,k)) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
          !$OMP CRITICAL
          write(infoline,'("Primitive2ConservativeCells: Unphysical Lorentz factor ",g15.6,&
            &" for data g = [",5(g15.6,","),g15.6,"] vel = [",&
            &2(g15.6,","),g15.6,"] Bvec = [",2(g15.6,","),g15.6,"], det=",g15.6," occured in iteration ",i6," at location [",&
            &2(g15.6,","),g15.6,"] on reflevel ",i2," weight ",f4.2)') &
            w_lorentz(i,j,k), &
            gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
            velx(i,j,k),vely(i,j,k),velz(i,j,k),localbx,localby,localbz,det,&
            cctk_iteration, x(i,j,k),y(i,j,k),z(i,j,k), whisky_reflevel, &
            Whisky_CarpetWeights(i,j,k)
          call CCTK_WARN(1,infoline)
          !$OMP END CRITICAL
        end if
        if ((tau(i,j,k).le.0 .or. tau(i,j,k).ne.tau(i,j,k)) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
          !$OMP CRITICAL
          write(infoline,'("Primitive2ConservativeCells: Unphysical tau = ",g15.6," found for data rho = ",&
            &g15.6,", eps = ",g15.6,", press = ",g15.6,", w = ",g15.6,&
            &" Bvec = [",2(g15.6,","),g15.6,"], det=",g15.6," occured in iteration ",i6," at location [",&
            &2(g15.6,","),g15.6,"] on reflevel ",i2," weight ",f4.2)') &
            tau(i,j,k), rho(i,j,k),eps(i,j,k),press(i,j,k),&
            w_lorentz(i,j,k),localbx,localby,localbz,det,&
            cctk_iteration, x(i,j,k),y(i,j,k),z(i,j,k), whisky_reflevel, &
            Whisky_CarpetWeights(i,j,k)
          call CCTK_WARN(1,infoline)
          !$OMP END CRITICAL
        end if

      end do
    end do
  end do
  !$OMP END PARALLEL DO

end subroutine Primitive2ConservativeCells


 /*@@
   @routine    Prim2ConservativePolytype
   @date       Tue Mar 19 22:52:21 2002
   @author     Ian Hawke
   @desc 
   Same as first routine, only for polytropes.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/


subroutine Prim2ConservativePolytype(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  integer :: i, j, k
  CCTK_REAL :: gxxl,gxyl,gxzl,gyyl,gyzl,gzzl,avg_detl,&
       gxxr,gxyr,gxzr,gyyr,gyzr,gzzr,avg_detr
  CCTK_REAL :: localbxplus, localbyplus, localbzplus
  CCTK_REAL :: localbxminus, localbyminus, localbzminus
  CHARACTER(len=600) :: infoline
  
  ! Intel 11.1 causes SEGFAULTs when declaring array whose size is unknown
  ! at compile time as FIRSTPRIVATE. We use SHARED insteadd.
  !$OMP PARALLEL DO PRIVATE (i,j,k, gxxl,gxyl,gxzl,gyyl,gyzl,gzzl,avg_detl,&
  !$OMP localbxplus,localbxminus,localbyplus,localbyminus,localbzplus,localbzminus,&
  !$OMP gxxr,gxyr,gxzr,gyyr,gyzr,gzzr,avg_detr,infoline)
  do k = whisky_stencil,cctk_lsh(3)-whisky_stencil+1
    do j = whisky_stencil,cctk_lsh(2)-whisky_stencil+1
      do i = whisky_stencil,cctk_lsh(1)-whisky_stencil+1
        
        gxxl = 0.5d0 * (gxx(i,j,k) + gxx(i-xoffset,j-yoffset,k-zoffset))
        gxyl = 0.5d0 * (gxy(i,j,k) + gxy(i-xoffset,j-yoffset,k-zoffset))
        gxzl = 0.5d0 * (gxz(i,j,k) + gxz(i-xoffset,j-yoffset,k-zoffset))
        gyyl = 0.5d0 * (gyy(i,j,k) + gyy(i-xoffset,j-yoffset,k-zoffset))
        gyzl = 0.5d0 * (gyz(i,j,k) + gyz(i-xoffset,j-yoffset,k-zoffset))
        gzzl = 0.5d0 * (gzz(i,j,k) + gzz(i-xoffset,j-yoffset,k-zoffset))
        gxxr = 0.5d0 * (gxx(i,j,k) + gxx(i+xoffset,j+yoffset,k+zoffset))
        gxyr = 0.5d0 * (gxy(i,j,k) + gxy(i+xoffset,j+yoffset,k+zoffset))
        gxzr = 0.5d0 * (gxz(i,j,k) + gxz(i+xoffset,j+yoffset,k+zoffset))
        gyyr = 0.5d0 * (gyy(i,j,k) + gyy(i+xoffset,j+yoffset,k+zoffset))
        gyzr = 0.5d0 * (gyz(i,j,k) + gyz(i+xoffset,j+yoffset,k+zoffset))
        gzzr = 0.5d0 * (gzz(i,j,k) + gzz(i+xoffset,j+yoffset,k+zoffset))

        avg_detl = SPATIAL_DET(gxxl,gxyl,gxzl,gyyl,gyzl,gzzl)
        avg_detr = SPATIAL_DET(gxxr,gxyr,gxzr,gyyr,gyzr,gzzr)

        if ( whisky_mhd_handle.gt.1 ) then
            call prim2conmhdpolytype(whisky_eos_handle, gxxl,gxyl,gxzl,&
                 gyyl,gyzl,gzzl, &
                 avg_detl,densminus(i,j,k),sxminus(i,j,k),&
                 syminus(i,j,k),szminus(i,j,k),tauminus(i,j,k),rhominus(i,j,k), &
                 velxminus(i,j,k),velyminus(i,j,k),velzminus(i,j,k),&
                 epsminus(i,j,k),pressminus(i,j,k),w_lorentzminus(i, j, k), &
                 bnxminus(i,j,k), bnyminus(i,j,k), bnzminus(i,j,k), &
                 bxminus(i,j,k), byminus(i,j,k), bzminus(i,j,k))
            call prim2conmhdpolytype(whisky_eos_handle, gxxr,gxyr,gxzr,&
                 gyyr,gyzr,gzzr, &
                 avg_detr,densplus(i,j,k),sxplus(i,j,k),&
                 syplus(i,j,k),szplus(i,j ,k),tauplus(i,j,k),&
                 rhoplus(i,j,k),velxplus(i,j,k),velyplus(i,j,k),&
                 velzplus(i,j,k),epsplus(i,j,k),pressplus(i,j,k),&
                 w_lorentzplus(i,j,k), bnxplus(i,j,k), bnyplus(i,j,k), &
                 bnzplus(i,j,k), bxplus(i,j,k), byplus(i,j,k), bzplus(i,j,k)) 
            localbxplus=bxplus(i,j,k)
            localbyplus=byplus(i,j,k)
            localbzplus=bzplus(i,j,k)
            localbxminus=bxminus(i,j,k)
            localbyminus=byminus(i,j,k)
            localbzminus=bzminus(i,j,k)
        else
            call prim2conpolytype(whisky_eos_handle, gxxl,gxyl,gxzl,&
                 gyyl,gyzl,gzzl, &
                 avg_detl,densminus(i,j,k),sxminus(i,j,k),&
                 syminus(i,j,k),szminus(i,j,k),tauminus(i,j,k),rhominus(i,j,k), &
                 velxminus(i,j,k),velyminus(i,j,k),velzminus(i,j,k),&
                 epsminus(i,j,k),pressminus(i,j,k),w_lorentzminus(i, j, k))
            call prim2conpolytype(whisky_eos_handle, gxxr,gxyr,gxzr,&
                 gyyr,gyzr,gzzr, &
                 avg_detr, densplus(i,j,k),sxplus(i,j,k),&
                 syplus(i,j,k),szplus(i,j ,k),tauplus(i,j,k),&
                 rhoplus(i,j,k),velxplus(i,j,k),velyplus(i,j,k),&
                 velzplus(i,j,k),epsplus(i,j,k),pressplus(i,j,k),&
                 w_lorentzplus(i,j,k)) 
            localbxplus=0
            localbyplus=0
            localbzplus=0
            localbxminus=0
            localbyminus=0
            localbzminus=0
        end if
        if ((w_lorentzminus(i, j, k).lt.1.d0 .or. w_lorentzminus(i,j,k).ne.w_lorentzminus(i,j,k)) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
          !$OMP CRITICAL
          write(infoline,'("Prim2ConservativePolytype: Unphysical Lorentz factor (l)",g15.6,&
            &" for data g = [",5(g15.6,","),g15.6,"] vel = [",&
            &2(g15.6,","),g15.6,"] Bvec=[",2(g15.6,","),g15.6,"], det=",g15.6," occured in iteration ",i6," at location [",&
            &2(g15.6,","),g15.6,"] on reflevel ",i2," weight ",f4.2)') &
            w_lorentzminus(i,j,k), gxxl,gxyl,gxzl,gyyl,gyzl,gzzl,&
            velxminus(i,j,k),velyminus(i,j,k),velzminus(i,j,k),&
            localbxminus, localbyminus, localbzminus,avg_detl,&
            cctk_iteration,&
            x(i,j,k)-0.5d0*xoffset*CCTK_DELTA_SPACE(1),&
            y(i,j,k)-0.5d0*yoffset*CCTK_DELTA_SPACE(2),z(i,j,k)-0.5d0*zoffset*CCTK_DELTA_SPACE(3),&
            whisky_reflevel, &
            Whisky_CarpetWeights(i,j,k)
          call CCTK_WARN(1,infoline)
          !$OMP END CRITICAL
        end if
        if ((tauminus(i,j,k).le.0 .or. tauminus(i,j,k).ne.tauminus(i,j,k)) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
          !$OMP CRITICAL
          write(infoline,'("Prim2ConservativePolytype: Unphysical tau (l) = ",g15.6," found for data rho = ",&
            &g15.6,", eps = ",g15.6,", press = ",g15.6,", w = ",g15.6,&
            &" Bvec=[",2(g15.6,","),g15.6,"], det=",g15.6," occured in iteration ",i6," at location [",&
            &2(g15.6,","),g15.6,"] on reflevel ",i2," weight ",f4.2)') &
            tauminus(i,j,k), rhominus(i,j,k),epsminus(i,j,k),pressminus(i,j,k),&
            w_lorentzminus(i,j,k),localbxminus,localbyminus,localbzminus,avg_detl,&
            cctk_iteration, &
            x(i,j,k)-0.5d0*xoffset*CCTK_DELTA_SPACE(1),&
            y(i,j,k)-0.5d0*yoffset*CCTK_DELTA_SPACE(2),z(i,j,k)-0.5d0*zoffset*CCTK_DELTA_SPACE(3),&
            whisky_reflevel, &
            Whisky_CarpetWeights(i,j,k)
          call CCTK_WARN(1,infoline)
          !$OMP END CRITICAL
        end if

        if ((w_lorentzplus(i, j, k).lt.1.d0 .or. w_lorentzplus(i,j,k).ne.w_lorentzplus(i,j,k)) &
            .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
          !$OMP CRITICAL
          write(infoline,'("Prim2ConservativePolytype: Unphysical Lorentz factor (r)",g15.6,&
            &" for data g = [",5(g15.6,","),g15.6,"] vel = [",&
            &2(g15.6,","),g15.6,"] Bvec=[",2(g15.6,","),g15.6,"], det=",g15.6," occured in iteration ",i6," at location [",&
            &2(g15.6,","),g15.6,"] on reflevel ",i2," weight ",f4.2)') &
            w_lorentzplus(i,j,k), gxxr,gxyr,gxzr,gyyr,gyzr,gzzr,&
            velxplus(i,j,k),velyplus(i,j,k),velzplus(i,j,k),&
            localbxplus, localbyplus, localbzplus,avg_detr,&
            cctk_iteration,&
            x(i,j,k)+0.5d0*xoffset*CCTK_DELTA_SPACE(1),&
            y(i,j,k)+0.5d0*yoffset*CCTK_DELTA_SPACE(2),z(i,j,k)+0.5d0*zoffset*CCTK_DELTA_SPACE(3),&
            whisky_reflevel, &
            Whisky_CarpetWeights(i,j,k)
          call CCTK_WARN(1,infoline)
          !$OMP END CRITICAL
        end if
        if ((tauplus(i,j,k).le.0 .or. tauplus(i,j,k).ne.tauplus(i,j,k)) &
            .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
          !$OMP CRITICAL
          write(infoline,'("Unphysical tau (r) = ",g15.6," found for data rho = ",&
            &g15.6,", eps = ",g15.6,", press = ",g15.6,", w = ",g15.6,&
            &" Bvec=[",2(g15.6,","),g15.6,"], det=",g15.6," occured in iteration ",i6," at location [",&
            &2(g15.6,","),g15.6,"] on reflevel ",i2," weight ",f4.2)') &
            tauplus(i,j,k), rhoplus(i,j,k),epsplus(i,j,k),pressplus(i,j,k),&
            w_lorentzplus(i,j,k), localbxplus, localbyplus, localbzplus,&
            avg_detr, cctk_iteration,&
            x(i,j,k)+0.5d0*xoffset*CCTK_DELTA_SPACE(1),&
            y(i,j,k)+0.5d0*yoffset*CCTK_DELTA_SPACE(2),z(i,j,k)+0.5d0*zoffset*CCTK_DELTA_SPACE(3),&
            whisky_reflevel, &
            Whisky_CarpetWeights(i,j,k)
          call CCTK_WARN(1,infoline)
          !$OMP END CRITICAL
        end if

      end do
    end do
  end do
  !$OMP END PARALLEL DO

end subroutine Prim2ConservativePolytype

 /*@@
   @routine    Primitive2ConservativePolyCells
   @date       Sun Mar 10 21:16:20 2002
   @author     
   @desc 
   Wrapper function that converts primitive to conservative at the 
     cell centres. 
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/


subroutine Primitive2ConservativePolyCells(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_INT :: i, j, k
  CCTK_REAL :: det, localbx,localby,localbz
  CHARACTER(len=600) :: infoline
  
  ! Intel 11.1 causes SEGFAULTs when declaring array whose size is unknown
  ! at compile time as FIRSTPRIVATE. We use SHARED insteadd.
  !$OMP PARALLEL DO PRIVATE (i,j,k,det,infoline,localbx,localby,localbz)
  do k = whisky_stencil,cctk_lsh(3)-whisky_stencil+1
    do j = whisky_stencil,cctk_lsh(2)-whisky_stencil+1
      do i = whisky_stencil,cctk_lsh(1)-whisky_stencil+1
        
        det = SPATIAL_DET(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))

        if (Whisky_CarpetWeights(i,j,k).eq.1.d0.and.(det .ne. det .or. det .lt. 0d0)) then
            !$OMP CRITICAL
            write(infoline,'("Metric determinant in iteration ",i6," at [",g15.6,",",g15.6,",",&
              &g15.6,"] on level ",i2," is : ",g15.6)') cctk_iteration,&
              x(i,j,k),y(i,j,k),z(i,j,k),whisky_reflevel,det
            call CCTK_WARN(1, infoline)
            !$OMP END CRITICAL
        end if
        if ( whisky_mhd_handle.gt.1 ) then
            call prim2conmhdpolytype(whisky_eos_handle,gxx(i,j,k),gxy(i,j,k),&
                 gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
                 det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
                 tau(i,j,k),rho(i,j,k),velx(i,j,k),vely(i,j,k),velz(i,j,k),&
                 eps(i,j,k),press(i,j,k),w_lorentz(i,j,k),bnx(i,j,k), &
                 bny(i,j,k),bnz(i,j,k),bvec(i,j,k,1),bvec(i,j,k,2),bvec(i,j,k,3))
            localbx = bvec(i,j,k,1)
            localby = bvec(i,j,k,2)
            localbz = bvec(i,j,k,3)
        else
            call prim2conpolytype(whisky_eos_handle,gxx(i,j,k),gxy(i,j,k),&
                 gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
                 det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
                 tau(i,j,k),rho(i,j,k),velx(i,j,k),vely(i,j,k),velz(i,j,k),&
                 eps(i,j,k),press(i,j,k),w_lorentz(i,j,k))
            localbx = 0
            localby = 0
            localbz = 0
        end if
        if ((w_lorentz(i, j, k).lt.1.d0 .or. w_lorentz(i,j,k).ne.w_lorentz(i,j,k)) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
          !$OMP CRITICAL
          write(infoline,'("Primitive2ConservativePolyCells: Unphysical Lorentz factor ",g15.6,&
            &" for data g = [",5(g15.6,","),g15.6,"] vel = [",&
            &2(g15.6,","),g15.6,"] Bvec=[",2(g15.6,","),g15.6,"] occured in iteration ",i6," at location [",&
            &2(g15.6,","),g15.6,"] on reflevel ",i2," weight ",f4.2)') &
            w_lorentz(i,j,k), &
            gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
            velx(i,j,k),vely(i,j,k),velz(i,j,k), localbx,localby,localbz,&
            cctk_iteration, x(i,j,k),y(i,j,k),z(i,j,k), whisky_reflevel, &
            Whisky_CarpetWeights(i,j,k)
          call CCTK_WARN(1,infoline)
          !$OMP END CRITICAL
        end if
        if ((tau(i,j,k).le.0 .or. tau(i,j,k).ne.tau(i,j,k)) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
          !$OMP CRITICAL
          write(infoline,'("Primitive2ConservativePolyCells: Unphysical tau = ",g15.6," found for data rho = ",&
            &g15.6,", eps = ",g15.6,", press = ",g15.6,", w = ",g15.6,&
            &" Bvec=[",2(g15.6,","),g15.6,"] occured in iteration ",i6," at location [",&
            &2(g15.6,","),g15.6,"] on reflevel ",i2," weight ",f4.2)') &
            tau(i,j,k), rho(i,j,k),eps(i,j,k),press(i,j,k),&
            w_lorentz(i,j,k),localbx,localby,localbz,&
            cctk_iteration, x(i,j,k),y(i,j,k),z(i,j,k), whisky_reflevel, &
            Whisky_CarpetWeights(i,j,k)
          call CCTK_WARN(1,infoline)
          !$OMP END CRITICAL
        end if

      end do
    end do
  end do
  !$OMP END PARALLEL DO

end subroutine Primitive2ConservativePolyCells

 /*@@
   @routine    Prim2ConservativeTracer
   @date       Mon Mar  8 13:32:32 2004
   @author     Ian Hawke
   @desc 
   Gets the conserved tracer variable from the primitive.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Prim2ConservativeTracer(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  integer :: i, j, k
  CCTK_REAL :: gxxl,gxyl,gxzl,gyyl,gyzl,gzzl,avg_detl,&
       gxxr,gxyr,gxzr,gyyr,gyzr,gzzr,avg_detr
  
  ! Intel 11.1 causes SEGFAULTs when declaring array whose size is unknown
  ! at compile time as FIRSTPRIVATE. We use SHARED insteadd.
  !$OMP PARALLEL DO PRIVATE (i,j,k, gxxl,gxyl,gxzl,gyyl,gyzl,gzzl,avg_detl,&
  !$OMP gxxr,gxyr,gxzr,gyyr,gyzr,gzzr,avg_detr)
  do k = whisky_stencil,cctk_lsh(3)-whisky_stencil+1
    do j = whisky_stencil,cctk_lsh(2)-whisky_stencil+1
      do i = whisky_stencil,cctk_lsh(1)-whisky_stencil+1
        
        gxxl = 0.5d0 * (gxx(i,j,k) + gxx(i-xoffset,j-yoffset,k-zoffset))
        gxyl = 0.5d0 * (gxy(i,j,k) + gxy(i-xoffset,j-yoffset,k-zoffset))
        gxzl = 0.5d0 * (gxz(i,j,k) + gxz(i-xoffset,j-yoffset,k-zoffset))
        gyyl = 0.5d0 * (gyy(i,j,k) + gyy(i-xoffset,j-yoffset,k-zoffset))
        gyzl = 0.5d0 * (gyz(i,j,k) + gyz(i-xoffset,j-yoffset,k-zoffset))
        gzzl = 0.5d0 * (gzz(i,j,k) + gzz(i-xoffset,j-yoffset,k-zoffset))
        gxxr = 0.5d0 * (gxx(i,j,k) + gxx(i+xoffset,j+yoffset,k+zoffset))
        gxyr = 0.5d0 * (gxy(i,j,k) + gxy(i+xoffset,j+yoffset,k+zoffset))
        gxzr = 0.5d0 * (gxz(i,j,k) + gxz(i+xoffset,j+yoffset,k+zoffset))
        gyyr = 0.5d0 * (gyy(i,j,k) + gyy(i+xoffset,j+yoffset,k+zoffset))
        gyzr = 0.5d0 * (gyz(i,j,k) + gyz(i+xoffset,j+yoffset,k+zoffset))
        gzzr = 0.5d0 * (gzz(i,j,k) + gzz(i+xoffset,j+yoffset,k+zoffset))

        avg_detl = SPATIAL_DET(gxxl,gxyl,gxzl,gyyl,gyzl,gzzl)
        avg_detr = SPATIAL_DET(gxxr,gxyr,gxzr,gyyr,gyzr,gzzr)

        cons_tracerplus(i,j,k,:) = tracerplus(i,j,k,:) * &
             sqrt(avg_detr) * rhoplus(i,j,k) / &
             sqrt(1.d0 - &
                   (gxxr * velxplus(i,j,k)**2 + &
                    gyyr * velyplus(i,j,k)**2 + &
                    gzzr * velzplus(i,j,k)**2 + &
                    2.d0 * (gxyr * velxplus(i,j,k) * velyplus(i,j,k) + &
                            gxzr * velxplus(i,j,k) * velzplus(i,j,k) + &
                            gyzr * velyplus(i,j,k) * velzplus(i,j,k) ) ) )
        cons_tracerminus(i,j,k,:) = tracerminus(i,j,k,:) * &
             sqrt(avg_detl) * rhominus(i,j,k) / &
             sqrt(1.d0 - &
                   (gxxl * velxminus(i,j,k)**2 + &
                    gyyl * velyminus(i,j,k)**2 + &
                    gzzl * velzminus(i,j,k)**2 + &
                    2.d0 * (gxyl * velxminus(i,j,k) * velyminus(i,j,k) + &
                            gxzl * velxminus(i,j,k) * velzminus(i,j,k) + &
                            gyzl * velyminus(i,j,k) * velzminus(i,j,k) ) ) )

      end do
    end do
  end do
  !$OMP END PARALLEL DO

end subroutine Prim2ConservativeTracer
