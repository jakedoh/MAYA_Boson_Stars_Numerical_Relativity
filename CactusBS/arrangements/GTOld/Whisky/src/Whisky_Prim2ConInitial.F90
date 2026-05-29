 /*@@
   @routine   Prim2ConInitial.F90 
   @date       Thu Aug  6 18:56:07 EDT 2009
   @author     Roland Haas
   @desc 
   Converts primitive to conserved variables for the initial data setup.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "Whisky_Utils.h"

 /*@@
   @routine    Prim2ConInitial
   @date       Thu Aug  6 18:56:07 EDT 2009
   @author     
   @desc 
   Wrapper function that converts primitive to conservative at the 
     cell centres everywhere
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Prim2ConInitial(CCTK_ARGUMENTS)

  implicit none
  
  INTERFACE
    subroutine prim2con(handle, gxx, gxy, gxz, gyy, gyz, gzz, det, ddens, &
      dsx, dsy, dsz, dtau , drho, dvelx, dvely, dvelz, deps, dpress, w) 
      implicit none
      CCTK_REAL :: gxx, gxy, gxz, gyy, gyz, gzz, det
      CCTK_REAL :: ddens, dsx, dsy, dsz, dtau, drho, dvelx, dvely, dvelz,&
         deps, dpress, w   
      CCTK_INT :: handle
    end subroutine
  END INTERFACE

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  CCTK_INT :: i, j, k
  CCTK_REAL :: det
  CHARACTER(len=400) :: infoline
  
  !$OMP PARALLEL DO PRIVATE (i,j,k,det,infoline)
  do k = 1,cctk_lsh(3)
    do j = 1,cctk_lsh(2)
      do i = 1,cctk_lsh(1)
        
        det = SPATIAL_DET(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))

        !! Actual conversion
        call prim2con(whisky_eos_handle,gxx(i,j,k),&
             gxy(i,j,k),gxz(i,j,k),&
             gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
             det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
             tau(i,j,k),rho(i,j,k),velx(i,j,k),vely(i,j,k),velz(i,j,k),&
             eps(i,j,k),press(i,j,k),w_lorentz(i,j,k))

        if ((w_lorentz(i, j, k).lt.1.d0 .or. w_lorentz(i,j,k).ne.w_lorentz(i,j,k)) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
          !$OMP CRITICAL
          write(infoline,'("Prim2ConInitial: Unphysical Lorentz factor ",g15.6,&
            &" for data g = [",5(g15.6,","),g15.6,"] vel = [",&
            &2(g15.6,","),g15.6,"] occured in iteration ",i6," at location [",&
            &2(g15.6,","),g15.6,"] on reflevel ",i2," weight ",f4.2)') &
            w_lorentz(i,j,k), &
            gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
            velx(i,j,k),vely(i,j,k),velz(i,j,k),&
            cctk_iteration, x(i,j,k),y(i,j,k),z(i,j,k), whisky_reflevel, &
            Whisky_CarpetWeights(i,j,k)
          call CCTK_WARN(1,infoline)
          !$OMP END CRITICAL
        end if
        if ((tau(i,j,k).le.0 .or. tau(i,j,k).ne.tau(i,j,k)) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
          !$OMP CRITICAL
          write(infoline,'("Prim2ConInitial: Unphysical tau = ",g15.6," found for data rho = ",&
            &g15.6,", eps = ",g15.6,", press = ",g15.6,", w = ",g15.6,&
            &" occured in iteration ",i6," at location [",&
            &2(g15.6,","),g15.6,"] on reflevel ",i2," weight ",f4.2)') &
            tau(i,j,k), rho(i,j,k),eps(i,j,k),press(i,j,k),&
            w_lorentz(i,j,k),&
            cctk_iteration, x(i,j,k),y(i,j,k),z(i,j,k), whisky_reflevel, &
            Whisky_CarpetWeights(i,j,k)
          call CCTK_WARN(1,infoline)
          !$OMP END CRITICAL
        end if

      end do
    end do
  end do
  !$OMP END PARALLEL DO

end subroutine Prim2ConInitial

subroutine Prim2ConMHDInitial(CCTK_ARGUMENTS)

  implicit none
  
  INTERFACE
    subroutine prim2conmhd(handle, gxx, gxy, gxz, gyy, gyz, gzz, det, ddens, &
      dsx, dsy, dsz, dtau, drho, dvelx, dvely, dvelz, deps, dpress, w, &
      bnx, bny, bnz, bvecx, bvecy, bvecz) 
      implicit none
      CCTK_REAL :: gxx, gxy, gxz, gyy, gyz, gzz, det
      CCTK_REAL :: ddens, dsx, dsy, dsz, dtau, drho, dvelx, dvely, dvelz,&
         deps, dpress, w, bnx, bny, bnz, bvecx, bvecy, bvecz 
      CCTK_INT :: handle
    end subroutine
  END INTERFACE

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  CCTK_INT :: i, j, k
  CCTK_REAL :: det
  CHARACTER(len=400) :: infoline
  
  !$OMP PARALLEL DO PRIVATE (i,j,k,det,infoline)
  do k = 1,cctk_lsh(3)
    do j = 1,cctk_lsh(2)
      do i = 1,cctk_lsh(1)
        
        det = SPATIAL_DET(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
        if (Whisky_CarpetWeights(i,j,k).eq.1.d0.and.(det .ne. det .or. det .lt. 0d0)) then
            !$OMP CRITICAL
            write(infoline,'("Metric determinant in iteration ",i6," at [",g15.6,",",g15.6,",",&
              &g15.6,"] on level ",i2," is : ",g15.6)') cctk_iteration,&
              x(i,j,k),y(i,j,k),z(i,j,k),whisky_reflevel,det
            call CCTK_WARN(1, infoline)
            !$OMP END CRITICAL
        end if

        !! Actual conversion
        call prim2conmhd(whisky_eos_handle,gxx(i,j,k),&
             gxy(i,j,k),gxz(i,j,k),&
             gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
             det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
             tau(i,j,k),rho(i,j,k),velx(i,j,k),vely(i,j,k),velz(i,j,k),&
             eps(i,j,k),press(i,j,k),w_lorentz(i,j,k),&
             bnx(i,j,k),bny(i,j,k),bnz(i,j,k),bvec(i,j,k,1), &
             bvec(i,j,k,2),bvec(i,j,k,3))

        if ((w_lorentz(i, j, k).lt.1.d0 .or. w_lorentz(i,j,k).ne.w_lorentz(i,j,k)) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
          !$OMP CRITICAL
          write(infoline,'("Prim2ConInitial: Unphysical Lorentz factor ",g15.6,&
            &" for data g = [",5(g15.6,","),g15.6,"] vel = [",&
            &2(g15.6,","),g15.6,"] occured in iteration ",i6," at location [",&
            &2(g15.6,","),g15.6,"] on reflevel ",i2," weight ",f4.2)') &
            w_lorentz(i,j,k), &
            gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
            velx(i,j,k),vely(i,j,k),velz(i,j,k),&
            cctk_iteration, x(i,j,k),y(i,j,k),z(i,j,k), whisky_reflevel, &
            Whisky_CarpetWeights(i,j,k)
          call CCTK_WARN(1,infoline)
          !$OMP END CRITICAL
        end if
        if ((tau(i,j,k).le.0 .or. tau(i,j,k).ne.tau(i,j,k)) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
          !$OMP CRITICAL
          write(infoline,'("Prim2ConInitial: Unphysical tau = ",g15.6," found for data rho = ",&
            &g15.6,", eps = ",g15.6,", press = ",g15.6,", w = ",g15.6,&
            &" occured in iteration ",i6," at location [",&
            &2(g15.6,","),g15.6,"] on reflevel ",i2," weight ",f4.2)') &
            tau(i,j,k), rho(i,j,k),eps(i,j,k),press(i,j,k),&
            w_lorentz(i,j,k),&
            cctk_iteration, x(i,j,k),y(i,j,k),z(i,j,k), whisky_reflevel, &
            Whisky_CarpetWeights(i,j,k)
          call CCTK_WARN(1,infoline)
          !$OMP END CRITICAL
        end if

      end do
    end do
  end do
  !$OMP END PARALLEL DO

end subroutine Prim2ConMHDInitial
