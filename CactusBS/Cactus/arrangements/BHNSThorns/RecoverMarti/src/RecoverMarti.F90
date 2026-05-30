/*@@
   @file      RecoverMarti.c
   @date      Sat Jan 26 01:06:01 2002
   @author    The Whisky Developers
   @desc 
   The routines for converting conservative to primitive variables.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "SpaceMask.h"

 /*@@
   @routine    Conservative2Primitive
   @date       Sat Jan 26 01:08:33 2002
   @author     Ian Hawke
   @desc 
   Wrapper routine that converts from conserved to primitive variables
   at every grid cell centre.
   @enddesc 
   @calls     
   @calledby   
   @history 
   Trimmed and altered from the GR3D routines, original author Mark Miller.
   2007?: Bruno escluded the points in the atmosphere and excision region from the computation.
   Aug. 2008: Luca added a check on whether a failure at a given point may be disregarded, because that point will then be restricted from a finer level.
   This should be completely safe only if *regridding happens at times when all levels are evolved.*
   @endhistory 

@@*/

subroutine Conservative2Primitive(CCTK_ARGUMENTS)
  
  use Recovery_Scalars

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
#ifdef _EOS_BASE_INC_
#undef _EOS_BASE_INC_
#endif
#include "EOS_Base.inc"

  integer :: i, j, k, itracer, nx, ny, nz
  CCTK_REAL :: uxx, uxy, uxz, uyy, uyz, uzz, det, psi4pt, pmin, enthalpy
  logical :: epsnegative
  character(len=400) :: warnline
  
  CCTK_INT :: type_bits, atmosphere
  CCTK_INT :: type2_bits, excised

  CCTK_REAL :: local_min_tracer


  call SpaceMask_GetTypeBits(type_bits, "Hydro_Atmosphere")
  call SpaceMask_GetStateBits(atmosphere, "Hydro_Atmosphere", "in_atmosphere")
    type2_bits = -1
    excised = -1

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)
  
  if (use_min_tracer .ne. 0) then
    local_min_tracer = min_tracer
  else
    local_min_tracer = 0d0
  end if
           
  pmin=EOS_Pressure(whisky_polytrope_handle, whisky_rho_min, 1.0d0)
  puncture_utau_min = whisky_rho_min * EOS_SpecificIntEnergy(whisky_polytrope_handle, whisky_rho_min, 1.0d0)
!!$  do k = whisky_stencil + 1, nz - whisky_stencil
!!$    do j = whisky_stencil + 1, ny - whisky_stencil
!!$      do i = whisky_stencil + 1, nx - whisky_stencil
  !$OMP PARALLEL DO PRIVATE (i, j, k, itracer, uxx, uxy, uxz, uyy, uyz, &
  !$OMP uzz, det, psi4pt, enthalpy, epsnegative, warnline)
  do k = 1, nz 
    do j = 1, ny 
      do i = 1, nx

         !do not compute if in atmosphere
         !or in an excised region
         !if (SpaceMask_CheckStateBitsF90(space_mask, i, j, k, type_bits, atmosphere)) cycle
         !do compute, to ensure consitent primitives (the metric might change
         !after all)
         
        epsnegative = .false.

        call SpatialDeterminant(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),&
             gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),det)
        call UpperMetric(uxx,uxy,uxz,uyy,uyz,uzz,det,&
             gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),&
             gyz(i,j,k),gzz(i,j,k))        
        if ((det .ne. det .or. det .lt. 0d0) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
            !$OMP CRITICAL
            write(warnline,'("Conservative2Primitive: Metric determinant in iteration ",i6," at ",g15.6,",",g15.6,",",&
              &g15.6," on level ",i2," is : ",g15.6, " from data g = [",&
              &5(g15.6,","),g15.6,"]")') cctk_iteration,&
              x(i,j,k),y(i,j,k),z(i,j,k),whisky_reflevel,det,&
              gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k)
            call CCTK_WARN(1, warnline)
            !$OMP END CRITICAL
        end if

        if (evolve_tracer .ne. 0) then
           do itracer=1,number_of_tracers
              call Con2Prim_ptTracer(cons_tracer(i,j,k,itracer), tracer(i,j,k,itracer), &
                   dens(i,j,k))

              if (use_min_tracer .ne. 0) then
                if (tracer(i,j,k,itracer) .le. local_min_tracer) then
                  tracer(i,j,k,itracer) = local_min_tracer
                end if
              end if

           enddo
           
         endif

        call Con2Prim_pt(whisky_eos_handle, dens(i,j,k),sx(i,j,k),sy(i,j,k), &
             sz(i,j,k),tau(i,j,k),rho(i,j,k),velx(i,j,k),vely(i,j,k), &
             velz(i,j,k),eps(i,j,k),press(i,j,k),w_lorentz(i,j,k), &
             uxx,uxy,uxz,uyy,uyz,uzz,det,x(i,j,k),y(i,j,k), &
             z(i,j,k),r(i,j,k),epsnegative,whisky_rho_min,pmin, whisky_reflevel, &
             Whisky_CarpetWeights(i,j,k), cctk_iteration)

        if (epsnegative) then
          if (store_recovery_fails.ne.0) then
             failed_con2prim_mask(i,j,k) = failed_con2prim_mask(i,j,k) + 1 
          end if 
          call Con2Prim_ptPolytype(whisky_polytrope_handle, &
               dens(i,j,k), sx(i,j,k), sy(i,j,k), sz(i,j,k), &
               tau(i,j,k), rho(i,j,k), velx(i,j,k), vely(i,j,k), &
               velz(i,j,k), eps(i,j,k), press(i,j,k), w_lorentz(i,j,k), &
               uxx, uxy, uxz, uyy, uyz, uzz, det, x(i,j,k), y(i,j,k), &
               z(i,j,k), r(i,j,k),whisky_rho_min, whisky_reflevel, &
               Whisky_CarpetWeights(i,j,k), cctk_iteration)
          if ((tau(i,j,k).le.0 .or. tau(i,j,k).ne.tau(i,j,k)) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
            !$OMP CRITICAL
            write(warnline,'("Conservative2Primitive: Unphysical tau = ",g15.6," found for data rho = ",&
              &g15.6,", eps = ",g15.6,", press = ",g15.6,", w = ",g15.6,&
              &" occured in iteration ",i6," at location [",&
              &2(g15.6,","),g15.6,"] on reflevel ",i2," weight ",f4.2)') &
              tau(i,j,k), rho(i,j,k),eps(i,j,k),press(i,j,k),&
              w_lorentz(i,j,k),&
              cctk_iteration, x(i,j,k),y(i,j,k),z(i,j,k), whisky_reflevel, &
              Whisky_CarpetWeights(i,j,k)
            call CCTK_WARN(1,warnline)
            !$OMP END CRITICAL
          end if
        end if
        if ((w_lorentz(i, j, k).lt.1.d0 .or. w_lorentz(i,j,k).ne.w_lorentz(i,j,k)) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
          !$OMP CRITICAL
          write(warnline,'("Conservative2Primitive: Unphysical Lorentz factor 1+(",g15.6,&
            &") for data g = [",5(g15.6,","),g15.6,"] vel = [",&
            &2(g15.6,","),g15.6,"] occured in iteration ",i6," at location [",&
            &2(g15.6,","),g15.6,"] on reflevel ",i2," weight ",f4.2)') &
            w_lorentz(i,j,k)-1.d0, gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
            velx(i,j,k),vely(i,j,k),velz(i,j,k),&
            cctk_iteration, x(i,j,k),y(i,j,k),z(i,j,k), whisky_reflevel, &
            Whisky_CarpetWeights(i,j,k)
          call CCTK_WARN(1,warnline)
          !$OMP END CRITICAL
        end if 

        if (eps(i,j,k) .lt. 0.0d0) then 
           !$OMP CRITICAL
           if (whisky_reflevel.ge.warn_from_reflevel .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
             call CCTK_WARN(1,'Con2Prim: stopping the code.')
             call CCTK_WARN(1, '   specific internal energy just went below 0! ')
             write(warnline,'(a28,i2)') 'on carpet reflevel: ',whisky_reflevel
             call CCTK_WARN(1,warnline)
             write(warnline,'(a20,3g16.7)') 'xyz location: ',&
                  x(i,j,k),y(i,j,k),z(i,j,k)
             call CCTK_WARN(1,warnline)
             write(warnline,'(a20,g16.7)') 'radius: ',r(i,j,k)
             call CCTK_WARN(1,warnline)
             write(warnline,'(a20,3g16.7)') 'velocities: ',&
                  velx(i,j,k),vely(i,j,k),velz(i,j,k)
             call CCTK_WARN(1,warnline)
             call CCTK_WARN(warnlevel, "Specific internal energy negative")
           else
             write(warnline,'(a60,i2)') 'Con2Prim: eps negative, but I was told to ignore level: ',whisky_reflevel
             call CCTK_WARN(1,warnline)
           endif
           !$OMP END CRITICAL
       endif
      end do
    end do
  end do
  !$OMP END PARALLEL DO

end subroutine Conservative2Primitive

 /*@@
   @routine    Conservative2PrimitiveBoundaries
   @date       Tue Mar 12 18:04:40 2002
   @author     The Whisky Developers    
   @desc 
        This routine is used only if the reconstruction is performed on the conserved variables. It computes the primitive variables on cell boundaries.
        Since reconstruction on conservative had not proved to be very successful, some of the improvements to the C2P routines (e.g. the check about 
        whether a failure happens in a point that will be restriced anyway) are not implemented here yet.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/


subroutine Conservative2PrimitiveBounds(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  integer :: i, j, k, itracer, nx, ny, nz
  CCTK_REAL :: uxxl, uxyl, uxzl, uyyl, uyzl, uzzl,&
       uxxr, uxyr, uxzr, uyyr, uyzr, uzzr
  CCTK_REAL :: gxxl,gxyl,gxzl,gyyl,gyzl,gzzl,avg_detl,&
       gxxr,gxyr,gxzr,gyyr,gyzr,gzzr,avg_detr,psi4r,psi4l
  logical :: epsnegative
  character(len=100) warnline
 
  CCTK_INT :: type_bits, atmosphere
  CCTK_INT :: type2_bits, excised

  CCTK_REAL :: local_min_tracer
                                                                               
  call SpaceMask_GetTypeBits(type_bits, "Hydro_Atmosphere")
  call SpaceMask_GetStateBits(atmosphere, "Hydro_Atmosphere", "in_atmosphere")
    type2_bits = -1
    excised = -1

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3) 
  
  if (use_min_tracer .ne. 0) then
    local_min_tracer = min_tracer
  else
    local_min_tracer = 0d0
  end if
  
  !$OMP PARALLEL DO PRIVATE (i, j, k, itracer, uxxl, uxyl, uxzl, uyyl, &
  !$OMP uyzl, uzzl, uxxr, uxyr, uxzr, uyyr, uyzr, uzzr, gxxl,gxyl,gxzl,gyyl, &
  !$OMP gyzl,gzzl,avg_detl, gxxr,gxyr,gxzr,gyyr,gyzr,gzzr,avg_detr,psi4r,psi4l, &
  !$OMP epsnegative, warnline)
  do k = whisky_stencil, nz - whisky_stencil + 1
    do j = whisky_stencil, ny - whisky_stencil + 1
      do i = whisky_stencil, nx - whisky_stencil + 1

        !do not compute if in atmosphere
        if (SpaceMask_CheckStateBitsF90(space_mask, i, j, k, type_bits, atmosphere)) cycle
         
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

        epsnegative = .false.

          call SpatialDeterminant(gxxl,gxyl,gxzl,gyyl, gyzl,gzzl,avg_detl)
          call SpatialDeterminant(gxxr,gxyr,gxzr,gyyr, gyzr,gzzr,avg_detr)
          call UpperMetric(uxxl,uxyl,uxzl,uyyl,uyzl,uzzl,avg_detl,&
               gxxl,gxyl,gxzl,gyyl,gyzl,gzzl)        
          call UpperMetric(uxxr,uxyr,uxzr,uyyr,uyzr,uzzr,avg_detr,&
               gxxr,gxyr,gxzr,gyyr,gyzr,gzzr)        

        if (evolve_tracer .ne. 0) then
           do itracer=1,number_of_tracers
              call Con2Prim_ptTracer(cons_tracer(i,j,k,itracer), &
                   tracer(i,j,k,itracer), dens(i,j,k))
           enddo

           if (use_min_tracer .ne. 0) then
             if (tracer(i,j,k,itracer) .le. local_min_tracer) then
               tracer(i,j,k,itracer) = local_min_tracer
             end if
           end if

        endif

        call Con2Prim_pt(whisky_eos_handle, densminus(i,j,k),&
             sxminus(i,j,k),syminus(i,j,k),szminus(i,j,k),&
             tauminus(i,j,k),rhominus(i,j,k),velxminus(i,j,k),&
             velyminus(i,j,k),velzminus(i,j,k),epsminus(i,j,k),&
             pressminus(i,j,k),w_lorentzminus(i,j,k),&
             uxxl,uxyl,uxzl,uyyl,uyzl,uzzl,avg_detl,&
             x(i,j,k)-0.5d0*CCTK_DELTA_SPACE(1),&
             y(i,j,k)-0.5d0*CCTK_DELTA_SPACE(2), &
             z(i,j,k)-0.5d0*CCTK_DELTA_SPACE(3),r(i,j,k),&
             epsnegative,whisky_rho_min, whisky_reflevel,&
             Whisky_CarpetWeights(i,j,k), cctk_iteration)
        if (epsnegative) then
          if (store_recovery_fails.ne.0) then
             failed_con2prim_mask(i,j,k) = failed_con2prim_mask(i,j,k)+10
          end if 
          call Con2Prim_ptPolytype(whisky_polytrope_handle, densminus(i,j,k),&
               sxminus(i,j,k),syminus(i,j,k),szminus(i,j,k),&
               tauminus(i,j,k),rhominus(i,j,k),velxminus(i,j,k),&
               velyminus(i,j,k),velzminus(i,j,k),epsminus(i,j,k),&
               pressminus(i,j,k),w_lorentzminus(i,j,k),&
               uxxl,uxyl,uxzl,uyyl,uyzl,uzzl,avg_detl,&
               x(i,j,k)-0.5d0*CCTK_DELTA_SPACE(1),&
               y(i,j,k)-0.5d0*CCTK_DELTA_SPACE(2), &
               z(i,j,k)-0.5d0*CCTK_DELTA_SPACE(3),r(i,j,k),whisky_rho_min, whisky_reflevel, &
               Whisky_CarpetWeights(i,j,k), cctk_iteration)
        end if

        if (epsminus(i,j,k) .lt. 0.0d0) then
          if (whisky_reflevel.ge.warn_from_reflevel .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
            !$OMP CRITICAL
            call CCTK_WARN(1,'Con2Prim: stopping the code.')
            call CCTK_WARN(1, '   specific internal energy just went below 0! ')
            write(warnline,'(a28,i2)') 'on carpet reflevel: ',whisky_reflevel
            call CCTK_WARN(1,warnline)
            write(warnline,'(a20,3g16.7)') 'xyz location: ',&
                 x(i,j,k),y(i,j,k),z(i,j,k)
            call CCTK_WARN(1,warnline)
            write(warnline,'(a20,g16.7)') 'radius: ',r(i,j,k)
            call CCTK_WARN(1,warnline)
            write(warnline,'(a20,3g16.7)') 'velocities: ',&
                 velxminus(i,j,k),velyminus(i,j,k),velzminus(i,j,k)
            call CCTK_WARN(1,warnline)
            call CCTK_WARN(warnlevel, "Specific internal energy negative")
            !$OMP END CRITICAL
            cycle
          endif
        endif

        epsnegative = .false.
        call Con2Prim_pt(whisky_eos_handle, densplus(i,j,k),&
             sxplus(i,j,k),syplus(i,j,k),szplus(i,j,k),&
             tauplus(i,j,k),rhoplus(i,j,k),velxplus(i,j,k),&
             velyplus(i,j,k),velzplus(i,j,k),epsplus(i,j,k),&
             pressplus(i,j,k),w_lorentzplus(i,j,k),&
             uxxr,uxyr,uxzr,uyyr,uyzr,uzzr,avg_detr,&
             x(i,j,k)+0.5d0*CCTK_DELTA_SPACE(1),&
             y(i,j,k)+0.5d0*CCTK_DELTA_SPACE(2), &
             z(i,j,k)+0.5d0*CCTK_DELTA_SPACE(3),r(i,j,k),&
             epsnegative,whisky_rho_min, whisky_reflevel,&
             Whisky_CarpetWeights(i,j,k), cctk_iteration)
        if (epsnegative) then
          if (store_recovery_fails.ne.0) then
             failed_con2prim_mask(i,j,k) = failed_con2prim_mask(i,j,k)+100
          end if 
          call Con2Prim_ptPolytype(whisky_polytrope_handle, densplus(i,j,k),&
               sxplus(i,j,k),syplus(i,j,k),szplus(i,j,k),&
               tauplus(i,j,k),rhoplus(i,j,k),velxplus(i,j,k),&
               velyplus(i,j,k),velzplus(i,j,k),epsplus(i,j,k),&
               pressplus(i,j,k),w_lorentzplus(i,j,k),&
               uxxr,uxyr,uxzr,uyyr,uyzr,uzzr,avg_detr,&
               x(i,j,k)+0.5d0*CCTK_DELTA_SPACE(1),&
               y(i,j,k)+0.5d0*CCTK_DELTA_SPACE(2), &
               z(i,j,k)+0.5d0*CCTK_DELTA_SPACE(3),r(i,j,k),whisky_rho_min, whisky_reflevel,&
               Whisky_CarpetWeights(i,j,k), cctk_iteration)
        end if

        if (epsplus(i,j,k) .lt. 0.0d0) then
          !$OMP CRITICAL
          if (whisky_reflevel.ge.warn_from_reflevel .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
            call CCTK_WARN(1,'Con2Prim: stopping the code.')
            call CCTK_WARN(1, '   specific internal energy just went below 0! ')
            write(warnline,'(a28,i2)') 'on carpet reflevel: ',whisky_reflevel
            call CCTK_WARN(1,warnline)
            write(warnline,'(a20,3g16.7)') 'xyz location: ',&
                 x(i,j,k),y(i,j,k),z(i,j,k)
            call CCTK_WARN(1,warnline)
            write(warnline,'(a20,g16.7)') 'radius: ',r(i,j,k)
            call CCTK_WARN(1,warnline)
            write(warnline,'(a20,3g16.7)') 'velocities: ',&
                 velxplus(i,j,k),velyplus(i,j,k),velzplus(i,j,k)
            call CCTK_WARN(1,warnline)
            call CCTK_WARN(warnlevel, "Specific internal energy negative")
          else
            write(warnline,'(a60,i2)') 'Con2Prim: eps negative, but I was told to ignore level: ',whisky_reflevel
            call CCTK_WARN(1,warnline)
          endif
          !$OMP END CRITICAL
        endif

     end do
    end do
  end do
  !$OMP END PARALLEL DO
  
end subroutine Conservative2PrimitiveBounds


 /*@@
   @routine    Con2PrimPolytype
   @date       Tue Mar 19 11:43:06 2002
   @author     Ian Hawke
   @desc 
   All routines below are identical to those above, just
   specialised from polytropic type EOS.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Conservative2PrimitivePolytype(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: i, j, k, itracer, nx, ny, nz
  CCTK_REAL :: uxx, uxy, uxz, uyy, uyz, uzz, det, psi4pt

  CCTK_INT :: type_bits, atmosphere
  CCTK_INT :: type2_bits, excised

  CCTK_REAL :: local_min_tracer

  character(len=400) :: warnline
                                                                               
  call SpaceMask_GetTypeBits(type_bits, "Hydro_Atmosphere")
  call SpaceMask_GetStateBits(atmosphere, "Hydro_Atmosphere", "in_atmosphere")
    type2_bits = -1
    excised = -1
 
  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)
  
  if (use_min_tracer .ne. 0) then
    local_min_tracer = min_tracer
  else
    local_min_tracer = 0d0
  end if
  
!!$  do k = whisky_stencil + 1, nz - whisky_stencil
!!$    do j = whisky_stencil + 1, ny - whisky_stencil
!!$      do i = whisky_stencil + 1, nx - whisky_stencil
  !$OMP PARALLEL DO PRIVATE (i, j, k, itracer, uxx, uxy, uxz, uyy, uyz, &
  !$OMP uzz, det, psi4pt, warnline)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx

         !do not compute if in atmosphere
         !or in an excised region
         !if (SpaceMask_CheckStateBitsF90(space_mask, i, j, k, type_bits, atmosphere)) cycle
         !do compute, to ensure consitent primitives (the metric might change
         !after all)

        call SpatialDeterminant(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),&
             gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),det)
        call UpperMetric(uxx,uxy,uxz,uyy,uyz,uzz,det,&
             gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),&
             gyz(i,j,k),gzz(i,j,k))        
        if ((det .ne. det .or. det .lt. 0d0) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
            !$OMP CRITICAL
            write(warnline,'("Conservative2PrimitivePolytype: Metric determinant in iteration ",i6," at ",g15.6,",",g15.6,",",&
              &g15.6," on level ",i2," is : ",g15.6, " from data g = [",&
              &5(g15.6,","),g15.6,"]")') cctk_iteration,&
              x(i,j,k),y(i,j,k),z(i,j,k),whisky_reflevel,det,&
              gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k)
            call CCTK_WARN(1, warnline)
            !$OMP END CRITICAL
        end if

        if (evolve_tracer .ne. 0) then
           do itracer=1,number_of_tracers
              call Con2Prim_ptTracer(cons_tracer(i,j,k,itracer), & 
                   tracer(i,j,k,itracer), dens(i,j,k))
           enddo

           if (use_min_tracer .ne. 0) then
             if (tracer(i,j,k,itracer) .le. local_min_tracer) then
               tracer(i,j,k,itracer) = local_min_tracer
             end if
           end if

        endif

        call Con2Prim_ptPolytype(whisky_eos_handle, dens(i,j,k),&
             sx(i,j,k),sy(i,j,k), &
             sz(i,j,k),tau(i,j,k),rho(i,j,k),velx(i,j,k),vely(i,j,k), &
             velz(i,j,k),eps(i,j,k),press(i,j,k),w_lorentz(i,j,k), &
             uxx,uxy,uxz,uyy,uyz,uzz,det,x(i,j,k),y(i,j,k), &
             z(i,j,k),r(i,j,k),whisky_rho_min, whisky_reflevel, &
             Whisky_CarpetWeights(i,j,k), cctk_iteration)
        if ((tau(i,j,k).le.0 .or. tau(i,j,k).ne.tau(i,j,k)) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
          !$OMP CRITICAL
          write(warnline,'("Conservative2PrimitivePolytype: Unphysical tau = ",g15.6," found for data rho = ",&
            &g15.6,", eps = ",g15.6,", press = ",g15.6,", w = ",g15.6,&
            &" occured in iteration ",i6," at location [",&
            &2(g15.6,","),g15.6,"] on reflevel ",i2," weight ",f4.2)') &
            tau(i,j,k), rho(i,j,k),eps(i,j,k),press(i,j,k),&
            w_lorentz(i,j,k),&
            cctk_iteration, x(i,j,k),y(i,j,k),z(i,j,k), whisky_reflevel, &
            Whisky_CarpetWeights(i,j,k)
          call CCTK_WARN(1,warnline)
          !$OMP END CRITICAL
        end if
        if ((w_lorentz(i, j, k).lt.1.d0 .or. w_lorentz(i,j,k).ne.w_lorentz(i,j,k)) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
          !$OMP CRITICAL
          write(warnline,'("Conservative2PrimitivePolytype: Unphysical Lorentz factor 1+(",g15.6,&
            &") for data g = [",5(g15.6,","),g15.6,"] vel = [",&
            &2(g15.6,","),g15.6,"] occured in iteration ",i6," at location [",&
            &2(g15.6,","),g15.6,"] on reflevel ",i2," weight ",f4.2)') &
            w_lorentz(i,j,k)-1.d0, gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
            velx(i,j,k),vely(i,j,k),velz(i,j,k),&
            cctk_iteration, x(i,j,k),y(i,j,k),z(i,j,k), whisky_reflevel, &
            Whisky_CarpetWeights(i,j,k)
          call CCTK_WARN(1,warnline)
          !$OMP END CRITICAL
        end if 
        if ( store_recovery_fails.ne.0 .and. ( tau(i,j,k).le.0 .or. tau(i,j,k).ne.tau(i,j,k) .or. w_lorentz(i,j,k).lt.1.d0 .or. w_lorentz(i,j,k).ne.w_lorentz(i,j,k) ) ) then
             failed_con2prim_mask(i,j,k) = failed_con2prim_mask(i,j,k) + 1 
        end if 
        
     end do
    end do
  end do
  !$OMP END PARALLEL DO
  
end subroutine Conservative2PrimitivePolytype


 /*@@
   @routine    Cons2PrimBoundsPolytype
   @date       Tue Mar 12 18:04:40 2002
   @author     The Whisky Developers
   @desc 
        This routine is used only if the reconstruction is performed on the conserved variables. It computes the primitive variables on cell boundaries.
        Since reconstruction on conservative had not proved to be very successful, some of the improvements to the C2P routines (e.g. the check about 
        whether a failure happens in a point that will be restriced anyway) are not implemented here yet.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Con2PrimBoundsPolytype(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  integer :: i, j, k, nx, ny, nz
  CCTK_REAL :: uxxl, uxyl, uxzl, uyyl, uyzl, uzzl,&
       uxxr, uxyr, uxzr, uyyr, uyzr, uzzr
  CCTK_REAL :: gxxl,gxyl,gxzl,gyyl,gyzl,gzzl,avg_detl,&
       gxxr,gxyr,gxzr,gyyr,gyzr,gzzr,avg_detr,psi4r,psi4l

  CCTK_INT :: type_bits, atmosphere
  CCTK_INT :: type2_bits, excised
  character(len=400) :: warnline
                                                                               
  call SpaceMask_GetTypeBits(type_bits, "Hydro_Atmosphere")
  call SpaceMask_GetStateBits(atmosphere, "Hydro_Atmosphere", "in_atmosphere")
    type2_bits = -1
    excised = -1

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  !$OMP PARALLEL DO PRIVATE (i, j, k, uxxl, uxyl, uxzl, uyyl, uyzl, uzzl, &
  !$OMP uxxr, uxyr, uxzr, uyyr, uyzr, uzzr, gxxl,gxyl,gxzl,gyyl,gyzl,gzzl, &
  !$OMP avg_detl, gxxr,gxyr,gxzr,gyyr,gyzr,gzzr,avg_detr,psi4r,psi4l, warnline)
  do k = whisky_stencil, nz - whisky_stencil + 1
    do j = whisky_stencil, ny - whisky_stencil + 1
      do i = whisky_stencil, nx - whisky_stencil + 1
        
        !do not compute if in atmosphere
         if (SpaceMask_CheckStateBitsF90(space_mask, i, j, k, type_bits, atmosphere)) cycle
        
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

        call SpatialDeterminant(gxxl,gxyl,gxzl,gyyl, gyzl,gzzl,avg_detl)
        call SpatialDeterminant(gxxr,gxyr,gxzr,gyyr, gyzr,gzzr,avg_detr)
        call UpperMetric(uxxl,uxyl,uxzl,uyyl,uyzl,uzzl,avg_detl,&
             gxxl,gxyl,gxzl,gyyl,gyzl,gzzl)        
        call UpperMetric(uxxr,uxyr,uxzr,uyyr,uyzr,uzzr,avg_detr,&
             gxxr,gxyr,gxzr,gyyr,gyzr,gzzr)        
        if (Whisky_CarpetWeights(i,j,k) .eq.1.d0 .and. (avg_detl .ne. avg_detl .or. avg_detl .lt. 0d0)) then
            !$OMP CRITICAL
            write(warnline,'("Metric determinant (l) in iteration ",i6," at [",g15.6,",",g15.6,",",&
              &g15.6,"] on level ",i2," is : ",g15.6)') cctk_iteration,&
              x(i,j,k)-0.5d0*CCTK_DELTA_SPACE(1),&
              y(i,j,k)-0.5d0*CCTK_DELTA_SPACE(2),z(i,j,k)-0.5d0*CCTK_DELTA_SPACE(3),&
              whisky_reflevel,avg_detl
            call CCTK_WARN(1, warnline)
            !$OMP END CRITICAL
        end if
        if (Whisky_CarpetWeights(i,j,k) .eq.1.d0 .and. (avg_detr .ne. avg_detr .or. avg_detr .lt. 0d0)) then
            !$OMP CRITICAL
            write(warnline,'("Metric determinant (r) in iteration ",i6," at [",g15.6,",",g15.6,",",&
              &g15.6,"] on level ",i2," is : ",g15.6)') cctk_iteration,&
              x(i,j,k)+0.5d0*CCTK_DELTA_SPACE(1),&
              y(i,j,k)+0.5d0*CCTK_DELTA_SPACE(2),z(i,j,k)+0.5d0*CCTK_DELTA_SPACE(3),&
              whisky_reflevel,avg_detr
            call CCTK_WARN(1, warnline)
            !$OMP END CRITICAL
        end if

        call Con2Prim_ptPolytype(whisky_eos_handle, densminus(i,j,k),&
             sxminus(i,j,k),syminus(i,j,k),szminus(i,j,k),&
             tauminus(i,j,k),rhominus(i,j,k),velxminus(i,j,k),&
             velyminus(i,j,k),velzminus(i,j,k),epsminus(i,j,k),&
             pressminus(i,j,k),w_lorentzminus(i,j,k),&
             uxxl,uxyl,uxzl,uyyl,uyzl,uzzl,avg_detl,&
             x(i,j,k)-0.5d0*CCTK_DELTA_SPACE(1),&
             y(i,j,k)-0.5d0*CCTK_DELTA_SPACE(2), &
             z(i,j,k)-0.5d0*CCTK_DELTA_SPACE(3),r(i,j,k),whisky_rho_min, &
             whisky_reflevel, Whisky_CarpetWeights(i,j,k), &
             cctk_iteration)
        call Con2Prim_ptPolytype(whisky_eos_handle, densplus(i,j,k),&
             sxplus(i,j,k),syplus(i,j,k),szplus(i,j,k),&
             tauplus(i,j,k),rhoplus(i,j,k),velxplus(i,j,k),&
             velyplus(i,j,k),velzplus(i,j,k),epsplus(i,j,k),&
             pressplus(i,j,k),w_lorentzplus(i,j,k),&
             uxxr,uxyr,uxzr,uyyr,uyzr,uzzr,avg_detr,&
             x(i,j,k)+0.5d0*CCTK_DELTA_SPACE(1),&
             y(i,j,k)+0.5d0*CCTK_DELTA_SPACE(2), &
             z(i,j,k)+0.5d0*CCTK_DELTA_SPACE(3),r(i,j,k),whisky_rho_min, &
             whisky_reflevel, Whisky_CarpetWeights(i,j,k),&
             cctk_iteration)
        
     end do
    end do
  end do
  !$OMP END PARALLEL DO
  
end subroutine Con2PrimBoundsPolytype


 /*@@
   @routine    Con2PrimBoundsTracer
   @date       Mon Mar  8 13:41:55 2004
   @author     Ian Hawke
   @desc 
   
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Con2PrimBoundsTracer(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  integer :: i, j, k, itracer, nx, ny, nz
  CCTK_REAL :: uxxl, uxyl, uxzl, uyyl, uyzl, uzzl,&
       uxxr, uxyr, uxzr, uyyr, uyzr, uzzr
  CCTK_REAL :: gxxl,gxyl,gxzl,gyyl,gyzl,gzzl,avg_detl,&
       gxxr,gxyr,gxzr,gyyr,gyzr,gzzr,avg_detr,psi4r,psi4l
 
  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)
  
  !$OMP PARALLEL DO PRIVATE (i, j, k, itracer, uxxl, uxyl, uxzl, uyyl, uyzl, uzzl,&
  !$OMP uxxr, uxyr, uxzr, uyyr, uyzr, uzzr, gxxl,gxyl,gxzl,gyyl,gyzl,gzzl,avg_detl,&
  !$OMP gxxr,gxyr,gxzr,gyyr,gyzr,gzzr,avg_detr,psi4r,psi4l)
  do k = whisky_stencil, nz - whisky_stencil + 1
    do j = whisky_stencil, ny - whisky_stencil + 1
      do i = whisky_stencil, nx - whisky_stencil + 1
        
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

        call SpatialDeterminant(gxxl,gxyl,gxzl,gyyl, gyzl,gzzl,avg_detl)
        call SpatialDeterminant(gxxr,gxyr,gxzr,gyyr, gyzr,gzzr,avg_detr)
        call UpperMetric(uxxl,uxyl,uxzl,uyyl,uyzl,uzzl,avg_detl,&
             gxxl,gxyl,gxzl,gyyl,gyzl,gzzl)        
        call UpperMetric(uxxr,uxyr,uxzr,uyyr,uyzr,uzzr,avg_detr,&
             gxxr,gxyr,gxzr,gyyr,gyzr,gzzr)        

        do itracer=1,number_of_tracers
           call Con2Prim_ptBoundsTracer(cons_tracerplus(i,j,k,itracer), &
                tracerplus(i,j,k,itracer), &
                rhoplus(i,j,k), velxplus(i,j,k), velyplus(i,j,k), &
                velzplus(i,j,k), 1.d0/sqrt(1.d0 - &
                (gxxr * velxplus(i,j,k)**2 + &
                gyyr * velyplus(i,j,k)**2 + &
                gzzr * velzplus(i,j,k)**2 + &
                2.d0 * (gxyr * velxplus(i,j,k) * velyplus(i,j,k) + &
                gxzr * velxplus(i,j,k) * velzplus(i,j,k) + &
                gyzr * velyplus(i,j,k) * velzplus(i,j,k) ) ) ), &
                avg_detr)
           call Con2Prim_ptBoundsTracer(cons_tracerminus(i,j,k,itracer), &
                tracerminus(i,j,k,itracer), &
                rhominus(i,j,k), velxminus(i,j,k), velyminus(i,j,k), &
                velzminus(i,j,k), 1.d0/sqrt(1.d0 - &
                (gxxr * velxminus(i,j,k)**2 + &
                gyyr * velyminus(i,j,k)**2 + &
                gzzr * velzminus(i,j,k)**2 + &
                2.d0 * (gxyr * velxminus(i,j,k) * velyminus(i,j,k) + &
                gxzr * velxminus(i,j,k) * velzminus(i,j,k) + &
                gyzr * velyminus(i,j,k) * velzminus(i,j,k) ) ) ), &
                avg_detl)
        enddo

!!$        tracerplus(i,j,k) = cons_tracerplus(i,j,k) / &
!!$             sqrt(avg_detr) / rhoplus(i,j,k) * &
!!$             sqrt(1.d0 - &
!!$                   (gxxr * velxplus(i,j,k)**2 + &
!!$                    gyyr * velyplus(i,j,k)**2 + &
!!$                    gzzr * velzplus(i,j,k)**2 + &
!!$                    2.d0 * (gxyr * velxplus(i,j,k) * velyplus(i,j,k) + &
!!$                            gxzr * velxplus(i,j,k) * velzplus(i,j,k) + &
!!$                            gyzr * velyplus(i,j,k) * velzplus(i,j,k) ) ) )
!!$        tracerminus(i,j,k) = cons_tracerminus(i,j,k) / &
!!$             sqrt(avg_detl) / rhominus(i,j,k) * &
!!$             sqrt(1.d0 - &
!!$                   (gxxl * velxminus(i,j,k)**2 + &
!!$                    gyyl * velyminus(i,j,k)**2 + &
!!$                    gzzl * velzminus(i,j,k)**2 + &
!!$                    2.d0 * (gxyl * velxminus(i,j,k) * velyminus(i,j,k) + &
!!$                            gxzl * velxminus(i,j,k) * velzminus(i,j,k) + &
!!$                            gyzl * velyminus(i,j,k) * velzminus(i,j,k) ) ) )
        
      end do
    end do
  end do
  !$OMP END PARALLEL DO
  
end subroutine Con2PrimBoundsTracer

 /*@@
   @routine    RecoverMarti_setCenterPars
   @date       Fri Dec 12 16:52:15 2008 
   @author     Tanja Bode
   @desc 
     Set the center and radii of where to impose the Faber limit.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine RecoverMarti_setCenterVars(CCTK_ARGUMENTS)

   use Recovery_Scalars
   implicit none

   DECLARE_CCTK_ARGUMENTS
   DECLARE_CCTK_PARAMETERS
   DECLARE_CCTK_FUNCTIONS

   integer :: pct, s_ind, sr_ind
   logical :: get_radius, get_center, skip
   character(len=200) warnline

   do pct = 1,number_of_punctures

     radii(pct) = -1d0;
     center_x(pct) = 0d0;
     center_y(pct) = 0d0;
     center_z(pct) = 0d0;

     get_center = CCTK_EQUALS(puncture_center_from(pct),"spherical surface")
     get_radius = CCTK_EQUALS(puncture_radius_from(pct),"spherical surface")

     !! Parameters only.
     if ( .not.get_center ) then
        center_x(pct) = puncture_center_x(pct)
        center_y(pct) = puncture_center_y(pct)
        center_z(pct) = puncture_center_z(pct)
     end if
     if ( .not.get_radius ) then
        radii(pct) = puncture_radius(pct)
     end if

     if ( get_center .or. get_radius ) then

        s_ind = 0
        sr_ind = 0
        skip = .FALSE.

        !! Handle any spherical surface indicies, default radius surface is center surface
        if ( get_center ) then
           s_ind = puncture_surface(pct)+1
        end if
        if ( get_radius .and. puncture_radius_surface(pct).lt.0 .and. get_center ) then
           sr_ind = s_ind
        else 
           sr_ind = puncture_radius_surface(pct)+1
        end if

        if ( get_center.and.(s_ind.lt.1 .or. sf_active(s_ind).lt.1) ) then
            write(warnline,'(a14,i2,a36,i2,a64)') "Surface index ", s_ind-1, " (center surface) for Faber puncture", pct, &
                 " limit is not a valid surface this timestep. Not applying limit."
            call CCTK_WARN(1,warnline)
            skip=.TRUE.
        else if ( sf_valid(s_ind).ge.0.d0 ) then 
            center_x(pct) = sf_centroid_x(s_ind)
            center_y(pct) = sf_centroid_y(s_ind)
            center_z(pct) = sf_centroid_z(s_ind)
        end if
            
        if ( (get_radius.and.sf_active(sr_ind).lt.1) ) then
            write(warnline,'(a14,i2,a36,i2,a64)') "Surface index ", sr_ind-1, " (radius surface) for Faber puncture", pct, &
                 " limit is not a valid surface this timestep. Not applying limit."
            call CCTK_WARN(1,warnline)
        else if ( sf_valid(s_ind).ge.0d0.and.(.not.skip) ) then
            radii(pct) =  puncture_radius_factor(pct)*sf_min_radius(sr_ind)
        end if
  
     end if

   end do

end subroutine RecoverMarti_setCenterVars

logical function check_interval(x_lower,x_upper,epsrel,epsabs)
  
  implicit none

  CCTK_REAL :: x_lower, x_upper, epsrel, epsabs
  CCTK_REAL :: abs_lower, abs_upper

  CCTK_REAL :: min_abs, tolerance;

  abs_upper = abs(x_upper)
  abs_lower = abs(x_lower)
  if ((x_lower > 0.0 .and. x_upper > 0.0) .or. (x_lower < 0.0 .and. x_upper < 0.0)) then
    min_abs = min(abs_lower, abs_upper)
  else
    min_abs = 0;
  end if

  tolerance = epsabs + epsrel * min_abs

  if (abs(x_upper - x_lower) <= tolerance) then
    check_interval = .true.
  else
    check_interval = .false.
  end if

end function check_interval

subroutine Check_Symmetry(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  integer :: i, j, k
  integer :: nx,ny,nz
  integer :: gx,gy,gz
  character(len=600) warnline
  logical,external :: check_interval

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  gx = cctk_nghostzones(1)
  gy = cctk_nghostzones(2)
  gz = cctk_nghostzones(3)

  if(cctk_iteration.eq.0) return

  if( (z(1,1,gz+1) .gt. CCTK_DELTA_SPACE(3)/2.d0 .or. z(1,1,gz+1) .lt. -CCTK_DELTA_SPACE(3)/2.d0)) then
    if (verbose .gt. 0) then
       write(warnline,'("z(gz+1) = ",g15.6," skipping.")') z(1,1,gz+1)
       call CCTK_WARN(1,warnline)
    end if
    return
  end if
  
  !$OMP PARALLEL DO PRIVATE (i,j,k,warnline)
  do k = 1, gz
    do j = gy+1, ny-gy
      do i = gx+1, nx-gx      

100     FORMAT ("Symmetry broken at index (",i6,",",i6,",",i6,") coordinates (",&
            g15.6,",",g15.6,",",g15.6,") for variable ",a," in iteration ",i6,g15.6," != ",g15.6,&
            " difference ",g15.6)
#define CHECK(X,Y,fact) \
        if(.not. check_interval(X(i,j,gz+1+k), fact*X(i,j,gz+1-k), zsymmetry_rel_err, zsymmetry_abs_err)             \
           .and. Whisky_CarpetWeights(i,j,gz+1+k).eq.1.d0) then                                                    &&\
          write(warnline,100) i,j,gz+1+k,x(i,j,gz+1+k),y(i,j,gz+1+k),z(i,j,gz+1+k),Y,cctk_iteration,X(i,j,gz+1+k),  &\
                              X(i,j,gz+1-k),X(i,j,gz+1+k)-(fact)*X(i,j,gz+1-k)                                     &&\
          call CCTK_WARN(1,warnline)                                                                               &&\
        end if
  
        CHECK(rho,"rho",1.0)
        CHECK(eps,"eps",1.0)
        CHECK(press,"press",1.0)
        CHECK(velx,"velx",1.0)
        CHECK(vely,"vely",1.0)
        CHECK(velz,"velz",-1.0)
        CHECK(w_lorentz,"lorentz",1.0)
        CHECK(dens,"dens",1.0)
        CHECK(tau,"tau",1.0)
        CHECK(sx,"sx",1.0)
        CHECK(sy,"sy",1.0)
        CHECK(sz,"sz",-1.0)
        CHECK(gxx,"gxx",1.0)
        CHECK(gxy,"gxy",1.0)
        CHECK(gxz,"gxz",-1.0)
        CHECK(gyy,"gyy",1.0)
        CHECK(gyz,"gyz",-1.0)
        CHECK(gzz,"gzz",1.0)
        CHECK(alp,"alp",1.0)
        CHECK(betax,"betax",1.0)
        CHECK(betay,"betay",1.0)
        CHECK(betaz,"betaz",-1.0)
      end do
    end do
  end do
  !$OMP END PARALLEL DO
end subroutine Check_Symmetry

 /*@@
   @routine    RecoverMarti_ResetFailMask
   @date       Fri Dec 12 16:52:15 2008 
   @author     Tanja Bode
   @desc 
     Reset the con2prim failure mask
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine RecoverMarti_ResetFailMask(CCTK_ARGUMENTS)

  implicit none
  DECLARE_CCTK_ARGUMENTS

  integer :: i, j, k
  integer :: nx,ny,nz

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  !$OMP PARALLEL DO PRIVATE (i,j,k)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
         failed_con2prim_mask(i,j,k)=0
      end do 
    end do 
  end do 
  !$OMP END PARALLEL DO

end subroutine RecoverMarti_ResetFailMask

 /*@@
   @routine    RecoverMarti_CheckFailMask
   @date       Fri Dec 12 16:52:15 2008 
   @author     Tanja Bode
   @desc 
     Check the con2prim failure mask
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine RecoverMarti_CheckFailMask(CCTK_ARGUMENTS)

  implicit none
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  integer :: i, j, k, c2p_fail_warnlevel
  integer :: nx,ny,nz
  character(len=600) :: warnline

  CCTK_INT :: type_bits, atmosphere
  CCTK_INT :: type2_bits, excised

  call SpaceMask_GetTypeBits(type_bits, "Hydro_Atmosphere")
  call SpaceMask_GetStateBits(atmosphere, "Hydro_Atmosphere", "in_atmosphere")
    type2_bits = -1
    excised = -1

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  if ( abort_on_fails.gt.0 ) then
     c2p_fail_warnlevel = 0
  else
     c2p_fail_warnlevel = 1
  end if

  !$OMP PARALLEL DO PRIVATE(i,j,k,warnline)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
         if ( failed_con2prim_mask(i,j,k).gt.0 .and. &
              .not.SpaceMask_CheckStateBitsF90(space_mask, i, j, k, type_bits, atmosphere) .and. &
              Whisky_CarpetWeights(i,j,k).gt.0.5 .and. alp(i,j,k).gt.0.25 .and. &
              rho(i,j,k) > 1.01*whisky_rho_min ) then

            !! Warn. Failure codes interpreted as follows:
            !!    1-9    : Failure on initial recovery attempt (accumulates over substeps)
            !!    10-20  : Failure on minus side of reconstructed variables 
            !!    100-110: Failure on plus side of reconstructed variables 

            !$OMP CRITICAL
            write(warnline,'("Con2Prim failure code ",g15.6," at ",g15.6,",",g15.6,",",&
              &g15.6," on level ",i2,", weight ",g15.6," from conservatives (dens,tau,sx,sy,sz) = (",&
              &4(g15.6,","),g15.6,") with relative rho/rho_atmo = ",g15.6)') failed_con2prim_mask(i,j,k),&
              x(i,j,k),y(i,j,k),z(i,j,k),whisky_reflevel,Whisky_CarpetWeights(i,j,k),&
              dens(i,j,k),tau(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
              rho(i,j,k)/whisky_rho_min
            
            call CCTK_WARN(c2p_fail_warnlevel, warnline)
            !$OMP END CRITICAL

        end if

      end do 
    end do 
  end do 
  !$OMP END PARALLEL DO

end subroutine RecoverMarti_CheckFailMask

