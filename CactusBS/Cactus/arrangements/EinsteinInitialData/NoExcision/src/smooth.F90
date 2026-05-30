! $Header$

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine NoExcision_Smooth (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL, parameter :: zero=0
  CCTK_REAL :: mask (cctk_lsh(1), cctk_lsh(2), cctk_lsh(3))
  CCTK_REAL :: cx, cy, cz, radx, rady, radz
  integer   :: n
  integer   :: iter
  integer   :: ierr
  
  do iter = 1, smoothing_iterations
     
     do n = 1, num_regions
        
        cx = centre_x(n)
        cy = centre_y(n)
        cz = centre_z(n)
        if (CCTK_EQUALS(region_shape(n), "sphere")) then
           radx = radius(n)
           rady = radius(n)
           radz = radius(n)
        else if (CCTK_EQUALS(region_shape(n), "ellipsoid")) then
           radx = radius_x(n)
           rady = radius_y(n)
           radz = radius_z(n)
        else
           call CCTK_WARN (0, "internal error")
        end if
        
        mask = 1 - sqrt (((x - cx) / radx)**2 + &
             &           ((y - cy) / rady)**2 + &
             &           ((z - cz) / radz)**2)
        where (mask <= 0)
           mask = 0             ! outside
        elsewhere (mask >= smoothing_zone_width(n))
           mask = 1             ! far inside
        elsewhere
           mask = mask / smoothing_zone_width(n) ! a bit inside
        end where
        
        if (overwrite_geometry(n) /= 0) then
           
           if (conformal_state >= 1) then
              call smooth (psi)
           end if
           if (conformal_state >= 2) then
              call smooth (psix)
              call smooth (psiy)
              call smooth (psiz)
           end if
           if (conformal_state >= 3) then
              call smooth (psixx)
              call smooth (psixy)
              call smooth (psixz)
              call smooth (psiyy)
              call smooth (psiyz)
              call smooth (psizz)
           end if
           
           call smooth (gxx)
           call smooth (gxy)
           call smooth (gxz)
           call smooth (gyy)
           call smooth (gyz)
           call smooth (gzz)
           call smooth (kxx)
           call smooth (kxy)
           call smooth (kxz)
           call smooth (kyy)
           call smooth (kyz)
           call smooth (kzz)
           
        end if
        
        if (overwrite_lapse(n) /= 0) then
           
           call smooth (alp)
           
        end if
        
        if (overwrite_shift(n) /= 0) then
           
           if (shift_state /= 0) then
              call smooth (betax)
              call smooth (betay)
              call smooth (betaz)
           end if
           
        end if
        
     end do
     
     if (overwrite_geometry(n) /= 0) then
        
        if (conformal_state >= 1) then
           call CCTK_SyncGroup (ierr, cctkGH, "StaticConformal::confac")
        end if
        if (conformal_state >= 2) then
           call CCTK_SyncGroup (ierr, cctkGH, "StaticConformal::confac_1derivs")
        end if
        if (conformal_state >= 3) then
           call CCTK_SyncGroup (ierr, cctkGH, "StaticConformal::confac_2derivs")
        end if
        
        call CCTK_SyncGroup (ierr, cctkGH, "ADMBase::metric")
        call CCTK_SyncGroup (ierr, cctkGH, "ADMBase::curv")
        
     end if
     
     if (overwrite_lapse(n) /= 0) then
        
        call CCTK_SyncGroup (ierr, cctkGH, "ADMBase::lapse")
        
     end if
     
     if (overwrite_shift(n) /= 0) then
        
        if (shift_state /= 0) then
           call CCTK_SyncGroup (ierr, cctkGH, "ADMBase::shift")
        end if
        
     end if
     
  end do
  
contains
  
  subroutine smooth (val)
    CCTK_REAL, intent(inout) :: val(:,:,:)
    
    where (mask > 0)
       val =  (1 - mask * smoothing_factor) * val &
            & + (mask * smoothing_factor / 6) &
            &   * (+ eoshift(val, +1, dim=1) + eoshift(val, -1, dim=1) &
            &      + eoshift(val, +1, dim=2) + eoshift(val, -1, dim=2) &
            &      + eoshift(val, +1, dim=3) + eoshift(val, -1, dim=3))
    end where
  end subroutine smooth
  
end subroutine NoExcision_Smooth
