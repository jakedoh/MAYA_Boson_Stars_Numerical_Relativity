! $Header$

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine NoExcision_OverwriteBSSN (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  CCTK_REAL, parameter :: zero=0
  CCTK_REAL :: dist2 (cctk_lsh(1), cctk_lsh(2), cctk_lsh(3))
  CCTK_REAL :: cx, cy, cz, radx, rady, radz, width, scale
  integer   :: timelevels
  integer   :: n
  
  integer, parameter :: sm_linear = 1
  integer, parameter :: sm_spline = 2
  integer :: sm_type
  
  do n = 1, num_regions
     
     if (cctk_iteration == iteration(n)) then
        
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

        if (CCTK_EQUALS (smoothing_function(n), "linear")) then
           sm_type = sm_linear
        else if (CCTK_EQUALS (smoothing_function(n), "spline")) then
           sm_type = sm_spline
        else
           call CCTK_WARN (0, "internal error")
        end if
        width = smoothing_zone_width(n)
        
        dist2 =   ((x - cx) / radx)**2 + ((y - cy) / rady)**2 &
             &  + ((z - cz) / radz)**2

        if (overwrite_geometry(n) /= 0) then
           
           scale = Minkowski_scale(n)
           
           if (conformal_state >= 1) then
              where (dist2 <= 1)
                 psi = 1
              end where
           end if
           if (conformal_state >= 2) then
              where (dist2 <= 1)
                 psix = 0
                 psiy = 0
                 psiz = 0
              end where
           end if
           if (conformal_state >= 3) then
              where (dist2 <= 1)
                 psixx = 0
                 psixy = 0
                 psixz = 0
                 psiyy = 0
                 psiyz = 0
                 psizz = 0
              end where
           end if

           where (dist2 <= 1)
              ADM_BS_phi = smooth (ADM_BS_phi, zero , dist2)
              ADM_BS_gxx = smooth (ADM_BS_gxx, scale, dist2)
              ADM_BS_gxy = smooth (ADM_BS_gxy, zero , dist2)
              ADM_BS_gxz = smooth (ADM_BS_gxz, zero , dist2)
              ADM_BS_gyy = smooth (ADM_BS_gyy, scale, dist2)
              ADM_BS_gyz = smooth (ADM_BS_gyz, zero , dist2)
              ADM_BS_gzz = smooth (ADM_BS_gzz, scale, dist2)
              ADM_BS_Gx  = smooth (ADM_BS_Gx , zero , dist2)
              ADM_BS_Gy  = smooth (ADM_BS_Gy , zero , dist2)
              ADM_BS_Gz  = smooth (ADM_BS_Gz , zero , dist2)
              ADM_BS_K   = smooth (ADM_BS_K  , zero , dist2)
              ADM_BS_Axx = smooth (ADM_BS_Axx, zero , dist2)
              ADM_BS_Axy = smooth (ADM_BS_Axy, zero , dist2)
              ADM_BS_Axz = smooth (ADM_BS_Axz, zero , dist2)
              ADM_BS_Ayy = smooth (ADM_BS_Ayy, zero , dist2)
              ADM_BS_Ayz = smooth (ADM_BS_Ayz, zero , dist2)
              ADM_BS_Azz = smooth (ADM_BS_Azz, zero , dist2)
           end where
           call CCTK_ActiveTimeLevelsVN &
                (timelevels, cctkGH, "ADM_BSSN::ADM_BS_phi")
           if (timelevels > 1) then
              where (dist2 <= 1)
                 ADM_BS_phi_p = smooth (ADM_BS_phi_p, zero , dist2)
                 ADM_BS_gxx_p = smooth (ADM_BS_gxx_p, scale, dist2)
                 ADM_BS_gxy_p = smooth (ADM_BS_gxy_p, zero , dist2)
                 ADM_BS_gxz_p = smooth (ADM_BS_gxz_p, zero , dist2)
                 ADM_BS_gyy_p = smooth (ADM_BS_gyy_p, scale, dist2)
                 ADM_BS_gyz_p = smooth (ADM_BS_gyz_p, zero , dist2)
                 ADM_BS_gzz_p = smooth (ADM_BS_gzz_p, scale, dist2)
                 ADM_BS_Gx_p  = smooth (ADM_BS_Gx_p , zero , dist2)
                 ADM_BS_Gy_p  = smooth (ADM_BS_Gy_p , zero , dist2)
                 ADM_BS_Gz_p  = smooth (ADM_BS_Gz_p , zero , dist2)
                 ADM_BS_K_p   = smooth (ADM_BS_K_p  , zero , dist2)
                 ADM_BS_Axx_p = smooth (ADM_BS_Axx_p, zero , dist2)
                 ADM_BS_Axy_p = smooth (ADM_BS_Axy_p, zero , dist2)
                 ADM_BS_Axz_p = smooth (ADM_BS_Axz_p, zero , dist2)
                 ADM_BS_Ayy_p = smooth (ADM_BS_Ayy_p, zero , dist2)
                 ADM_BS_Ayz_p = smooth (ADM_BS_Ayz_p, zero , dist2)
                 ADM_BS_Azz_p = smooth (ADM_BS_Azz_p, zero , dist2)
              end where
           end if
           if (timelevels > 2) then
              where (dist2 <= 1)
                 ADM_BS_phi_p_p = smooth (ADM_BS_phi_p_p, zero , dist2)
                 ADM_BS_gxx_p_p = smooth (ADM_BS_gxx_p_p, scale, dist2)
                 ADM_BS_gxy_p_p = smooth (ADM_BS_gxy_p_p, zero , dist2)
                 ADM_BS_gxz_p_p = smooth (ADM_BS_gxz_p_p, zero , dist2)
                 ADM_BS_gyy_p_p = smooth (ADM_BS_gyy_p_p, scale, dist2)
                 ADM_BS_gyz_p_p = smooth (ADM_BS_gyz_p_p, zero , dist2)
                 ADM_BS_gzz_p_p = smooth (ADM_BS_gzz_p_p, scale, dist2)
                 ADM_BS_Gx_p_p  = smooth (ADM_BS_Gx_p_p , zero , dist2)
                 ADM_BS_Gy_p_p  = smooth (ADM_BS_Gy_p_p , zero , dist2)
                 ADM_BS_Gz_p_p  = smooth (ADM_BS_Gz_p_p , zero , dist2)
                 ADM_BS_K_p_p   = smooth (ADM_BS_K_p_p  , zero , dist2)
                 ADM_BS_Axx_p_p = smooth (ADM_BS_Axx_p_p, zero , dist2)
                 ADM_BS_Axy_p_p = smooth (ADM_BS_Axy_p_p, zero , dist2)
                 ADM_BS_Axz_p_p = smooth (ADM_BS_Axz_p_p, zero , dist2)
                 ADM_BS_Ayy_p_p = smooth (ADM_BS_Ayy_p_p, zero , dist2)
                 ADM_BS_Ayz_p_p = smooth (ADM_BS_Ayz_p_p, zero , dist2)
                 ADM_BS_Azz_p_p = smooth (ADM_BS_Azz_p_p, zero , dist2)
              end where
           end if

        end if

        if (overwrite_lapse(n) /= 0) then
           
           scale = lapse_scale(n)
           
           where (dist2 <= 1)
              alp          = smooth (alp         , scale, dist2)
              ADM_BS_dtalp = smooth (ADM_BS_dtalp, zero , dist2)
           end where
           call CCTK_ActiveTimeLevelsVN &
                (timelevels, cctkGH, "ADMBase::alp")
           if (timelevels > 1) then
              where (dist2 <= 1)
                 alp_p          = smooth (alp_p         , scale, dist2)
                 ADM_BS_dtalp_p = smooth (ADM_BS_dtalp_p, zero , dist2)
              end where
           end if
           if (timelevels > 2) then
              where (dist2 <= 1)
                 alp_p_p          = smooth (alp_p_p         , scale, dist2)
                 ADM_BS_dtalp_p_p = smooth (ADM_BS_dtalp_p_p, zero , dist2)
              end where
           end if
           
        end if

        if (overwrite_shift(n) /= 0) then

           if (shift_state /= 0) then
              where (dist2 <= 1)
                 betax     = smooth (betax    , zero, dist2)
                 betay     = smooth (betay    , zero, dist2)
                 betaz     = smooth (betaz    , zero, dist2)
                 ADM_BS_Bx = smooth (ADM_BS_Bx, zero, dist2)
                 ADM_BS_By = smooth (ADM_BS_By, zero, dist2)
                 ADM_BS_Bz = smooth (ADM_BS_Bz, zero, dist2)
              end where
              call CCTK_ActiveTimeLevelsVN &
                   (timelevels, cctkGH, "ADMBase::betax")
              if (timelevels > 1) then
                 where (dist2 <= 1)
                    betax_p     = smooth (betax_p    , zero, dist2)
                    betay_p     = smooth (betay_p    , zero, dist2)
                    betaz_p     = smooth (betaz_p    , zero, dist2)
                    ADM_BS_Bx_p = smooth (ADM_BS_Bx_p, zero, dist2)
                    ADM_BS_By_p = smooth (ADM_BS_By_p, zero, dist2)
                    ADM_BS_Bz_p = smooth (ADM_BS_Bz_p, zero, dist2)
                 end where
              end if
              if (timelevels > 2) then
                 where (dist2 <= 1)
                    betax_p_p     = smooth (betax_p_p    , zero, dist2)
                    betay_p_p     = smooth (betay_p_p    , zero, dist2)
                    betaz_p_p     = smooth (betaz_p_p    , zero, dist2)
                    ADM_BS_Bx_p_p = smooth (ADM_BS_Bx_p_p, zero, dist2)
                    ADM_BS_By_p_p = smooth (ADM_BS_By_p_p, zero, dist2)
                    ADM_BS_Bz_p_p = smooth (ADM_BS_Bz_p_p, zero, dist2)
                 end where
              end if
           end if

        end if

     end if

  end do
  
contains
  
  elemental function smooth (oldval, goodval, dist2) result (newval)
    CCTK_REAL             :: newval
    CCTK_REAL, intent(in) :: oldval
    CCTK_REAL, intent(in) :: goodval
    CCTK_REAL, intent(in) :: dist2
    
    CCTK_REAL, parameter :: two=2, half=1/two
    CCTK_REAL :: dist
    CCTK_REAL :: x, f
    
    if (dist2 <= (1 - width)**2) then
       newval = goodval
    else if (width == 0 .or. dist2 >= 1) then
       newval = oldval
    else
       dist = sqrt(dist2)
       x = 1 - (1 - dist) / width
       select case (sm_type)
       case (sm_linear)
          f = x
       case (sm_spline)
          if (x < half) then
             f = 2 * x**2
          else
             f = 1 - 2 * (1-x)**2
          end if
       end select
       newval = oldval * f + goodval * (1 - f)
    end if
    
  end function smooth
  
end subroutine NoExcision_OverwriteBSSN
