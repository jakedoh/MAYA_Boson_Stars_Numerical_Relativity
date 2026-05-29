 /*@@
   @file      Whisky_Puncture.F90
   @date      Fri Dec 5 16:52:15 2008 
   @author    Tanja Bode
   @desc 
   Check physicallity of matter.  This is 
   useful around punctures.  See Faber 0708.2436
   Sec. IIE for discussion.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "SpaceMask.h"
#include "Whisky_Utils.h"

 /*@@
   @routine    Whisky_PunctureCheck
   @date       Fri Dec 5 16:52:15 2008 
   @author     Tanja Bode
   @desc 
   Impose limit at puncture.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Whisky_PunctureCheck(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL, dimension(3) :: radii, center_x, center_y, center_z 
  CCTK_REAL :: gxxl, gxyl, gxzl, gyyl, gyzl, gzzl, sxl, syl, szl, taul, &
               uxx, uxy, uxz, uyy, uyz, uzz, det 
  CCTK_REAL :: s2, maxs2, flimit, xl, yl, zl, rr
  CCTK_INT :: type_bits, not_atmosphere
  integer :: i,j,k, nx,ny,nz, h, r_in_h

  !! DEBUG !!
  character(len=150) infoline
  !! DEBUG !!

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  call SpaceMask_GetTypeBits(type_bits, "Hydro_Atmosphere")
  call SpaceMask_GetStateBits(not_atmosphere, "Hydro_Atmosphere", "not_in_atmosphere")  
  call WhiskyPuncture_SetCenterVars(CCTK_PASS_FTOF,radii,center_x, center_y,center_z)

  !$OMP PARALLEL DO PRIVATE (gxxl, gxyl, gxzl, gyyl, gyzl, gzzl, sxl, syl, szl, taul, &
  !$OMP uxx, uxy, uxz, uyy, uyz, uzz, det, s2, maxs2, flimit, xl, yl, zl, rr, &
  !$OMP i,j,k, h, r_in_h, infoline)
  do k = 1 + whisky_stencil, nz - whisky_stencil 
    do j = 1 + whisky_stencil, ny - whisky_stencil
      do i = 1 + whisky_stencil, nx - whisky_stencil

         xl = x(i,j,k)
         yl = y(i,j,k)
         zl = z(i,j,k)

         r_in_h = 0
         do h=1,number_of_punctures 
            if ( radii(h).ge.0 ) then
               rr = sqrt( (xl - center_x(h))**2 + (yl - center_y(h))**2 + (zl-center_z(h))**2 )
               if ( rr.le.radii(h) ) then
                  r_in_h = 1
                  exit 
               end if
            end if
         end do
         
!         if ( (alp(i,j,k).le.whisky_puncture_lapse).and.&
         if ( (r_in_h.ne.0).and.(SpaceMask_CheckStateBitsF90(space_mask, i, j, k, \
              type_bits, not_atmosphere)) ) then

            !! Load quantities 
            taul = tau(i,j,k)
            sxl = sx(i,j,k)
            syl = sy(i,j,k)
            szl = sz(i,j,k)

            gxxl = gxx(i,j,k)
            gxyl = gxy(i,j,k)
            gxzl = gxz(i,j,k)
            gyyl = gyy(i,j,k)
            gyzl = gyz(i,j,k)
            gzzl = gzz(i,j,k)

            det = SPATIAL_DET(gxxl, gxyl, gxzl, gyyl, gyzl, gzzl)

            call UpperMetric(uxx, uxy, uxz, uyy, uyz, uzz, det, gxxl,&
                 gxyl, gxzl, gyyl, gyzl, gzzl)

            !! Set a minimum tau --> Currently dealt with in atmosphere
            if ( taul.le.0d0 ) then
               taul = puncture_tau_min
               tau(i,j,k) = taul
               !$OMP CRITICAL
               call CCTK_INFO("Puncture Check: tau negative. This should have been checked as atmosphere!")
               !$OMP END CRITICAL
            end if

            !! Impose upper limit of |S|^2
            s2 = uxx*sxl**2 + uyy*syl**2 + uzz*szl**2 + 2.0d0*uxy*sxl*syl + &
                 2.0d0*uxz*sxl*szl + 2.0d0*uyz*syl*szl
            maxs2 = puncture_flimit * taul * (taul + 2.0d0*dens(i,j,k))

            if ( s2.gt.maxs2 ) then
                flimit = sqrt(maxs2/s2)
                if (flimit.ne.flimit .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
                   !$OMP CRITICAL
                   write(infoline, '(a61,4g15.7,a1)') 'Whisky_Puncture: Scale Factor is NaN! (tau, rho, s2, maxs2)=(',&
                         taul, rho(i,j,k), s2, maxs2, ')'
                   call CCTK_WARN(0,infoline)
                   !$OMP END CRITICAL
                end if
                sx(i,j,k) = flimit*sxl
                sy(i,j,k) = flimit*syl
                sz(i,j,k) = flimit*szl

                !! DEBUG !!
                !write(infoline, '(a20,g15.7,a17,g15.7)') 'S_i rescaled at r = ',&
                !      sqrt(xl**2+yl**2+zl**2),': scale factor = ', flimit
                !call CCTK_INFO(infoline)
                !! DEBUG !!

            end if

         end if

      end do
    end do
  end do 
  !$OMP END PARALLEL DO

end subroutine Whisky_PunctureCheck 

 /*@@
   @routine    WhiskyPuncture_setCenterPars
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

subroutine WhiskyPuncture_setCenterVars(CCTK_ARGUMENTS, radii, center_x, center_y, center_z)

   implicit none

   DECLARE_CCTK_ARGUMENTS
   DECLARE_CCTK_PARAMETERS
   DECLARE_CCTK_FUNCTIONS

   CCTK_REAL, dimension(3) :: radii, center_x, center_y, center_z 
   integer :: pct, s_ind
   character(len=200) warnline

   do pct = 1,number_of_punctures

     radii(pct) = -1d0;
     center_x(pct) = 0d0;
     center_y(pct) = 0d0;
     center_z(pct) = 0d0;

     if ( CCTK_EQUALS(puncture_center_from(pct),"parameter") .and. CCTK_EQUALS(puncture_radius_from(pct),"parameter") ) then

        radii(pct) = puncture_radius(pct)
        center_x(pct) = puncture_center_x(pct)
        center_y(pct) = puncture_center_y(pct)
        center_z(pct) = puncture_center_z(pct)

     else 

        s_ind = puncture_surface(pct)+1
        if ( sf_active(s_ind).ne.0 ) then

           if ( sf_valid(s_ind).lt.0d0 ) then
              write(warnline,'(a13,i2,a19,i2,a64)') "Surface index", s_ind, " for Faber puncture", pct, &
                    " limit is not a valid surface this timestep. Not applying limit."
              call CCTK_WARN(1,warnline)
           else 

              center_x(pct) = sf_centroid_x(s_ind)
              center_y(pct) = sf_centroid_y(s_ind)
              center_z(pct) = sf_centroid_z(s_ind)

              if ( CCTK_EQUALS(puncture_radius_from(pct),"spherical surface") ) then
                 radii(pct) = sf_min_radius(s_ind)
              else 
                 radii(pct) = puncture_radius(pct)
              end if

           end if

        else 
           write(warnline,'(a13,i2,a19,i2,a59)') "Surface index", s_ind, " for Faber puncture", pct, &
                 " limit is not active! Not applying limit for this puncture."
           call CCTK_WARN(1,warnline) 
        end if
        

     end if

   end do

end subroutine WhiskyPuncture_setCenterVars
