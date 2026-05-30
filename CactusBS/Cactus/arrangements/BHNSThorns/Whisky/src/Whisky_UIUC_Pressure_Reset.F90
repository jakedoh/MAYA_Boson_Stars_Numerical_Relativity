 /*@@
   @file      Whisky_UIUC_Pressure_Reset.F90
   @date      Thu Jan 28 14:05:27 EST 2010
   @author    Roland Haas
   @desc      
      Implements the pressure limit described in arXiv:0812.2245v2 ie keeps
      pressure in a certain range in low density regions and inside the
      AH
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "Whisky_Utils.h"
#include "SpaceMask.h"

 /*@@
   @routine    Whisky_UIUC_Pressure_Reset
   @date       Thu Jan 28 14:05:27 EST 2010
   @author     Roland Haas
   @desc 
      Implements the pressure limit described in arXiv:0812.2245v2 ie keeps
      pressure in a certain range in low density regions and inside the
      AH
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Whisky_UIUC_Pressure_Reset(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, k
  CCTK_REAL :: detg
  CCTK_REAL :: rho_max, press_max, press_min, detg_min, press_poly
  CHARACTER(len=400) :: infoline

  CCTK_INT :: type_bits, atmosphere, not_atmosphere
  
#ifdef _EOS_BASE_INC_
#undef _EOS_BASE_INC_
#endif
#include "EOS_Base.inc"

  call SpaceMask_GetTypeBits(type_bits, "Hydro_Atmosphere")
  call SpaceMask_GetStateBits(atmosphere, "Hydro_Atmosphere",&
                              "in_atmosphere")
  call SpaceMask_GetStateBits(not_atmosphere, "Hydro_Atmosphere",&
                              "not_in_atmosphere")

  ! compute the actual min/max values from the Cactus parameters
  rho_max = uiuc_rho_max_fact * whisky_rho_min
  detg_min = uiuc_psi6_min**2

  !$OMP PARALLEL DO &
  !$OMP PRIVATE (i, j, k, detg, press_max, press_min, press_poly, &
  !$OMP          infoline)
  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)
        
        ! skip points already marked as atmosphere (hopefully this is
        ! the right thing to do)
        if ( (atmosphere_mask(i, j, k) .eq. 1) &
             &.or. (SpaceMask_CheckStateBitsF90(space_mask,i, j, k, type_bits,\
                    atmosphere)) &
             &) then
             cycle
        end if

        ! is the pressure out of bounds? 
        press_poly = EOS_Pressure(whisky_polytrope_handle, rho(i,j,k), -1.d0)
        if (((press_poly .lt. 0d0) .or. (press_poly+1. .eq. press_poly )) &
            .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
          !$OMP CRITICAL
          write(infoline, '("Whisky_UIUC_Pressure_Reset: press_poly = ",g15.6,&
                         &" invalid for rho = ",g15.6)') press_poly, rho(i,j,k)
          call CCTK_WARN(1, infoline)
          !$OMP END CRITICAL
        end if
        press_min = uiuc_press_min_fact*press_poly
        press_max = uiuc_press_max_fact*press_poly
        if (press(i,j,k) .ge. press_min .and. press(i,j,k) .le. press_max) then
          cycle ! pressure is just fine
        end if

        ! if we get here, then the pressure is out of bounds, should we reset
        ! since we are in atmosphere or the AH?

        ! compute determinant of metric both to decide and for prim2con
        detg = SPATIAL_DET(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
        if ((detg .ne. detg .or. detg .le. 0d0) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
            !$OMP CRITICAL
            write(infoline,'("Whisky_UIUC_Pressure_Reset: Metric determinant in iteration ",i6," at [",g15.6,",",g15.6,",",&
              &g15.6,"] on level ",i2," is : ",g15.6)') cctk_iteration,&
              x(i,j,k),y(i,j,k),z(i,j,k),whisky_reflevel,detg
             call CCTK_WARN(1, infoline)
             !$OMP END CRITICAL
        end if

        ! are we inside a low density region or inside an AH?
        if (.not.(rho(i,j,k) .lt. rho_max .or. detg .gt. detg_min)) then 
          cycle
        end if

        ! if we get here we are in atmosphere or inside the AH and in either
        ! case the pressure is outside of the bounds
        if (press(i,j,k) .lt. press_min) then
          press(i,j,k) = press_min
        else if (press(i,j,k) .gt. press_max) then
          press(i,j,k) = press_max
        else
          !$OMP CRITICAL
          call CCTK_WARN(1,"Whisky_UIUC_Pressure_Reset: logic error")
          write(infoline,'("pressure ",g15.6, "is *not* out of bounds")') press(i,j,k)
          call CCTK_WARN(1,infoline)
          write(infoline,'("p_min ",g15.6, "p_max ",g15.6)') press_min, press_max
          call CCTK_WARN(1,infoline)
          call CCTK_WARN(0, "Aborting")
          !$OMP END CRITICAL
          cycle ! NOTREACHED
        end if

        ! recompute dependend quantities
        eps(i,j,k) = EOS_SpecificIntEnergy(whisky_eos_handle, rho(i,j,k), press(i,j,k))
        call prim2con(whisky_eos_handle, &
             gxx(i,j,k), gxy(i,j,k), gxz(i,j,k), &
             gyy(i,j,k), gyz(i,j,k), gzz(i,j,k), detg, &
             dens(i,j,k), sx(i,j,k), sy(i,j,k), sz(i,j,k), &
             tau(i,j,k), rho(i,j,k), velx(i,j,k), vely(i,j,k), &
             velz(i,j,k), eps(i,j,k), press(i,j,k), w_lorentz(i,j,k))

        if ((w_lorentz(i, j, k).lt.1.d0 .or. w_lorentz(i,j,k).ne.w_lorentz(i,j,k)) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
          !$OMP CRITICAL
          write(infoline,'("Whisky_UIUC_Pressure_Reset: Unphysical Lorentz factor ",g15.6,&
            &" for data g = [",5(g15.6,","),g15.6,"] vel = [",&
            &2(g15.6,","),g15.6,"] occured in iteration ",i6," at location [",&
            &2(g15.6,","),g15.6,"] on reflevel ",i2," weight ",f4.2)') &
            w_lorentz(i,j,k), gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
            velx(i,j,k),vely(i,j,k),velz(i,j,k),&
            cctk_iteration, x(i,j,k),y(i,j,k),z(i,j,k), whisky_reflevel, &
            Whisky_CarpetWeights(i,j,k)
          call CCTK_WARN(1,infoline)
          !$OMP END CRITICAL
        end if
        if ((tau(i,j,k).le.0 .or. tau(i,j,k).ne.tau(i,j,k)) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
          !$OMP CRITICAL
          write(infoline,'("Whisky_UIUC_Pressure_Reset: Unphysical tau = ",g15.6," found for data rho = ",&
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

end subroutine Whisky_UIUC_Pressure_Reset
