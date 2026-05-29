 /*@@
   @file      Whisky_UpdateMask.F90
   @date      Wed Mar 13 14:18:38 2002
   @author    
   @desc 
   Alter the update terms if inside the atmosphere or excision region
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "Whisky_Utils.h"
#include "SpaceMask.h"

#define Axrhs(i,j,k) Avecrhs(i,j,k,1)
#define Ayrhs(i,j,k) Avecrhs(i,j,k,2)
#define Azrhs(i,j,k) Avecrhs(i,j,k,3)

subroutine WhiskyKillIsolatedCells(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i,j,k
  CCTK_INT :: alive_bits, kill_bits, alive, dead, will_live, will_die

  call SpaceMask_GetTypeBits(alive_bits, "Hydro_AliveCells")
  call SpaceMask_GetStateBits(alive, "Hydro_AliveCells", "alive")
  call SpaceMask_GetStateBits(dead, "Hydro_AliveCells", "dead")
  call SpaceMask_GetTypeBits(kill_bits, "Hydro_KillCells")
  call SpaceMask_GetStateBits(will_live, "Hydro_KillCells", "will_live")
  call SpaceMask_GetStateBits(will_die, "Hydro_KillCells", "will_die")

  ! kill isolated cells or cells whose density is too low
  ! NB: this will leave a one cell thick zone of alive cells, which are required
  ! to have material flow into grid level
  !$OMP PARALLEL DO PRIVATE (i,j,k)
  do k = 2,cctk_lsh(3)-1
    do j = 2,cctk_lsh(2)-1
      do i = 2,cctk_lsh(1)-1

        ! the following means: dead if rho<rho_min+*(1+tol) or an inner points
        ! surrounded by only dead points
        if ( rho(i,j,k) .le. whisky_rho_min*(1.d0+atmo_tolerance) .or.  &
             ( i.gt.1           .and. SpaceMask_CheckStateBitsF90(space_mask, i-1,j  ,k  , alive_bits, dead) .and. &
               i.lt.cctk_lsh(1) .and. SpaceMask_CheckStateBitsF90(space_mask, i+1,j  ,k  , alive_bits, dead) .and. &
               j.gt.1           .and. SpaceMask_CheckStateBitsF90(space_mask, i  ,j-1,k  , alive_bits, dead) .and. &
               j.lt.cctk_lsh(2) .and. SpaceMask_CheckStateBitsF90(space_mask, i  ,j+1,k  , alive_bits, dead) .and. &
               k.gt.1           .and. SpaceMask_CheckStateBitsF90(space_mask, i  ,j  ,k-1, alive_bits, dead) .and. &
               k.lt.cctk_lsh(3) .and. SpaceMask_CheckStateBitsF90(space_mask, i  ,j  ,k+1, alive_bits, dead) ) &
            ) then
        
          SpaceMask_SetStateBitsF90(space_mask, i, j, k, kill_bits, will_die)

        else

          SpaceMask_SetStateBitsF90(space_mask, i, j, k, kill_bits, will_live)

        end if

      end do
    end do
  end do
  !$OMP END PARALLEL DO

end subroutine WhiskyKillIsolatedCells

 /*@@
   @routine    WhiskySelectMaskBCs
   @date       Tue Oct  6 17:20:19 EDT 2009
   @author     
   @desc 
   Calls the appropriate boundary routines for the space_mask
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine WhiskySelectMaskBCs(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  integer :: ierr

  CCTK_INT, parameter :: faces=CCTK_ALL_FACES
  CCTK_INT, parameter :: ione=1
  
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, whisky_stencil, -ione, &
                    "SpaceMask::space_mask_group", "none");

  if (ierr < 0) call CCTK_WARN(0, "problems with applying the chosen boundary condition")

end subroutine WhiskySelectMaskBCs

subroutine WhiskyBringCellsAlive(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i,j,k
  CCTK_INT :: alive_bits, kill_bits, alive, dead, will_live, will_die
  CHARACTER(len=400) :: infoline
  CCTK_REAL :: det

  call SpaceMask_GetTypeBits(alive_bits, "Hydro_AliveCells")
  call SpaceMask_GetStateBits(alive, "Hydro_AliveCells", "alive")
  call SpaceMask_GetStateBits(dead, "Hydro_AliveCells", "dead")
  call SpaceMask_GetTypeBits(kill_bits, "Hydro_KillCells")
  call SpaceMask_GetStateBits(will_live, "Hydro_KillCells", "will_live")
  call SpaceMask_GetStateBits(will_die, "Hydro_KillCells", "will_die")

  ! execute the verdicts
  !$OMP PARALLEL DO PRIVATE (i,j,k,det, infoline)
  do k = 1,cctk_lsh(3)
    do j = 1,cctk_lsh(2)
      do i = 1,cctk_lsh(1)

        if ( SpaceMask_CheckStateBitsF90(space_mask, i,j,k  , kill_bits, will_live) ) then
          
          SpaceMask_SetStateBitsF90(space_mask, i,j,k, alive_bits, alive)

        else
          
          SpaceMask_SetStateBitsF90(space_mask, i,j,k  , alive_bits, dead)

          rho(i,j,k) = whisky_rho_min
          velx(i,j,k) = 0.0d0
          vely(i,j,k) = 0.0d0
          velz(i,j,k) = 0.0d0
          det = SPATIAL_DET(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
          if ((det .ne. det .or. det .lt. 0d0) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
              !$OMP CRITICAL
              write(infoline,'("WhiskyBringCellsAlive: Metric determinant in iteration ",i6," at [",g15.6,",",g15.6,",",&
                &g15.6,"] on level ",i2," is : ",g15.6)') cctk_iteration,&
                x(i,j,k),y(i,j,k),z(i,j,k),whisky_reflevel,det
              call CCTK_WARN(1, infoline)
              !$OMP END CRITICAL
          end if
          call prim2conpolytype(whisky_polytrope_handle, &
               gxx(i,j,k), gxy(i,j,k), gxz(i,j,k), &
               gyy(i,j,k), gyz(i,j,k), gzz(i,j,k), det, &
               dens(i,j,k), sx(i,j,k), sy(i,j,k), sz(i,j,k), &
               tau(i,j,k), rho(i,j,k), velx(i,j,k), vely(i,j,k), &
               velz(i,j,k), eps(i,j,k), press(i,j,k), w_lorentz(i,j,k))
          if ((w_lorentz(i, j, k).lt.1.d0 .or. w_lorentz(i,j,k).ne.w_lorentz(i,j,k)) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
            !$OMP CRITICAL
            write(infoline,'("WhiskyBringCellsAlive: Unphysical Lorentz factor ",g15.6,&
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
            write(infoline,'("WhiskyBringCellsAlive: Unphysical tau = ",g15.6," found for data rho = ",&
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
        end if

      end do
    end do
  end do
  !$OMP END PARALLEL DO

end subroutine WhiskyBringCellsAlive

subroutine Whisky_ScanForDeadCells(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT :: type_bits, alive, dead, i, j, k
  CHARACTER(len=400) :: infoline
  CCTK_REAL :: det

  call SpaceMask_GetTypeBits(type_bits, "Hydro_AliveCells")
  call SpaceMask_GetStateBits(dead, "Hydro_AliveCells", "dead")
  call SpaceMask_GetStateBits(alive, "Hydro_AliveCells", "alive")

  !$OMP PARALLEL DO PRIVATE (i,j,k,det, infoline)
  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)

        if (rho(i,j,k) .le. whisky_rho_min*(1.d0+atmo_tolerance)) then
  
          SpaceMask_SetStateBitsF90(space_mask, i, j, k, type_bits, dead)

          rho(i,j,k) = whisky_rho_min
          velx(i,j,k) = 0.0d0
          vely(i,j,k) = 0.0d0
          velz(i,j,k) = 0.0d0
          det = SPATIAL_DET(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
          if ((det .ne. det .or. det .lt. 0d0) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
              !$OMP CRITICAL
              write(infoline,'("KillIsolatedCells: Metric determinant in iteration ",i6," at [",g15.6,",",g15.6,",",&
                &g15.6,"] on level ",i2," is : ",g15.6)') cctk_iteration,&
                x(i,j,k),y(i,j,k),z(i,j,k),whisky_reflevel,det
              call CCTK_WARN(1, infoline)
              !$OMP END CRITICAL
          end if
          call prim2conpolytype(whisky_polytrope_handle, &
               gxx(i,j,k), gxy(i,j,k), gxz(i,j,k), &
               gyy(i,j,k), gyz(i,j,k), gzz(i,j,k), det, &
               dens(i,j,k), sx(i,j,k), sy(i,j,k), sz(i,j,k), &
               tau(i,j,k), rho(i,j,k), velx(i,j,k), vely(i,j,k), &
               velz(i,j,k), eps(i,j,k), press(i,j,k), w_lorentz(i,j,k))
          if ((w_lorentz(i, j, k).lt.1.d0 .or. w_lorentz(i,j,k).ne.w_lorentz(i,j,k)) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
            !$OMP CRITICAL
            write(infoline,'("KillIsolatedCells: Unphysical Lorentz factor ",g15.6,&
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
            write(infoline,'("KillIsolatedCells: Unphysical tau = ",g15.6," found for data rho = ",&
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

        else

          SpaceMask_SetStateBitsF90(space_mask, i, j, k, type_bits, alive)

        end if

      end do
    end do
  end do
  !$OMP END PARALLEL DO
  
end subroutine Whisky_ScanForDeadCells
  

subroutine WhiskyUpdateAtmosphereMask(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i,j,k
  CCTK_REAL :: frac

  CCTK_INT :: type_bits, atmosphere

  call SpaceMask_GetTypeBits(type_bits, "Hydro_Atmosphere")
  call SpaceMask_GetStateBits(atmosphere, "Hydro_Atmosphere", "in_atmosphere")

  frac = CCTK_DELTA_TIME

  !$OMP PARALLEL DO PRIVATE (i,j,k)
  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)
        if ( (SpaceMask_CheckStateBitsF90(space_mask, i, j, k, \
                                          type_bits, atmosphere)) .or. &
             (tau(i,j,k) + frac * taurhs(i,j,k) .le. 0.d0).or.&
             (dens(i,j,k) + frac * densrhs(i,j,k) .le. 0.d0) ) then
          densrhs(i,j,k) = 0.0d0
          sxrhs(i,j,k)   = 0.0d0
          syrhs(i,j,k)   = 0.0d0
          szrhs(i,j,k)   = 0.0d0
          taurhs(i,j,k)  = 0.0d0
          atmosphere_mask(i,j,k) = 1

          SpaceMask_SetStateBitsF90(space_mask, i, j, k, type_bits, atmosphere)
        end if
      end do
    end do
  end do
  !$OMP END PARALLEL DO

end subroutine WhiskyUpdateAtmosphereMask

subroutine WhiskyDampAtmosphere(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i,j,k
  CCTK_REAL :: damping_factor

  !$OMP PARALLEL DO PRIVATE (i,j,k,damping_factor)
  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)
          if ( rho(i,j,k) >= whisky_rho_min ) then
             !! Damp only within atmosphere
             cycle
          end if

          damping_factor = exp( 1. - whisky_rho_min/rho(i,j,k) )
          densrhs(i,j,k) = densrhs(i,j,k)*damping_factor
          sxrhs(i,j,k)   = sxrhs(i,j,k)*damping_factor 
          syrhs(i,j,k)   = syrhs(i,j,k)*damping_factor
          szrhs(i,j,k)   = szrhs(i,j,k)*damping_factor
          taurhs(i,j,k)  = taurhs(i,j,k)*damping_factor

          if ( whisky_mhd_handle .eq. 2 ) then
             Bnxrhs(i,j,k) = Bnxrhs(i,j,k)*damping_factor
             Bnyrhs(i,j,k) = Bnyrhs(i,j,k)*damping_factor
             Bnzrhs(i,j,k) = Bnzrhs(i,j,k)*damping_factor
          end if
          if ( whisky_mhd_handle .ge. 3 ) then
             Axrhs(i,j,k) = Axrhs(i,j,k)*damping_factor
             Ayrhs(i,j,k) = Ayrhs(i,j,k)*damping_factor
             Azrhs(i,j,k) = Azrhs(i,j,k)*damping_factor
          end if
          if ( whisky_mhd_handle .eq. 4 ) then
             Aphirhs(i,j,k) = Aphirhs(i,j,k)*damping_factor
          end if

      end do
    end do
  end do
  !$OMP END PARALLEL DO

end subroutine WhiskyDampAtmosphere

 /*@@
   @routine    Whisky_CheckForAtmosphere
   @date       Wed Feb  4 08:31:16 EST 2009
   @author     Roland Haas
   @desc 
   Check newly updated grid functions for invalid data and reset to previous
   value of found. Set atmoshpere mask.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Whisky_CheckForAtmosphere(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i,j,k

  CCTK_INT :: type_bits, not_atmosphere, atmosphere
  CCTK_REAL :: det
  CHARACTER(len=400) :: infoline
  
  call SpaceMask_GetTypeBits(type_bits, "Hydro_Atmosphere")
  call SpaceMask_GetStateBits(atmosphere, "Hydro_Atmosphere", "in_atmosphere")
  call SpaceMask_GetStateBits(not_atmosphere, "Hydro_Atmosphere", "not_in_atmosphere")

  !$OMP PARALLEL DO PRIVATE (i,j,k,det, infoline)
  do k = 1+whisky_stencil, cctk_lsh(3)-whisky_stencil
    do j = 1+whisky_stencil, cctk_lsh(2)-whisky_stencil
      do i = 1+whisky_stencil, cctk_lsh(1)-whisky_stencil
        if ( (SpaceMask_CheckStateBitsF90(space_mask, i, j, k, \
                                          type_bits, not_atmosphere)) .and. ( &
!             (tau(i,j,k) .ne. tau(i,j,k)) .or. &
              (tau(i,j,k) .le. 0.d0) .or. &
!             (dens(i,j,k) .ne. dens(i,j,k)) .or. &
             (dens(i,j,k) .le. 0.d0) ) ) then


          rho(i,j,k) = whisky_rho_min
          velx(i,j,k) = 0.0d0
          vely(i,j,k) = 0.0d0
          velz(i,j,k) = 0.0d0
          det = SPATIAL_DET(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
          if ((det .ne. det .or. det .lt. 0d0) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
              !$OMP CRITICAL
              write(infoline,'("CheckForAtmosphere: Metric determinant in iteration ",i6," at [",g15.6,",",g15.6,",",&
                &g15.6,"] on level ",i2," is : ",g15.6)') cctk_iteration,&
                x(i,j,k),y(i,j,k),z(i,j,k),whisky_reflevel,det
              call CCTK_WARN(1, infoline)
              !$OMP END CRITICAL
          end if
          call prim2conpolytype(whisky_polytrope_handle, &
               gxx(i,j,k), gxy(i,j,k), gxz(i,j,k), &
               gyy(i,j,k), gyz(i,j,k), gzz(i,j,k), det, &
               dens(i,j,k), sx(i,j,k), sy(i,j,k), sz(i,j,k), &
               tau(i,j,k), rho(i,j,k), velx(i,j,k), vely(i,j,k), &
               velz(i,j,k), eps(i,j,k), press(i,j,k), w_lorentz(i,j,k))
          if ((w_lorentz(i, j, k).lt.1.d0 .or. w_lorentz(i,j,k).ne.w_lorentz(i,j,k)) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
            !$OMP CRITICAL
            write(infoline,'("CheckForAtmosphere: Unphysical Lorentz factor ",g15.6,&
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
            write(infoline,'("CheckForAtmosphere: Unphysical tau = ",g15.6," found for data rho = ",&
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

          atmosphere_mask(i,j,k) = 1
          SpaceMask_SetStateBitsF90(space_mask, i, j, k, type_bits, atmosphere)
        end if
      end do
    end do
  end do
  !$OMP END PARALLEL DO

end subroutine Whisky_CheckForAtmosphere

 /*@@
   @routine    Whisky_RemoveNaNs
   @date       Sun Nov  8 17:22:39 EST 2009
   @author     Roland Haas
   @desc 
   Check newly updated grid functions for NaNs data and reset to atmosphere
   value if found. Set atmoshpere mask.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Whisky_ClearNanMask(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT :: i,j,k

  !$OMP PARALLEL DO PRIVATE (i,j,k)
  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)
       whisky_nankind(i,j,k) = 0d0
      end do
    end do
  end do
  !$OMP END PARALLEL DO

end subroutine

subroutine Whisky_RemoveNaNsInRHS(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT :: i,j,k,ikind,ikind_seen,nancount
  CHARACTER(len=400) :: infoline

  ikind_seen = 0 ! shared
  nancount = 0 ! shared
  !$OMP PARALLEL DO PRIVATE (i,j,k,ikind)
  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)
        ikind = whisky_nankind(i,j,k)

        if(densrhs(i,j,k).ne.densrhs(i,j,k)) then
          ikind = ior(ikind, ishft(1,20))
        else if(sxrhs(i,j,k).ne.sxrhs(i,j,k)) then
          ikind = ior(ikind, ishft(1,21))
        else if(syrhs(i,j,k).ne.syrhs(i,j,k)) then
          ikind = ior(ikind, ishft(1,22))
        else if(szrhs(i,j,k).ne.szrhs(i,j,k)) then
          ikind = ior(ikind, ishft(1,23))
        else if(taurhs(i,j,k).ne.taurhs(i,j,k)) then
          ikind = ior(ikind, ishft(1,24))
        end if

        if (ikind .ne. 0d0 .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
            !$OMP CRITICAL
            if(iand(NOT(ikind_seen),ikind) .ne. 0) then 
              write(infoline,'("RemoveAllNaNsRHS: iteration ",i6," found nan of kind ",i11," at [",g15.6,",",g15.6,",",&
                &g15.6,"] on level ",i2)') cctk_iteration,ikind,&
                x(i,j,k),y(i,j,k),z(i,j,k),whisky_reflevel
              call CCTK_WARN(1, infoline)
              ikind_seen = ior(ikind_seen, ikind)
            end if
            nancount = nancount + 1
            !$OMP END CRITICAL
        end if

        if(ikind .gt. 0) then ! found a NaN in the hydro variables
          densrhs(i,j,k) = 0d0
          sxrhs(i,j,k) = 0d0
          syrhs(i,j,k) = 0d0
          szrhs(i,j,k) = 0d0
          taurhs(i,j,k) = 0d0
        end if

        whisky_nankind(i,j,k) = ikind
      end do
    end do
  end do
  !$OMP END PARALLEL DO
  if(nancount .gt. 0) then 
    write(infoline,'("RemoveAllNaNsRHS: iteration ",i6," found ",i6," nans")') &
      cctk_iteration, nancount
    call CCTK_WARN(1, infoline)
    ikind_seen = ior(ikind_seen, ikind)
  end if

end subroutine

subroutine Whisky_RemoveAllNaNs(CCTK_ARGUMENTS)

  implicit none
  
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i,j,k,ikind,ikind_seen,nancount
  CCTK_REAL :: det
  CHARACTER(len=400) :: infoline

  ikind_seen = 0 ! shared
  nancount = 0 ! shared
  !$OMP PARALLEL DO PRIVATE (i,j,k,det)
  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)
        ikind = whisky_nankind(i,j,k)

        if(dens(i,j,k).ne.dens(i,j,k)) then
          ikind = ior(ikind, ishft(1,0))
        else if(sx(i,j,k).ne.sx(i,j,k)) then
          ikind = ior(ikind, ishft(1,1))
        else if(sy(i,j,k).ne.sy(i,j,k)) then
          ikind = ior(ikind, ishft(1,2))
        else if(sz(i,j,k).ne.sz(i,j,k)) then
          ikind = ior(ikind, ishft(1,3))
        else if(tau(i,j,k).ne.tau(i,j,k)) then
          ikind = ior(ikind, ishft(1,4))
        else if(rho(i,j,k).ne.rho(i,j,k)) then
          ikind = ior(ikind, ishft(1,5))
        else if(velx(i,j,k).ne.velx(i,j,k)) then
          ikind = ior(ikind, ishft(1,6))
        else if(vely(i,j,k).ne.vely(i,j,k)) then
          ikind = ior(ikind, ishft(1,7))
        else if(velz(i,j,k).ne.velz(i,j,k)) then
          ikind = ior(ikind, ishft(1,8))
        else if(eps(i,j,k).ne.eps(i,j,k)) then
          ikind = ior(ikind, ishft(1,9))
        else if(press(i,j,k).ne.press(i,j,k)) then
          ikind = ior(ikind, ishft(1,10))
        else if(w_lorentz(i,j,k).ne.w_lorentz(i,j,k)) then
          ikind = ior(ikind, ishft(1,11))
        else if(gxx(i,j,k).ne.gxx(i,j,k)) then
          ikind = ior(ikind, ishft(1,12))
        else if(gxy(i,j,k).ne.gxy(i,j,k)) then
          ikind = ior(ikind, ishft(1,13))
        else if(gxz(i,j,k).ne.gxz(i,j,k)) then
          ikind = ior(ikind, ishft(1,14))
        else if(gyy(i,j,k).ne.gyy(i,j,k)) then
          ikind = ior(ikind, ishft(1,15))
        else if(gyz(i,j,k).ne.gyz(i,j,k)) then
          ikind = ior(ikind, ishft(1,16))
        else if(gzz(i,j,k).ne.gzz(i,j,k)) then
          ikind = ior(ikind, ishft(1,17))
        else if(alp(i,j,k).ne.alp(i,j,k)) then
          ikind = ior(ikind, ishft(1,16))
        else if(betax(i,j,k).ne.betax(i,j,k)) then
          ikind = ior(ikind, ishft(1,17))
        else if(betay(i,j,k).ne.betay(i,j,k)) then
          ikind = ior(ikind, ishft(1,18))
        else if(betaz(i,j,k).ne.betaz(i,j,k)) then
          ikind = ior(ikind, ishft(1,19))
        end if

        if (ikind .ne. 0 .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
            !$OMP CRITICAL
            if(iand(NOT(ikind_seen),ikind) .ne. 0) then 
              write(infoline,'("RemoveAllNaNs: iteration ",i6," found nan of kind ",i11," at [",g15.6,",",g15.6,",",&
                &g15.6,"] on level ",i2)') cctk_iteration,ikind,&
                x(i,j,k),y(i,j,k),z(i,j,k),whisky_reflevel
              call CCTK_WARN(1, infoline)
              ikind_seen = ior(ikind_seen, ikind)
            end if
            nancount = nancount + 1
            !$OMP END CRITICAL
        end if

        if(ikind .gt. 0 .and. ikind .lt. ishft(1,12)) then ! found a NaN in the hydro variables

          rho(i,j,k) = whisky_rho_min
          velx(i,j,k) = 0.0d0
          vely(i,j,k) = 0.0d0
          velz(i,j,k) = 0.0d0

          det = SPATIAL_DET(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
          call prim2conpolytype(whisky_polytrope_handle, &
               gxx(i,j,k), gxy(i,j,k), gxz(i,j,k), &
               gyy(i,j,k), gyz(i,j,k), gzz(i,j,k), det, &
               dens(i,j,k), sx(i,j,k), sy(i,j,k), sz(i,j,k), &
               tau(i,j,k), rho(i,j,k), velx(i,j,k), vely(i,j,k), &
               velz(i,j,k), eps(i,j,k), press(i,j,k), w_lorentz(i,j,k))

          if ((det .ne. det .or. det .lt. 0d0) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
              !$OMP CRITICAL
              write(infoline,'("RemoveNaNs: Metric determinant in iteration ",i6," at [",g15.6,",",g15.6,",",&
                &g15.6,"] on level ",i2," is : ",g15.6)') cctk_iteration,&
                x(i,j,k),y(i,j,k),z(i,j,k),whisky_reflevel,det
              call CCTK_WARN(1, infoline)
              !$OMP END CRITICAL
          end if
          if ((w_lorentz(i, j, k).lt.1.d0 .or. w_lorentz(i,j,k).ne.w_lorentz(i,j,k)) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
            !$OMP CRITICAL
            write(infoline,'("RemoveNaNs: Unphysical Lorentz factor ",g15.6,&
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
            write(infoline,'("RemoveNaNs: Unphysical tau = ",g15.6," found for data rho = ",&
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

        end if

        whisky_nankind(i,j,k) = ikind
      end do
    end do
  end do
  !$OMP END PARALLEL DO
  if(nancount .gt. 0) then 
    write(infoline,'("RemoveAllNaNs: iteration ",i6," found ",i6," nans")') &
      cctk_iteration, nancount
    call CCTK_WARN(1, infoline)
    ikind_seen = ior(ikind_seen, ikind)
  end if

end subroutine

 /*@@
   @routine    Whisky_SetupMask
   @date       Thu Jun 20 13:27:28 2002
   @author     Ian Hawke
   @desc 
   Initialize the mask to be zero.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Whisky_SetupMask(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT :: type_bits, not_atmosphere, i, j, k

  ! Initialize all rhs variables and the mask.
  ! The former vars need to be initialized since there is
  ! no rhs computation in CCTK_INITIAL or POSTINITIAL.

  densrhs = 0.0d0
  taurhs  = 0.0d0 
  sxrhs   = 0.0d0
  syrhs   = 0.0d0
  szrhs   = 0.0d0

  if (evolve_tracer .ne. 0) then
     cons_tracerrhs = 0.0d0
  endif

  atmosphere_mask = 0
  
  call SpaceMask_GetTypeBits(type_bits, "Hydro_Atmosphere")
  call SpaceMask_GetStateBits(not_atmosphere, &
       &"Hydro_Atmosphere", "not_in_atmosphere")

  !$OMP PARALLEL DO PRIVATE (i,j,k)
  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)
  
        SpaceMask_SetStateBitsF90(space_mask, i, j, k, type_bits, not_atmosphere)

      end do
    end do
  end do
  !$OMP END PARALLEL DO

  call CCTK_INFO("Setting up the atmosphere mask: all points are not_atmosphere")
  
end subroutine Whisky_SetupMask

 /*@@
   @routine    Whisky_InitAtmosMask
   @date       Thu Jun 20 13:27:28 2002
   @author     Ian Hawke
   @desc 
   Initialize the mask based on rho_min
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Whisky_InitAtmosMask(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT :: type_bits, atmosphere, i, j, k

  call SpaceMask_GetTypeBits(type_bits, "Hydro_Atmosphere")
  call SpaceMask_GetStateBits(atmosphere, &
       &"Hydro_Atmosphere", "in_atmosphere")

  !$OMP PARALLEL DO PRIVATE (i,j,k)
  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)

        if (rho(i,j,k) .le. whisky_rho_min*(1.d0+atmo_tolerance)) then
  
          SpaceMask_SetStateBitsF90(space_mask, i, j, k, type_bits, atmosphere)
          atmosphere_mask(i,j,k) = 1

        end if
        
      end do
    end do
  end do
  !$OMP END PARALLEL DO

  call CCTK_INFO("Setting up the atmosphere mask: points with rho<rho_min are set to be atmosphere")
  
end subroutine Whisky_InitAtmosMask

 /*@@
   @routine    Whisky_AtmosphereReset
   @date       Thu Jun 20 13:30:51 2002
   @author     Ian Hawke
   @desc 
   After MoL has evolved, if a point is supposed to be reset then do so.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Whisky_AtmosphereReset(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, k
  CCTK_REAL :: det
  CHARACTER(len=400) :: infoline

  CCTK_INT :: type_bits, atmosphere, not_atmosphere
  
  call SpaceMask_GetTypeBits(type_bits, "Hydro_Atmosphere")
  call SpaceMask_GetStateBits(atmosphere, "Hydro_Atmosphere",&
                              "in_atmosphere")
  call SpaceMask_GetStateBits(not_atmosphere, "Hydro_Atmosphere",&
                              "not_in_atmosphere")

  !$OMP PARALLEL DO PRIVATE (i,j,k,det,infoline)
  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)
        ! with synchronization of the primitives we only have to run in the
        ! interior, without we have to run everywhere
        
        if ( (atmosphere_mask(i, j, k) .eq. 1) &
             &.or. (SpaceMask_CheckStateBitsF90(space_mask,i, j, k, type_bits,\
                    atmosphere)) &
             &) then

!!$          write(*,*) 'Resetting at ',i,j,k, atmosphere_mask(i, j, k), &
!!$             & (SpaceMask_CheckStateBitsF90(space_mask,i, j, k, type_bits,\
!!$                    atmosphere)) 

          rho(i,j,k) = whisky_rho_min
          velx(i,j,k) = 0.0d0
          vely(i,j,k) = 0.0d0
          velz(i,j,k) = 0.0d0
          det = SPATIAL_DET(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
          if ((det .ne. det .or. det .lt. 0d0) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
              !$OMP CRITICAL
              write(infoline,'("AtmosphereReset: Metric determinant in iteration ",i6," at [",g15.6,",",g15.6,",",&
                &g15.6,"] on level ",i2," is : ",g15.6)') cctk_iteration,&
                x(i,j,k),y(i,j,k),z(i,j,k),whisky_reflevel,det
              call CCTK_WARN(1, infoline)
              !$OMP END CRITICAL
          end if
          call prim2conpolytype(whisky_polytrope_handle, &
               gxx(i,j,k), gxy(i,j,k), gxz(i,j,k), &
               gyy(i,j,k), gyz(i,j,k), gzz(i,j,k), det, &
               dens(i,j,k), sx(i,j,k), sy(i,j,k), sz(i,j,k), &
               tau(i,j,k), rho(i,j,k), velx(i,j,k), vely(i,j,k), &
               velz(i,j,k), eps(i,j,k), press(i,j,k), w_lorentz(i,j,k))
          if ((w_lorentz(i, j, k).lt.1.d0 .or. w_lorentz(i,j,k).ne.w_lorentz(i,j,k)) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
            !$OMP CRITICAL
            write(infoline,'("AtmosphereReset: Unphysical Lorentz factor ",g15.6,&
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
            write(infoline,'("AtmosphereReset: Unphysical tau = ",g15.6," found for data rho = ",&
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
          if (wk_atmosphere .eq. 0 .and. excision_atmosphere .eq. 0) then
            atmosphere_mask(i, j, k) = 0
            SpaceMask_SetStateBitsF90(space_mask, i, j, k, type_bits,\
                                      not_atmosphere)
          end if

        end if

      end do
    end do
  end do
  !$OMP END PARALLEL DO

!!$  call Whisky_Boundaries(CCTK_PASS_FTOF)

end subroutine Whisky_AtmosphereReset

 /*@@
   @routine    Whisky_ScanForAtmosphere
   @date       Mon Mar  2 15:22:05 EST 2009
   @author     Roland Haas
   @desc 
   Set atmosphere mask based on rho and rho_min
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Whisky_ScanForAtmosphere(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT :: type_bits, atmosphere, not_atmosphere, i, j, k

  if ( (scan_for_atmosphere_every .gt. 0 .and. &
        MOD(cctk_iteration,scan_for_atmosphere_every) .eq. 0) .or. &
       (excision_atmosphere .ne. 0) ) then

    call SpaceMask_GetTypeBits(type_bits, "Hydro_Atmosphere")
    call SpaceMask_GetStateBits(atmosphere, &
         &"Hydro_Atmosphere", "in_atmosphere")
    call SpaceMask_GetStateBits(not_atmosphere, &
         &"Hydro_Atmosphere", "not_in_atmosphere")

    !$OMP PARALLEL DO PRIVATE (i,j,k)
    do k = 1, cctk_lsh(3)
      do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)

          if (rho(i,j,k) .le. whisky_rho_min*(1.d0+atmo_tolerance)) then
    
            SpaceMask_SetStateBitsF90(space_mask, i, j, k, type_bits, atmosphere)
            atmosphere_mask(i,j,k) = 1

          else

            SpaceMask_SetStateBitsF90(space_mask, i, j, k, type_bits, not_atmosphere)
            atmosphere_mask(i,j,k) = 0

          end if
          
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    call CCTK_INFO("Scanning through the domain: points with rho<rho_min are set to be atmosphere")

  end if
  
end subroutine Whisky_ScanForAtmosphere

 /*@@
   @routine    Whisky_ExpandNonAtmosphere
   @date       Mon Mar  2 21:38:05 EST 2009
   @author     Roland Haas
   @desc 
   Surroudn the non-atmosphere region with a one cell wide region where
   atmosphere evolves
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Whisky_ExpandNonAtmosphere(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT :: type_bits, atmosphere, not_atmosphere, i, j, k 
!  CCTK_INT, DIMENSION(6,3) :: imin, imax
!  CCTK_INT :: face, dir, d

  call SpaceMask_GetTypeBits(type_bits, "Hydro_Atmosphere")
  call SpaceMask_GetStateBits(atmosphere, &
       &"Hydro_Atmosphere", "in_atmosphere")
  call SpaceMask_GetStateBits(not_atmosphere, &
       &"Hydro_Atmosphere", "not_in_atmosphere")

  ! this currently only creates a one cell wide growth region. If this is
  ! considered too small (it limits the speed with which matter can move into
  ! the atmosphere region to 1/MoLNumIntegratorSubsteps), then one might be able
  ! to define a new paramter excision_atmosphere_width and run the loop below
  ! several times, each time incrementing the mask type used (prevents unchecked
  ! undead cell growth)

!  ! finding in between cells is not totally local so we have to set the
!  ! atmosphere mask in the ghost zones as well (one layer only)
!  do face = 1,6 ! extends of a one layer larger cube
!    do dir = 1,3
!      imin(face,dir) = whisky_stencil
!      imax(face,dir) = cctk_lsh(dir) - whisky_stencil + 1
!    end do
!  end do
!  do dir = 1,3    ! fix extends of normal direction
!    imax(2*(dir-1)+1,dir) = whisky_stencil
!    imin(2*dir      ,dir) = cctk_lsh(dir) - whisky_stencil + 1
!  end do
!  do face = 1,6 
!    do k = imin(face,3), imax(face,3)
!      do j = imin(face,2), imax(face,2)
!        do i = imin(face,1), imax(face,1)
!
!          ! find atmospheric cells
!          if (rho(i,j,k) .le. whisky_rho_min*(1.d0+atmo_tolerance) ) then
!    
!            SpaceMask_SetStateBitsF90(space_mask, i, j, k, type_bits, atmosphere)
!            atmosphere_mask(i,j,k) = 1
!
!          end if
!        end do
!      end do
!    end do
!  end do

  ! find cells that are in between atmosphere (dead) and  non-atmosphere (live)
  ! cells 
  ! note that we leave the boundary region untreated
  !$OMP PARALLEL DO PRIVATE (i,j,k)
  do k = whisky_stencil + 1, cctk_lsh(3) - whisky_stencil
    do j = whisky_stencil + 1, cctk_lsh(2) - whisky_stencil
      do i = whisky_stencil + 1, cctk_lsh(1) - whisky_stencil

        if ( atmosphere_mask(i,j,k) .eq. 1 ) then
          
          ! hopefully the treatment above makes this all well defined

          if ( (atmosphere_mask(i-1,j  ,k  ) .eq. 0) .or. &
               (atmosphere_mask(i+1,j  ,k  ) .eq. 0) .or. &
               (atmosphere_mask(i  ,j-1,k  ) .eq. 0) .or. &
               (atmosphere_mask(i  ,j+1,k  ) .eq. 0) .or. &
               (atmosphere_mask(i  ,j  ,k-1) .eq. 0) .or. &
               (atmosphere_mask(i  ,j  ,k+1) .eq. 0) &
              ) then

            ! dead cell bordering a live cell => limbo
            SpaceMask_SetStateBitsF90(space_mask, i, j, k, type_bits, not_atmosphere)
            atmosphere_mask(i, j, k) = 2

          end if

       end if

      end do
    end do
  end do
  !$OMP END PARALLEL DO

end subroutine Whisky_ExpandNonAtmosphere

 /*@@
   @routine    Whisky_CheckExpandedNonAtmosphere
   @date       Mon Mar  3 21:38:05 EST 2009
   @author     Roland Haas
   @desc 
   Bring survisving undead cells to live
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Whisky_CheckExpandedNonAtmosphere(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT :: type_bits, atmosphere, not_atmosphere, i, j, k

  call SpaceMask_GetTypeBits(type_bits, "Hydro_Atmosphere")
  call SpaceMask_GetStateBits(atmosphere, "Hydro_Atmosphere", "in_atmosphere")
  call SpaceMask_GetStateBits(not_atmosphere, &
       &"Hydro_Atmosphere", "not_in_atmosphere")

  !$OMP PARALLEL DO PRIVATE (i,j,k)
  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)

        ! find atmospheric cells
        if (rho(i,j,k) .le. whisky_rho_min*(1.d0+atmo_tolerance)) then
  
          SpaceMask_SetStateBitsF90(space_mask, i, j, k, type_bits, atmosphere)
          atmosphere_mask(i,j,k) = 1

        end if

        ! promote surviving undead cell to live cells
        if ( atmosphere_mask(i,j,k) .eq. 2 ) then

          SpaceMask_SetStateBitsF90(space_mask, i, j, k, type_bits, not_atmosphere)
          atmosphere_mask(i, j, k) = 0

        end if

      end do
    end do
  end do
  !$OMP END PARALLEL DO

end subroutine Whisky_CheckExpandedNonAtmosphere

subroutine Whisky_SanityCheck(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, k
  CCTK_REAL :: v2
  CHARACTER(len=400) :: infoline

  !$OMP PARALLEL DO PRIVATE (i,j,k,v2,infoline)
  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)

        v2 =   gxx(i,j,k)*velx(i,j,k)*velx(i,j,k) +   gyy(i,j,k)*vely(i,j,k)*vely(i,j,k) + &
               gzz(i,j,k)*velz(i,j,k)*velz(i,j,k) + 2*gxy(i,j,k)*velx(i,j,k)*vely(i,j,k) + &
             2*gxz(i,j,k)*velx(i,j,k)*velz(i,j,k) + 2*gyz(i,j,k)*vely(i,j,k)*velz(i,j,k) 
        
        if ((w_lorentz(i, j, k).lt.1.d0 .or. w_lorentz(i,j,k).ne.w_lorentz(i,j,k)) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
          !$OMP CRITICAL
          write(infoline,'("SanityCheck: Unphysical Lorentz factor ",g15.6,&
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

        if ((v2.ge.1.d0 .or. v2.ne.v2) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
            !$OMP CRITICAL
          write(infoline,'("SanityCheck: Unphysical speed ",g15.6,&
            &" for data g = [",5(g15.6,","),g15.6,"] vel = [",&
            &2(g15.6,","),g15.6,"] occured in iteration ",i6," at location [",&
            &2(g15.6,","),g15.6,"] on reflevel ",i2," weight ",f4.2)') &
            v2, gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
            velx(i,j,k),vely(i,j,k),velz(i,j,k),&
            cctk_iteration, x(i,j,k),y(i,j,k),z(i,j,k), whisky_reflevel, &
            Whisky_CarpetWeights(i,j,k)
          call CCTK_WARN(1,infoline)
          !$OMP END CRITICAL
        end if

      end do
    end do
  end do
  !$OMP END PARALLEL DO

end subroutine Whisky_SanityCheck
