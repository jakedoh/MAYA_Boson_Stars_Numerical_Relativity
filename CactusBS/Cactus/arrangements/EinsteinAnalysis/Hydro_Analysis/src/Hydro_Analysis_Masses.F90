/*@@
  @file      Hydro_Analysis_Masses.F90
  @date      Mon Nov 5 16:44 2018
  @author    Narges Shahamat, Federico Cipolletta, Bruno Giacomazzo
             Based on the original code written by Luca Baiotti
  @desc
  Calculate the baryonic mass. Carpet compatible.
  @enddesc
@@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine Hydro_Analysis_Masses_Local(CCTK_ARGUMENTS)

 implicit none

 DECLARE_CCTK_ARGUMENTS
 DECLARE_CCTK_PARAMETERS
 DECLARE_CCTK_FUNCTIONS

 CCTK_INT :: i, j, k
 CCTK_INT :: lst, rst1, rst2, rst3, m
 CCTK_REAL :: dens

! initialize all the temporary grid-variables
 Hydro_Analysis_masses_temps = 0.d0

!we exclude ghost-points from calculation
 lst = 1
 rst1 = cctk_lsh(1)
 rst2 = cctk_lsh(2)
 rst3 = cctk_lsh(3)

!$OMP parallel do private(i,j,k,m, dens) collapse(4)
 do m = 1, restmass_masses_nr+1 !do all the masses-within-radii together; number 1 is the total rest mass.

   do k = lst, rst3
     do j = lst, rst2
       do i = lst, rst1

        if (rho(i,j,k) > restmass_rho_min) then !exclude points with densities below restmass_rho_min from mass computation 
                                                !(restmass_rho_min could be used to exclude atmosphere points)
           
           ! Call the function Hydro_Analysis_Dens(...) in order to compute dens conserved variable from primitive variables and
           ! metric components
           ! dens = w_lorentz * sqrt(det) * rho
           call Hydro_Analysis_Dens(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),rho(i,j,k),w_lorentz(i,j,k),dens)

           if (m==1) then

                  Hydro_Analysis_masses_temps(i,j,k,m) = dens

           else if (r(i,j,k) <= restmass_ref_radius_mass(m-1)) then

                Hydro_Analysis_masses_temps(i,j,k,m) = dens

           end if

         end if

       end do
     end do
   end do

 end do !end do on different radii
!$OMP end parallel do

 grid_spacing_product = product (cctk_delta_space) !must be computed in local mode and stored

end subroutine Hydro_Analysis_Masses_Local


subroutine Hydro_Analysis_Masses_Global(CCTK_ARGUMENTS)

 implicit none

 DECLARE_CCTK_ARGUMENTS
 DECLARE_CCTK_PARAMETERS

 CCTK_INT :: ierr, Reduction_Handle, VarIndex, m
 character(len=300) :: infoline

 call CCTK_ReductionHandle(Reduction_Handle, "sum")

 call CCTK_VarIndex(VarIndex, "Hydro_Analysis::Hydro_Analysis_masses_temps[0]")
 call CCTK_Reduce(ierr, cctkGH, -1, Reduction_Handle, &
      1, CCTK_VARIABLE_REAL, total_rest_mass, 1, VarIndex)
 if (ierr .ne. 0) then
   write (infoline, '("Failed to integrate Hydro_Analysis_masses_temps[0], ierr = ",i4)') ierr
   call CCTK_ERROR(infoline)
   STOP
 end if

 total_rest_mass = total_rest_mass * grid_spacing_product

 do m=1,restmass_masses_nr

   call CCTK_VarIndex(VarIndex, "Hydro_Analysis::Hydro_Analysis_masses_temps[0]")
   VarIndex = VarIndex + m
   call CCTK_Reduce(ierr, cctkGH, -1, Reduction_Handle, &
        1, CCTK_VARIABLE_REAL, Hydro_Analysis_masses(m), 1, VarIndex)
   if (ierr .ne. 0) then
     write (infoline, '("Failed to integrate Hydro_Analysis_masses_temps[",i4,"], ierr = ",i4)') m, ierr
     call CCTK_ERROR(infoline)
     STOP
   end if

   Hydro_Analysis_masses(m) = Hydro_Analysis_masses(m) * grid_spacing_product
   Hydro_Analysis_masses_fractions(m) = Hydro_Analysis_masses(m) / total_rest_mass

 end do

end subroutine Hydro_Analysis_Masses_Global


subroutine Hydro_Analysis_Dens(gxx,gxy,gxz,gyy,gyz,gzz,rho,w_lorentz,dens)

  implicit none

  CCTK_REAL :: dens
  CCTK_REAL :: gxx,gxy,gxz,gyy,gyz,gzz,rho,w_lorentz

  dens = rho * w_lorentz * sqrt(-gxz**2*gyy + 2*gxy*gxz*gyz - gxx*gyz**2 - gxy**2*gzz + gxx*gyy*gzz)

  return

end subroutine Hydro_Analysis_Dens
