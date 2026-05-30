/*@@
@file      Whisky_EOSResetHydro.F90
@date      Sat Jan 26 01:36:57 2002
@author     Ian Hawke
@desc 
   This routine will reset the specific internal energy using the polytropic
   EOS that will be used at evolution. This is necessary if the EoS changes
   between setting up the initial data and evolving
  @enddesc 
  @@*/
  
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "Whisky_Utils.h"
 
! the intel compiler do not like one of the OMP statements below, since this is really
! not time critical, leave it out for now (Feb 2011)

 /*@@
   @routine    Whisky_EOSResetHydro
   @date       Sat Jan 26 01:38:12 2002
   @author     Ian Hawke
   @desc 
   see above
   @enddesc 
   @calls     
   @calledby   
   @history 
   
   @endhistory 

@@*/


subroutine Whisky_EoSChangeGamma(CCTK_ARGUMENTS)

  USE EOS_Polytrope_Scalars

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, k
  CCTK_REAL :: det
  CHARACTER(len=400) :: infoline
  
  CCTK_REAL :: local_gamma
  
#include "EOS_Base.inc"

!!$  Set up the fluid constants

  local_Gamma = 1.d0 + EOS_Pressure(whisky_polytrope_handle, 1.d0, 1.d0) / &
       EOS_SpecificIntEnergy(whisky_polytrope_handle, 1.d0, 1.d0)

!!$  Change the pressure and specific internal energy

  press = p_geom_factor * eos_k_cgs * &
       (rho * rho_geom_factor_inv)**local_Gamma 
  eps = press / (rho * (local_Gamma - 1.d0))

!!$  Get the conserved variables. Hardwired polytrope EoS!!!
!!$  Note that this call also sets pressure and eps
    
  !!OMP PARALLEL DO PRIVATE (i,j,k,det)
  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)
        
        det = SPATIAL_DET(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
        if ((det .ne. det .or. det .lt. 0d0) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
            !!OMP CRITICAL
            write(infoline,'("EoSChangeGamma: Metric determinant in iteration ",i6," at [",g15.6,",",g15.6,",",&
              &g15.6,"] on level ",i2," is : ",g15.6)') cctk_iteration,&
              x(i,j,k),y(i,j,k),z(i,j,k),whisky_reflevel,det
            call CCTK_WARN(1, infoline)
            !!OMP END CRITICAL
        end if
        if ( whisky_mhd_handle.gt.1 ) then
           call prim2conmhdpolytype(whisky_polytrope_handle,gxx(i,j,k),gxy(i,j,k),&
                gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
                det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
                tau(i,j,k),rho(i,j,k),velx(i,j,k),vely(i,j,k),velz(i,j,k),&
                eps(i,j,k),press(i,j,k),w_lorentz(i,j,k),bnx(i,j,k),bny(i,j,k),&
                bnz(i,j,k),bvec(i,j,k,1),bvec(i,j,k,2),bvec(i,j,k,3))
        else
           call prim2conpolytype(whisky_polytrope_handle,gxx(i,j,k),gxy(i,j,k),&
                gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
                det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
                tau(i,j,k),rho(i,j,k),velx(i,j,k),vely(i,j,k),velz(i,j,k),&
                eps(i,j,k),press(i,j,k),w_lorentz(i,j,k))
        end if
        if ((w_lorentz(i, j, k).lt.1.d0 .or. w_lorentz(i,j,k).ne.w_lorentz(i,j,k)) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
          !!OMP CRITICAL
          write(infoline,'("EoSChangeGamma: Unphysical Lorentz factor ",g15.6,&
            &" for data g = [",5(g15.6,","),g15.6,"] vel = [",&
            &2(g15.6,","),g15.6,"] occured in iteration ",i6," at location [",&
            &2(g15.6,","),g15.6,"] on reflevel ",i2," weight ",f4.2)') &
            w_lorentz(i,j,k), gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
            velx(i,j,k),vely(i,j,k),velz(i,j,k),&
            cctk_iteration, x(i,j,k),y(i,j,k),z(i,j,k), whisky_reflevel, &
            Whisky_CarpetWeights(i,j,k)
          call CCTK_WARN(1,infoline)
          !!OMP END CRITICAL
        end if
        if ((tau(i,j,k).le.0 .or. tau(i,j,k).ne.tau(i,j,k)) .and. Whisky_CarpetWeights(i,j,k).eq.1.d0) then
            !!OMP CRITICAL
            write(infoline,'("EoSChangeGamma: Unphysical tau = ",g15.6," found for data rho = ",&
              &g15.6,", eps = ",g15.6,", press = ",g15.6,", w = ",g15.6,&
              &" occured in iteration ",i6," at location [",&
              &2(g15.6,","),g15.6,"] on reflevel ",i2," weight ",f4.2)') &
              tau(i,j,k), rho(i,j,k),eps(i,j,k),press(i,j,k),&
              w_lorentz(i,j,k),&
              cctk_iteration, x(i,j,k),y(i,j,k),z(i,j,k), whisky_reflevel, &
              Whisky_CarpetWeights(i,j,k)
            call CCTK_WARN(1,infoline)
            !!OMP END CRITICAL
        end if

      end do
    end do
  end do
  !!OMP END PARALLEL DO


end subroutine Whisky_EoSChangeGamma

 /*@@
   @routine    Whisky_EoSChangeK
   @date       Mon Oct 20 12:56:14 2003
   @author     Ian Hawke
   @desc 
   Reset the hydro variables when K is changed.
   Unlike the routine above, this actually gives a solution to
   the constraints.

   Only two cases are given as the general case is transcendental.
   We find this by holding rho * enthalpy fixed and assuming a
   polytropic EOS.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Whisky_EoSChangeK(CCTK_ARGUMENTS)

  USE EOS_Polytrope_Scalars

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT :: i, j, k
  CCTK_REAL :: det
  
  CCTK_REAL :: local_gamma, local_k
  
  CCTK_REAL, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)) :: Q

#include "EOS_Base.inc"

!!$  Set up the fluid constants

  local_Gamma = 1.d0 + EOS_Pressure(whisky_polytrope_handle, 1.d0, 1.d0) / &
       EOS_SpecificIntEnergy(whisky_polytrope_handle, 1.d0, 1.d0)

  local_K = EOS_Pressure(whisky_polytrope_handle, 1.d0, 1.d0)

  if (abs(local_Gamma - 2.d0) < 1.d-10) then

    rho = -0.5d0/local_k+sqrt(0.25d0/local_k**2+(rho+initial_k*rho**2)/local_k)

  else if (abs(local_Gamma - 3.d0) < 1.d-10) then

!!$ This code is probably just wrong
    Q = -9.d0 * local_k**2 * rho * (2.d0 + 3.d0 * initial_k * rho**2) + &
         sqrt(local_k**3 * (32.d0 + 81.d0 * local_k * rho**2 * &
         (2.d0 + 3.d0 * initial_k * rho**2)**2))
    
    Q = Q**(1.d0/3.d0)
    
    rho = (2**(7.d0/3.d0) * local_k - 2**(2.d0/3.d0) * Q**2) / &
         (6.d0 * local_k * Q)

  else
    call CCTK_WARN(0, "EoSChangeK only knows how to do Gamma=2 or 3!")
  end if
  
  press = local_k * rho**local_gamma
  eps = (local_gamma - 1.d0) * local_k * rho**local_gamma

  call Primitive2ConservativePolyCells(CCTK_ARGUMENTS)

end subroutine Whisky_EoSChangeK



 /*@@
   @routine    Whisky_EoSChangeGammaK_Shibata
   @date       Jan. 2005
   @author     Christian D. Ott
   @desc 
   Reset the hydro variables when K and Gamma are changed.

   This is according to Shibata in astro-ph/0412243 in
   which he switches K and Gamma after initial data setup,
   but keeps the internal energy constant.

   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Whisky_EoSChangeGammaK_Shibata(CCTK_ARGUMENTS)

  USE EOS_Polytrope_Scalars

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT :: i, j, k
  CCTK_REAL :: det
  
  CCTK_REAL :: local_Gamma, local_k, eos_k_initial_cgs, pressmax, epsmax

  CCTK_REAL, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)) :: Q

  character(len=100) infoline

#include "EOS_Base.inc"

!!$  Set up the fluid constants

  call CCTK_INFO("Pulling the rug...by changing K and Gamma")

  local_Gamma = 1.d0 + EOS_Pressure(whisky_polytrope_handle, 1.d0, 1.d0) / &
       EOS_SpecificIntEnergy(whisky_polytrope_handle, 1.d0, 1.d0)

  local_K = EOS_Pressure(whisky_polytrope_handle, 1.d0, 1.d0)

  eos_k_initial_cgs = initial_k * rho_geom_factor**initial_Gamma / p_geom_factor

  press = (local_Gamma - 1.d0) / (initial_Gamma - 1.0d0 ) * p_geom_factor * eos_k_initial_cgs * &
  	     (rho * rho_geom_factor_inv) ** initial_Gamma

  eps = p_geom_factor * eos_k_initial_cgs * & 
     (rho * rho_geom_factor_inv) ** initial_Gamma / &
     (rho * (initial_Gamma - 1.0d0))

  pressmax=0.0d0
  epsmax=0.0d0

  !!OMP PARALLEL DO PRIVATE (i,j,k) &
  !!OMP REDUCTION (MAX : pressmax, epsmax)
  do k = whisky_stencil,cctk_lsh(3)-whisky_stencil+1
    do j = whisky_stencil,cctk_lsh(2)-whisky_stencil+1
      do i = whisky_stencil,cctk_lsh(1)-whisky_stencil+1

	if(press(i,j,k).ge.pressmax) pressmax=press(i,j,k)
	if(eps(i,j,k).ge.epsmax) epsmax=eps(i,j,k)

      enddo
    enddo
  enddo
  !!OMP END PARALLEL DO

  write(infoline,'(a20,1P1E16.7)') "P_max=",pressmax
  call CCTK_INFO(infoline)

  write(infoline,'(a20,1P1E16.7)') "eps_max=",epsmax
  call CCTK_INFO(infoline)

  !!OMP PARALLEL DO PRIVATE (i,j,k,det)
  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)
        
        det = SPATIAL_DET(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
        if ( whisky_mhd_handle.gt.1 ) then
           call CCTK_WARN(0,'Did not extend Whisky_EoSChangeGammaK_Shibata to MHD.')
        end if
        call prim2con(whisky_eos_handle,gxx(i,j,k),gxy(i,j,k),&
             gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),&
             det, dens(i,j,k),sx(i,j,k),sy(i,j,k),sz(i,j,k),&
             tau(i,j,k),rho(i,j,k),velx(i,j,k),vely(i,j,k),velz(i,j,k),&
             eps(i,j,k),press(i,j,k),w_lorentz(i,j,k))

      end do
    end do
  end do
  !!OMP END PARALLEL DO


end subroutine Whisky_EoSChangeGammaK_Shibata
