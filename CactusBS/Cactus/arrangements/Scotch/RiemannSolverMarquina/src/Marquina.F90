  /*@@
   @file      Whisky_Marquina.f90
   @date      Thu Jan  11 11:03:32 2002
   @author    Pedro Montero, Toni Font    
   @desc 
   Routine to obtain the Marquina Fluxes. Note that this is the 
   MODIFIED Marquina formula as given by Aloy et.al. 
   (ApJ Supp 122 (1999) p.151) and not the full Marquina flux 
   of Donat and Marquina.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "SpaceMask.h"

 /*@@
   @routine    Whisky_Marquina.f90 
   @date       Wed Feb 13 11:03:32 2002
   @author     Pedro Montero, Toni Font
   @desc 
   Routine to obtain the Marquina Fluxes
   @enddesc 
   @calls     
   @calledby   
   @history 
   Based on routines by Toni Font
   @endhistory 

@@*/


subroutine Whisky_Marquina(CCTK_ARGUMENTS)
    
    implicit none

!!$ #ifdef _EOS_BASE_INC_
!!$ #undef _EOS_BASE_INC_
!!$ #endif
!!$ #include "EOS_Base.inc"

    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    DECLARE_CCTK_FUNCTIONS
    CCTK_REAL, dimension(5) :: marquinaflux, &
         consp,consm_i,fplus,fminus,f_marquina,primp,primm_i
    CCTK_REAL :: avg_alp,avg_beta,gxxh,gxyh,gxzh,gyyh,gyzh,gzzh, &
         pressp,pressm_i, &
         w_lorentzp,w_lorentzm_i
    CCTK_REAL :: alpp,betap,gxxp,gxyp,gxzp,gyyp,gyzp,gzzp, &
         detp,uxxp,uxyp,uxzp,uyyp,uyzp,uzzp,usendp
    CCTK_REAL :: alpm_i,betam_i,gxxm_i,gxym_i,gxzm_i,gyym_i,gyzm_i,gzzm_i, &
         detm_i,uxxm_i,uxym_i,uxzm_i,uyym_i,uyzm_i,uzzm_i,usendm_i
    CCTK_REAL :: v2, weight
    integer :: m
    integer :: i,j,k
    CCTK_INT :: eos_handle
    character(len=400) :: warnline
    
    CCTK_INT :: type_bits, trivial

    CCTK_INT :: alive_bits, alive, dead

    call SpaceMask_GetTypeBits(alive_bits, "Hydro_AliveCells")
    call SpaceMask_GetStateBits(alive, "Hydro_AliveCells", "alive")
    call SpaceMask_GetStateBits(dead, "Hydro_AliveCells", "dead")

    if (flux_direction == 1) then
      call SpaceMask_GetTypeBits(type_bits, "Hydro_RiemannProblemX")
      call SpaceMask_GetStateBits(trivial, "Hydro_RiemannProblemX", &
           &"trivial")
    else if (flux_direction == 2) then
      call SpaceMask_GetTypeBits(type_bits, "Hydro_RiemannProblemY")
      call SpaceMask_GetStateBits(trivial, "Hydro_RiemannProblemY", &
           &"trivial")
    else if (flux_direction == 3) then
      call SpaceMask_GetTypeBits(type_bits, "Hydro_RiemannProblemZ")
      call SpaceMask_GetStateBits(trivial, "Hydro_RiemannProblemZ", &
           &"trivial")
    else
      call CCTK_WARN(0, "Flux direction not x,y,z")
    end if
    
    !$OMP PARALLEL DO PRIVATE (marquinaflux, &
    !$OMP consp,consm_i,fplus,fminus,f_marquina,primp,primm_i, &
    !$OMP avg_alp,avg_beta,gxxh,gxyh,gxzh,gyyh,gyzh,gzzh, &
    !$OMP pressp,pressm_i, w_lorentzp,w_lorentzm_i, &
    !$OMP alpp,betap,gxxp,gxyp,gxzp,gyyp,gyzp,gzzp, &
    !$OMP detp,uxxp,uxyp,uxzp,uyyp,uyzp,uzzp,usendp, &
    !$OMP alpm_i,betam_i,gxxm_i,gxym_i,gxzm_i,gyym_i,gyzm_i,gzzm_i, &
    !$OMP detm_i,uxxm_i,uxym_i,uxzm_i,uyym_i,uyzm_i,uzzm_i,usendm_i, &
    !$OMP v2, weight, m,i,j,k, warnline, eos_handle)
    do k = whisky_stencil, cctk_lsh(3) - whisky_stencil
      do j = whisky_stencil, cctk_lsh(2) - whisky_stencil
        do i = whisky_stencil, cctk_lsh(1) - whisky_stencil

!!$       weight is used only to decide when to print a warning
          weight = max(Whisky_CarpetWeights(i,j,k),Whisky_CarpetWeights(i+xoffset,j+yoffset,k+zoffset));

!!$       Set the left (p for plus) and right (m_i for minus, i+1) states

          consp(1)   = densplus(i,j,k) 
          consp(2)   = sxplus(i,j,k)
          consp(3)   = syplus(i,j,k)
          consp(4)   = szplus(i,j,k)
          consp(5)   = tauplus(i,j,k)
          
          consm_i(1) = densminus(i+xoffset,j+yoffset,k+zoffset)
          consm_i(2) = sxminus(i+xoffset,j+yoffset,k+zoffset)
          consm_i(3) = syminus(i+xoffset,j+yoffset,k+zoffset)
          consm_i(4) = szminus(i+xoffset,j+yoffset,k+zoffset)
          consm_i(5) = tauminus(i+xoffset,j+yoffset,k+zoffset) 
          
          primp(1)   = rhoplus(i,j,k) 
          primp(2)   = velxplus(i,j,k)
          primp(3)   = velyplus(i,j,k) 
          primp(4)   = velzplus(i,j,k)
          primp(5)   = epsplus(i,j,k)
          
          primm_i(1) = rhominus(i+xoffset,j+yoffset,k+zoffset)
          primm_i(2) = velxminus(i+xoffset,j+yoffset,k+zoffset)
          primm_i(3) = velyminus(i+xoffset,j+yoffset,k+zoffset)
          primm_i(4) = velzminus(i+xoffset,j+yoffset,k+zoffset)
          primm_i(5) = epsminus(i+xoffset,j+yoffset,k+zoffset) 

          marquinaflux = 0.d0
        
!!$        Set metric terms at interface
          
          if (flux_direction == 1) then
            avg_beta = 0.5d0 * (betax(i+xoffset,j+yoffset,k+zoffset) + &
                 betax(i,j,k))
          else if (flux_direction == 2) then
            avg_beta = 0.5d0 * (betay(i+xoffset,j+yoffset,k+zoffset) + &
                 betay(i,j,k))
          else if (flux_direction == 3) then
            avg_beta = 0.5d0 * (betaz(i+xoffset,j+yoffset,k+zoffset) + &
                 betaz(i,j,k))
          else
            !$OMP CRITICAL
            call CCTK_WARN(0, "Flux direction not x,y,z")
            !$OMP END CRITICAL
            cycle ! NOTREACHED
          end if

          avg_alp = 0.5 * (alp(i,j,k) + alp(i+xoffset,j+yoffset,k+zoffset))

          gxxh = 0.5d0 * (gxx(i+xoffset,j+yoffset,k+zoffset) + &
               gxx(i,j,k))
          gxyh = 0.5d0 * (gxy(i+xoffset,j+yoffset,k+zoffset) + &
               gxy(i,j,k))
          gxzh = 0.5d0 * (gxz(i+xoffset,j+yoffset,k+zoffset) + &
               gxz(i,j,k))
          gyyh = 0.5d0 * (gyy(i+xoffset,j+yoffset,k+zoffset) + &
               gyy(i,j,k))
          gyzh = 0.5d0 * (gyz(i+xoffset,j+yoffset,k+zoffset) + &
               gyz(i,j,k))
          gzzh = 0.5d0 * (gzz(i+xoffset,j+yoffset,k+zoffset) + &
               gzz(i,j,k))

          ! copy smooth values into left and right hand variables,
          ! these are to be used if everything is ok (default behaviour)
          gxxp = gxxh;gxxm_i = gxxh
          gxyp = gxyh;gxym_i = gxyh
          gxzp = gxzh;gxzm_i = gxzh
          gyyp = gyyh;gyym_i = gyyh
          gyzp = gyzh;gyzm_i = gyzh
          gzzp = gzzh;gzzm_i = gzzh
          betap = avg_beta;betam_i = avg_beta
          alpp = avg_alp;alpm_i = avg_alp

          ! check that left state is physical
          v2 = gxxp*primp(2)*primp(2) + &
                 gyyp*primp(3)*primp(3) + gzzp*primp(4)*primp(4) + &
                 2*gxyp*primp(2)*primp(3) + 2*gxzp*primp(2) *primp(4) + &
                 2*gyzp*primp(3)*primp(4)
          if (speed_limiter.ne.0 .and. v2.ge.speed_limit) then
            consp(1)   = dens(i,j,k) 
            consp(2)   = sx(i,j,k)
            consp(3)   = sy(i,j,k)
            consp(4)   = sz(i,j,k)
            consp(5)   = tau(i,j,k)
            
            primp(1)   = rho(i,j,k) 
            primp(2)   = velx(i,j,k)
            primp(3)   = vely(i,j,k) 
            primp(4)   = velz(i,j,k)
            primp(5)   = eps(i,j,k)

            
            if (flux_direction == 1) then
              betap = betax(i,j,k)
            else if (flux_direction == 2) then
              betap = betay(i,j,k)
            else if (flux_direction == 3) then
              betap = betaz(i,j,k)
            else
              !$OMP CRITICAL
              call CCTK_WARN(0, "Flux direction not x,y,z")
              !$OMP END CRITICAL
              cycle !NOTREACHED
            end if

            alpp = alp(i,j,k)

            gxxp = gxx(i,j,k)
            gxyp = gxy(i,j,k)
            gxzp = gxz(i,j,k)
            gyyp = gyy(i,j,k)
            gyzp = gyz(i,j,k)
            gzzp = gzz(i,j,k)

          end if 
          v2 = gxxp*primp(2)*primp(2) + &
                 gyyp*primp(3)*primp(3) + gzzp*primp(4)*primp(4) + &
                 2*gxyp*primp(2)*primp(3) + 2*gxzp*primp(2) *primp(4) + &
                 2*gyzp*primp(3)*primp(4)
          if (speed_limiter.ne.0 .and. v2.ge.1.d0 .and. weight.eq.1.d0) then
            !$OMP CRITICAL
            call CCTK_WARN(1, "v2p is larger than unity even after limiting")
            write (warnline,*) 'speed: ',v2
            call CCTK_WARN(1, warnline)
            write (warnline,*) 'metric: ',gxxp,gxyp,gxzp,gyyp,gyzp,gzzp
            call CCTK_WARN(1, warnline)
            write (warnline,*) 'velocity: ',primp(2),primp(3),primp(4)
            call CCTK_WARN(1, warnline)
            write(warnline,'(a28,i2)') 'on carpet reflevel: ',whisky_reflevel
            call CCTK_WARN(1,warnline)
            write(warnline,'(a60,g15.6)') 'weight function (from CarpetReduce): ',weight
            call CCTK_WARN(1,warnline)      
            write (warnline,*) 'x, y, z: ',0.5d0*(x(i+xoffset,j+yoffset,k+zoffset)+x(i,j,k)),&
              0.5d0*(y(i+yoffset,j+yoffset,k+zoffset)+y(i,j,k)),&
              0.5d0*(z(i+zoffset,j+yoffset,k+zoffset)+z(i,j,k))
            call CCTK_WARN(1, warnline)
            !$OMP END CRITICAL
          end if

          ! check that right state is physical
          v2 = gxxm_i*primm_i(2)*primm_i(2) + &
                 gyym_i*primm_i(3)*primm_i(3) + gzzm_i*primm_i(4)*primm_i(4) + &
                 2*gxym_i*primm_i(2)*primm_i(3) + 2*gxzm_i*primm_i(2) *primm_i(4) + &
                 2*gyzm_i*primm_i(3)*primm_i(4)
          if (speed_limiter.ne.0 .and. v2.ge.speed_limit) then
            consm_i(1)   = dens(i+xoffset,j+yoffset,k+zoffset) 
            consm_i(2)   = sx(i+xoffset,j+yoffset,k+zoffset)
            consm_i(3)   = sy(i+xoffset,j+yoffset,k+zoffset)
            consm_i(4)   = sz(i+xoffset,j+yoffset,k+zoffset)
            consm_i(5)   = tau(i+xoffset,j+yoffset,k+zoffset)
            
            primm_i(1)   = rho(i+xoffset,j+yoffset,k+zoffset) 
            primm_i(2)   = velx(i+xoffset,j+yoffset,k+zoffset)
            primm_i(3)   = vely(i+xoffset,j+yoffset,k+zoffset) 
            primm_i(4)   = velz(i+xoffset,j+yoffset,k+zoffset)
            primm_i(5)   = eps(i+xoffset,j+yoffset,k+zoffset)

            
            if (flux_direction == 1) then
              betam_i = betax(i+xoffset,j+yoffset,k+zoffset)
            else if (flux_direction == 2) then
              betam_i = betay(i+xoffset,j+yoffset,k+zoffset)
            else if (flux_direction == 3) then
              betam_i = betaz(i+xoffset,j+yoffset,k+zoffset)
            else
              !$OMP CRITICAL
              call CCTK_WARN(0, "Flux direction not x,y,z")
              !$OMP END CRITICAL
              cycle !NOTREACHED
            end if

            alpm_i = alp(i+xoffset,j+yoffset,k+zoffset)

            gxxm_i = gxx(i+xoffset,j+yoffset,k+zoffset)
            gxym_i = gxy(i+xoffset,j+yoffset,k+zoffset)
            gxzm_i = gxz(i+xoffset,j+yoffset,k+zoffset)
            gyym_i = gyy(i+xoffset,j+yoffset,k+zoffset)
            gyzm_i = gyz(i+xoffset,j+yoffset,k+zoffset)
            gzzm_i = gzz(i+xoffset,j+yoffset,k+zoffset)

          end if 
          v2 = gxxm_i*primm_i(2)*primm_i(2) + &
                 gyym_i*primm_i(3)*primm_i(3) + gzzm_i*primm_i(4)*primm_i(4) + &
                 2*gxym_i*primm_i(2)*primm_i(3) + 2*gxzm_i*primm_i(2) *primm_i(4) + &
                 2*gyzm_i*primm_i(3)*primm_i(4)
          if (speed_limiter.ne.0 .and. v2.ge.1.d0 .and. weight.eq.1.d0) then
            !$OMP CRITICAL
            call CCTK_WARN(1, "v2m_i is larger than unity even after limiting")
            write (warnline,*) 'speed: ',v2
            call CCTK_WARN(1, warnline)
            write (warnline,*) 'metric: ',gxxm_i,gxym_i,gxzm_i,gyym_i,gyzm_i,gzzm_i
            call CCTK_WARN(1, warnline)
            write (warnline,*) 'velocity: ',primm_i(2),primm_i(3),primm_i(4)
            call CCTK_WARN(1, warnline)
            write(warnline,'(a28,i2)') 'on carpet reflevel: ',whisky_reflevel
            call CCTK_WARN(1,warnline)
            write(warnline,'(a60,g15.6)') 'weight function (from CarpetReduce): ',weight
            call CCTK_WARN(1,warnline)      
            write (warnline,*) 'x, y, z: ',0.5d0*(x(i+xoffset,j+yoffset,k+zoffset)+x(i,j,k)),&
              0.5d0*(y(i+yoffset,j+yoffset,k+zoffset)+y(i,j,k)),&
              0.5d0*(z(i+zoffset,j+yoffset,k+zoffset)+z(i,j,k))
            call CCTK_WARN(1, warnline)
            !$OMP END CRITICAL
          end if

          ! trash the *h variables to make sure no one uses them accidentally
          ! DEBUG
          gxxh = -1d238
          gxyh = -1d238
          gxzh = -1d238
          gyyh = -1d238
          gyzh = -1d238
          gzzh = -1d238
          avg_beta = -1d238
          avg_alp = -1d238

!!$ routine to calculate the determinant of the metric

          call SpatialDeterminant(gxxp,gxyp,gxzp,gyyp,gyzp,gzzp,detp)
          if ((detp .ne. detp .or. detp .lt. 0d0) .and. weight.eq.1.d0) then
              !$OMP CRITICAL
              write(warnline,'("Marquina: Metric determinant (p) in iteration ",i6," at ",g15.6,",",g15.6,",",&
                &g15.6," on level ",i2," is : ",g15.6, " from data g = [",&
                &5(g15.6,","),g15.6,"]")') cctk_iteration,&
                0.5d0*(x(i+xoffset,j+yoffset,k+zoffset)+x(i,j,k)),&
                0.5d0*(y(i+yoffset,j+yoffset,k+zoffset)+y(i,j,k)),&
                0.5d0*(z(i+zoffset,j+yoffset,k+zoffset)+z(i,j,k)),whisky_reflevel,detp,&
                gxxp,gxyp,gxzp,gyyp,gyzp,gzzp
              call CCTK_WARN(1, warnline)
              !$OMP END CRITICAL
          end if
          call SpatialDeterminant(gxxm_i,gxym_i,gxzm_i,gyym_i,gyzm_i,gzzm_i,detm_i)
          if ((detm_i .ne. detm_i .or. detm_i .lt. 0d0) .and. weight.eq.1.d0) then
              !$OMP CRITICAL
              write(warnline,'("Marquina: Metric determinant (m_i) in iteration ",i6," at ",g15.6,",",g15.6,",",&
                &g15.6," on level ",i2," is : ",g15.6, " from data g = [",&
                &5(g15.6,","),g15.6,"]")') cctk_iteration,&
                0.5d0*(x(i+xoffset,j+yoffset,k+zoffset)+x(i,j,k)),&
                0.5d0*(y(i+yoffset,j+yoffset,k+zoffset)+y(i,j,k)),&
                0.5d0*(z(i+zoffset,j+yoffset,k+zoffset)+z(i,j,k)),whisky_reflevel,detm_i,&
                gxxm_i,gxym_i,gxzm_i,gyym_i,gyzm_i,gzzm_i
              call CCTK_WARN(1, warnline)
              !$OMP END CRITICAL
          end if
          
    
!!$ check if left and right states are both dead or below atmosphere treshold
            if (enforce_no_atmosphere_flux .ne. 0 .and. &
                ((primm_i(1) .le. whisky_rho_min*(1.d0+atmo_tolerance) .and. &
                   primp(1) .le. whisky_rho_min*(1.d0+atmo_tolerance)) .or. &
                (SpaceMask_CheckStateBitsF90(space_mask, i,j,k,                         alive_bits, dead) .and. &
                 SpaceMask_CheckStateBitsF90(space_mask, i+zoffset,j+yoffset,k+zoffset, alive_bits, dead)))) then
              f_marquina = 0.d0
            else if (SpaceMask_CheckStateBitsF90(space_mask, i, j, k, type_bits, trivial)) then
!!$ If the Riemann problem is trivial, just calculate the fluxes from the 
!!$ left state and skip to the next cell

            if (flux_direction == 1) then
              call num_x_flux(consp(1),consp(2),consp(3),consp(4),consp(5),&
                   f_marquina(1),f_marquina(2),f_marquina(3),&
                   f_marquina(4),f_marquina(5),&
                   velxplus(i,j,k),pressplus(i,j,k),&
                   detp,alpp,betap)
            else if (flux_direction == 2) then
              call num_x_flux(consp(1),consp(3),consp(4),consp(2),consp(5),&
                   f_marquina(1),f_marquina(3),f_marquina(4),&
                   f_marquina(2),f_marquina(5),&
                   velyplus(i,j,k),pressplus(i,j,k),&
                   detp,alpp,betap)
            else if (flux_direction == 3) then
              call num_x_flux(consp(1),consp(4),consp(2),consp(3),consp(5),&
                   f_marquina(1),f_marquina(4),f_marquina(2),&
                   f_marquina(3),f_marquina(5),&
                   velzplus(i,j,k),pressplus(i,j,k),&
                   detp,alpp,betap)
            else
              !$OMP CRITICAL
              call CCTK_WARN(0, "Flux direction not x,y,z")
              !$OMP END CRITICAL
              cycle !NOTREACHED
            end if
            
          else !!! The end of this branch is right at the bottom of the routine
            
            call UpperMetric(uxxp, uxyp, uxzp, uyyp, uyzp, uzzp, &
                 detp,gxxp, gxyp, gxzp, gyyp, gyzp, gzzp)
            call UpperMetric(uxxm_i, uxym_i, uxzm_i, uyym_i, uyzm_i, uzzm_i, &
                 detm_i,gxxm_i, gxym_i, gxzm_i, gyym_i, gyzm_i, gzzm_i)
            
            if (flux_direction == 1) then
              usendp = uxxp
              usendm_i = uxxm_i
            else if (flux_direction == 2) then
              usendp = uyyp
              usendm_i = uyym_i
            else if (flux_direction == 3) then
              usendp = uzzp
              usendm_i = uzzm_i
            else
              !$OMP CRITICAL
              call CCTK_WARN(0, "Flux direction not x,y,z")
              !$OMP END CRITICAL
              cycle !NOTREACHED
            end if

!!$left state

            w_lorentzp = 1.d0 / sqrt(1.d0 - (gxxp*primp(2)*primp(2) + &
                 gyyp*primp(3)*primp(3) + gzzp*primp(4)*primp(4) + &
                 2*gxyp*primp(2)*primp(3) + 2*gxzp*primp(2) *primp(4) + &
                 2*gyzp*primp(3)*primp(4)))  
          if (speed_limiter.ne.0 .and.  w_lorentzp.ne.w_lorentzp .and. weight.eq.1.d0) then
            !$OMP CRITICAL
            call CCTK_WARN(1, "w_lorentzp is NaN, this should not happen")
            write (warnline,*) 'metric: ',gxxp,gxyp,gxzp,gyyp,gyzp,gzzp
            call CCTK_WARN(1, warnline)
            write (warnline,*) 'velocity: ',primp(2),primp(3),primp(4)
            call CCTK_WARN(1, warnline)
            write(warnline,'(a28,i2)') 'on carpet reflevel: ',whisky_reflevel
            call CCTK_WARN(1,warnline)
            write(warnline,'(a60,g15.6)') 'weight function (from CarpetReduce): ',weight
            call CCTK_WARN(1,warnline)      
            write (warnline,*) 'x, y, z: ',0.5d0*(x(i+xoffset,j+yoffset,k+zoffset)+x(i,j,k)),&
              0.5d0*(y(i+yoffset,j+yoffset,k+zoffset)+y(i,j,k)),&
              0.5d0*(z(i+zoffset,j+yoffset,k+zoffset)+z(i,j,k))
            call CCTK_WARN(1, warnline)
            !$OMP END CRITICAL
          end if
                        
!!$            pressp = EOS_Pressure(whisky_eos_handle,primp(1),primp(5))
        
!!$right state

            w_lorentzm_i = 1.d0 / sqrt(1.d0 - (gxxm_i*primm_i(2)*primm_i(2) + &
                 gyym_i*primm_i(3)*primm_i(3) + gzzm_i*primm_i(4)*primm_i(4) + &
                 2*gxym_i*primm_i(2)*primm_i(3) + &
                 2*gxzm_i*primm_i(2) *primm_i(4)+ &
                 2*gyzm_i*primm_i(3)*primm_i(4)))  
          if (speed_limiter.ne.0 .and.  w_lorentzm_i.ne.w_lorentzm_i .and. weight.eq.1.d0) then
            !$OMP CRITICAL
            call CCTK_WARN(1, "w_lorentzm_i is NaN, this should not happen")
            write (warnline,*) 'metric: ',gxxm_i,gxym_i,gxzm_i,gyym_i,gyzm_i,gzzm_i
            call CCTK_WARN(1, warnline)
            write (warnline,*) 'velocity: ',primm_i(2),primm_i(3),primm_i(4)
            call CCTK_WARN(1, warnline)
            write(warnline,'(a28,i2)') 'on carpet reflevel: ',whisky_reflevel
            call CCTK_WARN(1,warnline)
            write(warnline,'(a60,g15.6)') 'weight function (from CarpetReduce): ',weight
            call CCTK_WARN(1,warnline)      
            write (warnline,*) 'x, y, z: ',0.5d0*(x(i+xoffset,j+yoffset,k+zoffset)+x(i,j,k)),&
              0.5d0*(y(i+yoffset,j+yoffset,k+zoffset)+y(i,j,k)),&
              0.5d0*(z(i+zoffset,j+yoffset,k+zoffset)+z(i,j,k))
            call CCTK_WARN(1, warnline)
            !$OMP END CRITICAL
          end if
            
!!$            pressm_i = EOS_Pressure(whisky_eos_handle,primm_i(1),primm_i(5))
            
!!$eigenvalues and right eigenvectors, calculate the fluxes
            fplus = 0.d0
            fminus = 0.d0
            
            !!$ if both interfaces are within atmosphere, then evolve using the
            !!$ polytropic EOS used for atmosphere reset
            eos_handle = whisky_eos_handle
            if (polytropic_atmosphere .ne. 0) then
              if (primm_i(1) .le. whisky_rho_min*(1.d0+atmo_tolerance) .and. &
                  primp(1)   .le. whisky_rho_min*(1.d0+atmo_tolerance)) then
                eos_handle = whisky_polytrope_handle
              end if
            end if
            if (flux_direction == 1) then
              
              call eigenproblem_marquina(eos_handle,&
                   primm_i(1),primm_i(2), & 
                   primm_i(3),primm_i(4),primm_i(5),primp(1), &
                   primp(2),primp(3),primp(4),primp(5), &
                   gxxm_i,gxym_i,gxzm_i,gyym_i,gyzm_i,gzzm_i, &
                   usendm_i,detm_i,alpm_i,betam_i,&
                   gxxp,gxyp,gxzp,gyyp,gyzp,gzzp, &
                   usendp,detp,alpp,betap,consp(1),consp(2),&
                   consp(3), consp(4), consp(5),consm_i(1),consm_i(2), &
                   consm_i(3),consm_i(4),consm_i(5),marquinaflux(1), &
                   marquinaflux(2),marquinaflux(3),marquinaflux(4), &
                   marquinaflux(5), weight)
              
              call num_x_flux(consp(1),consp(2),consp(3),consp(4),consp(5), &
                   fplus(1),fplus(2),fplus(3),fplus(4), &
                   fplus(5),velxplus(i,j,k),pressplus(i,j,k), &
                   detp,alpp,betap)
              
              call num_x_flux(consm_i(1),consm_i(2),consm_i(3), &
                   consm_i(4),consm_i(5),fminus(1),fminus(2),fminus(3), &
                   fminus(4), fminus(5), &
                   velxminus(i+xoffset,j+yoffset,k+zoffset), &
                   pressminus(i+xoffset,j+yoffset,k+zoffset), &
                   detm_i,alpm_i,betam_i)
              
            else if (flux_direction == 2) then
              
              call eigenproblem_marquina(eos_handle,&
                   primm_i(1),primm_i(3), & 
                   primm_i(4),primm_i(2),primm_i(5),primp(1), &
                   primp(3),primp(4),primp(2),primp(5), &
                   gyym_i,gyzm_i,gxym_i,gzzm_i,gxzm_i,gxxm_i, &
                   usendm_i,detm_i,alpm_i,betam_i, &
                   gyyp,gyzp,gxyp,gzzp,gxzp,gxxp, &
                   usendp,detp,alpp,betap,consp(1),consp(3),&
                   consp(4), consp(2), consp(5),consm_i(1),consm_i(3), &
                   consm_i(4),consm_i(2),consm_i(5),marquinaflux(1), &
                   marquinaflux(3),marquinaflux(4),marquinaflux(2), &
                   marquinaflux(5), weight)
              
              call num_x_flux(consp(1),consp(3),consp(4),consp(2),consp(5), &
                   fplus(1),fplus(3),fplus(4),fplus(2), &
                   fplus(5),velyplus(i,j,k),pressplus(i,j,k), &
                   detp,alpp,betap)
              
              call num_x_flux(consm_i(1),consm_i(3),consm_i(4), &
                   consm_i(2),consm_i(5),fminus(1),fminus(3),fminus(4), &
                   fminus(2), fminus(5), &
                   velyminus(i+xoffset,j+yoffset,k+zoffset), &
                   pressminus(i+xoffset,j+yoffset,k+zoffset), &
                   detm_i,alpm_i,betam_i)
              
            else if (flux_direction == 3) then
              
              call eigenproblem_marquina(eos_handle,&
                   primm_i(1),primm_i(4), & 
                   primm_i(2),primm_i(3),primm_i(5),primp(1), &
                   primp(4),primp(2),primp(3),primp(5), &
                   gzzm_i,gxzm_i,gyzm_i,gxxm_i,gxym_i,gyym_i, &
                   usendm_i,detm_i,alpm_i,betam_i, &
                   gzzp,gxzp,gyzp,gxxp,gxyp,gyyp, &
                   usendp,detp,alpp,betap,consp(1),consp(4),&
                   consp(2), consp(3), consp(5),consm_i(1),consm_i(4), &
                   consm_i(2),consm_i(3),consm_i(5),marquinaflux(1), &
                   marquinaflux(4),marquinaflux(2),marquinaflux(3), &
                   marquinaflux(5), weight)

              call num_x_flux(consp(1),consp(4),consp(2),consp(3),consp(5), &
                   fplus(1),fplus(4),fplus(2),fplus(3), &
                   fplus(5),velzplus(i,j,k),pressplus(i,j,k),detp, &
                   alpp,betap)
              
              call num_x_flux(consm_i(1),consm_i(4),consm_i(2), &
                   consm_i(3),consm_i(5),fminus(1),fminus(4),fminus(2), &
                   fminus(3), fminus(5), &
                   velzminus(i+xoffset,j+yoffset,k+zoffset), &
                   pressminus(i+xoffset,j+yoffset,k+zoffset), &
                   detm_i,alpm_i,betam_i)
              
            else
              !$OMP CRITICAL
              call CCTK_WARN(0, "Flux direction not x,y,z")
              !$OMP END CRITICAL
              cycle !NOTREACHED
            end if
            
!!$ Marquina flux
            
            do m = 1,5
              
              f_marquina(m) = 0.5d0 * (fplus(m) + fminus(m) - marquinaflux(m))
              
            end do

          end if !!! The end of the SpaceMask check for a trivial RP.

          densflux(i,j,k) = f_marquina(1)
          sxflux(i,j,k)   = f_marquina(2)
          syflux(i,j,k)   = f_marquina(3)
          szflux(i,j,k)   = f_marquina(4)
          tauflux(i,j,k)  = f_marquina(5)
          
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    if (evolve_tracer .ne. 0) then

      !$OMP PARALLEL DO PRIVATE (i,j,k)
      do k = whisky_stencil, cctk_lsh(3) - whisky_stencil
        do j = whisky_stencil, cctk_lsh(2) - whisky_stencil
          do i = whisky_stencil, cctk_lsh(1) - whisky_stencil

            if (densflux(i, j, k) > 0.d0) then

              cons_tracerflux(i, j, k,:) = &
                   tracerplus(i, j, k,:) * &
                   densflux(i, j, k)

            else

              cons_tracerflux(i, j, k,:) = &
                   tracerminus(i + xoffset, j + yoffset, k + zoffset,:) * &
                   densflux(i, j, k)

            end if

          end do
        end do
      end do
      !$OMP END PARALLEL DO
      
    end if
    
    return
end subroutine Whisky_Marquina

