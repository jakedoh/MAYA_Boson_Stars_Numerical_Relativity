/*@@
   @file      Pointwise.c
   @date      Sat Jan 26 01:06:01 2002
   @author    The Whisky Developers
   @desc 
   The routines for converting conservative to primitive variables at specific points.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "SpaceMask.h"

/*@@
@routine    Con2Prim_pt
@date       Sat Jan 26 01:09:39 2002
@author     The Whisky Develpoers    
@desc 
Given all the appropriate data, recover the primitive variables at
a single point.
@enddesc 
 @calls     
@calledby   
@history 
@endhistory 

@@*/

subroutine Con2Prim_pt(handle, dens, sx, sy, sz, tau, rho, velx, vely, &
     velz, epsilon, press, w_lorentz, uxx, uxy, uxz, uyy, &
     uyz, uzz, det, x, y, z, r, epsnegative,whisky_rho_min,pmin, &
     whisky_reflevel, weight, cctk_iteration)
  
  use Recovery_Scalars

  implicit none
  
  DECLARE_CCTK_PARAMETERS
  
  CCTK_REAL dens, sx, sy, sz, tau, rho, velx, vely, velz, epsilon, &
       press, uxx, uxy, uxz, uyy, uyz, uzz, det, w_lorentz, x, &
       y, z, r, whisky_rho_min
  
  CCTK_REAL s2, c0, c1, c2, c3, c4, f, df, ftol, v2, w, vlowx, vlowy, &
       vlowz 
  CCTK_INT count, i, handle, whisky_reflevel, cctk_iteration
  CCTK_REAL udens, usx, usy, usz, utau, pold, pnew, epsold, epsnew, w2, &
       w2mhalf, temp1,drhobydpress,depsbydpress,dpressbydeps,dpressbydrho,&
       pmin, weight, rhoold, plow
  CCTK_REAL flimit,maxs2,rr
  logical do_impose_s2_limit
  character(len=400) warnline
  logical epsnegative

#ifdef _EOS_BASE_INC_
#undef _EOS_BASE_INC_
#endif
#include "EOS_Base.inc"
  
  do_impose_s2_limit = impose_s2_limit.ne.0

!!$ are we inside an AH and are asked to limit s2 there?
  if (puncture_trick.ne.0) then
    do i=1,number_of_punctures 
      if ( radii(i).ge.0 ) then
        rr = (x - center_x(i))**2 + (y - center_y(i))**2 + (z-center_z(i))**2
        if ( rr.le.radii(i)**2 ) then
          do_impose_s2_limit = .true.
          exit 
        end if
      end if
    end do
  end if

!!$  Undensitize the variables 

  udens = dens /sqrt(det)
  usx = sx /sqrt(det)
  usy = sy /sqrt(det)
  usz = sz /sqrt(det)
  utau = tau /sqrt(det)
  s2 = usx*usx*uxx + usy*usy*uyy + usz*usz*uzz + 2.*usx*usy*uxy + &
       2.*usx*usz*uxz + 2.*usy*usz*uyz

  ! used to restore rho should the recovery fail
  rhoold = rho 

  if (do_impose_s2_limit) then
    if ( tau.le.0d0 .and. puncture_tau_min.ne.-2d0 ) then ! is tau nonsensical, do we care?
      !if(weight.eq.1.d0) then ! only print warning if on finest level
      !  write(warnline, '("Unphysical tau encountered=",g15.7," at (",3g15.7,")")')  &
      !       utau, x,y,z
      !  call CCTK_WARN(1,warnline)
      !end if
    
      if (puncture_tau_min.eq.-1d0) then
        tau = sqrt(det) * puncture_utau_min ! utau_min in Recovery_Scalars 
        utau = puncture_utau_min
      else
        tau = puncture_tau_min
        utau = puncture_tau_min /sqrt(det)
      end if
    end if

    maxs2 = recovery_flimit * utau * (utau + 2.0d0*udens)
    if ( s2.gt.maxs2 ) then

        ! limit appears twice
        flimit = sqrt(recovery_flimit * maxs2/s2)
        if (flimit.ne.flimit .and. weight.eq.1.d0 .and. udens.gt.whisky_rho_min) then
           write(warnline, '("Scale factor is NaN! (utau, udens, s2, maxs2)=(",4g15.7,") at (",3g15.7,")")')  &
                 utau, udens, s2, maxs2, x,y,z
           call CCTK_WARN(1,warnline)
           
        end if
        sx = flimit*sx
        sy = flimit*sy
        sz = flimit*sz
        usx = flimit*usx
        usy = flimit*usy
        usz = flimit*usz
        s2 = flimit**2*s2

    end if
  end if

!!$  Set initial guess for pressure:

  pold = max((1.d0+press_tolerance)*pmin,EOS_Pressure(handle, rho, epsilon)) ! start a little bit above 

!!$  Check that the variables have a chance of being physical

  if( (utau + pold + udens)**2 - s2 .le. 0.0d0) then

    ! Frankly I don't know why they call it pressure_reset, given that it
    ! resets only the initial guess (Roland)
    if (reset_pressure .ne. 0) then
      pold = pmin + (sqrt(s2) - utau - udens)
    else 
      if (weight .eq. 1d0) then
        call CCTK_WARN(1,'Con2Prim_pt: physical variables unphysical!')
        call CCTK_WARN(1,'Even with mesh refinement, this point is not restricted from a finer level, so this is really an error')
        Write(warnline,'(a14,i8)') 'in iteration: ',cctk_iteration
        call CCTK_WARN(1,warnline)
        Write(warnline,'(a28,i2)') 'on carpet reflevel: ',whisky_reflevel
        call CCTK_WARN(1,warnline)
        write(warnline,'(a60,g15.6)') 'weight function (from CarpetReduce): ',weight
        call CCTK_WARN(1,warnline)      
        write(warnline,'(a25,g15.6)') 'undensitized energy: ',utau
        call CCTK_WARN(1,warnline)
        write(warnline,'(a25,g15.6)') 'undensitized density: ',udens
        call CCTK_WARN(1,warnline)
        write(warnline,'(a25,g15.6)') 'undensitized s2: ',s2
        call CCTK_WARN(1,warnline)
        write(warnline,'(a25,g15.6)') 'pressure guess: ',pold
        call CCTK_WARN(1,warnline)
        write(warnline,'(a25,4g15.6)') 'coordinates: x,y,z,r:',x,y,z,r
        call CCTK_WARN(1,warnline)
        call CCTK_WARN(warnlevel, "Unphysical variables")
      endif

      rho = rhoold         ! not really needed here yet since rho.eq.rhoold anyway
      epsnegative = .true. ! Conservative2Primitive takes this as an indication
                           ! to use the polytropic con2prim 
      return
    end if
  endif

!!$ smallest pressure we accept (physical & larger than atmosphere) 
  plow = max(pmin, sqrt(s2) - utau - udens)

!!$  Calculate rho and epsilon 

  rho = udens * sqrt( (utau + pold + udens)**2 - s2)/(utau + pold + udens)
  w_lorentz = 1.d0 / sqrt( 1.d0 - s2/((utau + pold + udens)**2) )
  epsilon = (sqrt( (utau + pold + udens)**2 - s2) - pold * w_lorentz - &
       udens)/udens
  
!!$  Calculate the function

  f = pold - EOS_Pressure(handle, rho, epsilon)

!!$Find the root
  
  count = 0
  pnew = pold
  ! reset what must be atmosphere (since w_lorentz .ge. 1) to atmosphere
  ! this has to be after pnew is set because pnew is used below
  if (udens .le. whisky_rho_min*(1.d0+atmo_tolerance)) then
    rho = whisky_rho_min ! triggers atmosphere reset at label 50
    goto 50
  end if
  do while ( ((abs(pnew - pold)/abs(pnew) .gt. whisky_perc_ptol) .and. &
       (abs(pnew - pold) .gt. whisky_del_ptol))  .or. &
       (count .lt. countmin))
    count = count + 1
    if (count > countmax) then
      ! print warnings only on a level for which this is required
      if (whisky_reflevel.ge.warn_from_reflevel .and. weight.eq.1.d0) then
        call CCTK_WARN(1, 'Con2Prim_pt(press): error: did not converge in ')
        write(warnline,'(a20,g20.7,a10)') '              ',countmax,' steps'
        call CCTK_WARN(1,warnline)
        call CCTK_WARN(1,'Even with mesh refinement, this point is not restricted from a finer level, so this is really an error')
        Write(warnline,'(a14,i8)') 'in iteration: ',cctk_iteration
        call CCTK_WARN(1,warnline)
        write(warnline,'(a28,i2)') 'on carpet reflevel: ',whisky_reflevel
        call CCTK_WARN(1,warnline)
        write(warnline,'(a60,g15.6)') 'weight function (from CarpetReduce): ',weight
        call CCTK_WARN(1,warnline)      
        write(warnline,'(a20,g20.7)') 'de is ',dens
        call CCTK_WARN(1,warnline)
        write(warnline,'(a20,g20.7)') 'sx is ',sx
        call CCTK_WARN(1,warnline)
        write(warnline,'(a20,g20.7)') 'sy is ',sy
        call CCTK_WARN(1,warnline)
        write(warnline,'(a20,g20.7)') 'sz is ',sz
        call CCTK_WARN(1,warnline)
        write(warnline,'(a20,g20.7)') 'en is ',tau
        call CCTK_WARN(1,warnline)
        write(warnline,'(a20,g20.7)') 'press func: ',f
        call CCTK_WARN(1,warnline)
        write(warnline,'(a20,g20.7)') 'deriv of press func: ',df
        call CCTK_WARN(1,warnline)
        write(warnline,'(a20,g20.7)') 'pnew: ',pnew
        call CCTK_WARN(1,warnline)
        write(warnline,'(a20,g20.7)') 'pold: ',pold
        call CCTK_WARN(1,warnline)
        write(warnline,'(a20,g20.7)') 'rho: ',rho
        call CCTK_WARN(1,warnline)
        write(warnline,'(a20,g20.7)') 'epsilon: ',epsilon
        call CCTK_WARN(1,warnline)
        write(warnline,'(a20,g20.7)') 'w_lorentz: ',w_lorentz
        call CCTK_WARN(1,warnline)        
        write(warnline,'(a20,3g20.7)') 'xyz location: ',x,y,z
        call CCTK_WARN(1,warnline)
        write(warnline,'(a20,g20.7)') 'radius: ',r
        call CCTK_WARN(1,warnline)
        call CCTK_WARN(1, "Did not converge for General EOS")
      end if 

      epsnegative = .true. ! Conservative2Primitive takes this as an indication
                           ! to use the polytropic con2prim 
      rho = rhoold ! restore rho to previous value so that the polytropic EoS has
                   ! good starting point, press and epsilon are still nonsensical
      return

    endif

    dpressbydrho = EOS_DPressByDRho(handle,rho,epsilon)
    dpressbydeps = EOS_DPressByDEps(handle,rho,epsilon)
    temp1 = (utau+udens+pnew)**2 - s2
    drhobydpress = udens * s2 / (sqrt(temp1)*(udens+utau+pnew)**2)
    depsbydpress = pnew * s2 / (rho * (udens + utau + pnew) * temp1)
    df = 1.0d0 - dpressbydrho*drhobydpress - &
         dpressbydeps*depsbydpress

    pold = pnew
!!$    pnew = max(pold - f/df, pmin)
    pnew = pold - f/df
    
    ! if the new value would be negtative, do a fake bisection step with 
    ! plow as the left-hand side
    temp1 = (utau+udens+pnew)**2 - s2
    if (pnew .le. plow .or. temp1 .le. 0d0) then
      pnew = (pold+plow) / 2.d0

      ! check that this did bring us into the allowed region
      ! if not fails try exiting the loop with the old value and hope for the best
      temp1 = (utau+udens+pnew)**2 - s2
      if (pnew .le. plow .or. temp1 .lt. 0d0) then
        if (weight.eq.1.d0 .and. pnew.ne.plow) then !! Only a problem if haven't bottomed out at atmosphere 
          write(warnline,'("Con2Prim_pt: pressure invalid even after pressure reset &
            & in iteration ",i8,&
            &" (numerical accuracy exhausted) temp1 = ",g15.6,&
            &" for data utau = ",g15.6,", pnew = ",g15.6,", udens = ",g15.6,&
            &", s2 = ",g15.6," detg = ",g15.6,", pmin = ",g15.6,&
            &", plow = ",g15.6," pold = ",g15.6," at (",3g15.6,")")') &
            cctk_iteration, temp1, utau, pnew, udens, s2, sqrt(det), pmin,plow,pold, x,y,z
          call CCTK_WARN(1,warnline)
        end if

        ! give up and go with whatever we currently have found
        pnew = pold ! this causes us to exit the loop

        ! check for disasters
        temp1 = (utau+udens+pnew)**2 - s2
        if (temp1 .le. 0d0) then
          write(warnline,'("Con2Prim_pt: could not find *any* good value for pnew &
            & in iteration ",i8,&
            &" (numerical accuracy exhausted) ",g15.6,&
            &" for data utau = ",g15.6,", pnew = ",g15.6,", udens = ",g15.6,&
            &", s2 = ",g15.6," detg = ",g15.6,", pmin = ",g15.6,&
            &", plow = ",g15.6," pold = ",g15.6," at (",3g15.6,")")') &
            cctk_iteration, temp1, utau, pnew, udens, s2, sqrt(det), pmin,plow,pold, x,y,z
          call CCTK_WARN(0,warnline)
        end if
      end if
    end if

!!$    Recalculate primitive variables and function (let's hope these are never
!!$    negative or otherwise unphysical)
       
    temp1 = (utau+udens+pnew)**2 - s2
    rho = udens * sqrt(temp1)/(utau + pnew + udens)
    w_lorentz = 1.d0/sqrt( 1.d0 - s2/((utau + pnew + udens)**2))
    epsilon = (sqrt(temp1) - pnew * w_lorentz - udens)/udens

    f = pnew - EOS_Pressure(handle, rho, epsilon)

  enddo
  
  !!$  Polish the root? 

  do i=1,whisky_polish

    dpressbydrho = EOS_DPressByDRho(handle,rho,epsilon)
    dpressbydeps = EOS_DPressByDEps(handle,rho,epsilon)
    temp1 = (utau+udens+pnew)**2 - s2
    drhobydpress = udens * s2 / (sqrt(temp1)*(udens+utau+pnew)**2)
    depsbydpress = pnew * s2 / (rho * (udens + utau + pnew) * temp1)
    df = 1.0d0 - dpressbydrho*drhobydpress - &
         dpressbydeps*depsbydpress
    pold = pnew
    pnew = pold - f/df
    
!!$    Recalculate primitive variables and function

    rho = udens * sqrt( (utau + pnew + udens)**2 - s2)/(utau + pnew + udens)
    w_lorentz = 1.d0 / sqrt( 1.d0 - s2/(utau + pnew + udens)**2 )
    epsilon = (sqrt( (utau + pnew + udens)**2 - s2) - pnew * w_lorentz - &
         udens)/udens

    f = pold - EOS_Pressure(handle, rho, epsilon) ! this is pnew in the iteration above

  enddo
      
!!$  Calculate primitive variables from root

50 continue
  if (rho .le. whisky_rho_min*(1.d0+atmo_tolerance) ) then
    if (use_polytropic_atmosphere .ne. 0) then
      ! set minimum density as for polytype
      ! since we do not have the polytype-handle available we rely on:
      ! Gamma_general = Gamma_polytype and pmin = pmin_polytype
      rho = whisky_rho_min
      pnew = pmin
      sx = 0.d0
      sy = 0.d0
      sz = 0.d0
      s2 = 0.d0
      usx = 0.d0
      usy = 0.d0
      usz = 0.d0
      w_lorentz = 1.d0
      epsilon = EOS_SpecificIntEnergy(handle, rho, pnew)
      !tau = sqrt(det) * ((rho * enthalpy) * w_lorentz**2 - press) - dens
      ! hand-crafted expression to involve only positive terms
      tau = sqrt(det) * rho*epsilon
      utau = rho*epsilon
      ! It is a bad thing to have a consistent dens, apparently (try anyway)
      dens = sqrt(det) * rho
      udens = rho
    else
      ! NB: does not reset tau or p (Roland)
      rho = whisky_rho_min
      udens = rho
      dens = sqrt(det) * rho
      epsilon = (sqrt( (utau + pnew + udens)**2) - pnew - &
           udens)/udens ! = utau/udens
      sx = 0.d0
      sy = 0.d0
      sz = 0.d0
      s2 = 0.d0
      usx = 0.d0
      usy = 0.d0
      usz = 0.d0
      w_lorentz = 1.d0
    end if
  end if

  press = pnew
  vlowx = usx / ( (rho + rho*epsilon + press) * w_lorentz**2)
  vlowy = usy / ( (rho + rho*epsilon + press) * w_lorentz**2)
  vlowz = usz / ( (rho + rho*epsilon + press) * w_lorentz**2)
  velx = uxx * vlowx + uxy * vlowy + uxz * vlowz
  vely = uxy * vlowx + uyy * vlowy + uyz * vlowz
  velz = uxz * vlowx + uyz * vlowy + uzz * vlowz

!!$If all else fails, use the polytropic EoS

  if(epsilon .lt. 0.0d0 .or. epsilon.ne.epsilon) then
    epsnegative = .true.
    rho = rhoold ! restore rho to previous value so that the polytropic EoS has
                 ! good starting point, press and epsilon are still nonsensical
    return
  endif

!!$ Checks physicality of *incoming* rho ...
  if ((rho.lt.0.d0 .or. rho.ne.rho) .and. weight.eq.1.d0) then
    write(warnline,'("Con2Prim_pt: Unphysical density ",g15.6,&
      &" in iteration ",i8,&
      &" for data utau = ",g15.6,", pnew = ",g15.6,", udens = ",g15.6,&
      &", s2 = ",g15.6," detg = ",g15.6,", pmin = ",g15.6,&
      &" at (",3g15.6,")")') &
      rho, cctk_iteration, utau, pnew, udens, s2, sqrt(det), pmin, x,y,z
    call CCTK_WARN(1,warnline)
  end if 
  if ((w_lorentz.lt.1.d0 .or. w_lorentz.ne.w_lorentz) .and. weight.eq.1.d0) then
    write(warnline,'("Con2Prim_pt: Unphysical Lorentz factor 1+(",g15.6,&
      &") in iteration ",i8,&
      &" for data utau = ",g15.6,", pnew = ",g15.6,", udens = ",g15.6,&
      &", s2 = ",g15.6," detg = ",g15.6,", pmin = ",g15.6,&
      &" at (",3g15.6,")")') &
      w_lorentz-1.d0, cctk_iteration, utau, pnew, udens, s2, sqrt(det), pmin, x,y,z
    call CCTK_WARN(1,warnline)
  end if 
  if ((epsilon.ne.epsilon) .and. weight.eq.1.d0) then ! epsilon.lt.0.d0 is handled by the polytrope
    write(warnline,'("Con2Prim_pt: Unphysical internal energy ",g15.6,&
      &" in iteration ",i8,&
      &" for data utau = ",g15.6,", pnew = ",g15.6,", udens = ",g15.6,&
      &", s2 = ",g15.6," w_lorentz = ",g15.6," detg = ",g15.6,&
      &", pmin = ",g15.6," at (",3g15.6,")")') &
      epsilon, cctk_iteration, utau, pnew, udens, s2, w_lorentz, sqrt(det), pmin, x,y,z
    call CCTK_WARN(1,warnline)
  end if 

  ! this should never happen
  if(velx*vlowx+vely*vlowy+velz*vlowz.gt.warn_if_speed_above .and. weight.eq.1.d0) then
    write(warnline,'(a60)') 'Caught speed above RecoverMarti::warn_if_speed_above'
    call CCTK_WARN(1,warnline)
    write(warnline,'(a14,i8)') 'in iteration: ',cctk_iteration
    call CCTK_WARN(1,warnline)
    write(warnline,'(a60,g15.6)') 'weight function (from CarpetReduce): ',weight
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,3g16.7)') 'xyz location: ',x,y,z
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,g16.7)') 'speed: ',velx*vlowx+vely*vlowy+velz*vlowz
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,g16.7)') 'w_lorentz: ',w_lorentz
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,g16.7)') 'rho: ',rho
    call CCTK_WARN(1,warnline)
  end if

  if (velx*vlowx+vely*vlowy+velz*vlowz .gt. 1.0d0 .and. weight.eq.1.d0) then
    call CCTK_WARN(1,'Con2prim_pt: stopping the code.')
    call CCTK_WARN(1, '   velocity is larger than c! ')
    call CCTK_WARN(1,'Even with mesh refinement, this point is not restricted from a finer level, so this is really an error')
    write(warnline,'(a14,i8)') 'in iteration: ',cctk_iteration
    call CCTK_WARN(1,warnline)
    write(warnline,'(a28,i2)') 'on carpet reflevel: ',whisky_reflevel
    call CCTK_WARN(1,warnline)
    write(warnline,'(a60,g15.6)') 'weight function (from CarpetReduce): ',weight
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,3g16.7)') 'xyz location: ',x,y,z
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,g16.7)') 'radius: ',r
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,3g16.7)') 'velocities: ',velx,vely,velz
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,g16.7)') 'speed: ',velx*vlowx+vely*vlowy+velz*vlowz
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,g16.7)') 'rho: ',rho
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,g16.7)') 'epsilon: ',epsilon
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,g16.7)') 'pressure: ',press
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,g16.7)') 'w_lorentz: ',w_lorentz
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,g16.7)') 'df: ',df
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,g16.7)') 'det: ',det
    call CCTK_WARN(1,warnline)
!!$      write(warnline,'(a60,3g16.7)') 'uxx, uxy, uxz: ',uxx, uxy, uxz
!!$      call CCTK_WARN(1,warnline)
!!$      write(warnline,'(a60,3g16.7)') 'uyy, uyz, uzz: ',uyy, uyz, uzz
!!$      call CCTK_WARN(1,warnline)
    call CCTK_WARN(warnlevel, "speed larger than c")
    return
  endif
  
end subroutine Con2Prim_pt

/*@@
@routine    Con2Prim_ptPolytype
@date       Sat Jan 26 01:09:39 2002
@author     The Whisky Developers   
@desc 
Given all the appropriate data, recover the primitive variables at
a single point.
@enddesc 
 @calls     
@calledby   
@history 

@endhistory 

@@*/

subroutine Con2Prim_ptPolytype(handle, dens, sx, sy, sz, tau, rho, &
     velx, vely, velz, epsilon, press, w_lorentz, uxx, uxy, uxz, uyy, &
     uyz, uzz, det, x, y, z, r, whisky_rho_min, whisky_reflevel, weight, &
     cctk_iteration)
  
  use Recovery_Scalars

  implicit none
  
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL dens, sx, sy, sz, tau, rho, velx, vely, velz, epsilon, &
       press, uxx, uxy, uxz, uyy, uyz, uzz, det, w_lorentz, x, &
       y, z, r, whisky_rho_min
  
  CCTK_REAL s2, f, df, ftol, vlowx, vlowy, vlowz, weight
  CCTK_INT count, i, handle, whisky_reflevel, cctk_iteration
  CCTK_REAL udens, usx, usy, usz, rhoold, rhonew, epsold, epsnew, w2, &
       w2mhalf, enthalpy, denthalpy
  CCTK_REAL utau,flimit,maxs2,rr
  logical s2_limited, do_impose_s2_limit
  character(len=400) warnline

#ifdef _EOS_BASE_INC_
#undef _EOS_BASE_INC_
#endif
#include "EOS_Base.inc"

  s2_limited = .false.
  do_impose_s2_limit = impose_s2_limit.ne.0

!!$ are we inside an AH and are asked to limit s2 there?
  if (puncture_trick.ne.0) then
    do i=1,number_of_punctures 
      if ( radii(i).ge.0 ) then
        rr = (x - center_x(i))**2 + (y - center_y(i))**2 + (z-center_z(i))**2
        if ( rr.le.radii(i)**2 ) then
          do_impose_s2_limit = .true.
          exit 
        end if
      end if
    end do
  end if

!!$  Undensitize the variables 


10 udens = dens /sqrt(det)
  usx = sx /sqrt(det)
  usy = sy /sqrt(det)
  usz = sz /sqrt(det)
  s2 = usx*usx*uxx + usy*usy*uyy + usz*usz*uzz + 2.d0*usx*usy*uxy + &
       2.d0*usx*usz*uxz + 2.d0*usy*usz*uyz

!!$  Set initial guess for rho:

  rhoold = max(whisky_rho_min,rho)
  
!!$  Calculate w_lorentz etc. 

  enthalpy = 1.d0 + EOS_SpecificIntEnergy(handle, rhoold, press) + &
       EOS_Pressure(handle, rhoold, epsilon) / rhoold
  w_lorentz = sqrt(1.d0 + s2 / ( (udens*enthalpy)**2 ))
  press = EOS_Pressure(handle, rhoold, epsilon)

!!$  Calculate the function

  f = rhoold * w_lorentz - udens

!!$Find the root
  
  count = 0
  rhonew = rhoold
  do while ( ((abs(rhonew - rhoold)/abs(rhonew) .gt. whisky_perc_ptol) .and. &
       (abs(rhonew - rhoold) .gt. whisky_del_ptol))  .or. &
       (count .lt. countmin))
    count = count + 1
    if (count > countmax) then
      ! print warnings only on a level for which this is required
      if (whisky_reflevel.ge.warn_from_reflevel .and. weight.eq.1.d0) then
        call CCTK_WARN(1, 'Con2Prim_ptPolytype(press): error: did not converge in ')
        write(warnline,'(a14,i4,a6)') '              ',countmax,' steps'
        call CCTK_WARN(1,warnline)
        call CCTK_WARN(1,'Even with mesh refinement, this point is not restricted from a finer level, so this is really an error')
        write(warnline,'(a14,i8)') 'in iteration: ',cctk_iteration
        call CCTK_WARN(1,warnline)
        write(warnline,'(a28,i2)') 'on carpet reflevel: ',whisky_reflevel
        call CCTK_WARN(1,warnline)
        write(warnline,'(a60,g15.6)') 'weight function (from CarpetReduce): ',weight
        call CCTK_WARN(1,warnline)      
        write(warnline,'(a6,es20.12)') 'de is ',dens
        call CCTK_WARN(1,warnline)
        write(warnline,'(a6,es20.12)') 'sx is ',sx
        call CCTK_WARN(1,warnline)
        write(warnline,'(a6,es20.12)') 'sy is ',sy
        call CCTK_WARN(1,warnline)
        write(warnline,'(a6,es20.12)') 'sz is ',sz
        call CCTK_WARN(1,warnline)
        write(warnline,'(a10,es20.12)') 'rho func: ',f
        call CCTK_WARN(1,warnline)
        write(warnline,'(a19,es20.12)') 'deriv of rho func: ',df
        call CCTK_WARN(1,warnline)
        write(warnline,'(a8,es20.12)') 'rhonew: ',rhonew
        call CCTK_WARN(1,warnline)
        write(warnline,'(a8,es20.12)') 'rhoold: ',rhoold
        call CCTK_WARN(1,warnline)
        write(warnline,'(a7,es20.12)') 'press: ',press
        call CCTK_WARN(1,warnline)
        write(warnline,'(a9,es20.12)') 'epsilon: ',epsilon
        call CCTK_WARN(1,warnline)
        write(warnline,'(a11,es20.12)') 'w_lorentz: ',w_lorentz
        call CCTK_WARN(1,warnline)
        write(warnline,'(a14,3es20.12)') 'xyz location: ',x,y,z
        call CCTK_WARN(1,warnline)
        write(warnline,'(a8,es20.12)') 'radius: ',r
        call CCTK_WARN(1,warnline)
        call CCTK_WARN(warnlevel, "Did not converge")
        return
      else 
        write(warnline,'(a60,i2)') 'Con2Prim: eps negative, but I was told to ignore level: ',whisky_reflevel
        call CCTK_WARN(1,warnline)
        return
      endif
    endif

!!$    I have the feeling that this is an error in general. But for the
!!$    2D_Polytrope it gives the right answers.

    denthalpy = EOS_DPressByDrho(handle, rhonew, epsilon) / rhonew

    df = w_lorentz - rhonew * s2 * denthalpy / &
         (w_lorentz*(udens**2)*(enthalpy**3))

    rhoold = rhonew
    rhonew = rhoold - f/df

    !!$ impose positivity on rho
    if (rhonew .le. 0.d0) then
      ! if the new value would be negtative, do a fake bisection step with 
      ! 0d0 as the left-hand side
      rhonew = rhoold / 2.d0
    end if

!!$    Recalculate primitive variables and function
       
    enthalpy = 1.d0 + EOS_SpecificIntEnergy(handle, rhonew, press) + &
         EOS_Pressure(handle, rhonew, epsilon) / rhonew
    w_lorentz = sqrt(1.d0 + s2 / ( (udens*enthalpy)**2 ))
    press = EOS_Pressure(handle, rhonew, epsilon)
    !tau = sqrt(det) * ((rho * enthalpy) * w_lorentz**2 - press) - dens
    
    f = rhonew * w_lorentz - udens

  enddo

!!$  Calculate primitive variables from root

  if (rhonew .le. whisky_rho_min*(1.d0+atmo_tolerance) ) then
    rhonew = whisky_rho_min
    sx = 0.d0
    sy = 0.d0
    sz = 0.d0
    s2 = 0.d0
    usx = 0.d0
    usy = 0.d0
    usz = 0.d0
!!$    It is a bad thing to have a consistent dens, apparently (try anyway)
    dens  = sqrt(det) * rhonew
    udens = rhonew
  
  end if 

  rho = rhonew

  enthalpy = 1.d0 + EOS_SpecificIntEnergy(handle, rhonew, press) + &
       EOS_Pressure(handle, rhonew, epsilon) / rhonew
  w_lorentz = sqrt(1.d0 + s2 / ( (udens*enthalpy)**2 ))
  press = EOS_Pressure(handle, rhonew, epsilon)
  epsilon = EOS_SpecificIntEnergy(handle, rhonew, press)
  !tau = sqrt(det) * ((rho * enthalpy) * w_lorentz**2 - press) - dens
  ! hand-crafted expression to involve only positive terms
  tau = sqrt(det) * (rho*(w_lorentz-1.d0)*w_lorentz &
    + rho*epsilon*w_lorentz**2 + press*(w_lorentz**2-1.d0))
  vlowx = usx / ( (rho + rho*epsilon + press) * w_lorentz**2)
  vlowy = usy / ( (rho + rho*epsilon + press) * w_lorentz**2)
  vlowz = usz / ( (rho + rho*epsilon + press) * w_lorentz**2)
  velx = uxx * vlowx + uxy * vlowy + uxz * vlowz
  vely = uxy * vlowx + uyy * vlowy + uyz * vlowz
  velz = uxz * vlowx + uyz * vlowy + uzz * vlowz
  
  if (do_impose_s2_limit) then

    utau = (rhonew*(w_lorentz-1.d0)*w_lorentz &
      + rhonew*EOS_SpecificIntEnergy(handle, rhonew, 0.d0)*w_lorentz**2 + press*(w_lorentz**2-1.d0))
    if (.not.s2_limited) then 
      maxs2 = recovery_flimit * utau * (utau + 2.0d0*udens)
    else
      maxs2 =                   utau * (utau + 2.0d0*udens)
    end if
    if ( s2.gt.maxs2 ) then
        if(s2_limited .and. weight.eq.1.d0) then 

          ! we have been here before and it obviously did not help
          call CCTK_WARN(1,'Con2prim_ptPolytype: s2 limit could not be imposed,
          continuing anyway.')
          call CCTK_WARN(1, '   s2 is too large! ')
          call CCTK_WARN(1,'Even with mesh refinement, this point is not restricted from a finer level, so this is really an error')
          write(warnline,'(a14,i8)') 'in iteration: ',cctk_iteration
          call CCTK_WARN(1,warnline)
          write(warnline,'(a20,i2)') 'on carpet reflevel: ',whisky_reflevel
          call CCTK_WARN(1,warnline)
          write(warnline,'(a20,g15.6)') 'weight function: ',weight
          call CCTK_WARN(1,warnline)
          write(warnline,'(a20,3g16.7)') 'xyz location: ',x,y,z
          call CCTK_WARN(1,warnline)
          write(warnline,'(a20,g16.7)') 'radius: ',r
          call CCTK_WARN(1,warnline)
          write(warnline,'(a20,3g16.7)') 'currents: ',sx,sy,sz
          call CCTK_WARN(1,warnline)
          write(warnline,'(a20,g16.7)') 's2: ',s2
          call CCTK_WARN(1,warnline)
          write(warnline,'(a20,g16.7)') 'maxs2: ',maxs2
          call CCTK_WARN(1,warnline)
          write(warnline,'(a20,g16.7)') 'tau: ',sqrt(det)*utau
          call CCTK_WARN(1,warnline)
          write(warnline,'(a20,g16.7)') 'dens: ',dens
          call CCTK_WARN(1,warnline)
          write(warnline,'(a20,g16.7)') 'det: ',det
          call CCTK_WARN(1,warnline)
          write(warnline,'(a20,g16.7)') 'rhonew: ',rhonew
          call CCTK_WARN(1,warnline)
          write(warnline,'(a20,g16.7)') 'rhoold: ',rhoold
          call CCTK_WARN(1,warnline)
          call CCTK_WARN(1, "s2 is too large")
          !return

        else 

          ! limit appears twice
          flimit = sqrt(recovery_flimit * maxs2/s2)
          if (flimit.ne.flimit .and. weight.eq.1.d0) then
             write(warnline, '(a62,6g15.7,a6,3g15.7,a1)') 'Scale Factor is NaN! '&
             &'(tau, rho, udens, w_lorentz, s2, maxs2)=(',&
                   utau, rhonew, udens, w_lorentz, s2, maxs2, ') at (',x,y,z,')'
             call CCTK_WARN(0,warnline)
          end if
          sx = flimit*sx
          sy = flimit*sy
          sz = flimit*sz

          s2_limited = .true.
          if (recalc_rho_after_limit.ne.0) then
            goto 10 ! sue me...
          end if

        end if

    end if
  end if

  ! Note the "goto 10" above, we only get here we pass the s2 limit

  if(velx*vlowx+vely*vlowy+velz*vlowz.gt.warn_if_speed_above .and. weight.eq.1.d0) then
    write(warnline,'(a60)') 'Caught speed above RecoverMarti::warn_if_speed_above'
    call CCTK_WARN(1,warnline)
    write(warnline,'(a14,i8)') 'in iteration: ',cctk_iteration
    call CCTK_WARN(1,warnline)
    write(warnline,'(a60,g15.6)') 'weight function (from CarpetReduce): ',weight
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,3g16.7)') 'xyz location: ',x,y,z
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,g16.7)') 'speed: ',velx*vlowx+vely*vlowy+velz*vlowz
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,g16.7)') 'w_lorentz: ',w_lorentz
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,g16.7)') 'rho: ',rho
    call CCTK_WARN(1,warnline)
  end if

  ! this should never happen (only happens if rho*w_lorentz != udens)
  if (velx*vlowx+vely*vlowy+velz*vlowz .gt. 1.0d0 .and. weight.eq.1.d0) then
    call CCTK_WARN(1,'Con2prim_ptPolytype: stopping the code.')
    call CCTK_WARN(1, '   velocity is larger than c! ')
    call CCTK_WARN(1,'Even with mesh refinement, this point is not restricted from a finer level, so this is really an error')
    write(warnline,'(a14,i8)') 'in iteration: ',cctk_iteration
    call CCTK_WARN(1,warnline)
    write(warnline,'(a28,i2)') 'on carpet reflevel: ',whisky_reflevel
    call CCTK_WARN(1,warnline)
    write(warnline,'(a60,g15.6)') 'weight function (from CarpetReduce): ',weight
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,3g16.7)') 'xyz location: ',x,y,z
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,g16.7)') 'radius: ',r
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,3g16.7)') 'velocities: ',velx,vely,velz
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,g16.7)') 'speed: ',velx*vlowx+vely*vlowy+velz*vlowz
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,g16.7)') 'rho: ',rho
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,g16.7)') 'epsilon: ',epsilon
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,g16.7)') 'pressure: ',press
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,g16.7)') 'w_lorentz: ',w_lorentz
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,g16.7)') 'udens: ',udens
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,g16.7)') 'df: ',df
    call CCTK_WARN(1,warnline)
    write(warnline,'(a20,g16.7)') 'det: ',det
    call CCTK_WARN(1,warnline)
!!$      write(warnline,'(a60,3g16.7)') 'uxx, uxy, uxz: ',uxx, uxy, uxz
!!$      call CCTK_WARN(1,warnline)
!!$      write(warnline,'(a60,3g16.7)') 'uyy, uyz, uzz: ',uyy, uyz, uzz
!!$      call CCTK_WARN(1,warnline)
    call CCTK_WARN(warnlevel, "speed larger than c")
    return
  endif

  if (epsilon .lt. 0.0d0) then
    if (whisky_reflevel.ge.warn_from_reflevel .and. weight.eq.1d0) then
      call CCTK_WARN(1,'Con2prim_ptPolytype: stopping the code.')
      call CCTK_WARN(1, '   specific internal energy just went below 0! ')
      call CCTK_WARN(1,'Even with mesh refinement, this point is not restricted from a finer level, so this is really an error')
      write(warnline,'(a14,i8)') 'in iteration: ',cctk_iteration
      call CCTK_WARN(1,warnline)
      write(warnline,'(a28,i2)') 'on carpet reflevel: ',whisky_reflevel
      call CCTK_WARN(1,warnline)
      write(warnline,'(a60,g15.6)') 'weight function (from CarpetReduce): ',weight
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,3g16.7)') 'xyz location: ',x,y,z
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,g16.7)') 'radius: ',r
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,3g16.7)') 'velocities: ',velx,vely,velz
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,g16.7)') 'rho: ',rho
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,g16.7)') 'epsilon: ',epsilon
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,g16.7)') 'pressure: ',press
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,g16.7)') 'w_lorentz: ',w_lorentz
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,g16.7)') 'df: ',df
      call CCTK_WARN(1,warnline)
      write(warnline,'(a20,g16.7)') 'det: ',det
      call CCTK_WARN(1,warnline)
  !!$      write(warnline,'(a60,3g16.7)') 'uxx, uxy, uxz: ',uxx, uxy, uxz
  !!$      call CCTK_WARN(1,warnline)
  !!$      write(warnline,'(a60,3g16.7)') 'uyy, uyz, uzz: ',uyy, uyz, uzz
  !!$      call CCTK_WARN(1,warnline)
      call CCTK_WARN(warnlevel, "Specific internal energy negative")
      return
    else
      write(warnline,'(a60,i2)') 'Con2Prim: eps negative, but I was told to ignore level: ',whisky_reflevel
      call CCTK_WARN(1,warnline)
      write(warnline,'(a14,i8)') 'in iteration: ',cctk_iteration
      call CCTK_WARN(1,warnline)
      return            
    endif
  endif

  return

end subroutine Con2Prim_ptPolytype

 /*@@
   @routine    Con2Prim_ptTracer
   @date       Mon Mar  8 14:26:29 2004
   @author     Ian Hawke
   @desc 
   
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Con2Prim_ptTracer(cons_tracer, tracer, dens) 
  
  implicit none
  
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL cons_tracer, tracer, dens

  tracer = cons_tracer / dens

end subroutine Con2Prim_ptTracer


 /*@@
   @routine    Con2Prim_ptTracerBounds
   @date       Mon Mar  8 14:26:29 2004
   @author     Ian Hawke
   @desc 
   
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Con2Prim_ptBoundsTracer(cons_tracer, tracer, rho, &
     velx, vely, velz, w_lorentz, det)
  
  implicit none
  
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL cons_tracer, tracer, rho, &
     velx, vely, velz, w_lorentz, det

  tracer = cons_tracer / sqrt(det) / rho / w_lorentz

end subroutine Con2Prim_ptBoundsTracer


