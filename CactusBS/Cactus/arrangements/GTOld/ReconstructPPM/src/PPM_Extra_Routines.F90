 /*@@
   @file      PPM.F90
   @date      Sun Feb 10 16:53:29 2002
   @author    Ian Hawke, Toni Font, Luca Baiotti, Frank Loeffler
   @desc 
   Routines to do PPM reconstruction.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

module PPM_Extra_Routines

 CONTAINS
  subroutine PPM_TVD(origm, orig, origp, bextm, bextp)
    CCTK_REAL :: origm, orig, origp, bextm, bextp
    CCTK_REAL :: dloc, dupw, delta

    dupw = orig - origm
    dloc = origp - orig
    if (dupw*dloc < 0.d0) then
      delta=0.d0
    else if (abs(dupw) < abs(dloc)) then
      delta=dupw
    else
      delta=dloc
    end if
    bextm = orig - 0.5d0 * delta
    bextp = orig + 0.5d0 * delta
  end subroutine PPM_TVD

  subroutine SimplePPM_tracer_1d(nx,dx,rho,velx,vely,velz, &
       tracer,tracerminus,tracerplus,press)
  
    USE PPM_Scalars
  
    implicit none
  
    DECLARE_CCTK_PARAMETERS
  
    CCTK_INT :: nx
    CCTK_REAL :: dx
    CCTK_REAL, dimension(nx) :: rho,velx,vely,velz
    CCTK_REAL, dimension(nx,number_of_tracers) :: tracer,tracerminus,tracerplus
    CCTK_REAL :: tracerflatomega
  
  
    CCTK_INT :: i,s,itracer
    CCTK_REAL, dimension(nx) :: press,dpress,tilde_flatten
    CCTK_REAL, dimension(nx,number_of_tracers) :: dmtracer, dtracer, tracerflat!, d2tracer
    CCTK_REAL :: dpress2,w,flatten,dvel
    CCTK_REAL :: eta, etatilde
  
  !!$  Average slopes delta_m(a). See (1.7) of Colella and Woodward, p.178
  !!$  This is the expression for an even grid.
  
  
    do i = 2, nx - 1
      dpress(i) = press(i+1) - press(i-1)
    end do
  
    do itracer=1,number_of_tracers
       do i = 2, nx - 1
          dtracer(i,itracer) = 0.5d0 * (tracer(i+1,itracer) - tracer(i-1,itracer))
  !        d2tracer(i,itracer) = (tracer(i+1) - 2.d0 * tracer(i) + tracer(i-1))! / 6.d0 / dx / dx
  !    ! since we use d2tracer only for the condition d2tracer(i+1)*d2tracer(i-1)<0 
  !    ! the denominator is not necessary
       enddo
    enddo
  
  !!$  Steepened slope. See (1.8) of Colella and Woodward, p.178
  
    do itracer=1,number_of_tracers
       do i = 2, nx - 1
          if( (tracer(i+1,itracer) - tracer(i,itracer)) * &
               (tracer(i,itracer) - tracer(i-1,itracer)) > 0.0d0 ) then
             dmtracer(i,itracer) = sign(1.0d0,dtracer(i,itracer)) * &
                  min(abs(dtracer(i,itracer)), 2.0d0 * &
                  abs(tracer(i,itracer) - tracer(i-1,itracer)), &
                  2.0d0 * abs(tracer(i+1,itracer) - tracer(i,itracer)))
          else
             dmtracer(i,itracer) = 0.0d0
          endif
       end do
    enddo
  
  !!$  Initial boundary states. See (1.9) of Colella and Woodward, p.178
  
      do itracer=1,number_of_tracers
         do i = 2, nx - 2
            tracerplus(i,itracer) = 0.5d0 * (tracer(i,itracer) + tracer(i+1,itracer)) + &
                 (dmtracer(i,itracer) - dmtracer(i+1,itracer)) / 6.d0
            tracerminus(i+1,itracer) = tracerplus(i,itracer)
         enddo
      enddo
      
  
  !!$Discontinuity steepening. See (1.14-17) of C&W.
  !!$This is the detect routine which mat be activated with the ppm_detect parameter
  !!$Note that this part really also depends on the grid being even. 
  !!$Note also that we do not have access to the gas constant gamma.
  !!$So this is just dropped from eq. (3.2) of C&W.
  !!$We can get around this by just rescaling the constant k0 (ppm_k0 here).
  
  !!! We might play around with this for the tracer. CURRENTLY TURNED OFF
  
#if 0
  if (ppm_detect .eq. 1000) then
     do itracer=1,number_of_tracers
        
        do i = 3, nx - 2
           if ( (dtracer(i+1,itracer)*dtracer(i-1,itracer) > 0.d0) & !make sure this is not an extremum
           .and.(abs(tracer(i+1,itracer)-tracer(i-1,itracer)) - & !this is to prevent steepening
                !of relatively small composition jumps
                ppm_epsilon_shock * min(tracer(i+1,itracer), tracer(i-1,itracer)) > 0.d0 )  & 
                .and. & ! the actual criterion from Plewa & Mueller
                 ((tracer(i+1,itracer)-tracer(i-1,itracer)) / &
                 (tracer(i+2,itracer)-tracer(i-2,itracer)) > ppm_omega1 ) ) then

           etatilde = (tracer(i-2,itracer) - tracer(i+2,itracer) + & 
                4.d0 * dtracer(i,itracer)) / (dtracer(i,itracer) * 12.d0)

           write(*,*) "Additional Steepening in Zone: ",i

        else
           etatilde = 0.d0
        end if
        eta = max(0.d0, min(1.d0, ppm_eta1 * (etatilde - ppm_eta2)))
        if (ppm_k0 * abs(dtracer(i,itracer)) * min(press(i-1),press(i+1)) < &
             abs(dpress(i)) * min(tracer(i-1,itracer), tracer(i+1,itracer))) then
           eta = 0.d0
        end if
        tracerminus(i,itracer) = tracerminus(i,itracer) * (1.d0 - eta) + &
             (tracer(i-1,itracer) + 0.5d0 * dmtracer(i-1,itracer)) * eta
        tracerplus(i,itracer) = tracerplus(i,itracer) * (1.d0 - eta) + &
             (tracer(i+1,itracer) - 0.5d0 * dmtracer(i+1,itracer)) * eta
     end do

  enddo

  end if
#endif

  !!$  Zone flattening. See appendix of C&W, p. 197-8.
  
    do i = 3, nx - 2
      dpress2 = press(i+2) - press(i-2)
      dvel = velx(i+1) - velx(i-1)
      if ( (abs(dpress(i)) >  ppm_epsilon * min(press(i-1),press(i+1))) .and. &
           (dvel < 0.d0) ) then
        w = 1.d0
      else
        w = 0.d0
      end if
      if (abs(dpress2) < ppm_small) then
        tilde_flatten(i) = 1.d0
      else
        tilde_flatten(i) = max(0.d0, 1.d0 - w * max(0.d0, ppm_omega2 * &
             (dpress(i) / dpress2 - ppm_omega1)))
      end if
    end do
  
    if (PPM3) then
       do itracer=1,number_of_tracers
          do i = 3, nx - 2
             flatten = tilde_flatten(i)
             tracerplus(i,itracer) = flatten * tracerplus(i,itracer) & 
                  + (1.d0 - flatten) * tracer(i,itracer)
             tracerminus(i,itracer) = flatten * tracerminus(i,itracer) & 
                  + (1.d0 - flatten) * tracer(i,itracer)
          end do
       enddo
    else  !!$ Really implement C&W, page 197; which requires stencil 4.
       do itracer=1,number_of_tracers
          do i = 4, nx - 3
             s=sign(1.d0, -dpress(i))
             flatten = max(tilde_flatten(i), tilde_flatten(i+s))  
             tracerplus(i,itracer) = flatten * tracerplus(i,itracer) + &
                  (1.d0 - flatten) * tracer(i,itracer)
             tracerminus(i,itracer) = flatten * tracerminus(i,itracer) & 
                  + (1.d0 - flatten) * tracer(i,itracer)
          end do
       enddo
    end if
  
  
  !! Additional flattening a la Plewa & Mueller                                                                 
  
#if 1
    do itracer=1,number_of_tracers
       do i = 2, nx - 1
          if ( ( tracer(i+1,itracer) - tracer(i,itracer) ) * &
             ( tracer(i,itracer) - tracer(i-1,itracer) ) < 0.0d0 ) then
             tracerflat(i,itracer) = 1.0d0
          else
             tracerflat(i,itracer) = 0.0d0
          endif
       enddo
    enddo
  
    do itracer=1,number_of_tracers
       do i = 3, nx -2
  
          tracerflatomega = 0.5d0 * max(tracerflat(i-1,itracer),2.0d0*tracerflat(i,itracer), &
               tracerflat(i+1,itracer)) * ppm_omega_tracer
  
          tracerplus(i,itracer) = tracerflatomega*tracer(i,itracer) + &
               (1.0d0 - tracerflatomega)*tracerplus(i,itracer)
  
          tracerminus(i,itracer) = tracerflatomega*tracer(i,itracer) + &
               (1.0d0 - tracerflatomega)*tracerminus(i,itracer)
  
       enddo
    enddo
#endif

  !!$ Monotonicity. See (1.10) of C&W.                                                                          
  
  
    do itracer=1,number_of_tracers
       do i = whisky_stencil, nx - whisky_stencil + 1
          if (((tracerplus(i,itracer)-tracer(i,itracer))*      &
             (tracer(i,itracer)-tracerminus(i,itracer)) .le. 0.d0)) then
             tracerminus(i,itracer) = tracer(i,itracer)
             tracerplus(i,itracer) = tracer(i,itracer)
          else if ((tracerplus(i,itracer) - tracerminus(i,itracer)) * (tracer(i,itracer) - 0.5d0 * &
               (tracerplus(i,itracer) + tracerminus(i,itracer))) > &
             (tracerplus(i,itracer) - tracerminus(i,itracer))**2 / 6.d0) then
             tracerminus(i,itracer) = 3.d0 * tracer(i,itracer) - 2.d0 * tracerplus(i,itracer)
          else if ((tracerplus(i,itracer) - tracerminus(i,itracer)) * (tracer(i,itracer) - 0.5d0 * &
               (tracerplus(i,itracer) + tracerminus(i,itracer))) <  &
             -(tracerplus(i,itracer) - tracerminus(i,itracer))**2 / 6.d0 ) then
             tracerplus(i,itracer) = 3.d0 * tracer(i,itracer) - 2.d0 * tracerminus(i,itracer)
          end if
       end do
    enddo
  
  
  
  end subroutine SimplePPM_tracer_1d

end module PPM_Extra_Routines

