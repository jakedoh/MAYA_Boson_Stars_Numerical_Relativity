! Maya
! analytic: Sets fields variables to analytic solution Teukolsky Wave
! $Header: /usr/center/raid1/hinder/IanCVS/arrangements/Ian/TeukolskyADMBase/src/teukolsky.F90,v 1.1 2006/01/31 02:00:25 ich Exp $

! Since the Teukolsky wave has 'so many degrees of freedom', more features
! can be added to this thorn. If you plan to do so watch out, however, for
! the irregular behaviour at the origin and on the z-axis. In particular
! the divisions by $r^5$ require very careful treatment and will normally
! exceed the limits of machine (!) precision.

#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

module TeukolskyWaveTools
  implicit none
  public   calcMetric
  public   calcMetricOrigin_Eppley
  public   calcShape_Eppley
  public   calcABC_Eppley

contains
!
!==============================================================================
!
subroutine calcMetric(gg, kk, coa, cob, coc, col, cok, coadot, cobdot, cocdot, &
                      coldot, cokdot, xx, mm, even)

  implicit none
  CCTK_REAL, dimension(3,3), intent(out) :: gg, kk
  CCTK_REAL, intent(in)                  :: coa, cob, coc, coadot, cobdot,   &
                                            cocdot
  CCTK_REAL, intent(in)                  :: col, cok, coldot, cokdot   
  CCTK_REAL, dimension(3), intent(in)    :: xx
  CCTK_INT,  intent(in)                  :: mm, even
  CCTK_REAL                                 x, y, z, x2, y2, z2
  CCTK_REAL                                 rad2, rad4, rho2, rad, rad3

  x    = xx(1)
  y    = xx(2)
  z    = xx(3)
  x2   = x * x
  y2   = y * y
  z2   = z * z
  rad2 = x2 + y2 + z2
  rad  = sqrt(rad2)
  rad3 = rad2 * rad
  rad4 = rad2 * rad2
  rho2 = x2 + y2

  select case(even)
    case(1)
      select case(mm)
        case(-2)
             gg(1,1) = 1 + 2 * y * x / rad4 * (x2 * coa + (y2 + z2 - x2)     &
                           * cob - (y2 + z2) * coc)
             gg(1,2) =     (-z2 / rad2 + 2 * x2 * y2 / rad4) * coa           &
                           + (1 - z2 / rad2 - 4 * x2 * y2 / rad4) * cob      &
                           + 2 * (z2 / rad2 + x2 * y2 / rad4) * coc
             gg(1,3) =     y * z / rad4 * ((rad2 + 2 * x2) * coa             &
                           + (rad2 - 4 * x2) * cob - 2 * (y2 + z2) * coc)
             gg(2,1) = gg(1,2)
             gg(2,2) = 1 + 2 * x * y / rad4 * (y2 * coa + (x2 + z2 - y2)     &
                           * cob - (x2 + z2) * coc)
             gg(2,3) =     x * z / rad4 * ((rad2 + 2 * y2) * coa + (rad2     &
                           - 4 * y2) * cob - 2 * (x2 + z2) * coc)
             gg(3,1) = gg(1,3)
             gg(3,2) = gg(2,3)
             gg(3,3) = 1 + 2 * x * y / rad4 * (-(x2 + y2) * coa - 2 * z2     &
                           * cob + (rad2 + z2) * coc)
                           
             kk(1,1) = 2 * y * x / rad4 * (x2 * coadot + (y2 + z2 - x2)      &
                       * cobdot - (y2 + z2) * cocdot)                        
             kk(1,2) = (-z2 / rad2 + 2 * x2 * y2 / rad4) * coadot            &
                       + (1 - z2 / rad2 - 4 * x2 * y2 / rad4) * cobdot       &
                       + 2 * (z2 / rad2 + x2 * y2 / rad4) * cocdot           
             kk(1,3) = y * z / rad4 * ((rad2 + 2 * x2) * coadot              &
                       + (rad2 - 4 * x2) * cobdot - 2 * (y2 + z2) * cocdot)
             kk(2,1) = kk(1,2)
             kk(2,2) = 2 * x * y / rad4 * (y2 * coadot + (x2 + z2 - y2)      &
                       * cobdot - (x2 + z2) * cocdot)
             kk(2,3) = x * z / rad4 * ((rad2 + 2 * y2) * coadot + (rad2      &
                       - 4 * y2) * cobdot - 2 * (x2 + z2) * cocdot)
             kk(3,1) = kk(1,3)
             kk(3,2) = kk(2,3)
             kk(3,3) = 2 * x * y / rad4 * (-(x2 + y2) * coadot - 2 * z2      &
                       * cobdot + (rad2 + z2) * cocdot)
             !----------------------------------------------------------------
         case(-1)
             gg(1,1) = 1 + 2 * y * z / rad4 * (-(y2 + z2) * coa - 2 * x2     &
                           * cob + (rad2 + x2) * coc)
             gg(1,2) =     z * x / rad4 * ((rad2 + 2 * y2) * coa + (rad2     &
                           - 4 * y2) * cob - 2 * (x2 + z2) * coc)
             gg(1,3) =     x * y / rad4 * ((rad2 + 2 * z2) * coa + (rad2     &
                           - 4 * z2) * cob - 2 * (x2 + y2) * coc)
             gg(2,1) = gg(1,2)
             gg(2,2) = 1 + 2 * y * z / rad4 * (y2 * coa + (rad2 - 2 * y2)    &
                           * cob - (x2 + z2) * coc)
             gg(2,3) =     (-x2 / rad2 + 2 * y2 * z2 / rad4) * coa           &
                           + (1 - x2 / rad2 - 4 * y2 * z2 / rad4) * cob      &
                           + 2 * (x2 / rad2 + y2 * z2 / rad4) * coc
             gg(3,1) = gg(1,3)
             gg(3,2) = gg(2,3)
             gg(3,3) = 1 + 2 * y * z / rad4 * (z2 * coa + (rad2 - 2 * z2)    &
                           * cob - (x2 + y2) * coc)
                           
             kk(1,1) = 2 * y * z / rad4 * (-(y2 + z2) * coadot - 2 * x2      &
                       * cobdot + (rad2 + x2) * cocdot)
             kk(1,2) = z * x / rad4 * ((rad2 + 2 * y2) * coadot + (rad2      &
                       - 4 * y2) * cobdot - 2 * (x2 + z2) * cocdot)
             kk(1,3) = x * y / rad4 * ((rad2 + 2 * z2) * coadot + (rad2      &
                       - 4 * z2) * cobdot - 2 * (x2 + y2) * cocdot)
             kk(2,1) = kk(1,2)
             kk(2,2) = 2 * y * z / rad4 * (y2 * coadot + (rad2 - 2 * y2)     &
                       * cobdot - (x2 + z2) * cocdot)
             kk(2,3) = (-x2 / rad2 + 2 * y2 * z2 / rad4) * coadot            &
                       + (1 - x2 / rad2 - 4 * y2 * z2 / rad4) * cobdot       &
                       + 2 * (x2 / rad2 + y2 * z2 / rad4) * cocdot
             kk(3,1) = kk(1,3)
             kk(3,2) = kk(2,3)
             kk(3,3) = 2 * y * z / rad4 * (z2 * coadot + (rad2 - 2 * z2)     &
                       * cobdot - (x2 + y2) * cocdot)
             !----------------------------------------------------------------
         case(0) 
             gg(1,1) = 1 + (-1 + 3 * y2 / rad2 + 3 * x2 * z2 / rad4) * coa   &
                           - 6 * z2 * x2 / rad4 * cob                        &
                           + 3 * (-y2 / rad2 + x2 * z2 / rad4) * coc
             gg(1,2) =     3 * x * y / rad4 * (-rho2 * coa - 2 * z2 * cob    &
                           + (rad2 + z2) * coc)
             gg(1,3) =     3 * x * z / rad4 * (z2 * coa + (rho2 - z2) * cob  &
                           - rho2 * coc)
             gg(2,1) = gg(1,2)
             gg(2,2) = 1 + (-1 + 3 * x2 / rad2 + 3 * y2 * z2 / rad4) * coa   &
                           - 6 * y2 * z2 / rad4 * cob                        &
                           + 3 * (-x2 / rad2 + y2 * z2 / rad4) * coc
             gg(2,3) =     3 * y * z / rad4 * (z2 * coa + (rho2 - z2) * cob  &
                           - rho2 * coc)
             gg(3,1) = gg(1,3)
             gg(3,2) = gg(2,3)
             gg(3,3) = 1 + (-1 + 3 * z2*z2 / rad4) * coa + 6 * z2 * rho2     &
                           / rad4 * cob + 3 * rho2*rho2 / rad4 * coc

             kk(1,1) = (-1 + 3 * y2 / rad2 + 3 * x2 * z2 / rad4) * coadot    &
                       - 6 * z2 * x2 / rad4 * cobdot                         &
                       + 3 * (-y2 / rad2 + x2 * z2 / rad4) * cocdot
             kk(1,2) = 3 * x * y / rad4 * (-rho2 * coadot - 2 * z2 * cobdot  &
                       + (rad2 + z2) * cocdot)
             kk(1,3) = 3 * x * z / rad4 * (z2 * coadot                       &
                       + (rho2 - z2) * cobdot - rho2 * cocdot)
             kk(2,1) = kk(1,2)
             kk(2,2) = (-1 + 3 * x2 / rad2 + 3 * y2 * z2 / rad4) * coadot    &
                       - 6 * y2 * z2 / rad4 * cobdot                         &
                       + 3 * (-x2 / rad2 + y2 * z2 / rad4) * cocdot
             kk(2,3) = 3 * y * z / rad4 * (z2 * coadot                       &
                       + (rho2 - z2) * cobdot - rho2 * cocdot)
             kk(3,1) = kk(1,3)
             kk(3,2) = kk(2,3)
             kk(3,3) = (-1 + 3 * z2*z2 / rad4) * coadot + 6 * z2 * rho2      &
                       / rad4 * cobdot + 3 * rho2*rho2 / rad4 * cocdot
             !----------------------------------------------------------------
         case(1) 
             gg(1,1) = 1 + 2 * x * z / rad4 * (x2 * coa + (rad2 - 2 * x2)    &
                           * cob - (y2 + z2) * coc)
             gg(1,2) = y * z / rad4 * ((rad2 + 2 * x2) * coa + (rad2         &
                           - 4 * x2) * cob - 2 * (y2 + z2) * coc)
             gg(1,3) =     (-y2 / rad2 + 2 * x2 * z2 / rad4) * coa           &
                           + (1 - y2 / rad2 - 4 * x2 * z2 / rad4) * cob      &
                           + 2 * (y2 / rad2 + x2 * z2 / rad4) * coc
             gg(2,1) = gg(1,2)
             gg(2,2) = 1 + 2 * x * z / rad4 * (-(x2 + z2) * coa - 2 * y2     &
                           * cob + (rad2 + y2) * coc)
             gg(2,3) =     x * y / rad4 * ((rad2 + 2 * z2) * coa + (rad2     &
                           - 4 * z2) * cob - 2 * (x2 + y2) * coc)
             gg(3,1) = gg(1,3)
             gg(3,2) = gg(2,3)
             gg(3,3) = 1 + 2 * x * z / rad4 * (z2 * coa + (rad2 - 2 * z2)    &
                           * cob - (x2 + y2) * coc)

             kk(1,1) = 2 * x * z / rad4 * (x2 * coadot + (rad2 - 2 * x2)     &
                       * cobdot - (y2 + z2) * cocdot)
             kk(1,2) = y * z / rad4 * ((rad2 + 2 * x2) * coadot + (rad2      &
                       - 4 * x2) * cobdot - 2 * (y2 + z2) * cocdot)
             kk(1,3) = (-y2 / rad2 + 2 * x2 * z2 / rad4) * coadot            &
                       + (1 - y2 / rad2 - 4 * x2 * z2 / rad4) * cobdot       &
                       + 2 * (y2 / rad2 + x2 * z2 / rad4) * cocdot
             kk(2,1) = kk(1,2)
             kk(2,2) = 2 * x * z / rad4 * (-(x2 + z2) * coadot - 2 * y2      &
                       * cobdot + (rad2 + y2) * cocdot)
             kk(2,3) = x * y / rad4 * ((rad2 + 2 * z2) * coadot + (rad2      &
                       - 4 * z2) * cobdot - 2 * (x2 + y2) * cocdot)
             kk(3,1) = kk(1,3)
             kk(3,2) = kk(2,3)
             kk(3,3) = 2 * x * z / rad4 * (z2 * coadot + (rad2 - 2 * z2)     &
                       * cobdot - (x2 + y2) * cocdot)
             !----------------------------------------------------------------
         case(2) 
             gg(1,1) = 1 + ((x2 - z2) / rad2 - x2 * (z2 + 2 * y2) / rad4)    &
                           * coa + 2 * x2 * (z2 + 2 * y2) / rad4 * cob       &
                           + ((y2 + 2 * z2) / rad2 - x2 * (z2 + 2 * y2)      &
                           / rad4) * coc
             gg(1,2) =     y * x * (x2 - y2) / rad4 * (coa - 2 * cob + coc)
             gg(1,3) =     x * z / rad4 * ((2 * x2 + z2) * coa + (-x2        &
                           + 3 * y2 + z2) * cob - (x2 + 2 * z2 + 3 * y2)     &
                           * coc)
             gg(2,1) = gg(1,2)
             gg(2,2) = 1 + ((z2 - y2) / rad2 + y2 * (2 * x2 + z2) / rad4)    &
                           * coa - 2 * y2 * (2 * x2 + z2) / rad4 * cob       &
                           + (-(x2 + 2 * z2) / rad2 + y2 * (2 * x2 + z2)     &
                           / rad4) * coc
             gg(2,3) =     y * z / rad4 * (-(z2 + 2 * y2) * coa - (3 * x2    &
                           - y2 + z2) * cob + (3 * x2 + 2 * z2 + y2) * coc)
             gg(3,1) = gg(1,3)
             gg(3,2) = gg(2,3)
             gg(3,3) = 1 + (y2*y2 - x2*x2) / rad4 *coa - 2 * z2 * (x2 - y2)  &
                           / rad4 * cob + (x2 - y2) * (rad2 + z2) / rad4     &
                           * coc
                           
             kk(1,1) = ((x2 - z2) / rad2 - x2 * (z2 + 2 * y2) / rad4)        &
                       * coadot + 2 * x2 * (z2 + 2 * y2) / rad4 * cobdot     &
                       + ((y2 + 2 * z2) / rad2 - x2 * (z2 + 2 * y2)          &
                       / rad4) * cocdot
             kk(1,2) = y * x * (x2 - y2) / rad4 * (coadot - 2 * cobdot       &
                       + cocdot)
             kk(1,3) = x * z / rad4 * ((2 * x2 + z2) * coadot + (-x2         &
                       + 3 * y2 + z2) * cobdot - (x2 + 2 * z2 + 3 * y2)      &
                       * cocdot)
             kk(2,1) = kk(1,2)
             kk(2,2) = ((z2 - y2) / rad2 + y2 * (2 * x2 + z2) / rad4)        &
                       * coadot - 2 * y2 * (2 * x2 + z2) / rad4 * cobdot     &
                       + (-(x2 + 2 * z2) / rad2 + y2 * (2 * x2 + z2)         &
                       / rad4) * cocdot
             kk(2,3) = y * z / rad4 * (-(z2 + 2 * y2) * coadot - (3 * x2     &
                       - y2 + z2) * cobdot + (3 * x2 + 2 * z2 + y2) * cocdot)
             kk(3,1) = kk(1,3)
             kk(3,2) = kk(2,3)
             kk(3,3) = (y2*y2 - x2*x2) / rad4 *coadot - 2 * z2 * (x2 - y2)   &
                       / rad4 * cobdot + (x2 - y2) * (rad2 + z2) / rad4      &
                       * cocdot
             !----------------------------------------------------------------
      end select
    case(0)
      select case(mm)
         case(0)
             gg(1,1) = 1 + x * y * z * (8 * cok + 2 * col) / rad3 
             gg(1,2) = - z * (x2 - y2) * (4 * cok + col) / rad3
             gg(1,3) = - y * (col * rho2 - 4 * cok * z2) / rad3
             gg(2,1) = gg(1,2)
             gg(2,2) = 1 - x * y * z * (8 * cok + 2 * col) / rad3
             gg(2,3) = x * (col * rho2 - 4 * cok * z2) / rad3
             gg(3,1) = gg(1,3)
             gg(3,2) = gg(2,3)
             gg(3,3) = 1

             kk(1,1) = x * y * z * (8 * cokdot + 2 * coldot) / rad3
             kk(1,2) = - z * (x2 - y2) * (4 * cokdot + coldot) / rad3
             kk(1,3) = - y * (coldot * rho2 - 4 * cokdot * z2) / rad3
             kk(2,1) = kk(1,2)
             kk(2,2) = - x * y * z * (8 * cokdot + 2 * coldot) / rad3
             kk(2,3) = x * (coldot * rho2 - 4 * cokdot * z2) / rad3
             kk(3,1) = kk(1,3)
             kk(3,2) = kk(2,3)
             kk(3,3) = 0

      end select
  end select
end subroutine calcMetric
  
!
!==============================================================================
!
subroutine calcMetricOrigin_Eppley(gg, kk, time, amp, sigma, mm, even, zbef, npow, norm_packet)

  implicit none
  CCTK_REAL, dimension(3,3), intent(out) :: gg, kk
  CCTK_REAL, intent(in)                  :: time
  CCTK_REAL, intent(in)                  :: amp, sigma
  CCTK_INT,  intent(in)                  :: mm, even, npow, zbef, norm_packet
  CCTK_REAL                                 tos, coa, coadot, sigma2, fac
  CCTK_REAL                                 cob, coc, cobdot, cocdot, ex2s2

  gg     = 0
  kk     = 0
  tos    = time / sigma
  sigma2 = sigma * sigma
  ex2s2  = amp * exp(-tos*tos) 
  fac = ex2s2 * sigma**(npow-5)
  if (norm_packet .eq. 1) then
     ex2s2 = ex2s2 * (2./npow)**(npow/2.) * exp(npow/2.) / (sigma**(npow-1))
     fac = fac * (2./npow)**(npow/2.) * exp(npow/2.) / (sigma**(npow-1))
  end if

  if (zbef .eq. 1) then
     gg=0
     kk=0
  else

  ! A, B, C are all equal; we only use A, Adot
    coa = -( 20*npow*(1+2*npow**2)*tos**(npow-1)                               &
          - 40*(3+2*npow*(npow+2))*tos**(npow+1)                              &
          + 80*(npow+2)*tos**(npow+3) - 32*tos**(npow+5) ) * 2 * fac/5.          
    if (npow >= 3) then
        coa = coa + 10*npow*(npow-2)*(npow-1)**2*tos**(npow-3) * 2 * fac/5.
        if (npow >= 5) then
           coa = coa - npow*(npow-1)*(npow-2)*(npow-3)*(npow-4)*tos**(npow-5) * 2 * fac/5.
        end if
    end if
    cob = coa
    coc = coa

  ! Absorb this annoying factor of -2 in the definition of the EK here
    coadot = - ( 40*(2*npow+1)*(3+2*npow*(npow+1))*tos**npow                       &
             - 240*(3+npow*(npow+3))*tos**(npow+2) + 96*(5+2*npow)*tos**(npow+4) &
             - 64*tos**(npow+6) ) * fac/5. / sigma
    if (npow >= 2) then
       coadot = coadot + 60*npow*(npow-1)*(1+npow*(npow-1))*tos**(npow-2) * fac/(5.*sigma)
       if (npow >= 4) then
          coadot = coadot - 6*npow*(npow-1)*(npow-2)*(npow-3)*(2*npow-3)*tos**(npow-4) * fac/(5*sigma)
          if (npow >= 6) then
             coadot = coadot + npow*(npow-1)*(npow-2)*(npow-3)*(npow-4)*(npow-5)*tos**(npow-6)*fac/(5*sigma)
          end if
       end if
    end if

    cobdot = coadot
    cocdot = coadot

  select case(even)
    case(1)         
      select case(mm)
         case(-2)
             gg(1,1) = 1
             gg(1,2) = cob
             gg(2,1) = cob
             gg(2,2) = 1
             gg(3,3) = 1
             kk(1,2) = cobdot
             kk(2,1) = cobdot
         case(-1)
             gg(1,1) = 1
             gg(2,2) = 1
             gg(2,3) = 2*coc-coa
             gg(3,2) = 2*coc-coa 
             gg(3,3) = 1
             kk(2,3) = 2*cocdot-coadot
             kk(3,2) = 2*cocdot-coadot
         case(0)
             gg(1,1) = 1 - coa
             gg(2,2) = 1 + 2*coa - 3*coc
             gg(3,3) = 1 - coa + 3*coc
             kk(1,1) = -coadot
             kk(2,2) = 2*coadot - 3*cocdot
             kk(3,3) = -coadot - 3*cocdot
         case(1)
             gg(1,1) = 1
             gg(1,3) = cob
             gg(2,2) = 1
             gg(3,1) = cob
             gg(3,3) = 1
             kk(1,3) = cobdot
             kk(3,1) = cobdot
         case(2)
             gg(1,1) = 1 + coa
             gg(2,2) = 1 - coc
             gg(3,3) = 1 + coc - coa
             kk(1,1) = coadot
             kk(2,2) = -cocdot
             kk(3,3) = cocdot - coadot
      end select
    case(0)
      select case(mm)
         case(0)
             gg(1,1) = 1
             gg(2,2) = 1
             gg(3,3) = 1 
      end select
  end select

  end if
end subroutine
!
!==============================================================================
!
subroutine calcShape_Eppley( dff, u, amp, sigma, zbef, npow, norm_packet)

  implicit none
  CCTK_REAL, dimension(0:5), intent(out) :: dff
  CCTK_REAL, intent(in)                  :: u
  CCTK_REAL, intent(in)                  :: amp, sigma
  CCTK_REAL                                 uos, uos2, ex2s2, uos4
  CCTK_INT, intent(in)                   :: zbef, npow, norm_packet

  uos   = u / sigma
  uos2  = uos * uos
  uos4  = uos2 * uos2
  ex2s2 = amp * exp(-uos2)
  if (norm_packet .eq. 1) then
     ex2s2 = ex2s2 * (2./npow)**(npow/2.) * exp(npow/2.) / sigma**(npow-1)
  end if

  if( (zbef .eq. 1) .AND. (u .le. 0.0) ) then
     dff(0) = 0
     dff(1) = 0
     dff(2) = 0
     dff(3) = 0
     dff(4) = 0
     dff(5) = 0
  else
     if (u .eq. 0.0) then
        dff(0) = 0
        dff(1) = 0
        dff(2) = 0
        dff(3) = 0
        dff(4) = 0
        dff(5) = 0
        if (npow .eq. 1) dff(1)=amp
        if (npow .eq. 2) dff(2)=2*amp
        if (npow .eq. 3) dff(3)=6*amp
        if (npow .eq. 4) dff(4)=24*amp
        if (npow .eq. 5) dff(4)=120*amp
     else
        dff(0) = u**npow * ex2s2
        dff(1) = (npow-2.*uos2) * u**(npow-1.) * ex2s2
        dff(2) = ( (npow-1)*npow - 2*(1+2*npow)*uos2 + 4*uos4) *  u**(npow-2) * ex2s2
        dff(3) = ( npow*(npow-1)*(npow-2) - 6*npow**2*uos2 + 12*(npow+1)*uos4      &
                 - 8*uos4*uos2  ) * u**(npow-3) * ex2s2
        dff(4) = ( npow*(npow-1)*(npow-2)*(npow-3) - 4*npow*(npow-1)*(2*npow-1)*uos2 &
                 + 12*(1+2*npow*(npow+1))*uos4 - 16*(2*npow+3)*uos4*uos2           &
                 + 16*uos4*uos4 ) * u**(npow-4) * ex2s2
        dff(5) = ( npow*(npow-1)*(npow-2)*(npow-3)*(npow-4) - 10*npow*(npow-2)*(npow-1)**2*uos2 &
                 + 20*npow*(2*npow**2+1)*uos4 - 40*(3+2*npow*(npow+2))*uos4*uos2   &
                 + 80*(npow+2)*uos4*uos4 - 32*uos4*uos4*uos2) * u**(npow-5) * ex2s2
      end if
! Backup: n=5 
!     dff(0) = u * uos4 * ex2s2
!     dff(1) = (5 - 2 * uos2) * uos4 * ex2s2
!     dff(2) = uos2 / sigma * uos * (20 - 22*uos2 + 4*uos4) * ex2s2
!     dff(3) = uos2 / (sigma*sigma) * (60 - 150*uos2 + 72*uos4 - 8*uos4*uos2) * ex2s2
!     dff(4) = uos / (sigma*sigma*sigma) * 4 * (30 - 180*uos2 + 183*uos4 & 
!              - 52*uos4*uos2+4*uos4*uos4) * ex2s2
!     dff(5) = 4 / (sigma*sigma*sigma*sigma) * (30 - 600 * uos2                  &
!              + 1275 * uos4 - 730 * uos4*uos2 + 140*uos4*uos4 - 8*uos4*uos4*uos2) &
!              * ex2s2

! Backup: n=1
!  dff(0) = u * ex2s2
!  dff(1) = (1 - 2 * uos2) * ex2s2
!  dff(2) = uos / sigma * (-6 + 4 * uos2) * ex2s2
!  dff(3) = 1 / (sigma*sigma) * (-6 + 24 * uos2 - 8 * uos2*uos2) * ex2s2
!  dff(4) = uos / (sigma*sigma*sigma) * (60 - 80 * uos2 + 16 * uos2*uos2)     &
!           * ex2s2
!  dff(5) = 1 / (sigma*sigma*sigma*sigma) * (60 - 360 * uos2                  &
!           + 240 * uos2*uos2 - 32 * uos2*uos2*uos2) * ex2s2
  end if

end subroutine calcShape_Eppley
!
!==============================================================================
!
subroutine calcABC_Eppley( coa, cob, coc, col, cok, coadot, cobdot, cocdot, coldot, &
                           cokdot, time, xx)

  implicit none
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  CCTK_REAL, intent(out)              :: coa, cob, coc, coadot, cobdot, cocdot
  CCTK_REAL, intent(out)              :: cok, col, cokdot, coldot
  CCTK_REAL, intent(in)               :: time
  CCTK_REAL, dimension(3), intent(in) :: xx
  CCTK_REAL, dimension(0:5)           :: dff
  CCTK_REAL                              rad, rad2, rad3, rad4, rad5
  CCTK_REAL                              u, tos, tos2
  CCTK_REAL                              r0p2, r0p4, r0p6, r0p8, sigp2, sigp4, sigp6, sigp8 

  rad  = sqrt( sum( xx*xx ) )
  rad2 = rad * rad
  rad3 = rad2 * rad
  rad4 = rad3 * rad
  rad5 = rad4 * rad

  if( rad <= texprad .and. CCTK_EQUALS( wave_dir, "timesym" ) ) then

       r0p2=rad0*rad0
       r0p4=r0p2*r0p2
       r0p6=r0p4*r0p2
       r0p8=r0p6*r0p2

       sigp2=sigma*sigma
       sigp4=sigp2*sigp2
       sigp6=sigp4*sigp2
       sigp8=sigp6*sigp2

       ! Evaluate the Taylor Expansions of the coefficients to 8th order
       coa = -((8*amp*(-8*r0p6 + 60*r0p4*sigp2 - 90*r0p2*sigp4 + 15*sigp6))/(5.*exp(r0p2/sigp2)*sigp8*sigp2) - &
             (8*amp*rad2*(16*r0p8 - 224*r0p6*sigp2 + 840*r0p4*sigp4 - 840*r0p2*sigp6 + 105*sigp8))/ &
             (35.*exp(r0p2/sigp2)*sigp8*sigp6) +                                                              &
             (4*amp*rad4*(-32*r0p6*r0p4 + 720*r0p8*sigp2 - 5040*r0p6*sigp4 + 12600*r0p4*sigp6 -               &
             9450*r0p2*sigp8 + 945*sigp6*sigp4))/(315.*exp(r0p2/sigp2)*sigp8*sigp8*sigp2) -                   &
             (4*amp*rad3*rad3*(64*r0p6*r0p6 - 2112*r0p8*r0p2*sigp2 + 23760*r0p8*sigp4 - 110880*r0p6*sigp6 +   &
             207900*r0p4*sigp8 - 124740*r0p2*sigp6*sigp4 + 10395*sigp6*sigp6))/(10395.*exp(r0p2/sigp2)*sigp8*sigp8*sigp6) + &
             (amp*rad4*rad4*(-128*r0p8*r0p6 + 5824*r0p6*r0p6*sigp2 - 96096*r0p8*r0p2*sigp4 + 720720*r0p8*sigp6 - &
             2522520*r0p6*sigp8 + 3783780*r0p4*sigp6*sigp4 - 1891890*r0p2*sigp6*sigp6 + 135135*sigp8*sigp6))/ &
             (135135.*exp(r0p2/sigp2)*sigp8*sigp8*sigp8*sigp2))
       cob = -((8*amp*(-8*r0p6 + 60*r0p4*sigp2 - 90*r0p2*sigp4 + 15*sigp6))/(5.*exp(r0p2/sigp2)*sigp8*sigp2) - &
             (8*amp*rad2*(16*r0p8 - 224*r0p6*sigp2 + 840*r0p4*sigp4 - 840*r0p2*sigp6 + 105*sigp8))/ &
             (21.*exp(r0p2/sigp2)*sigp8*sigp6) +                                                              &
             (4*amp*rad4*(-32*r0p6*r0p4 + 720*r0p8*sigp2 - 5040*r0p6*sigp4 + 12600*r0p4*sigp6 -               &
             9450*r0p2*sigp8 + 945*sigp6*sigp4))/(135.*exp(r0p2/sigp2)*sigp8*sigp8*sigp2) -                   &
             (4*amp*rad3*rad3*(64*r0p6*r0p6 - 2112*r0p8*r0p2*sigp2 + 23760*r0p8*sigp4 - 110880*r0p6*sigp6 +   &
             207900*r0p4*sigp8 - 124740*r0p2*sigp6*sigp4 + 10395*sigp6*sigp6))/(3465.*exp(r0p2/sigp2)*sigp8*sigp8*sigp6) + &
             (amp*rad4*rad4*(-128*r0p8*r0p6 + 5824*r0p6*r0p6*sigp2 - 96096*r0p8*r0p2*sigp4 + 720720*r0p8*sigp6 - &
             2522520*r0p6*sigp8 + 3783780*r0p4*sigp6*sigp4 - 1891890*r0p2*sigp6*sigp6 + 135135*sigp8*sigp6))/ &
             (36855.*exp(r0p2/sigp2)*sigp8*sigp8*sigp8*sigp2))
       coc = -((8*amp*(-8*r0p6 + 60*r0p4*sigp2 - 90*r0p2*sigp4 + 15*sigp6))/(5.*exp(r0p2/sigp2)*sigp8*sigp2) - &
             (8*amp*rad2*(16*r0p8 - 224*r0p6*sigp2 + 840*r0p4*sigp4 - 840*r0p2*sigp6 + 105*sigp8))/ &
             (15.*exp(r0p2/sigp2)*sigp8*sigp6) + &
             (52*amp*rad4*(-32*r0p6*r0p4 + 720*r0p8*sigp2 - 5040*r0p6*sigp4 + 12600*r0p4*sigp6 - &
             9450*r0p2*sigp8 + 945*sigp6*sigp4))/(945.*exp(r0p2/sigp2)*sigp8*sigp8*sigp2) - &
             (4*amp*rad3*rad3*(64*r0p6*r0p6 - 2112*r0p8*r0p2*sigp2 + 23760*r0p8*sigp4 - 110880*r0p6*sigp6 + &
             207900*r0p4*sigp8 - 124740*r0p2*sigp6*sigp4 + 10395*sigp6*sigp6))/(1485.*exp(r0p2/sigp2)*sigp8*sigp8*sigp6) + &
             (31.*amp*rad4*rad4*(-128*r0p8*r0p6 + 5824*r0p6*r0p6*sigp2 - 96096*r0p8*r0p2*sigp4 + 720720*r0p8*sigp6 - &
             2522520*r0p6*sigp8 + 3783780*r0p4*sigp6*sigp4 - 1891890*r0p2*sigp6*sigp6 + 135135*sigp8*sigp6))/ &
             (405405.*exp(r0p2/sigp2)*sigp8*sigp8*sigp8*sigp2))
       col = - (32*amp*rad*(8*r0p6 - 60*r0p4*sigp2 + 90*r0p2*sigp4 - 15*sigp6))/(15.*exp(r0p2/sigp2)*sigp6*sigp4) + &
             (16*amp*rad3*(16*r0p8 - 224*r0p6*sigp2 + 840*r0p4*sigp4 - 840*r0p2*sigp6 + 105*sigp8))/ &
             (35.*exp(r0p2/sigp2)*sigp8*sigp6) + &
             (32*amp*rad4*rad*(32*r0p8*r0p2 - 720*r0p8*sigp2 + 5040*r0p6*sigp4 - 12600*r0p4*sigp6 + 9450*r0p2*sigp8 - &
             945*sigp4*sigp6))/(945.*exp(r0p2/sigp2)*sigp8*sigp8*sigp2) + &
             (8*amp*rad3*rad4*(64*r0p8*r0p4 - 2112*r0p6*r0p4*sigp2 + 23760*r0p8*sigp4 - 110880*r0p6*sigp6 + &
             207900*r0p4*sigp8 - 124740*r0p2*sigp6*sigp4 + 10395*sigp8*sigp4))/(6237.*exp(r0p2/sigp2)*sigp8*sigp8*sigp6)
       cok = rad*coa/3.

       ! Evaluating the time derivatives
       coadot = (8*amp*rad0*(-8*r0p6 + 84*r0p4*sigp2 - 210*r0p2*sigp4 + 105*sigp6))/(5.*exp(r0p2/sigp2)*sigp8*sigp2) - &
             (8*amp*rad0*rad2*(16*r0p8 - 288*r0p6*sigp2 + 1512*r0p4*sigp4 - 2520*r0p2*sigp6 + 945*sigp8))/ &
             (35.*exp(r0p2/sigp2)*sigp8*sigp6) +                                                              &
             (4*amp*rad0*rad4*(-32*r0p6*r0p4 + 880*r0p8*sigp2 - 7920*r0p6*sigp4 + 27720*r0p4*sigp6 -               &
             34650*r0p2*sigp8 + 10395*sigp6*sigp4))/(315.*exp(r0p2/sigp2)*sigp8*sigp8*sigp2) -                   &
             (4*amp*rad3*rad3*(64*r0p6*r0p6 - 2496*r0p8*r0p2*sigp2 + 34320*r0p8*sigp4 - 205920*r0p6*sigp6 +   &
             540540*r0p4*sigp8 - 540540*r0p2*sigp6*sigp4 + 135135*sigp6*sigp6))/(10395.*exp(r0p2/sigp2)*sigp8*sigp8*sigp6) + &
             (amp*rad0*rad4*rad4*(-128*r0p8*r0p6 + 6720*r0p6*r0p6*sigp2 - 131040*r0p8*r0p2*sigp4 + 1201200*r0p8*sigp6 - &
             5405400*r0p6*sigp8 + 11351340*r0p4*sigp6*sigp4 - 9459450*r0p2*sigp6*sigp6 + 2027025*sigp8*sigp6))/ &
             (135135.*exp(r0p2/sigp2)*sigp8*sigp8*sigp8*sigp2)
       cobdot = (8*amp*rad0*(-8*r0p6 + 84*r0p4*sigp2 - 210*r0p2*sigp4 + 105*sigp6))/(5.*exp(r0p2/sigp2)*sigp8*sigp2) - &
             (8*amp*rad0*rad2*(16*r0p8 - 288*r0p6*sigp2 + 1512*r0p4*sigp4 - 2520*r0p2*sigp6 + 945*sigp8))/ &
             (21.*exp(r0p2/sigp2)*sigp8*sigp6) +                                                              &
             (4*amp*rad0*rad4*(-32*r0p6*r0p4 + 880*r0p8*sigp2 - 7920*r0p6*sigp4 + 27720*r0p4*sigp6 -               &
             34650*r0p2*sigp8 + 10395*sigp6*sigp4))/(135.*exp(r0p2/sigp2)*sigp8*sigp8*sigp2) -                   &
             (4*amp*rad3*rad3*(64*r0p6*r0p6 - 2496*r0p8*r0p2*sigp2 + 34320*r0p8*sigp4 - 205920*r0p6*sigp6 +   &
             540540*r0p4*sigp8 - 540540*r0p2*sigp6*sigp4 + 135135*sigp6*sigp6))/(3465.*exp(r0p2/sigp2)*sigp8*sigp8*sigp6) + &
             (amp*rad0*rad4*rad4*(-128*r0p8*r0p6 + 6720*r0p6*r0p6*sigp2 - 131040*r0p8*r0p2*sigp4 + 1201200*r0p8*sigp6 - &
             5405400*r0p6*sigp8 + 11351340*r0p4*sigp6*sigp4 - 9459450*r0p2*sigp6*sigp6 + 2027025*sigp8*sigp6))/ &
             (36855.*exp(r0p2/sigp2)*sigp8*sigp8*sigp8*sigp2)
       cocdot = (8*amp*rad0*(-8*r0p6 + 84*r0p4*sigp2 - 210*r0p2*sigp4 + 105*sigp6))/(5.*exp(r0p2/sigp2)*sigp8*sigp2) - &
             (8*amp*rad0*rad2*(16*r0p8 - 288*r0p6*sigp2 + 1512*r0p4*sigp4 - 2520*r0p2*sigp6 + 945*sigp8))/ &
             (15.*exp(r0p2/sigp2)*sigp8*sigp6) +                                                              &
             (4*amp*rad0*rad4*(-32*r0p6*r0p4 + 880*r0p8*sigp2 - 7920*r0p6*sigp4 + 27720*r0p4*sigp6 -               &
             34650*r0p2*sigp8 + 10395*sigp6*sigp4))/(945.*exp(r0p2/sigp2)*sigp8*sigp8*sigp2) -                   &
             (4*amp*rad3*rad3*(64*r0p6*r0p6 - 2496*r0p8*r0p2*sigp2 + 34320*r0p8*sigp4 - 205920*r0p6*sigp6 +   &
             540540*r0p4*sigp8 - 540540*r0p2*sigp6*sigp4 + 135135*sigp6*sigp6))/(1485.*exp(r0p2/sigp2)*sigp8*sigp8*sigp6) + &
             (amp*rad0*rad4*rad4*(-128*r0p8*r0p6 + 6720*r0p6*r0p6*sigp2 - 131040*r0p8*r0p2*sigp4 + 1201200*r0p8*sigp6 - &
             5405400*r0p6*sigp8 + 11351340*r0p4*sigp6*sigp4 - 9459450*r0p2*sigp6*sigp6 + 2027025*sigp8*sigp6))/ &
             (405405.*exp(r0p2/sigp2)*sigp8*sigp8*sigp8*sigp2)
       coldot = rad0 * (32*amp*rad*(8*r0p6 - 84*r0p4*sigp2 + 210*r0p2*sigp4 - 105*sigp6))/(15.*exp(r0p2/sigp2)*sigp6*sigp4) + &
             (16*amp*rad3*(16*r0p8 - 288*r0p6*sigp2 + 1512*r0p4*sigp4 - 2520*r0p2*sigp6 + 945*sigp8))/ &
             (35.*exp(r0p2/sigp2)*sigp8*sigp6) + &
             (32*amp*rad4*rad*(32*r0p8*r0p2 - 880*r0p8*sigp2 + 7920*r0p6*sigp4 - 27720*r0p4*sigp6 + 34650*r0p2*sigp8 - &
             10395*sigp4*sigp6))/(945.*exp(r0p2/sigp2)*sigp8*sigp8*sigp2) + &
             (8*amp*rad3*rad4*(64*r0p8*r0p4 - 2496*r0p6*r0p4*sigp2 + 34320*r0p8*sigp4 - 205920*r0p6*sigp6 + &
             540540*r0p4*sigp8 - 540540*r0p2*sigp6*sigp4 + 135135*sigp8*sigp4))/(6237.*exp(r0p2/sigp2)*sigp8*sigp8*sigp6)
       cokdot = rad*coadot/3.


  else if( CCTK_EQUALS( wave_dir, "outgoing" ) ) then
       ! Outgoing part
       !---------------------------------------------------------------------------
       u = time - (rad - rad0)
       call calcShape_Eppley( dff, u, amp, sigma, zbef, npow, norm_packet )
       coa    = 3 * (dff(2) / rad3 + 3 * dff(1) / rad4 + 3 * dff(0) / rad5)
       cob    = -(dff(3) / rad2 + 3 * dff(2) / rad3 + 6 * dff(1) / rad4           &
                + 6 * dff(0) / rad5)
       coc    = (dff(4) / rad + 2 * dff(3) / rad2 + 9 * dff(2) / rad3             &
                + 21 * dff(1) / rad4 + 21 * dff(0) / rad5) / 4
       col    = dff(3)/rad + 2 * dff(2) / rad2 + 3 * dff(1)/rad3 + 3 * dff(0) / rad4
       cok    = coa * rad / 3

       ! Absorb this annoying factor of -1/2 in the EK in the coefficients here
       coadot = -3 * (dff(3) / rad3 + 3 * dff(2) / rad4 + 3 * dff(1) / rad5) / 2
       cobdot = (dff(4) / rad2 + 3 * dff(3) / rad3 + 6 * dff(2) / rad4            &
                + 6 * dff(1) / rad5) / 2
       cocdot = -(dff(5) / rad + 2 * dff(4) / rad2 + 9 * dff(3) / rad3            &
                + 21 * dff(2) / rad4 + 21 * dff(1) / rad5) / 8
       coldot = - (dff(4)/rad + 2 * dff(3) / rad2 + 3 * dff(2)/rad3               & 
                + 3 * dff(1) / rad4)
       cokdot = coadot * rad / 3
       !---------------------------------------------------------------------------

   else if( CCTK_EQUALS( wave_dir, "ingoing" ) ) then

       ! Subtract ingoing part
       !---------------------------------------------------------------------------
       u = time + (rad - rad0)
       call calcShape_Eppley( dff, u, amp, sigma, zbef, npow, norm_packet )
       coa    = 3 * (dff(2) / rad3 - 3 * dff(1) / rad4 + 3 * dff(0) / rad5)
       cob    = - (-dff(3) / rad2 + 3 * dff(2) / rad3 - 6 * dff(1) / rad4     &
                + 6 * dff(0) / rad5)
       coc    = (dff(4) / rad - 2 * dff(3) / rad2 + 9 * dff(2) / rad3       &
                - 21 * dff(1) / rad4 + 21 * dff(0) / rad5) / 4
       cok = rad * coa / 3
       col = (-dff(3)/rad + 2 * dff(2) / rad2 - 3 * dff(1)/rad3 + 3 * dff(0) / rad4)

       ! Absorb this annoying factor of -1/2 in the EK in the coefficients here
       coadot = - 3 * (dff(3) / rad3 - 3 * dff(2) / rad4 + 3 * dff(1)      &
                / rad5) / 2
       cobdot = (-dff(4) / rad2 + 3 * dff(3) / rad3 - 6 * dff(2) / rad4  &
                + 6 * dff(1) / rad5) / 2
       cocdot = - (dff(5) / rad - 2 * dff(4) / rad2 + 9 * dff(3) / rad3    &
                - 21 * dff(2) / rad4 + 21 * dff(1) / rad5) / 8
       cokdot = rad * coadot / 3
       coldot = dff(4)/rad - 2 * dff(3) / rad2 + 3 * dff(2)/rad3 - 3 * dff(1) / rad4
       !---------------------------------------------------------------------------

    else if( CCTK_EQUALS( wave_dir, "timesym" ) ) then

       ! Outgoing part
       !---------------------------------------------------------------------------
       u = time - (rad - rad0)
       call calcShape_Eppley( dff, u, amp, sigma, zbef, npow, norm_packet )
       coa    = 3 * (dff(2) / rad3 + 3 * dff(1) / rad4 + 3 * dff(0) / rad5)
       cob    = -(dff(3) / rad2 + 3 * dff(2) / rad3 + 6 * dff(1) / rad4           &
                + 6 * dff(0) / rad5)
       coc    = (dff(4) / rad + 2 * dff(3) / rad2 + 9 * dff(2) / rad3             &
                + 21 * dff(1) / rad4 + 21 * dff(0) / rad5) / 4
       col    = dff(3)/rad + 2 * dff(2) / rad2 + 3 * dff(1)/rad3 + 3 * dff(0) / rad4

       ! Absorb this annoying factor of -1/2 in the EK in the coefficients here
       coadot = -3 * (dff(3) / rad3 + 3 * dff(2) / rad4 + 3 * dff(1) / rad5) / 2
       cobdot = (dff(4) / rad2 + 3 * dff(3) / rad3 + 6 * dff(2) / rad4            &
                + 6 * dff(1) / rad5) / 2
       cocdot = -(dff(5) / rad + 2 * dff(4) / rad2 + 9 * dff(3) / rad3            &
                + 21 * dff(2) / rad4 + 21 * dff(1) / rad5) / 8
       coldot = - (dff(4)/rad + 2 * dff(3) / rad2 + 3 * dff(2)/rad3               & 
                + 3 * dff(1) / rad4)
       !---------------------------------------------------------------------------


       ! Subtract ingoing part
       !---------------------------------------------------------------------------
       u = time + (rad - rad0)
       call calcShape_Eppley( dff, u, amp, sigma, zbef, npow, norm_packet)
       coa    = coa - 3 * (dff(2) / rad3 - 3 * dff(1) / rad4 + 3 * dff(0) / rad5)
       cob    = cob + (-dff(3) / rad2 + 3 * dff(2) / rad3 - 6 * dff(1) / rad4     &
                + 6 * dff(0) / rad5)
       coc    = coc - (dff(4) / rad - 2 * dff(3) / rad2 + 9 * dff(2) / rad3       &
                - 21 * dff(1) / rad4 + 21 * dff(0) / rad5) / 4
       cok = rad * coa / 3
       col = col - (-dff(3)/rad + 2 * dff(2) / rad2 - 3 * dff(1)/rad3 + 3 * dff(0) / rad4)

       ! Absorb this annoying factor of -1/2 in the EK in the coefficients here
       coadot = coadot + 3 * (dff(3) / rad3 - 3 * dff(2) / rad4 + 3 * dff(1)      &
                / rad5) / 2
       cobdot = cobdot - (-dff(4) / rad2 + 3 * dff(3) / rad3 - 6 * dff(2) / rad4  &
                + 6 * dff(1) / rad5) / 2
       cocdot = cocdot + (dff(5) / rad - 2 * dff(4) / rad2 + 9 * dff(3) / rad3    &
                - 21 * dff(2) / rad4 + 21 * dff(1) / rad5) / 8
       cokdot = rad * coadot / 3
       coldot = coldot - dff(4)/rad + 2 * dff(3) / rad2 - 3 * dff(2)/rad3 + 3 * dff(1) / rad4
       !---------------------------------------------------------------------------
  else
       call CCTK_WARN (0, "Parameter wave_dir doesn't have a legit value!") 
  end if


end subroutine calcABC_Eppley
!
!==============================================================================
!
end module TeukolskyWaveTools
!
!==============================================================================
!
subroutine calcTeukolskyWave( time, xx, gg, kk, alph, beta )

  use TeukolskyWaveTools
  implicit none
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  CCTK_REAL, intent(in)                   :: time
  CCTK_REAL, intent(in),  dimension(3)    :: xx
  CCTK_REAL, intent(out), dimension(3,3)  :: gg, kk
  CCTK_REAL, intent(out)                  :: alph
  CCTK_REAL, intent(out), dimension(3)    :: beta

  CCTK_REAL                                  rad, efftime
  CCTK_REAL                                  coa, cob, coc, coadot, cobdot,  &
                                             cocdot, col, cok, coldot, cokdot
  efftime = time - rad0
    
  ! Eppley packet
  if( CCTK_EQUALS( ffshape, "eppley" ) ) then
    if( even == 0 .and. mm /= 0) then
      call CCTK_WARN (0, "Sorry, only m=0 implemented for odd Teukolsky wave ...")
    else if( even == 0 ) then
      call CCTK_WARN (1, "Issues at the origin with mm=0, odd TW!")
    else
      ! r=0 requires special treatment
      rad = sqrt( sum( xx*xx ) )
      if( rad < 1.0e-10 ) then
        call calcMetricOrigin_Eppley(gg, kk, efftime, amp, sigma, mm, even, zbef, npow, norm_packet)
      else
        call calcABC_Eppley( coa, cob, coc, col, cok, coadot, cobdot, cocdot, &
                             coldot, cokdot, time, xx)
        call calcMetric( gg, kk, coa, cob, coc, col, cok, coadot, cobdot, cocdot, &
                         coldot, cokdot, xx, mm, even )
      end if
    end if
  end if

  !if (xx(2) > 0.55 .and. xx(2) < 0.65) then
  !write(*,*) 'xyz = ', xx(1), xx(2), xx(3)
  !write(*,*) 'g(1,:) = ', gg(1,:)
  !write(*,*) 'g(2,:) = ', gg(2,:)
  !write(*,*) 'g(3,:) = ', gg(3,:)
  !end if

  beta = 0
  alph = 1
  
end subroutine calcTeukolskyWave
