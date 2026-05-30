#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

!----------------------------------------------------------------------------
!
! Defines the background metric used in LWT. 
! Needs to initialize the following GF:
!
!                 + 3-metric, indices down. The components need to be named 
!                   gxx,gxy,gxz,gyy,gyz,gzz. The other (symmetric) components 
!                   are hopefully defined elsewhere in terms of these. 
!
!                 + lapse, needs to be named alpha
!
!                 + shift components, indices up. Need to be named
!                   betax, betay, betaz
!
!
!-------------------------------------------------------------------------------


subroutine LWT_init_metric (CCTK_ARGUMENTS)
  use constants
  use lapack
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  integer   :: i, j, k

  !----------------------------------------------------------
  ! local, pointwise variables
  ! convention:
  !           gg: 4-metric, indices down
  !           gama: 3-metric, indices down
  !           betal: shift, indices down
  !           beta: shift, indices up
  !           beta2: beta^i beta_i
  !           alfa: lapse 
  !-----------------------------------------------------------

  CCTK_REAL :: xx, yy, zz, rh, rr
  CCTK_REAL :: hh, ll(0:3), gg(0:3,0:3), ggtmp(0:3,0:3), jac(0:3,0:3)
  CCTK_REAL :: gama(3,3), betal(3), beta(3), beta2, alfa

  CCTK_REAL :: mat(3,3), rhs(3,1), masstmp
  integer   :: ipiv(3), info

  integer   :: a, b, c, d

  character :: msg*100



  ! Metric
  if (CCTK_EQUALS (metric, "Minkowski")) then

     gxx = 1
     gxy = 0
     gxz = 0
     gyy = 1
     gyz = 0
     gzz = 1

     alpha = lapse

     betax = shift(1) - shift_omega * y
     betay = shift(2) + shift_omega * x
     betaz = shift(3)

  else if (CCTK_EQUALS (metric, "Kerr-Schild")) then

     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)

              ! The event horizon is at r = M + sqrt (M^2 - a^2)
              ! with  r^4 - r^2 (rho^2 - a^2) - a^2 z^2 = 0
              ! or    rho^2 = r^2 + a^2 (1 - z^2 / r^2)
              ! where rho^2 = x^2 + y^2 + z^2

              xx = x(i,j,k)
              yy = y(i,j,k)
              zz = z(i,j,k)
              rh = r(i,j,k)

              rr = sqrt(+ (rh**2 - spin**2) / 2.0 &
                   &    + sqrt(+ (rh**2 - spin**2)**2 / 4.0 &
                   &           + spin**2 * zz**2))

              ! masstmp = mass * (1.0-exp(-rr**4/0.7**4))      
              masstmp = mass
              hh = masstmp * rr**3 / (rr**4 + spin**2 * zz**2)

              ll(0) = 1
              ll(1) = (rr*xx + spin*yy) / (rr**2 + spin**2)
              ll(2) = (rr*yy - spin*xx) / (rr**2 + spin**2)
              ll(3) = zz / rr

              do a=0,3
                 do b=0,3   ! eta is defined in TAT/TGRtensor
										if (rr<0.00001) then
										  gg(a,b) = eta4(a,b)
										else
                      gg(a,b) = eta4(a,b) + 2.0 * hh * ll(a) * ll(b)
									  end if 
                 end do
              end do

              do a=1,3
                 do b=1,3
                    gama(a,b) = gg(a,b)
                 end do
              end do

              do a=1,3
                 betal(a) = gg(0,a)
              end do

              ! getting the shift with indices up, take a look at TAT/TGRtensor
              do a=1,3
                 do b=1,3
                    mat(a,b) = gama(a,b)
                 end do
                 rhs(a,1) = betal(a)
              end do
              call gesv (3, 1, mat, 3, ipiv, rhs, 3, info)
              if (info /= 0) then
                 write (msg, '("Error in call to DGESV, info=",i4)') info
                 call CCTK_WARN (0, msg)
              end if
              do a=1,3
                 beta(a) = rhs(a,1)
              end do

              ! beta contracted with itself
              beta2 = 0.0 
              do a=1,3
                 beta2 = beta2 + betal(a) * beta(a)
              end do

              ! write pointwise vars into GF
              gxx(i,j,k) = gama(1,1)
              gxy(i,j,k) = gama(1,2)
              gxz(i,j,k) = gama(1,3)
              gyy(i,j,k) = gama(2,2)
              gyz(i,j,k) = gama(2,3)
              gzz(i,j,k) = gama(3,3)

              betax(i,j,k) = beta(1)
              betay(i,j,k) = beta(2)
              betaz(i,j,k) = beta(3)

              ! g_00 = - alpha^2 + beta^2
              ! alpha^2 = beta^2 - g_00
              alfa = sqrt(beta2 - gg(0,0))
              alpha(i,j,k) = alfa

           end do
        end do
     end do

  else if (CCTK_EQUALS (metric, "Kerr")) then

     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)

              ! The event horizon is at r = M + sqrt (M^2 - a^2)
              ! with  r^4 - r^2 (rho^2 - a^2) - a^2 z^2 = 0
              ! or    rho^2 = r^2 + a^2 (1 - z^2 / r^2)
              ! where rho^2 = x^2 + y^2 + z^2

              xx = x(i,j,k) - spin*y(i,j,k)/r(i,j,k)
              yy = y(i,j,k) + spin*x(i,j,k)/r(i,j,k)
              zz = z(i,j,k)
              rh = r(i,j,k)
              ! rh = xx**2 + yy**2 + zz**2

							rr = rh
              ! rr = sqrt(+ (rh**2 - spin**2) / 2.0 &
              !     &    + sqrt(+ (rh**2 - spin**2)**2 / 4.0 &
              !     &           + spin**2 * zz**2))

              hh = mass * rr**3 / (rr**4 + spin**2 * z(i,j,k)**2)
              ll(0) = 1
              ll(1) = (rr*xx + spin*yy) / (rr**2 + spin**2)
              ll(2) = (rr*yy - spin*xx) / (rr**2 + spin**2)
              ll(3) = zz / rr

              do a=0,3
                 do b=0,3   ! eta is defined in TAT/TGRtensor
                    gg(a,b) = eta4(a,b) + 2.0 * hh * ll(a) * ll(b)
                 end do
              end do

			 !=============================================================== 
			 != Transform Metric to Kerr Cartesian Coordinates
			 !=============================================================== 

        jac(0,0) = 1
			  jac(0,1) = 0
			  jac(0,2) = 0
			  jac(0,3) = 0
        
				jac(1,0) = 0 
			  jac(1,1) = 1 + (spin*x(i,j,k)*y(i,j,k)) / (rh**3)
			  jac(1,2) = - spin * (x(i,j,k)**2+z(i,j,k)**2) / (rh**3)
			  jac(1,3) = spin*y(i,j,k)*z(i,j,k) / (rh**3)
        
				jac(2,0) = 0
			  jac(2,1) = spin * (y(i,j,k)**2+z(i,j,k)**2) / (rh**3)
			  jac(2,2) = 1 - spin*x(i,j,k)*y(i,j,k) / (rh**3)
			  jac(2,3) = - spin*x(i,j,k)*z(i,j,k) / (rh**3)
        
				jac(3,0) = 0
			  jac(3,1) = 0
			  jac(3,2) = 0
			  jac(3,3) = 1
			 
            

			  do a=0,3
			  	do b=0,3
					ggtmp(a,b) = gg(a,b)
					gg(a,b) = 0
				end do
			  end do
			
			  do a=0,3
			    do b=0,3
				  do c=0,3
				    do d=0,3

						gg(a,b) = gg(a,b) + ggtmp(c,d)*jac(c,a)*jac(d,b) 					  
					
					end do
				  end do
				end do
			  end do

              do a=1,3
                 do b=1,3
                    gama(a,b) = gg(a,b)
                 end do
              end do

              do a=1,3
                 betal(a) = gg(0,a)
              end do

              ! getting the shift with indices up, take a look at TAT/TGRtensor
              do a=1,3
                 do b=1,3
                    mat(a,b) = gama(a,b)
                 end do
                 rhs(a,1) = betal(a)
              end do
              call gesv (3, 1, mat, 3, ipiv, rhs, 3, info)
              if (info /= 0) then
                 write (msg, '("Error in call to DGESV, info=",i4)') info
                 call CCTK_WARN (0, msg)
              end if
              do a=1,3
                 beta(a) = rhs(a,1)
              end do

              ! beta contracted with itself
              beta2 = 0.0 
              do a=1,3
                 beta2 = beta2 + betal(a) * beta(a)
              end do

              ! write pointwise vars into GF
              gxx(i,j,k) = gama(1,1)
              gxy(i,j,k) = gama(1,2)
              gxz(i,j,k) = gama(1,3)
              gyy(i,j,k) = gama(2,2)
              gyz(i,j,k) = gama(2,3)
              gzz(i,j,k) = gama(3,3)

              betax(i,j,k) = beta(1)
              betay(i,j,k) = beta(2)
              betaz(i,j,k) = beta(3)

              ! g_00 = - alpha^2 + beta^2
              ! alpha^2 = beta^2 - g_00
              alfa = sqrt(beta2 - gg(0,0))
              alpha(i,j,k) = alfa

           end do
        end do
     end do



  else

     call CCTK_WARN (0, "internal error")

  end if                        ! if metric
  
end subroutine LWT_init_metric
