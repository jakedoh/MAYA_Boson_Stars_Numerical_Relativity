#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine exactsolution (u_pt, rho_pt, x_pt, y_pt, z_pt, t_pt) 
	implicit none
	DECLARE_CCTK_FUNCTIONS
	DECLARE_CCTK_PARAMETERS
        CCTK_REAL, intent(out) :: u_pt, rho_pt
	CCTK_REAL, intent(in) :: x_pt, y_pt, z_pt, t_pt

        CCTK_REAL :: theta_pt, phi_pt
        CCTK_REAL, parameter :: one = 1, half = one / 2
	CCTK_REAL :: nx, ny, nz, nt
	CCTK_REAL :: pi
	CCTK_REAL :: omega
	CCTK_REAL :: w
        CCTK_REAL :: r1, rpt, rmt
        CCTK_REAL :: asmallnumber
        CCTK_REAL :: dummy
        CCTK_REAL :: ylm, tmptmp

        asmallnumber = 0.0000001

	if (CCTK_EQUALS (initial_data, "linear")) then
              u_pt   = x_pt + t_pt * y_pt
              rho_pt = y_pt

        else if (CCTK_EQUALS (initial_data, "plane")) then
              pi = 4 * atan(one)
              omega = sqrt(sum(wave_number**2))
              
              nx = wave_number(1) * (x_pt - space_offset(1))
              ny = wave_number(2) * (y_pt - space_offset(2))
              nz = wave_number(3) * (z_pt - space_offset(3))
              nt = omega          * (t_pt        - time_offset   )
              w = 2*pi * (nx + ny + nz + nt)
              
              u_pt   =   amplitude * cos(w)
              
              rho_pt = - amplitude * sin(w) * 2*pi * omega

        else if (CCTK_EQUALS (initial_data, "Gaussian")) then
              r1 = sqrt(x_pt**2 + y_pt**2 + z_pt**2 )
              rpt = r1 + t_pt
              rmt = r1 - t_pt
            
              if (r1 < asmallnumber) then
                      u_pt   =  Gaussian(t_pt)*(1.0-2.0*t_pt**2/width**2)
                      rho_pt = DGaussian(t_pt) * (1.0-2.0*t_pt**2/width**2) - 4.0 * Gaussian(t_pt) * t_pt / width**2
              else

                      u_pt = half * (rpt/r1 * Gaussian(rpt) + rmt/r1 * Gaussian(rmt))
              
                      rho_pt = half * (+ 1/r1 * Gaussian(rpt) + rpt/r1 * DGaussian(rpt) &
                             &         - 1/r1 * Gaussian(rmt) - rmt/r1 * DGaussian(rmt))
              endif
            else if (CCTK_EQUALS (initial_data, "GeneralMultipole")) then
              r1 = sqrt(x_pt**2 + y_pt**2 + z_pt**2 )
              rpt = r1 + t_pt
              rmt = r1 - t_pt
            
              if (r1 < asmallnumber) then
                      u_pt   =  Gaussian(t_pt)*t_pt**2 * (-2.0*t_pt**2+3.0*width**2)/width**2
                      rho_pt =  2.0*Gaussian(t_pt)*t_pt*(2.0*t_pt**4-7.0*t_pt**2*width**2+3*width**4)/width**4
              else

                      u_pt = half * (rpt**3/r1 * Gaussian(rpt) + rmt**3/r1 * Gaussian(rmt))
              
                      rho_pt = half/r1 * (+ 3.0*rpt**2 * Gaussian(rpt) - 2.0*rpt**4/width**2 * Gaussian(rpt) &
                             &         - 3.0*rmt**2 * Gaussian(rmt) + 2.0*rmt**4/width**2 * Gaussian(rmt))
              endif
              if (r1 < asmallnumber) then
                theta_pt = 0.0
              else
                theta_pt = acos(z_pt/r1) 
              endif
              if (sqrt(x_pt**2+y_pt**2) < asmallnumber) then
                phi_pt = 0.0
              else
                phi_pt = asin(y_pt/sqrt(x_pt**2 + y_pt**2))
              endif
              call sylm_re(multipole_s, multipole_l, multipole_m, theta_pt, phi_pt, ylm) 
              u_pt = u_pt * ylm
              rho_pt = rho_pt * ylm
! 	else
! 	  call CCTK_WARN(0, "Exact Solution for chosen initial data type does not exist")
	endif 

contains
  
  CCTK_REAL function Gaussian (rr)
    CCTK_REAL, intent(in) :: rr
    Gaussian = amplitude * exp(- ((rr - radius) / width)**2)
  end function Gaussian
  
  CCTK_REAL function DGaussian (rr)
    CCTK_REAL, intent(in) :: rr
    DGaussian = -2 * Gaussian(rr) * (rr - radius) / width**2
  end function DGaussian
 
end subroutine exactsolution

subroutine LWT_init (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  CCTK_REAL, parameter :: one = 1, half = one / 2
  CCTK_REAL :: pi
  CCTK_REAL :: t
  
  integer :: i, j, k
  
  CCTK_REAL :: omega, nx, ny, nz, nt, w
  CCTK_REAL :: r1, rpt, rmt, u0, rho0, vr0, v0(3), vr, ur, rxy
  CCTK_REAL :: Y02, cosTheta, sinTheta, cosTheta2
  CCTK_REAL :: dThetadx, dThetady, dThetadz, vTheta
  CCTK_REAL :: Y22, sinPhi, cosPhi, cos2Phi, cos3Phi, cos4Phi, sin2Phi, sin3Phi, sin4Phi 
  CCTK_REAL :: Y21, Y10, Y11, Y40, Y41, Y42, Y43, Y44  
  
  pi = 4 * atan(one)
  
  t = cctk_time
  
  
  
  ! Initial data
  if (CCTK_EQUALS (initial_data, "linear")) then
     
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              
               call exactsolution(u(i,j,k), rho(i,j,k), x(i,j,k), y(i,j,k), z(i,j,k), t)
             
           end do
        end do
     end do
     
  else if (CCTK_EQUALS (initial_data, "plane")) then
     
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              
               call exactsolution(u(i,j,k), rho(i,j,k), x(i,j,k), y(i,j,k), z(i,j,k), t)
             
           end do
        end do
     end do

  else if (CCTK_EQUALS (initial_data, "Debug")) then
     
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)

              u(i,j,k) = MultiPatch_GetMap(cctkGH)
              rho(i,j,k) = i-1 + cctk_gsh(1) * (j-1 + cctk_gsh(2) * (k-1))
              
           end do
        end do
     end do
      
  else if (CCTK_EQUALS (initial_data, "Gaussian")) then

      do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              
               call exactsolution(u(i,j,k), rho(i,j,k), x(i,j,k), y(i,j,k), z(i,j,k), t)
             
           end do
        end do
     end do

   else if (CCTK_EQUALS (initial_data, "GeneralMultipole")) then

      do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              
               call exactsolution(u(i,j,k), rho(i,j,k), x(i,j,k), y(i,j,k), z(i,j,k), t)
             
           end do
        end do
     end do
    
    
  else if (CCTK_EQUALS (initial_data, "GaussianNonLinear")) then
     
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              
              r1 = sqrt(epsx*x(i,j,k)**2 + epsy*y(i,j,k)**2 + z(i,j,k)**2 + eps**2)
              
              vr = ANL * exp(-((r1-RNL)/deltaNL)**2)
              
              u(i,j,k) = vr
              
              ur = -2*(r1-RNL)/deltaNL**2 * vr
              
              rho(i,j,k) = muNL*ur + omeNL*y(i,j,k)*x(i,j,k)*ur*(epsx-epsy)/r1
              
              vx(i,j,k) = x(i,j,k)/r1 * ur * epsx
              vy(i,j,k) = y(i,j,k)/r1 * ur * epsy
              vz(i,j,k) = z(i,j,k)/r1 * ur
              
           end do
        end do
     end do
     
  else if (CCTK_EQUALS (initial_data, "multipole")) then
     
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              
              r1 = sqrt(r(i,j,k)**2 + eps**2)
              rpt = r1 + t
              rmt = r1 - t
              
              u0 = half * (rpt/r1 * Gaussian(rpt) + rmt/r1 * Gaussian(rmt))
              
              rho0 = half * (+ 1/r1 * Gaussian(rpt) + rpt/r1 * DGaussian(rpt) &
                   &         - 1/r1 * Gaussian(rmt) - rmt/r1 * DGaussian(rmt))
              
              vr0 = half &
                   * (- t/r1**2 * Gaussian(rpt) + rpt/r1 * DGaussian(rpt) &
                   &  + t/r1**2 * Gaussian(rmt) + rmt/r1 * DGaussian(rmt))
              
              u(i,j,k) = u0 * z(i,j,k)/r1
              
              rho(i,j,k) = rho0 * z(i,j,k)/r1
              
              vr = vr0 * z(i,j,k)/r1 - u(i,j,k)/r1
              
              vx(i,j,k) = x(i,j,k)/r1 * vr
              vy(i,j,k) = y(i,j,k)/r1 * vr
              vz(i,j,k) = z(i,j,k)/r1 * vr + u0/r1
              
           end do
        end do
     end do
     
  else if (CCTK_EQUALS (initial_data, "multipole l=2")) then
     
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              
              r1 = sqrt(r(i,j,k)**2 + eps**2)
              rpt = r1 + t
              rmt = r1 - t
              
              u0 = half * (rpt/r1 * Gaussian(rpt) + rmt/r1 * Gaussian(rmt))
              
              rho0 = half * (+ 1/r1 * Gaussian(rpt) + rpt/r1 * DGaussian(rpt) &
                   &         - 1/r1 * Gaussian(rmt) - rmt/r1 * DGaussian(rmt))
              
              vr0 = half &
                   * (- t/r1**2 * Gaussian(rpt) + rpt/r1 * DGaussian(rpt) &
                   &  + t/r1**2 * Gaussian(rmt) + rmt/r1 * DGaussian(rmt))
              
              cosTheta = z(i,j,k)/r1
              sinTheta = sqrt(x(i,j,k)**2 + y(i,j,k)**2)/r1
              cosTheta2 = cosTheta**2
             
              Y02 = sqrt(5/(4*pi)) * half * (3*cosTheta2-1)

              u(i,j,k) = u0 * Y02
              
              rho(i,j,k) = rho0 * Y02
              
             
              dthetadx = z(i,j,k)*x(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetady = z(i,j,k)*y(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetadz = z(i,j,k)**2 &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2)) &
                   - 1/sqrt(x(i,j,k)**2+y(i,j,k)**2)
              
              vr = vr0 * Y02
              vtheta = - u0 * sqrt(5/(4*pi))*3*cosTheta*sinTheta
              
              vx(i,j,k) = x(i,j,k)/r1*vr + dthetadx*vtheta;
              vy(i,j,k) = y(i,j,k)/r1*vr + dthetady*vtheta;
              vz(i,j,k) = z(i,j,k)/r1*vr + dthetadz*vtheta; 
              
           end do
        end do
     end do

  else if (CCTK_EQUALS (initial_data, "multipole l=2, u=0")) then
     
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              
              r1 = sqrt(r(i,j,k)**2 + eps**2)
              rpt = r1 + t
              rmt = r1 - t
              
              u0 = 0.0
              
              rho0 = half * (rpt/r1 * Gaussian(rpt) + rmt/r1 * Gaussian(rmt))
              
              vr0 = 0.0
              
              cosTheta = z(i,j,k)/r1
              sinTheta = sqrt(x(i,j,k)**2 + y(i,j,k)**2)/r1
              cosTheta2 = cosTheta**2
             
              Y02 = sqrt(5/(4*pi)) * half * (3*cosTheta2-1)

              u(i,j,k) = u0 * Y02
              
              rho(i,j,k) = rho0 * Y02
              
             
              dthetadx = z(i,j,k)*x(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetady = z(i,j,k)*y(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetadz = z(i,j,k)**2 &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2)) &
                   - 1/sqrt(x(i,j,k)**2+y(i,j,k)**2)
              
              vr = vr0 * Y02
              vtheta = - u0 * sqrt(5/(4*pi))*3*cosTheta*sinTheta
              
              vx(i,j,k) = x(i,j,k)/r1*vr + dthetadx*vtheta;
              vy(i,j,k) = y(i,j,k)/r1*vr + dthetady*vtheta;
              vz(i,j,k) = z(i,j,k)/r1*vr + dthetadz*vtheta; 
              
           end do
        end do
     end do

     else if (CCTK_EQUALS (initial_data, "multipole l=1, m=0")) then
     
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              
              r1 = sqrt(r(i,j,k)**2 + eps**2)
              rpt = r1 + t
              rmt = r1 - t
              
              u0 = 0.0
              
              rho0 = half * (rpt/r1 * Gaussian(rpt) + rmt/r1 * Gaussian(rmt))
              
              vr0 = 0.0
              
              cosTheta = z(i,j,k)/r1
              sinTheta = sqrt(x(i,j,k)**2 + y(i,j,k)**2)/r1
              cosTheta2 = cosTheta**2
			  sinPhi = y(i,j,k) / sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k))
			  cosPhi = x(i,j,k) / sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k))
             
              Y10 =  sqrt(3./(4.*pi))*cosTheta

              u(i,j,k) = u0 * Y10
              
              rho(i,j,k) = rho0 * Y10
              
             
              dthetadx = z(i,j,k)*x(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetady = z(i,j,k)*y(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetadz = z(i,j,k)**2 &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2)) &
                   - 1/sqrt(x(i,j,k)**2+y(i,j,k)**2)
              
              vr = vr0 * Y10
              vtheta = - u0 * sqrt(5/(4*pi))*3*cosTheta*sinTheta
              
              vx(i,j,k) = x(i,j,k)/r1*vr + dthetadx*vtheta;
              vy(i,j,k) = y(i,j,k)/r1*vr + dthetady*vtheta;
              vz(i,j,k) = z(i,j,k)/r1*vr + dthetadz*vtheta; 
              
           end do
        end do
     end do

     else if (CCTK_EQUALS (initial_data, "multipole l=1, m=1")) then
     
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              
              r1 = sqrt(r(i,j,k)**2 + eps**2)
              rpt = r1 + t
              rmt = r1 - t
              
              u0 = 0.0
              
              rho0 = half * (rpt/r1 * Gaussian(rpt) + rmt/r1 * Gaussian(rmt))
              
              vr0 = 0.0
              
              cosTheta = z(i,j,k)/r1
              sinTheta = sqrt(x(i,j,k)**2 + y(i,j,k)**2)/r1
              cosTheta2 = cosTheta**2
			  sinPhi = y(i,j,k) / sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k))
			  cosPhi = x(i,j,k) / sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k))
             
              Y11 =  -0.5*sqrt(3./(2.*pi))*sinTheta*cosPhi

              u(i,j,k) = u0 * Y11
              
              rho(i,j,k) = rho0 * Y11
              
             
              dthetadx = z(i,j,k)*x(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetady = z(i,j,k)*y(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetadz = z(i,j,k)**2 &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2)) &
                   - 1/sqrt(x(i,j,k)**2+y(i,j,k)**2)
              
              vr = vr0 * Y11
              vtheta = - u0 * sqrt(5/(4*pi))*3*cosTheta*sinTheta
              
              vx(i,j,k) = x(i,j,k)/r1*vr + dthetadx*vtheta;
              vy(i,j,k) = y(i,j,k)/r1*vr + dthetady*vtheta;
              vz(i,j,k) = z(i,j,k)/r1*vr + dthetadz*vtheta; 
              
           end do
        end do
     end do


     else if (CCTK_EQUALS (initial_data, "multipole l=1, m=-1")) then
     
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              
              r1 = sqrt(r(i,j,k)**2 + eps**2)
              rpt = r1 + t
              rmt = r1 - t
              
              u0 = 0.0
              
              rho0 = half * (rpt/r1 * Gaussian(rpt) + rmt/r1 * Gaussian(rmt))
              
              vr0 = 0.0
              
              cosTheta = z(i,j,k)/r1
              sinTheta = sqrt(x(i,j,k)**2 + y(i,j,k)**2)/r1
              cosTheta2 = cosTheta**2
			  sinPhi = y(i,j,k) / sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k))
			  cosPhi = x(i,j,k) / sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k))
             
              Y11 =  -0.5*sqrt(3./(2.*pi))*sinTheta*sinPhi

              u(i,j,k) = u0 * Y11
              
              rho(i,j,k) = rho0 * Y11
              
             
              dthetadx = z(i,j,k)*x(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetady = z(i,j,k)*y(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetadz = z(i,j,k)**2 &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2)) &
                   - 1/sqrt(x(i,j,k)**2+y(i,j,k)**2)
              
              vr = vr0 * Y11
              vtheta = - u0 * sqrt(5/(4*pi))*3*cosTheta*sinTheta
              
              vx(i,j,k) = x(i,j,k)/r1*vr + dthetadx*vtheta;
              vy(i,j,k) = y(i,j,k)/r1*vr + dthetady*vtheta;
              vz(i,j,k) = z(i,j,k)/r1*vr + dthetadz*vtheta; 
              
           end do
        end do
     end do



     else if (CCTK_EQUALS (initial_data, "multipole l=4, m=0")) then
     
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              
              r1 = sqrt(r(i,j,k)**2 + eps**2)
              rpt = r1 + t
              rmt = r1 - t
              
              u0 = 0.0
              
              rho0 = half * (rpt/r1 * Gaussian(rpt) + rmt/r1 * Gaussian(rmt))
              
              vr0 = 0.0
              
              cosTheta = z(i,j,k)/r1
              sinTheta = sqrt(x(i,j,k)**2 + y(i,j,k)**2)/r1
              cosTheta2 = cosTheta**2
			  sinPhi = y(i,j,k) / sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k))
			  cosPhi = x(i,j,k) / sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k))
             
              Y40 =  1./sqrt(pi)*(9./16. - 45.*cosTheta*cosTheta/8. + 105.*cosTheta**4/16.)

              u(i,j,k) = u0 * Y40
              
              rho(i,j,k) = rho0 * Y40
              
             
              dthetadx = z(i,j,k)*x(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetady = z(i,j,k)*y(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetadz = z(i,j,k)**2 &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2)) &
                   - 1/sqrt(x(i,j,k)**2+y(i,j,k)**2)
              
              vr = vr0 * Y40
              vtheta = - u0 * sqrt(5/(4*pi))*3*cosTheta*sinTheta
              
              vx(i,j,k) = x(i,j,k)/r1*vr + dthetadx*vtheta;
              vy(i,j,k) = y(i,j,k)/r1*vr + dthetady*vtheta;
              vz(i,j,k) = z(i,j,k)/r1*vr + dthetadz*vtheta; 
              
           end do
        end do
     end do


     else if (CCTK_EQUALS (initial_data, "multipole l=4, m=1")) then
     
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              
              r1 = sqrt(r(i,j,k)**2 + eps**2)
              rpt = r1 + t
              rmt = r1 - t
              
              u0 = 0.0
              
              rho0 = half * (rpt/r1 * Gaussian(rpt) + rmt/r1 * Gaussian(rmt)) 
              
              vr0 = 0.0
              
              cosTheta = z(i,j,k)/r1
              sinTheta = sqrt(x(i,j,k)**2 + y(i,j,k)**2)/r1
              cosTheta2 = cosTheta**2
			  sinPhi = y(i,j,k) / sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k))
			  cosPhi = x(i,j,k) / sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k))
             
              Y41 =  9./8.*sqrt(5./pi)*cosTheta*cosPhi*sinTheta - 21./8.*sqrt(5./pi)*cosTheta**3*cosPhi*sinTheta

              u(i,j,k) = u0 * Y41
              
              rho(i,j,k) = rho0 * Y41
              
             
              dthetadx = z(i,j,k)*x(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetady = z(i,j,k)*y(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetadz = z(i,j,k)**2 &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2)) &
                   - 1/sqrt(x(i,j,k)**2+y(i,j,k)**2)
              
              vr = vr0 * Y41
              vtheta = - u0 * sqrt(5/(4*pi))*3*cosTheta*sinTheta
              
              vx(i,j,k) = x(i,j,k)/r1*vr + dthetadx*vtheta;
              vy(i,j,k) = y(i,j,k)/r1*vr + dthetady*vtheta;
              vz(i,j,k) = z(i,j,k)/r1*vr + dthetadz*vtheta; 
              
           end do
        end do
     end do

    else if (CCTK_EQUALS (initial_data, "multipole l=4, m=-1")) then
     
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              
              r1 = sqrt(r(i,j,k)**2 + eps**2)
              rpt = r1 + t
              rmt = r1 - t
              
              u0 = 0.0
              
              rho0 = half * (rpt/r1 * Gaussian(rpt) + rmt/r1 * Gaussian(rmt))
              
              vr0 = 0.0
              
              cosTheta = z(i,j,k)/r1
              sinTheta = sqrt(x(i,j,k)**2 + y(i,j,k)**2)/r1
              cosTheta2 = cosTheta**2
			  sinPhi = y(i,j,k) / sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k))
			  cosPhi = x(i,j,k) / sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k))
             
              Y41 =  9./8.*sqrt(5./pi)*cosTheta*sinPhi*sinTheta - 21./8.*sqrt(5./pi)*cosTheta**3*sinPhi*sinTheta

              u(i,j,k) = u0 * Y41
              
              rho(i,j,k) = rho0 * Y41
              
             
              dthetadx = z(i,j,k)*x(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetady = z(i,j,k)*y(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetadz = z(i,j,k)**2 &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2)) &
                   - 1/sqrt(x(i,j,k)**2+y(i,j,k)**2)
              
              vr = vr0 * Y41
              vtheta = - u0 * sqrt(5/(4*pi))*3*cosTheta*sinTheta
              
              vx(i,j,k) = x(i,j,k)/r1*vr + dthetadx*vtheta;
              vy(i,j,k) = y(i,j,k)/r1*vr + dthetady*vtheta;
              vz(i,j,k) = z(i,j,k)/r1*vr + dthetadz*vtheta; 
              
           end do
        end do
     end do

    else if (CCTK_EQUALS (initial_data, "multipole l=4, m=2")) then
     
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              
              r1 = sqrt(r(i,j,k)**2 + eps**2)
              rpt = r1 + t
              rmt = r1 - t
              
              u0 = 0.0
              
              rho0 = half * (rpt/r1 * Gaussian(rpt) + rmt/r1 * Gaussian(rmt))
              
              vr0 = 0.0
              
              cosTheta = z(i,j,k)/r1
              sinTheta = sqrt(x(i,j,k)**2 + y(i,j,k)**2)/r1
              cosTheta2 = cosTheta**2
			  sinPhi = y(i,j,k) / sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k))
			  cosPhi = x(i,j,k) / sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k))
			  cos2Phi = cosPhi*cosPhi - sinPhi*sinPhi
             
              Y42 =  -3./8.*sqrt(5./(2.*pi))*cos2Phi*sinTheta*sinTheta + 21./8.*sqrt(5./(2.*pi))*cosTheta*cosTheta*cos2Phi*sinTheta*sinTheta

              u(i,j,k) = u0 * Y42
              
              rho(i,j,k) = rho0 * Y42
              
             
              dthetadx = z(i,j,k)*x(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetady = z(i,j,k)*y(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetadz = z(i,j,k)**2 &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2)) &
                   - 1/sqrt(x(i,j,k)**2+y(i,j,k)**2)
              
              vr = vr0 * Y42
              vtheta = - u0 * sqrt(5/(4*pi))*3*cosTheta*sinTheta
              
              vx(i,j,k) = x(i,j,k)/r1*vr + dthetadx*vtheta;
              vy(i,j,k) = y(i,j,k)/r1*vr + dthetady*vtheta;
              vz(i,j,k) = z(i,j,k)/r1*vr + dthetadz*vtheta; 
              
           end do
        end do
     end do

    else if (CCTK_EQUALS (initial_data, "multipole l=4, m=-2")) then
     
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              
              r1 = sqrt(r(i,j,k)**2 + eps**2)
              rpt = r1 + t
              rmt = r1 - t
              
              u0 = 0.0
              
              rho0 = half * (rpt/r1 * Gaussian(rpt) + rmt/r1 * Gaussian(rmt))
              
              vr0 = 0.0
              
              cosTheta = z(i,j,k)/r1
              sinTheta = sqrt(x(i,j,k)**2 + y(i,j,k)**2)/r1
              cosTheta2 = cosTheta**2
			  sinPhi = y(i,j,k) / sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k))
			  cosPhi = x(i,j,k) / sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k))
			  cos2Phi = cosPhi*cosPhi - sinPhi*sinPhi
             
              Y42 =  -3./8.*sqrt(5./(2.*pi))*sin2Phi*sinTheta*sinTheta + 21./8.*sqrt(5./(2.*pi))*cosTheta*cosTheta*sin2Phi*sinTheta*sinTheta

              u(i,j,k) = u0 * Y42
              
              rho(i,j,k) = rho0 * Y42
              
             
              dthetadx = z(i,j,k)*x(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetady = z(i,j,k)*y(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetadz = z(i,j,k)**2 &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2)) &
                   - 1/sqrt(x(i,j,k)**2+y(i,j,k)**2)
              
              vr = vr0 * Y42
              vtheta = - u0 * sqrt(5/(4*pi))*3*cosTheta*sinTheta
              
              vx(i,j,k) = x(i,j,k)/r1*vr + dthetadx*vtheta;
              vy(i,j,k) = y(i,j,k)/r1*vr + dthetady*vtheta;
              vz(i,j,k) = z(i,j,k)/r1*vr + dthetadz*vtheta; 
              
           end do
        end do
     end do

    else if (CCTK_EQUALS (initial_data, "multipole l=4, m=3")) then
     
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              
              r1 = sqrt(r(i,j,k)**2 + eps**2)
              rpt = r1 + t
              rmt = r1 - t
              
              u0 = 0.0
              
              rho0 = half * (rpt/r1 * Gaussian(rpt) + rmt/r1 * Gaussian(rmt))
              
              vr0 = 0.0
              
              cosTheta = z(i,j,k)/r1
              sinTheta = sqrt(x(i,j,k)**2 + y(i,j,k)**2)/r1
              cosTheta2 = cosTheta**2
			  sinPhi = y(i,j,k) / sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k))
			  cosPhi = x(i,j,k) / sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k))
			  cos2Phi = cosPhi*cosPhi - sinPhi*sinPhi
			  sin2Phi = 2.*sinPhi*cosPhi
              cos3Phi = cosPhi*cos2Phi - sinPhi*sin2Phi 
 
			
              Y43 =  -3./8.*sqrt(35./pi)*cosTheta*cos3Phi*sinTheta**3

              u(i,j,k) = u0 * Y43
              
              rho(i,j,k) = rho0 * Y43
              
             
              dthetadx = z(i,j,k)*x(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetady = z(i,j,k)*y(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetadz = z(i,j,k)**2 &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2)) &
                   - 1/sqrt(x(i,j,k)**2+y(i,j,k)**2)
              
              vr = vr0 * Y43
              vtheta = - u0 * sqrt(5/(4*pi))*3*cosTheta*sinTheta
              
              vx(i,j,k) = x(i,j,k)/r1*vr + dthetadx*vtheta;
              vy(i,j,k) = y(i,j,k)/r1*vr + dthetady*vtheta;
              vz(i,j,k) = z(i,j,k)/r1*vr + dthetadz*vtheta; 
              
           end do
        end do
     end do

    else if (CCTK_EQUALS (initial_data, "multipole l=4, m=-3")) then
     
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              
              r1 = sqrt(r(i,j,k)**2 + eps**2)
              rpt = r1 + t
              rmt = r1 - t
              
              u0 = 0.0
              
              rho0 = half * (rpt/r1 * Gaussian(rpt) + rmt/r1 * Gaussian(rmt))
              
              vr0 = 0.0
              
              cosTheta = z(i,j,k)/r1
              sinTheta = sqrt(x(i,j,k)**2 + y(i,j,k)**2)/r1
              cosTheta2 = cosTheta**2
			  sinPhi = y(i,j,k) / sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k))
			  cosPhi = x(i,j,k) / sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k))
			  cos2Phi = cosPhi*cosPhi - sinPhi*sinPhi
			  sin2Phi = 2.*sinPhi*cosPhi
              cos3Phi = cosPhi*cos2Phi - sinPhi*sin2Phi 

              sin3Phi = sinPhi*cos2Phi + sin2Phi*cosPhi 
			
              Y43 =  -3./8.*sqrt(35./pi)*cosTheta*sin3Phi*sinTheta**3

              u(i,j,k) = u0 * Y43
              
              rho(i,j,k) = rho0 * Y43
              
             
              dthetadx = z(i,j,k)*x(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetady = z(i,j,k)*y(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetadz = z(i,j,k)**2 &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2)) &
                   - 1/sqrt(x(i,j,k)**2+y(i,j,k)**2)
              
              vr = vr0 * Y43
              vtheta = - u0 * sqrt(5/(4*pi))*3*cosTheta*sinTheta
              
              vx(i,j,k) = x(i,j,k)/r1*vr + dthetadx*vtheta;
              vy(i,j,k) = y(i,j,k)/r1*vr + dthetady*vtheta;
              vz(i,j,k) = z(i,j,k)/r1*vr + dthetadz*vtheta; 
              
           end do
        end do
     end do

    else if (CCTK_EQUALS (initial_data, "multipole l=4, m=4")) then
     
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              
              r1 = sqrt(r(i,j,k)**2 + eps**2)
              rpt = r1 + t
              rmt = r1 - t
              
              u0 = 0.0
              
              rho0 = half * (rpt/r1 * Gaussian(rpt) + rmt/r1 * Gaussian(rmt))
              
              vr0 = 0.0
              
              cosTheta = z(i,j,k)/r1
              sinTheta = sqrt(x(i,j,k)**2 + y(i,j,k)**2)/r1
              cosTheta2 = cosTheta**2
			  sinPhi = y(i,j,k) / sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k))
			  cosPhi = x(i,j,k) / sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k))
			  cos2Phi = cosPhi*cosPhi - sinPhi*sinPhi
			  sin2Phi = 2.*sinPhi*cosPhi
              cos3Phi = cosPhi*cos2Phi - sinPhi*sin2Phi 
              cos4Phi = cos2Phi*cos2Phi - sin2Phi*sin2Phi 
			
              Y44 =  3./16.*sqrt(35./(2.*pi))*cos4Phi*sinTheta**4

              u(i,j,k) = u0 * Y44
              
              rho(i,j,k) = rho0 * Y44
              
             
              dthetadx = z(i,j,k)*x(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetady = z(i,j,k)*y(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetadz = z(i,j,k)**2 &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2)) &
                   - 1/sqrt(x(i,j,k)**2+y(i,j,k)**2)
              
              vr = vr0 * Y44
              vtheta = - u0 * sqrt(5/(4*pi))*3*cosTheta*sinTheta
              
              vx(i,j,k) = x(i,j,k)/r1*vr + dthetadx*vtheta;
              vy(i,j,k) = y(i,j,k)/r1*vr + dthetady*vtheta;
              vz(i,j,k) = z(i,j,k)/r1*vr + dthetadz*vtheta; 
              
           end do
        end do
     end do

    else if (CCTK_EQUALS (initial_data, "multipole l=4, m=-4")) then
     
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              
              r1 = sqrt(r(i,j,k)**2 + eps**2)
              rpt = r1 + t
              rmt = r1 - t
              
              u0 = 0.0
              
              rho0 = half * (rpt/r1 * Gaussian(rpt) + rmt/r1 * Gaussian(rmt))
              
              vr0 = 0.0
              
              cosTheta = z(i,j,k)/r1
              sinTheta = sqrt(x(i,j,k)**2 + y(i,j,k)**2)/r1
              cosTheta2 = cosTheta**2
			  sinPhi = y(i,j,k) / sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k))
			  cosPhi = x(i,j,k) / sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k))
			  cos2Phi = cosPhi*cosPhi - sinPhi*sinPhi
			  sin2Phi = 2.*sinPhi*cosPhi
              cos3Phi = cosPhi*cos2Phi - sinPhi*sin2Phi 
              cos4Phi = cos2Phi*cos2Phi - sin2Phi*sin2Phi 
			
              Y44 =  3./16.*sqrt(35./(2.*pi))*sin4Phi*sinTheta**4

              u(i,j,k) = u0 * Y44
              
              rho(i,j,k) = rho0 * Y44
              
             
              dthetadx = z(i,j,k)*x(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetady = z(i,j,k)*y(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetadz = z(i,j,k)**2 &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2)) &
                   - 1/sqrt(x(i,j,k)**2+y(i,j,k)**2)
              
              vr = vr0 * Y44
              vtheta = - u0 * sqrt(5/(4*pi))*3*cosTheta*sinTheta
              
              vx(i,j,k) = x(i,j,k)/r1*vr + dthetadx*vtheta;
              vy(i,j,k) = y(i,j,k)/r1*vr + dthetady*vtheta;
              vz(i,j,k) = z(i,j,k)/r1*vr + dthetadz*vtheta; 
              
           end do
        end do
     end do




     else if (CCTK_EQUALS (initial_data, "multipole l=2, m=2")) then
     
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              
              r1 = sqrt(r(i,j,k)**2 + eps**2)
              rpt = r1 + t
              rmt = r1 - t
              rxy = sqrt(x(i,j,k)**2+y(i,j,k)**2)             
             
              if (rxy<1e-12) then
              	rxy = 1e-12
	      end if
 
              u0 = 0.0
              
              rho0 = half * (rpt/r1 * Gaussian(rpt) + rmt/r1 * Gaussian(rmt))
              
              vr0 = 0.0
              
              cosTheta = z(i,j,k)/r1
              sinTheta = rxy/r1
              cosTheta2 = cosTheta**2
	      sinPhi = y(i,j,k) / rxy
	      cosPhi = x(i,j,k) / rxy
             
              Y22 = 0.25*sqrt(15./(2.*pi))*sinTheta*sinTheta*(cosPhi*cosPhi-sinPhi*sinPhi)

              u(i,j,k) = u0 * Y22
              
              rho(i,j,k) = rho0 * Y22
              
             
              dthetadx = z(i,j,k)*x(i,j,k) &
                   / (r1**2*rxy)
              dthetady = z(i,j,k)*y(i,j,k) &
                   / (r1**2*rxy)
              dthetadz = z(i,j,k)**2 &
                   / (r1**2*rxy) &
                   - 1/rxy
              
              vr = vr0 * Y22
              vtheta = - u0 * sqrt(5/(4*pi))*3*cosTheta*sinTheta
              
              vx(i,j,k) = x(i,j,k)/r1*vr + dthetadx*vtheta;
              vy(i,j,k) = y(i,j,k)/r1*vr + dthetady*vtheta;
              vz(i,j,k) = z(i,j,k)/r1*vr + dthetadz*vtheta; 
              
           end do
        end do
     end do

     else if (CCTK_EQUALS (initial_data, "multipole l=2, m=-2")) then
     
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              
              r1 = sqrt(r(i,j,k)**2 + eps**2)
              rpt = r1 + t
              rmt = r1 - t
              
              u0 = 0.0
              
              rho0 = half * (rpt/r1 * Gaussian(rpt) + rmt/r1 * Gaussian(rmt))
              
              vr0 = 0.0
              
              cosTheta = z(i,j,k)/r1
              sinTheta = sqrt(x(i,j,k)**2 + y(i,j,k)**2)/r1
              cosTheta2 = cosTheta**2
			  sinPhi = y(i,j,k) / sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k))
			  cosPhi = x(i,j,k) / sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k))
             
              Y22 = 0.25*sqrt(15./(2.*pi))*sinTheta*sinTheta*(2.0*cosPhi*sinPhi)

              u(i,j,k) = u0 * Y22
              
              rho(i,j,k) = rho0 * Y22
              
             
              dthetadx = z(i,j,k)*x(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetady = z(i,j,k)*y(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetadz = z(i,j,k)**2 &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2)) &
                   - 1/sqrt(x(i,j,k)**2+y(i,j,k)**2)
              
              vr = vr0 * Y22
              vtheta = - u0 * sqrt(5/(4*pi))*3*cosTheta*sinTheta
              
              vx(i,j,k) = x(i,j,k)/r1*vr + dthetadx*vtheta;
              vy(i,j,k) = y(i,j,k)/r1*vr + dthetady*vtheta;
              vz(i,j,k) = z(i,j,k)/r1*vr + dthetadz*vtheta; 
              
           end do
        end do
     end do
 
   else if (CCTK_EQUALS (initial_data, "multipole l=2, m=1")) then
     
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              
              r1 = sqrt(r(i,j,k)**2 + eps**2)
              rpt = r1 + t
              rmt = r1 - t
              
              u0 = 0.0
              
              rho0 = half * (rpt/r1 * Gaussian(rpt) + rmt/r1 * Gaussian(rmt))
              
              vr0 = 0.0
              
              cosTheta = z(i,j,k)/r1
              sinTheta = sqrt(x(i,j,k)**2 + y(i,j,k)**2)/r1
              cosTheta2 = cosTheta**2
			  sinPhi = y(i,j,k) / sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k))
			  cosPhi = x(i,j,k) / sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k))
             
              Y21 = - 0.5*sqrt(15./(2.*pi))*cosTheta*sinTheta*cosPhi

              u(i,j,k) = u0 * Y21
              
              rho(i,j,k) = rho0 * Y21
              
             
              dthetadx = z(i,j,k)*x(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetady = z(i,j,k)*y(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetadz = z(i,j,k)**2 &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2)) &
                   - 1/sqrt(x(i,j,k)**2+y(i,j,k)**2)
              
              vr = vr0 * Y21
              vtheta = - u0 * sqrt(5/(4*pi))*3*cosTheta*sinTheta
              
              vx(i,j,k) = x(i,j,k)/r1*vr + dthetadx*vtheta;
              vy(i,j,k) = y(i,j,k)/r1*vr + dthetady*vtheta;
              vz(i,j,k) = z(i,j,k)/r1*vr + dthetadz*vtheta; 
              
           end do
        end do
     end do

     else if (CCTK_EQUALS (initial_data, "multipole l=2, m=-1")) then
     
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              
              r1 = sqrt(r(i,j,k)**2 + eps**2)
              rpt = r1 + t
              rmt = r1 - t
              
              u0 = 0.0
              
              rho0 = half * (rpt/r1 * Gaussian(rpt) + rmt/r1 * Gaussian(rmt))
              
              vr0 = 0.0
              
              cosTheta = z(i,j,k)/r1
              sinTheta = sqrt(x(i,j,k)**2 + y(i,j,k)**2)/r1
              cosTheta2 = cosTheta**2
			  sinPhi = y(i,j,k) / sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k))
			  cosPhi = x(i,j,k) / sqrt(x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k))
             
              Y21 = -0.5*sqrt(15./(2.*pi))*cosTheta*sinTheta*sinPhi

              u(i,j,k) = u0 * Y21
              
              rho(i,j,k) = rho0 * Y21
              
             
              dthetadx = z(i,j,k)*x(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetady = z(i,j,k)*y(i,j,k) &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2))
              dthetadz = z(i,j,k)**2 &
                   / (r1**2*sqrt(x(i,j,k)**2+y(i,j,k)**2)) &
                   - 1/sqrt(x(i,j,k)**2+y(i,j,k)**2)
              
              vr = vr0 * Y21
              vtheta = - u0 * sqrt(5/(4*pi))*3*cosTheta*sinTheta
              
              vx(i,j,k) = x(i,j,k)/r1*vr + dthetadx*vtheta;
              vy(i,j,k) = y(i,j,k)/r1*vr + dthetady*vtheta;
              vz(i,j,k) = z(i,j,k)/r1*vr + dthetadz*vtheta; 
              
           end do
        end do
     end do
    
  else if (CCTK_EQUALS (initial_data, "noise")) then
     
     do k=1,cctk_lsh(3)
        do j=1,cctk_lsh(2)
           do i=1,cctk_lsh(1)
              
              call random_number (u0)
              call random_number (rho0)
              call random_number (v0)
              
              u(i,j,k) = amplitude * (2*u0-1)
              rho(i,j,k) = amplitude * (2*rho0-1)
              vx(i,j,k) = amplitude * (2*v0(1)-1)
              vy(i,j,k) = amplitude * (2*v0(2)-1)
              vz(i,j,k) = amplitude * (2*v0(3)-1)
              
           end do
        end do
     end do
     
  else 
     
     call CCTK_WARN (0, "internal error")
     
  end if                        ! if initial_data
  
 
contains
  
  CCTK_REAL function Gaussian (rr)
    CCTK_REAL, intent(in) :: rr
    Gaussian = amplitude * exp(- ((rr - radius) / width)**2)
  end function Gaussian
  
  CCTK_REAL function DGaussian (rr)
    CCTK_REAL, intent(in) :: rr
    DGaussian = -2 * Gaussian(rr) * (rr - radius) / width**2
  end function DGaussian
  
end subroutine LWT_init



subroutine LWT_init_derivs (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  ! Calculate v_i numerically
  write(*,*) "WARNING: coputing vx,vy,vz numerically. This should NEVER happen"
  call Diff_gv (cctkGH, 0, u, vx, -1)
  call Diff_gv (cctkGH, 1, u, vy, -1)
  call Diff_gv (cctkGH, 2, u, vz, -1) 
  
end subroutine LWT_init_derivs
