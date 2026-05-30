#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
  subroutine initial_tov(CCTK_ARGUMENTS)

! --->   INITIAL DATA   <---

! This routine calculates the initial data for
! the numerical evolution.

  use global_numbers
  use ode
  implicit none
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS
!>>>>>>>>>>>>>>>>>>>>>>>>>Defining ARRAYS and Numbers>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  real(kind=8), allocatable, dimension (:) :: x_my
  real(kind=8), allocatable, dimension (:) :: a, alpha, m
  real(kind=8), allocatable, dimension (:) :: p,rho0


  CCTK_REAL dt_temp, grr,xp,yp,zp,rp
  integer i,j,k,l
  logical rk3
  integer num_steps
  CCTK_REAL aux,aux0
  integer antes,despues
  CCTK_REAL idx,p_c

! ->  THE RK4 stuff
    integer, parameter            :: neq = 3 ! # of eqs.
    real(kind=8), dimension(neq)  :: w       ! solution vector
    real(kind=8), dimension(neq)  :: rhs     ! rhs vector
    integer                       :: step
    logical                       :: done = .false.
    character(len=42)             :: method

    CCTK_REAL  rr, drr, r_final
! <-
    CCTK_REAL Total_mass,Rest_mass,rho0_here,mu_here
    integer i_surface

! Some numbers

   zero  = 0.0D0
   third = 1.0D0/3.0D0
   half  = 0.5D0
   one   = 1.0D0
   two   = 2.0D0
   pii = 4.0d0*atan(1.0d0)

    print *, 'pi=',pii
    print *, 'xmin=', myxmin
    print *, 'xmax=', myxmax
    print *, 'Nx=', Nx
    print *, 'courant=', courant
!>>>>>>>>>>>>>>>>>>>>>>>>>end defining arrays>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! --->   ALLOCATE ARRAYS

  allocate(a(0:Nx))
  allocate(x_my(0:Nx))
  allocate(alpha(0:Nx),m(0:Nx))
  allocate(p(0:Nx),rho0(0:Nx))


!===============================================

! -------------------------
!  Convention
!  u(1,:) a * D
!  u(2,:) a * S_r
!  u(3,:) a * tau
! -------------------------

! --->   INITIALIZE TIME

  t = zero


! Set up the radial 1D grid and some messgaes to the screen

     dx = (myxmax - myxmin)/dble(Nx)
     do i=0,Nx
       x_my(i)  =myxmin + dble(i) * dx
     end do

     dt = courant * dx
     print *, '---------------------'
     print *, 'Numerical grid'
     print *, 'xmin=',myxmin
     print *, 'xmax=',myxmax
     print *, 'dx=',dx
     print *, 'dt=',dt
     print *, 'courant=',dt/dx


! Initializing arrays


     print *,'----------------------------'
     print *,'|  Time step  |    Time    |'
     print *,'----------------------------'

     write(*,"(A5,I6,A6,ES9.2,A3)") ' |   ',0,'    | ',t,'  |'

  print *,'-----------------------------'


! --->   END   <---
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Finish initial>>>>>>>>>>>>>>>
! --->   NUMBERS   <---

  idx = 1.0D0/dx
  Total_mass = zero
  Rest_mass  = zero

! Numerical integration of 
!	Mass function m(r)
!	Pressure P(r)
!	Grav potential Phi(r) ->  alpha(r)
! according to eqs (1.77-1.79) in Baumgarte's book p15.

    p_c = K_poly * rho0_c ** ( mygamma )

    w(1) = zero		! m(r)
    w(2) = p_c		! p(r)
    w(3) = alpha_c	! alpha(r)

    step = 0
    rr    = x_my(0)
    r_final = x_my(Nx)
    method = 'rk4'

    print *, '+++++++++++++++++++'
    print *, 'rho0_c=',rho0_c
    print *, 'gamma=',mygamma
    print *, 'alpha_c=',alpha_c
    print *, 'p_c=',p_c
    print *, 'r=',rr
    print *, 'r_final=',r_final
    print *, '+++++++++++++++++++'


!    do while (.not.done)
    do i=0,Nx

!       print *, rr, w(1), w(2), w(3)
       m(i)     = w(1)
       p(i)     = w(2)
       alpha(i) = w(3)
       call odestep(w,rr, dx, method)


      if (w(2).gt.0.0d0) then
        rho0_here = ( w(2) / K_poly )**(1./mygamma)
      else
        rho0_here = 0.0d0
      end if
      mu_here = rho0_here + w(2) / ( mygamma - one)

      Total_mass = m(i)
      Rest_mass  = Rest_mass + 4.0d0*pii*rr**2*rho0_here*dx/sqrt(one-two*w(1)/rr)

      if (rho0_here.le.myfloor) then
        Star_radius     = rr
        Star_radius_ini = rr
!        print *, 'Reached the star surface at R=',rr
        i_surface = i
        Star_radius_etiqueta_ini = i
!        exit
      end if

    enddo

  print *, '-----------------------------------'
  print *, 'Star radius=',Star_radius
  print *, 'Total mass=',Total_mass
  print *, 'Rest  mass=',Rest_mass
  print *, '-----------------------------------'

! Define the metric function a
  a = sqrt(one/(one-two*m/x_my))
  a(0) = one

!  Match the solution with Schwarzschild's solution
  do i=i_surface+1,Nx
    a(i) = sqrt( one / ( one - two * Total_mass / x_my(i) ) )
    alpha(i) = one / a(i)
  end do

! Rescale the lapse inside the star, such that alpha matches Schwarzschild
! alpha_new = alpha_old / alpha_old(surface) / a(surface)
  do i=0,i_surface
    alpha(i) = alpha(i)/alpha(i_surface)/a(i_surface)
  end do


  rho0 = ( p / K_poly )**(1./mygamma)
  rho0 = max(myfloor,rho0)
  p = max(myfloor,p)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!======================================HERE COMES THE INTERPOLATION=========================
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>        
!=====Interpolation of the metric and curvature components from 1D to the 3D cactus grid=====
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!The convention from ADMBase is:             
! gij(:,:,:)= metric componentes (with ij=xx,yy,zz,xy,xz,yz)
! alp(:,:,:)= the lapse                           
! betai(:,:,:)= the shift components (with i =x,y,z)                   
! kij(:,:,:)= x extrinsic curvature components (with ij=xx,yy,zz,xy,xz,yz)    
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



   do k = 1, cctk_lsh(3)
   do j = 1, cctk_lsh(2)
   do i = 1, cctk_lsh(1)


        
        antes=0
!        locate  points before and after r(i,j,k)
           do l=0,Nx-1
             if (x_my(l).le.r(i,j,k).and.x(i,j,k).ge.0.0) then
               antes = l
             end if
           end do
        despues = antes+1

               alp(i,j,k) = (alpha(despues)-alpha(antes))/dx&
                          *(r(i,j,k)-x_my(antes)) + alpha(antes)
                       
               grr=one/(alp(i,j,k)**2)
               xp = x(i,j,k)
               yp = y(i,j,k)
               zp = z(i,j,k)
               rp = r(i,j,k)

               if (rp.eq.0) then

                  gxx(i,j,k) = 1.
                  gyy(i,j,k) = 1.
                  gzz(i,j,k) = 1.
                  gxy(i,j,k) = 0.
                  gxz(i,j,k) = 0.
                  gyz(i,j,k) = 0.

               else

                  gxx(i,j,k) = (xp/rp)**2 * grr + &
                      (zp**2 + yp**2)/(rp**2) + 1.0e-30
                  gyy(i,j,k) = (yp/rp)**2 * grr + &
                      (xp**2 + zp**2)/(rp**2) + 1.0e-30
                  gzz(i,j,k) = (zp/rp)**2 * grr + &
                      (xp**2 + yp**2)/(rp**2) + 1.0e-30
     
             if ((xp.eq.0) .and. (yp.eq.0)) then

                     gxy(i,j,k) = 0.
                     gyz(i,j,k) = 0.
                     gxz(i,j,k) = 0.

                  else

                     gxy(i,j,k) = (xp * yp/(rp**2) * grr&
                    + (xp*yp*zp**2)/(rp**2*(xp**2+yp**2))&
                    - (yp*xp)/(xp**2+yp**2))

                     gxz(i,j,k) = (xp * zp/(rp**2)&
                    * grr - (xp * zp/(rp**2)))

                     gyz(i,j,k) = (yp * zp/(rp**2)&
                    * grr - (yp * zp/(rp**2)))

                  endif

               endif
                   kxx(i,j,k)=zero
                   kyy(i,j,k)=zero
                   kzz(i,j,k)=zero
                   kxy(i,j,k)=zero
                   kxz(i,j,k)=zero
                   kyz(i,j,k)=zero
                   betax(i,j,k)=zero
                   betay(i,j,k)=zero
                   betaz(i,j,k)=zero
        enddo
        enddo
        enddo

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>        
!>>>>>>>>>>>>>>>>>>>>End defining metric and curvature>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>        
!>>>>>>>>>>>>>>>>>>>>Interpolation of the hydro variables from 1D to the 3D cactus grid>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>HydroBase press:=pressure, rho=density vel(i,j,k,1)=velx>
!The convention from Hydrobase is:             |
! press(:,:,:)=pressure                        |
! rho(:,:,:)=density                           |
! eps(:,:,:)=internal energy                   |
! vel(:,:,:,1)= x component of the velocity    |
! vel(:,:,:,2)= y component of the velocity    |
! vel(:,:,:,3)= z component of the velocity    | 
! w_lorentz(:,:,:)=lorentz factor              |
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


   do k = 1, cctk_lsh(3)
   do j = 1, cctk_lsh(2)
   do i = 1, cctk_lsh(1)
        antes=0
           do l=0,Nx-1
             if (x_my(l).le.r(i,j,k)) then
               antes = l
             end if
           end do
        despues = antes+1
 
          press(i,j,k) = (p(despues)-p(antes))/dx*(r(i,j,k)-x_my(antes)) + p(antes)
          rho(i,j,k) =  (press(i,j,k) / K_poly )**(1.0/mygamma)

          eps(i,j,k)= press(i,j,k) / ( mygamma - one ) / rho(i,j,k)
          vel(i,j,k,1)= zero
          vel(i,j,k,2)= zero
          vel(i,j,k,3)= zero

          w_lorentz(i,j,k) = 1/sqrt(1.0-(&
                           + gxx(i,j,k) * vel(i,j,k,1) * vel(i,j,k,1)&
                           + gyy(i,j,k) * vel(i,j,k,2) * vel(i,j,k,2)&
                           + gzz(i,j,k) * vel(i,j,k,3) * vel(i,j,k,3)&
                           + 2*gxy(i,j,k) * vel(i,j,k,1) * vel(i,j,k,2)&
                           + 2*gxz(i,j,k) * vel(i,j,k,1) * vel(i,j,k,3)&
                           + 2*gyz(i,j,k) * vel(i,j,k,2) * vel(i,j,k,3)))  



   end do
   end do
   end do

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>        
!>>>>>>>>>>>>>>>>>>>>End defining hydro variables>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>        
!>>>>>>>>>>>>>>>>>>>>Free the memory from the 1D arrays>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! --->   END   <---
  deallocate(a)
  deallocate(x_my)
  deallocate(alpha,m)
  deallocate(p,rho0)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  print *,
  print *, 'PROGRAM TOV HAS FINISHED AND WORKS'

  end subroutine initial_tov

subroutine calcrhs(w, rr,rhs)

!  use arrays
  use global_numbers
  implicit none
       DECLARE_CCTK_PARAMETERS
!       DECLARE_CCTK_ARGUMENTS


  real(kind=8), dimension(:),        intent(in)  :: w
  real(kind=8),                      intent(in)  :: rr
  real(kind=8), dimension(size(w)),  intent(out) :: rhs

  real(kind=8) rho0_here,mu_here
    

!!  mu = rho0 ( 1 + epsilon ) = rho0 + p/(mygamma-1)
!! because 
!! p = rho0 * epsilon( gamma - 1 ) => rho0 * epsilon = p / ( mygamma - 1 )

   if (w(2).gt.0.0d0) then
     rho0_here = ( w(2) / K_poly )**(1./mygamma)
   else
     rho0_here = 0.0d0
   end if
   mu_here = rho0_here + w(2) / ( mygamma - one)

     rhs(1) = 4.0d0 * pii * rr**2 * mu_here
     if (rr.lt.dx) then 
       rhs(2) = - ( mu_here + w(2) ) * 4.0d0*pii* rr * (mu_here/3.0d0 + w(2) ) / &
                         (one - two*4.0d0*pii * mu_here * rr**2 / 3.0d0)
       rhs(3) = w(3) * ( 4.0d0 * pii * mu_here * rr / 3.0d0  + 4.0d0 * pii * rr * w(2) ) / &
                ( one - 8.0d0 * pii * mu_here * rr**2 / 3.0d0 )
     else
       rhs(2) = - ( mu_here + w(2) ) * ( w(1)  + 4.0d0 * pii * rr**3 * w(2) ) / &
                         ( rr**2 - two * rr * w(1) ) 
       rhs(3) = w(3) * (w(1) + 4.0d0 * pii * rr**3 * w(2) ) / ( rr**2 - two * w(1) * rr )
     end if


end subroutine calcrhs

