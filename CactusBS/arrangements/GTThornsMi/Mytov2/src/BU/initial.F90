#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
  subroutine initial_tov(CCTK_ARGUMENTS)

! --->   INITIAL DATA   <---

! This routine calculates the initial data for
! the numerical evolution.

!  use arrays
!  use global_numbers
  use ode
  implicit none
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS
!>>>>>>>>>>>>>>>>>>>>>>>>>ARRAYS

  real(kind=8), allocatable, dimension (:) :: x_my
  real(kind=8), allocatable, dimension (:) :: a, alpha, m
  real(kind=8), allocatable, dimension (:) :: myvel,WW,p,rho0
  real(kind=8), allocatable, dimension (:) :: hh,epsilon,sound
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>GLOBAL NUMBERS

    real(kind=8) zero,third,half,one,two,pii
    real(kind=8)  t, dt,dx
    real(kind=8) diss_coef

    integer Nt,kx,ky,kz

! TOV stuff

    integer Star_radius_etiqueta,Star_radius_etiqueta_ini

    real(kind=8) Star_radius, Star_radius_ini



!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>AQUI VA LO QUE HABIA EN MAIN
  CCTK_REAL dt_temp
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


! --->   EVOLUTION   <---



!>>>>>>>>>>>>>>>>>>>>>>>>>>AQUI LO QUE IBA EN EVOLVE>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



! --->   ALLOCATE ARRAYS

  allocate(a(0:Nx))
  allocate(x_my(0:Nx))


  allocate(alpha(0:Nx),m(0:Nx))
  allocate(myvel(0:Nx),WW(0:Nx),p(0:Nx),rho0(0:Nx))
  allocate(hh(0:Nx),epsilon(0:Nx),sound(0:Nx))

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



!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



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
        print *, 'Reached the star surface at R=',rr
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

! ================================================
! =======     Conservative variables     =========
! ================================================

! Construction of intermediate and helpful variables

     myvel = zero		! for equilibrium initial configs
     WW = one / sqrt(one - a**2 * myvel**2)

! For the ideal gas

     epsilon = p / ( mygamma - one ) / rho0
     hh = one + epsilon + p / rho0
     sound = sqrt (   p * mygamma * ( mygamma - one) / &
             ( p * mygamma + rho0 * ( mygamma - one ) )   )



!>>>>>>>>>>>>>>>>>>>>>>DEFINE THE METRIC AND THE HYDRO VARIABLES  
  do k=1,cctk_lsh(3)
   do j=1,cctk_lsh(2)
    do i=1,cctk_lsh(1)
     gxx(i,j,k)= 1.0
     gyy(i,j,k)= 1.0
     gzz(i,j,k)= 1.0
     gxy(i,j,k)=0.0
     gxz(i,j,k)=0.0
     gyz(i,j,k) = 0.0
     kxx(i,j,k) = 0.0
     kyy(i,j,k) = 0.0
     kzz(i,j,k) = 0.0
     kxy(i,j,k) = 0.0
     kxz(i,j,k) = 0.0
     kyz(i,j,k) = 0.0
     alp(i,j,k) =0.0
     betax(i,j,k) = 0.0
     betay(i,j,k) = 0.0
     betaz(i,j,k) = 0.0
     end do
    end do
  end do
 rho =0.
 press=0.
 eps = 0.
 vel = 0.0
 w_lorentz =0.0
      print *,'gxx=', gxx(4,4,4),'kxx=',kxx(1,2,3),'rho=',rho(1,2,3),'press=',press(1,2,3),'eps=',eps(1,2,3)
      print *,'WW=',w_lorentz(1,2,3),'betaz=',betaz(1,2,3)

! --->   END   <---


! --->   END   <---

  print *,
  print *, 'PROGRAM TOV HAS FINISHED AND WORKS'

  end subroutine initial_tov

subroutine calcrhs(w, rr,rhs)

!  use arrays
!  use global_numbers
  implicit none
       DECLARE_CCTK_PARAMETERS
!       DECLARE_CCTK_ARGUMENTS


  real(kind=8), dimension(:),        intent(in)  :: w
  real(kind=8),                      intent(in)  :: rr
  real(kind=8), dimension(size(w)),  intent(out) :: rhs

  real(kind=8) rho0_here,mu_here,one,pii,two,dx
  one=1.0d0
  two=2.0d0
  pii = 4.0d0*atan(1.0d0)
    

!  mu = rho0 ( 1 + epsilon ) = rho0 + p/(mygamma-1)
! because 
! p = rho0 * epsilon( gamma - 1 ) => rho0 * epsilon = p / ( mygamma - 1 )

   if (w(2).gt.0.0d0) then
     rho0_here = ( w(2) / K_poly )**(1./mygamma)
   else
     rho0_here = 0.0d0
   end if
   mu_here = rho0_here + w(2) / ( mygamma - one)

! DEBUG
!     print *, 'r=',rr
!     print *, 'dx=',dx
!     print *, 'w(1)=',w(1)
!     print *, 'w(2)=',w(2)
!     print *, 'w(3)=',w(3)
!     print *, 'log(K_poly)=',log(K_poly)
!     print *, 'rho0_here=', rho0_here
!     print *, 'gamma=',mygamma
!     print *, 'K_poly=',K_poly

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

