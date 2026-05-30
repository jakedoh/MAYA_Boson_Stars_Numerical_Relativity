#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
  subroutine initial_tov(CCTK_ARGUMENTS)

! --->   INITIAL DATA   <---

! This routine calculates the initial data for
! the numerical evolution.

!  use arrays
  use global_numbers
  use ode
  implicit none
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS
!>>>>>>>>>>>>>>>>>>>>>>>>>ARRAYS

  real(kind=8), allocatable, dimension (:) :: r_my
  real(kind=8), allocatable, dimension (:) :: a, alpha, m
  real(kind=8), allocatable, dimension (:) :: myvel,WW,p,rho0,phi,log_r
  real(kind=8), allocatable, dimension (:) :: hh,epsilon,sound

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>MAIN
  CCTK_REAL dt_temp, grr,xp,yp,zp,rp
  integer i,j,k,l
  logical rk3
  integer num_steps
  CCTK_REAL aux,aux0
  integer antes,despues
  CCTK_REAL idr,p_c

! ->  THE RK4 stuff
    integer, parameter            :: neq = 4 ! # of eqs.
    real(kind=8), dimension(neq)  :: w       ! solution vector
    real(kind=8), dimension(neq)  :: rhs     ! rhs vector
    integer                       :: step
    logical                       :: done = .false.
    character(len=42)             :: method

    CCTK_REAL  rr, drr, r_final
! <-
    CCTK_REAL Total_mass,Rest_mass,rho0_here,mu_here,M_bh,mysigma
    integer i_surface

! Some numbers

   zero  = 0.0D0
   third = 1.0D0/3.0D0
   half  = 0.5D0
   one   = 1.0D0
   two   = 2.0D0
   pii = 4.0d0*atan(1.0d0)
   mysigma = one
    print *, 'pi=',pii
    print *, 'rmin=', rmin
    print *, 'rmax=', rmax
    print *, 'Nr=', Nr
    print *, 'courant=', courant



! --->   ALLOCATE ARRAYS

  allocate(a(0:Nr))
  allocate(r_my(0:Nr))

  
  allocate(alpha(0:Nr),m(0:Nr))
  allocate(myvel(0:Nr),WW(0:Nr),p(0:Nr),rho0(0:Nr))
  allocate(hh(0:Nr),epsilon(0:Nr),sound(0:Nr),phi(0:Nr),log_r(0:Nr))

!===============================================

! -------------------------


  t = zero


! Set up the radial 1D grid and some messgaes to the screen

     dr = -(rmax - rmin)/dble(Nr)
     do i=0,Nr
       r_my(i)  =rmin + dble(i) * dr
     end do

     dt = courant * dr
     print *, '---------------------'
     print *, 'Numerical grid'
     print *, 'rmin=',rmin
     print *, 'rmax=',rmax
     print *, 'dr=',dr
     print *, 'dt=',dt
     print *, 'courant=',dt/dr
!stop

! Initializing arrays


     print *,'----------------------------'
     print *,'|  Time step  |    Time    |'
     print *,'----------------------------'

     write(*,"(A5,I6,A6,ES9.2,A3)") ' |   ',0,'    | ',t,'  |'

  print *,'-----------------------------'


! --->   END   <---



!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



! --->   NUMBERS   <---

  idr = 1.0D0/dr
  Total_mass = zero
  Rest_mass  = zero
  M_bh = 4.0d0*12.65d0/5.0d0
! Numerical integration of 
!	Mass function m(r)
!	Pressure P(r)
!	Grav potential Phi(r) ->  alpha(r)
! according to eqs (1.77-1.79) in Baumgarte's book p15.
  !  rr = (3.0d0*(12.65 - M_bh)*(8.0d0*pii*M_bh)**4/(4.0d0*pii*mysigma))**(1.0d0/3.0d0)
    rr = (3.0d0*(12.65 - M_bh)*(8.0d0*pii*M_bh)**4/(4.0d0*pii*mysigma))**(1.0d0/3.0d0)
  !  p_c = K_poly * rho0_c ** ( MyGamma )
p_c = (MyGamma - 1.0d0)*mysigma*(8.0d0*pii*M_bh)**(-MyGamma/(MyGamma - 1.0d0))
    w(1) = 12.65d0		! m(r)
    w(2) = p_c		! p(r)
    w(3) =  0.5*log(1.0d0 - 2.0d0*w(1)/rr)	! alpha(r)
    w(4) = 0.0  ! log(r_hat/r)
     step = 0

  !  r_final = r_my(0)
    method = 'rk4'

    print *, '+++++++++++++++++++'
    print *, 'rho0_c=',rho0_c
    print *, 'gamma=',MyGamma
    print *, 'M_bh=',M_bh
    print *, 'p_c=',p_c
    print *, 'r=',rr
    print *, 'r_final=',r_final
    print *, '+++++++++++++++++++'

!    open(1,file='TovBack.dat',status='replace')
!            write(1,*) r_my,m,p,rho0_here,phi
!    do while (.not.done)
    do i=0,Nr

!       print *, rr, w(1), w(2), w(3)
       m(i)     = w(1)
       p(i)     = w(2)
       phi(i)   = w(3)
       log_r(i) = w(4)
       call odestep(w,rr, -dr, method)


!      if (w(2).gt.0.0d0) then
        rho0_here =  w(2)/(MyGamma-1.0d0)
!      else
!        rho0_here = 0.0d0
!      end if
      mu_here = rho0_here

      Total_mass = m(i)
      Rest_mass  = Rest_mass + 4.0d0*pii*rr**2*rho0_here*dr/sqrt(one-two*w(1)/rr)

!      if (rho0_here.le.MyAtmosphere) then
!        Star_radius     = rr
!        Star_radius_ini = rr
!        print *, 'Reached the star surface at R=',rr
!        i_surface = i
!        Star_radius_etiqueta_ini = i
!        exit
!      end if
 !           write(1,*) r_my,m,p,rho0_here,phi
    enddo

  print *, '-----------------------------------'
  print *, 'Star radius=',Star_radius
  print *, 'Total mass=',Total_mass
  print *, 'Rest  mass=',Rest_mass
  print *, '-----------------------------------'

 !  close(1)

!stop


close(1)
! Define the metric function a
!  a = sqrt(one/(one-two*m/r_my))
!  a(0) = one

!  Match the solution with Schwarzschild's solution
!  do i=i_surface+1,Nr
!    a(i) = sqrt( one / ( one - two * Total_mass / r_my(i) ) )
!    phi(i) = half*log(one - two * Total_mass / r_my(i))
!    alpha(i) = exp(phi(i))
!  end do

! Rescale the lapse inside the star, such that alpha matches Schwarzschild
! alpha_new = alpha_old / alpha_old(surface) / a(surface)
!  do i=0,i_surface
!    phi(i) = phi(i)/phi(i_surface)
!    alpha(i) = exp(phi(i))
!  end do


 ! rho0 = ( p / K_poly )**(1./MyGamma)
  rho0 = p /(MyGamma-1.0d0)
 ! rho0 = max(MyAtmosphere,rho0)
  p = max(MyAtmosphere,p)
   alpha = exp(phi)
! ================================================
! ================================================

! Construction of intermediate and helpful variables

     myvel = zero		! for equilibrium initial configs
     WW = one / sqrt(one - a**2 * myvel**2)

! For the ideal gas

     epsilon = p / ( MyGamma - one ) / rho0
     hh = one + epsilon + p / rho0
     sound = sqrt (   p * MyGamma * ( MyGamma - one) / &
             ( p * MyGamma + rho0 * ( MyGamma - one ) )   )
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>HERE COMES THE INTERPOLATION>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.>>>
    do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
       do i = 1, cctk_lsh(1)
              
         antes=1
 
!        locate  points before and after r(i,j,k)
           do l=1,Nr-1
             if (r_my(l).le.r(i,j,k)) then
               antes = l
             end if
           end do	
        despues = antes+1


      print *, antes,despues
!      print *, r(antes,antes,antes),r(1,1,1),r(1,1,2)
!     if((r_my(antes).lt.r(i,j,k)))) then
          rho(i,j,k) = (rho0(despues)-rho0(antes))/dr*(r(i,j,k)-r_my(antes)) + rho0(antes)

!         print*,rho(1,1,1),rho0(antes),r(1,1,1), '= rho'
!stop
  
!         else
!            rho(i,j,k)=MyAtmosphere

!          print*,rho(i,j,k),'=rho'
!       end if
!                  rho(i,j,k)= 10*exp(- (x(i,j,k) )**2 / 0.1**2 &
!                             - (y(i,j,k))**2 /0.1**2&
!                              - (z(i,j,k))**2 / 0.1**2 )
                           
!                  print*,rho(i,j,k),cctk_lsh(1),cctk_lsh(2),cctk_lsh(3),'.......======'!

          alp(i,j,k) = (alpha(despues)-alpha(antes))/dr*(r(i,j,k)-r_my(antes)) + alpha(antes)
 !                 print*,alp(i,j,k)
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
                      (zp**2 + yp**2)/(rp**2)
                  gyy(i,j,k) = (yp/rp)**2 * grr + &
                      (xp**2 + zp**2)/(rp**2)
                  gzz(i,j,k) = (zp/rp)**2 * grr + &
                      (xp**2 + yp**2)/(rp**2)
     
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
                   press(i,j,k)= K_poly*rho(i,j,k)**(MyGamma)
                   eps(i,j,k)= press(i,j,k) / ( MyGamma - one ) / rho(i,j,k)
!                   velx(i,j,k)=0.0
!         	   vely(i,j,k)=0.0
!                   velz(i,j,k)=0.0
                   w_lorentz(i,j,k) = 1/sqrt(1.0-(&
                                  gxx(i,j,k) * vel(0,i,j,k) * vel(0,i,j,k)+&
                                  gyy(i,j,k) * vel(1,i,j,k) * vel(1,i,j,k)+&
                                  gzz(i,j,k) * vel(2,i,j,k) * vel(2,i,j,k)+&
                                2*gxy(i,j,k) * vel(0,i,j,k) * vel(1,i,j,k)+&
                                2*gxz(i,j,k) * vel(0,i,j,k) * vel(2,i,j,k)+&
                                2*gyz(i,j,k) * vel(1,i,j,k) * vel(2,i,j,k)))  
!                    vel(i,j,k)=zero
!                    w_lorentz(i,j,k) = one
            enddo
          enddo
        enddo


!LUEGO DEFINES el resto de las variables hydrodinamicas
! --->   END   <---
  deallocate(a)
  deallocate(r_my)
  deallocate(alpha,m,phi)
  deallocate(myvel,WW,p,rho0)
  deallocate(hh,epsilon,sound)
!  print*,rho(1,1,1),rho(2,2,2),rho(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)),'>>>>>>>>>>>>>>>>>>CHECKING THE INTERPOLATION'
!stop

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
    

!!  mu = rho0 ( 1 + epsilon ) = rho0 + p/(MyGamma-1)
!! because 
!! p = rho0 * epsilon( gamma - 1 ) => rho0 * epsilon = p / ( mygamma - 1 )

!   if (w(2).gt.0.0d0) then
     rho0_here =  w(2) / (MyGamma-1.0d0)
!      print *, K_poly, "WE ARE HERE"
   

   mu_here = rho0_here 

     rhs(1) = 4.0d0 * pii * rr**2 * mu_here
!     if (rr.lt.dr) then 
!       rhs(2) = - ( mu_here + w(2) ) * 4.0d0*pii* rr * (mu_here/3.0d0 + w(2) ) / &
 !                        (one - two*4.0d0*pii * mu_here * rr**2 / 3.0d0)
 !      rhs(3) = w(3) * ( 4.0d0 * pii * mu_here * rr / 3.0d0  + 4.0d0 * pii * rr * w(2) ) / &
 !               ( one - 8.0d0 * pii * mu_here * rr**2 / 3.0d0 )
!     else
       rhs(2) = - ( mu_here + w(2) ) * ( w(1)  + 4.0d0 * pii * rr**3 * w(2) ) / &
                         ( rr**2 - 2.0d0 * rr * w(1) ) 
!       rhs(3) = w(3) * (w(1) + 4.0d0 * pii * rr**3 * w(2) ) / ( rr**2 - two * w(1) * rr )
        

       rhs(3) = w(1) + 4.0D0*pii* rr**3 * w(2)/( rr**2 - 2.0d0 * rr * w(1) )
       
       rhs(4) = (sqrt(rr) - sqrt(rr - 2.0D0*w(1)))/(rr*(sqrt(rr - 2.0D0 * w(1) ) ) )


!     end if


end subroutine calcrhs

