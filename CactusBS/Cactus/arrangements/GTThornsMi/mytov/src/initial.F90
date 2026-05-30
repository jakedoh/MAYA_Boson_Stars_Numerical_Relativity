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

  real(kind=8), allocatable, dimension (:) :: x_my,phi,xx,pp,rho00
  real(kind=8), allocatable, dimension (:) :: a, alpha, m
  real(kind=8), allocatable, dimension (:) :: myvel,WW,p,rho0
  real(kind=8), allocatable, dimension (:) :: hh,epsilon,sound
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>GLOBAL NUMBERS

!    real(kind=8) zero,third,half,one,two,pii,dr
!    real(kind=8)  t, dt,dx
!    real(kind=8) diss_coef

!    integer Nt,kx,ky,kz

! TOV stuff

!    integer Star_radius_etiqueta,Star_radius_etiqueta_ini

!    real(kind=8) Star_radius, Star_radius_ini



!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>AQUI VA LO QUE HABIA EN MAIN
  CCTK_REAL dt_temp, grr,xp,yp,zp,rp,dr
  integer i,j,k,l,ii
  logical rk3
  integer num_steps
  CCTK_REAL aux,aux0,x1,x2
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
  allocate(hh(0:Nx),epsilon(0:Nx),sound(0:Nx),phi(0:Num_datos),xx(0:Num_datos),pp(0:Num_datos),rho00(0:Num_datos))

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



  if(CCTK_EQUALS(initial_profile,"ZurekPage")) then

     
     
        open(unit = 78,file='Zurek.txt',status='old')
    print *, 'Leyendo datos de la función'
    print *
    do ii=1,Num_datos
    read(78,*) xx(ii),rho0(ii),p(ii),m(ii),phi(ii)  
        alpha(ii) = exp(phi(ii))  
        rho00(ii) = ( (1.0d0/3.0d0)*rho0(ii) )**(3.0d0/4.0d0)    
    end do
    print *, 'dr obtenido del archivo =',xx(2)-xx(1)
    dx = xx(2)-xx(1)
    print *, 'Lo ha leido!! ',xx(1),rho0(1),pp(1),m(1),phi(1),alpha(1)  

    ! stop



!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>HERE COMES THE INTERPOLATION>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.>>>
    do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
       do i = 1, cctk_lsh(1)

!         dx=xx(2)-xx(1)
         antes=1
             
!        locate  points before and after r(i,j,k)
           do l=1,Num_datos
             if (xx(l).le.r(i,j,k)) then
               antes = l
             end if
           end do	
        despues = antes+1

         x1 = xx(antes)
         x2 = xx(despues+1)

    !      print*, "dx=",dx


   !   print *, antes,despues
!      print *, r(antes,antes,antes),r(1,1,1),r(1,1,2)
!     if((x_my(antes).lt.r(i,j,k)))) then
!          rho(i,j,k) = (rho0(despues)-rho0(antes))/dx*(r(i,j,k)-xx(antes)) + rho0(antes)
          rho(i,j,k) = (rho00(despues)-rho00(antes))/dx*(r(i,j,k)-xx(antes)) + rho00(antes)
          press(i,j,k)= mygamma*rho(i,j,k)
 !         rho(i,j,k) = 1e-8
!          press(i,j,k) = (pp(despues)-pp(antes))/dx*(r(i,j,k)-xx(antes)) + pp(antes)
 !         press(i,j,k) = 1e-8
!         print*,rho(1,1,1),rho0(antes),r(1,1,1), '= rho'
!stop
  
!         else
!            rho(i,j,k)=myfloor

!          print*,rho(i,j,k),'=rho'
!       end if
!                  rho(i,j,k)= 10*exp(- (x(i,j,k) )**2 / 0.1**2 &
!                             - (y(i,j,k))**2 /0.1**2&
!                              - (z(i,j,k))**2 / 0.1**2 )
                           
!                  print*,rho(i,j,k),cctk_lsh(1),cctk_lsh(2),cctk_lsh(3),'.......======'!

          alp(i,j,k) = (alpha(despues)-alpha(antes))/dx*(r(i,j,k)-xx(antes)) + alpha(antes)
 !                 print*,alp(i,j,k)
               grr=one/(alp(i,j,k)**2)
               xp = x(i,j,k)
               yp = y(i,j,k)
               zp = z(i,j,k)
               rp = r(i,j,k)

               if (rp.eq.0) then

                  gxx(i,j,k) = 10.
                  gyy(i,j,k) = 10.
                  gzz(i,j,k) = 10.
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
!                   press(i,j,k)= mygamma*rho(i,j,k)
                   
                   eps(i,j,k)= press(i,j,k) / ( mygamma - one ) / rho(i,j,k)
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

        close(unit=78)

            deallocate(a,phi,pp,xx)
            deallocate(x_my)
            deallocate(alpha,m)
            deallocate(myvel,WW,p,rho0,rho00)
            deallocate(hh,epsilon,sound)

  else if (CCTK_EQUALS(initial_profile,"MyTov")) then

             print*, 'Are we here?? REALLY!!! that is not correct'

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

!stop

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
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>HERE COMES THE INTERPOLATION>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.>>>
    do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
       do i = 1, cctk_lsh(1)
              
         antes=1
 
!        locate  points before and after r(i,j,k)
           do l=1,Nx-1
             if (x_my(l).le.r(i,j,k)) then
               antes = l
             end if
           end do	
        despues = antes+1


      print *, antes,despues
!      print *, r(antes,antes,antes),r(1,1,1),r(1,1,2)
!     if((x_my(antes).lt.r(i,j,k)))) then
          rho(i,j,k) = (rho0(despues)-rho0(antes))/dx*(r(i,j,k)-x_my(antes)) + rho0(antes)

!         print*,rho(1,1,1),rho0(antes),r(1,1,1), '= rho'
!stop
  
!         else
!            rho(i,j,k)=myfloor

!          print*,rho(i,j,k),'=rho'
!       end if
!                  rho(i,j,k)= 10*exp(- (x(i,j,k) )**2 / 0.1**2 &
!                             - (y(i,j,k))**2 /0.1**2&
!                              - (z(i,j,k))**2 / 0.1**2 )
                           
!                  print*,rho(i,j,k),cctk_lsh(1),cctk_lsh(2),cctk_lsh(3),'.......======'!

          alp(i,j,k) = (alpha(despues)-alpha(antes))/dx*(r(i,j,k)-x_my(antes)) + alpha(antes)
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
!                   press(i,j,k)= K_poly*rho(i,j,k)**(mygamma)
!                   eps(i,j,k)= press(i,j,k) / ( mygamma - one ) / rho(i,j,k)
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
  deallocate(a,phi,pp,xx)
  deallocate(x_my)
  deallocate(alpha,m)
  deallocate(myvel,WW,p,rho0,rho00)
  deallocate(hh,epsilon,sound)
!  print*,rho(1,1,1),rho(2,2,2),rho(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)),'>>>>>>>>>>>>>>>>>>CHECKING THE INTERPOLATION'
!stop
  end if
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

! DEBUG
!     print *, 'r=',rr
!     print *, 'dx=',dx
!     print *, 'w(1)=',w(1)
!     print *, 'w(2)=',w(2)
!1     print *, 'w(3)=',w(3)
!1     print *, 'log(K_poly)=',log(K_poly)
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

