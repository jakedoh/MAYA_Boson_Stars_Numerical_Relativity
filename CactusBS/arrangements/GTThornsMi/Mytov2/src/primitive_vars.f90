
  subroutine primitive_vars( )

! --->  Right hand sides   <---


  use arrays
  use global_numbers
!  use derivatives

  implicit none

  integer i,j

  real(kind=8) idxh,idx
  real(kind=8) pressK,pressK_p,fun_of_p,funprime_of_p,tol,press_guess

! --->   Numbers

  idxh = half/dx
  idx  = 1.0D0/dx

     ! ============================================
     !   Reconstruct p
     ! ============================================

     if (EoS.eq.2) then     ! IDEAL GAS

! DEBUG
!     do i=1,Nx-1   OPTION 1
!     do i=1,Star_radius_etiqueta+50*2**(res_num-1)    !OPTION 2 works ok
     do i=1,Star_radius_etiqueta    ! running

!       press_guess = 0.5d0*(p(i) + p(i+1))
       press_guess = p(i)
       pressK = press_guess
       do j=1,1000000
         pressK_p = pressK

       fun_of_p = pressK_p-rho0(i)*(gamma-one)*(u(3,i)+ &
                 u(1,i)*(one-one/sqrt(one-u(2,i)**2/(u(3,i)+ &
                 a(i)*pressK_p+u(1,i))**2/a(i)**2 &
                 ))+a(i)*pressK_p*(one-one/(one-u(2,i)**2/ &
                 (u(3,i)+a(i)*pressK_p+u(1,i))**2/a(i)**2)))/ &
                 u(1,i)*sqrt(one-u(2,i)**2/(u(3,i)+a(i)* &
                 pressK_p+u(1,i))**2/a(i)**2)


       funprime_of_p = one-rho0(i)*(gamma-one)*(u(1,i)/ &
                 sqrt(one-u(2,i)**2/(u(3,i)+a(i)*pressK_p+ &
                 u(1,i))**2/a(i)**2)**3*u(2,i)**2/(u(3,i)+ &
                 a(i)*pressK_p+u(1,i))**3/a(i)+a(i)*(one-one/ &
                 (one-u(2,i)**2/(u(3,i)+a(i)*pressK_p+u(1,i))**2 &
                 /a(i)**2))+2*pressK_p/(one-u(2,i)**2/(u(3,i)+ &
                 a(i)*pressK_p+u(1,i))**2/a(i)**2)**2*u(2,i)**2/ &
                 (u(3,i)+a(i)*pressK_p+u(1,i))**3)/u(1,i)* &
                 sqrt(one-u(2,i)**2/(u(3,i)+a(i)*pressK_p+ &
                 u(1,i))**2/a(i)**2)-rho0(i)*(gamma-one)*(u(3,i)+ &
                 u(1,i)*(one-one/sqrt(one-u(2,i)**2/(u(3,i)+ &
                 a(i)*pressK_p+u(1,i))**2/a(i)**2))+a(i)*pressK_p &
                 *(one-one/(one-u(2,i)**2/(u(3,i)+a(i)*pressK_p+ &
                 u(1,i))**2/a(i)**2)))/u(1,i)/sqrt(one-u(2,i)**2/ &
                 (u(3,i)+a(i)*pressK_p+u(1,i))**2/a(i)**2)*u(2,i)**2 &
                 /(u(3,i)+a(i)*pressK_p+u(1,i))**3/a(i)

! --------------------------------------------------


         pressK = pressK_p - fun_of_p/funprime_of_p
         tol = 2.0d0*(pressK - pressK_p)/(pressK + pressK_p)

!         THE ORIGINAL tol is 1.e-6
         if (tol.lt.1.e-8) then
           p(i) = pressK
    !      print *, 'stops at iteration=',i
! DEBUG
! print *, 'p(',i,')=',p(i)
           exit
         else
         end if
         if (j.gt.500000) then
           print *, 'The Newton_Rapson routine did not converge after 10000 iterations'
           print *, 'At point x(i)=', x(i)
           stop
         end if

       end do

     end do

     p(0)  = 3.0D0*p(1) - 3.0D0*p(2) + p(3)
     p(Nx) = 3.0D0*p(Nx-1) - 3.0D0*p(Nx-2) + p(Nx-3)

     else if (EoS.eq.1) then
       p=zero
     end if

     ! ----->     Reconstruct vel,WW,rho0     <-----

     p = max(floor,p)

     vel = u(2,:) / a**3 / (u(3,:)/a + u(1,:)/a + p)
     WW  = one / sqrt(one - a**2 * vel**2)

     rho0 = max( floor, u(1,:) / a / WW )


!DEBUG
do i=Star_radius_etiqueta_ini,Star_radius_etiqueta   ! seems to work fine
  if (vel(i).lt.zero) vel(i)=floor
end do

! DEBUG
!     do i=Star_radius_etiqueta+51*2**(res_num-1),Nx-1	!works
     do i=Star_radius_etiqueta+1,Nx	!running
       rho0(i) = floor
       p(i) = floor
       epsilon(i) = floor ! -> new NOTHING happens
       vel(i) = floor    ! TESTING
     end do
!     do i=1,Nx
!       if (x(i).ge.Star_radius) vel(i)=floor 			!OPTION works
!     end do

     epsilon = p/(gamma - one)/rho0
     hh = one + epsilon + p / rho0
     sound = sqrt (   p * gamma * ( gamma - one) / &
             ( p * gamma + rho0 * ( gamma - one ) )   )

!     Some diagnostics

     do i=0,Nx
       if (rho0(i).lt.floor) then
          print *, 'Densidad menor que floor en x=', x(i)
          stop
       end if
       if ((p(i).le.0.0).and.(Eos.eq.2)) then
          print *, 'Presion no positiva en x=', x(i), 'p=',p(i)
!          stop
       end if
       if (p(i).lt.0.0) then
          print *, 'Presion negativa en x=', x(i), 'p=',p(i)
!          stop
       end if
     end do

     ! ============================================

! --->   END   <---

  end subroutine primitive_vars


