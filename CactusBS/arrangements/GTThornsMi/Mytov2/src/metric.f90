
  subroutine metric( )

! --->   INITIAL DATA   <---

! This routine calculates the initial data for
! the numerical evolution.

  use arrays
  use derivatives
  use global_numbers
  use ode_fly

  implicit none

  integer i,j

  real(kind=8) idx,aux,p_c

! ->  THE RK4 stuff
    integer, parameter            :: neq = 1 ! # of eqs.
    real(kind=8), dimension(neq)  :: w       ! solution vector
    real(kind=8), dimension(neq)  :: rhs     ! rhs vector
    integer                       :: step
    logical                       :: done = .false.
    character(len=42)             :: method

    real(kind=8) :: r, dr, r_final
! <-

! ----->     Convention of grid functions     <-----

! w(1): alpha

! --->   NUMBERS   <---

    w(1) = alpha(0)

    step = 0
    r    = x(0)
    r_final = x(Nx)
    method = 'rk4'

    do i=0,Nx

       alpha(i) = w(1)

       call odestep_fly(w, r, dx, method)

    enddo

! Rescale the lapse such that a=1/alpha at the boundary

    alpha = alpha/(a(Nx)*alpha(Nx))

! --->   END   <---

  end subroutine metric


! ===========================================================================
! ===========================================================================
! ===========================     Ode rhss    ===============================
! ===========================================================================
! ===========================================================================

subroutine calcrhs_fly(w, r, rhs)

  use arrays
  use global_numbers

  implicit none

  real(kind=8), dimension(:),        intent(in)  :: w
  real(kind=8),                      intent(in)  :: r
  real(kind=8), dimension(size(w)),  intent(out) :: rhs
  real(kind=8) rho0_here,mu_here,u2_here,vel_here,a_here

  integer etiqueta

  etiqueta = r/dx

!  print *, 'r=',r,'etiqueta=',etiqueta,'x(etiqueta)=',x(etiqueta), 'r/dx=',r/dx


      rho0_here = ( p(etiqueta) / K_poly )**(1./gamma)
      mu_here = rho0_here + p(etiqueta) / ( gamma - one)
      u2_here = u(2,etiqueta)
      vel_here = vel(etiqueta)
      a_here = a(etiqueta)


!     if (r.lt.dx) then
!       rhs(1) = w(1) * ( 4.0d0 * pii * mu_here * r / 3.0d0  &
!                + 4.0d0 * pii * r * p(etiqueta) ) / &
!                ( one - 8.0d0 * pii * mu_here * r**2 / 3.0d0 )
!     else
!       rhs(1) = w(1) * ( half*(one-one/a(etiqueta)**2)*r + 4.0d0 * &
!                pii * r**3 * p(etiqueta) ) &
!                / ( r**2 - two * (half*(one-one/a(etiqueta)**2)*r) * r )
!     end if

     if (r.lt.dx) then
       rhs(1) = w(1) * ( 4.0d0 * pii * mu_here * r / 3.0d0  &
                + 4.0d0 * pii * r * ( p(etiqueta) + u2_here * vel_here /a_here) ) / &
                ( one - 8.0d0 * pii * mu_here * r**2 / 3.0d0 )
     else
       rhs(1) = w(1) * ( half*(one-one/a(etiqueta)**2)*r + 4.0d0 * &
                pii * r**3 * ( p(etiqueta) + u2_here * vel_here / a_here ) ) &
                / ( r**2 - two * (half*(one-one/a(etiqueta)**2)*r) * r )
     end if

end subroutine calcrhs_fly

