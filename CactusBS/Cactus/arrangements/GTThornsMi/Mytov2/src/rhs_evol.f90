
  subroutine rhs_evol( )

! --->  Right hand sides   <---

! Calculates rhs for (phi,psi,pi).

  use arrays
  use global_numbers

  implicit none

  integer i,j

  real(kind=8) idxh,idx
  real(kind=8) grr_int,alpha_int,beta_int

  integer before_geo,after_geo

! ---> numbers

   idx = one/dx

! --->   Primitive vars and metric

   call primitive_vars( )

   call metric( )

! --->   Calculating the rhs

   if (num_method.eq.1) then

     flux1 = u(1,:) * vel / a
     flux2 = u(2,:) * vel / a 
     flux3 = ( u(3,:) / a + p ) * vel
     flux4 = p

     source1 = zero
     source2 = -alpha * a**3 * ( half * x * (one-one/a**2) ) * (u(2,:) * vel/a + u(3,:)/a + p + u(1,:)/a ) &
                /x**2
     source3 = - alpha * ( half * x * (one-one/a**2) ) * u(2,:) / x**2

     do i=1,Nx-1
       rhs_u(1,i) = - 3.0d0 * ( alpha(i+1) * a(i+1) * x(i+1)**2 * flux1(i+1) &
                             -  alpha(i-1) * a(i-1) * x(i-1)**2 * flux1(i-1)) &
                              / (x(i+1)**3 - x(i-1)**3) &
                    + source1(i) &
                    + diss_coef*(u(1,i+1) - two*u(1,i) + u(1,i-1))/dx**2

       rhs_u(2,i) = - 3.0d0 * ( alpha(i+1) * a(i+1) * x(i+1)**2 * flux2(i+1) &
                             -  alpha(i-1) * a(i-1) * x(i-1)**2 * flux2(i-1)) &
                              / (x(i+1)**3 - x(i-1)**3) &
                    - half * ( alpha(i+1) * a(i+1) * flux4(i+1) - alpha(i-1) * a(i-1) * flux4(i-1) )/dx &
                    + source2(i) &
                    + diss_coef*(u(2,i+1) - two*u(2,i) + u(2,i-1))/dx**2

       rhs_u(3,i) = - 3.0d0 * ( alpha(i+1) * a(i+1) * x(i+1)**2 * flux3(i+1) &
                             -  alpha(i-1) * a(i-1) * x(i-1)**2 * flux3(i-1)) &
                              / (x(i+1)**3 - x(i-1)**3) &
                    + source3(i) &
                    + diss_coef*(u(3,i+1) - two*u(3,i) + u(3,i-1))/dx**2
     end do

     rhs_u(:,0)  = 3.0D0*rhs_u(:,1) - 3.0D0*rhs_u(:,2) + rhs_u(:,3)
     rhs_u(:,Nx) = 3.0D0*rhs_u(:,Nx-1) - 3.0D0*rhs_u(:,Nx-2) + rhs_u(:,Nx-3)

   else if (num_method.eq.2) then

     call fluxes_hll( )

     source1 = zero
     source2 = -alpha * a**3 * ( half * x * (one-one/a**2) ) * (u(2,:) * vel/a + u(3,:)/a + p + u(1,:)/a ) &
                /x**2
     source3 = - alpha * ( half * x * (one-one/a**2) ) * u(2,:) / x**2

     do i=1,Nx-1
       rhs_u(1,i) = - 3.0d0*(alpha(i  )*a(i  )*x(i  )**2*flux1(i  ) &
                           - alpha(i-1)*a(i-1)*x(i-1)**2*flux1(i-1))/(x(i)**3 - x(i-1)**3) &
                    + source1(i)
       rhs_u(2,i) = - 3.0d0*(alpha(i  )*a(i  )*x(i  )**2*flux2(i  ) &
                           - alpha(i-1)*a(i-1)*x(i-1)**2*flux2(i-1))/(x(i)**3 - x(i-1)**3) &
                    - (       alpha(i  )*a(i  )*flux4(i  ) &
                            - alpha(i-1)*a(i-1)*flux4(i-1))/dx &
                    + source2(i)
       rhs_u(3,i) = - 3.0d0*(alpha(i  )*a(i  )*x(i  )**2*flux3(i  ) &
                           - alpha(i-1)*a(i-1)*x(i-1)**2*flux3(i-1))/(x(i)**3 - x(i-1)**3) &
                    + source3(i)

     end do

   end if

! ===========================================================
! =============     The metric function a     ===============
! ===========================================================

   rhs_a = - 4.0d0 * pii * x * alpha * u(2,:)

! ======================================================
! =============     GEODESICS     ======================
! ======================================================


   rhs_u(:,Nx) = rhs_u(:,Nx-1)
!   rhs_u(:,Nx) = 3.0D0*rhs_u(:,Nx-1) - 3.0D0*rhs_u(:,Nx-2) + rhs_u(:,Nx-3)

  before_geo = 0
  after_geo  = 0

   do j=1,num_geos
      do i=0,Nx
        if (xgeoout(j).lt.x(i)) then
          before_geo = i-1
          after_geo = i
          exit
        end if
      end do
      alpha_int = (alpha(after_geo)-alpha(before_geo))* &
                  (xgeoout(j)-x(before_geo))/(x(after_geo)-x(before_geo)) &
                  + alpha(before_geo)
      beta_int = zero
      grr_int = (a(after_geo)**2-a(before_geo)**2)* &
                (xgeoout(j)-x(before_geo))/(x(after_geo)-x(before_geo)) &
                + a(before_geo)**2

!      rhs_xgeoout(j) = ( - grr_int * beta_int + sqrt(grr_int) * alpha_int ) / &
!                       ( - alpha_int**2 + grr_int * beta_int**2 )

      rhs_xgeoout(j) = - beta_int + alpha_int / sqrt( grr_int )
   end do



! --->   END   <---

  end subroutine rhs_evol

