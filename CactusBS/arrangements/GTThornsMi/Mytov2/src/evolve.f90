
  subroutine evolve( )

! --->   NUMERICAL EVOLUTION   <---

  use arrays
  use global_numbers

  implicit none

  logical rk3

  integer i,j,k,l
  integer num_steps

  real(kind=8) dt_temp
  real(kind=8) aux,aux0

! --->   ALLOCATE ARRAYS

  allocate(x(0:Nx))
  allocate(u(1:nvars,-g_cells:Nx+g_cells),u_p(1:nvars,-g_cells:Nx+g_cells))
  allocate(rhs_u(1:nvars,-g_cells:Nx+g_cells))
  allocate(a(0:Nx),da(0:Nx),a_p(0:Nx),rhs_a(0:Nx))

  allocate(alpha(0:Nx),m(0:Nx))
  allocate(flux1(0:Nx),flux2(0:Nx),flux3(0:Nx),flux4(0:Nx))
  allocate(source1(0:Nx),source2(0:Nx),source3(0:Nx))

  allocate(vel(0:Nx),WW(0:Nx),p(0:Nx),rho0(0:Nx))
  allocate(hh(0:Nx),epsilon(0:Nx),sound(0:Nx))

  allocate(Ham(0:Nx))
!  allocate(detgamma_sqrt(0:Nx),Ham(0:Nx))

  allocate(xgeoout(1:num_geos),xgeoout_p(1:num_geos),rhs_xgeoout(1:num_geos),xgeoout_i(1:num_geos))

! -------------------------
!  Convention
!  u(1,:) a * D
!  u(2,:) a * S_r
!  u(3,:) a * tau
! -------------------------

! --->   INITIALIZE TIME

  t = zero


! Set up the grid and some messgaes to the screen

     dx = (xmax - xmin)/dble(Nx)
     do i=0,Nx
       x(i)  = xmin + dble(i) * dx
     end do
     dt = courant * dx

     print *, '---------------------'
     print *, 'Numerical grid'
     print *, 'xmin=',xmin
     print *, 'xmax=',xmax
     print *, 'dx=',dx
     print *, 'dt=',dt
     print *, 'courant=',dt/dx
     print *, '---------------------'

! Initializing arrays

     u_p = zero
     a_p = zero

     print *,'----------------------------'
     print *,'|  Time step  |    Time    |'
     print *,'----------------------------'
     
     write(*,"(A5,I6,A6,ES9.2,A3)") ' |   ',0,'    | ',t,'  |'

     call initial()
     call diagnostics()

! --->   SAVE THE INITIAL DATA


     call save0Ddata(Nx,res_num,t,dx,a,'a',0)
     call save0Ddata(Nx,res_num,t,dx,m,'m',0)
     call save0Ddata(Nx,res_num,t,dx,alpha,'alpha',0)
     call save0Ddata(Nx,res_num,t,dx,rho0,'rho0',0)
     call save0Ddata(Nx,res_num,t,dx,Ham,'Ham',0)

     call save1Ddata(Nx,res_num,t,x,a,'a',0)
     call save1Ddata(Nx,res_num,t,x,alpha,'alpha',0)
     call save1Ddata(Nx,res_num,t,x,m,'m',0)
     call save1Ddata(Nx,res_num,t,x,vel,'vel',0)
     call save1Ddata(Nx,res_num,t,x,p,'p',0)
     call save1Ddata(Nx,res_num,t,x,Ham,'Ham',0)

     call save1Ddata(Nx,res_num,t,x,u(1,:),'u1',0)
     call save1Ddata(Nx,res_num,t,x,u(2,:),'u2',0)
     call save1Ddata(Nx,res_num,t,x,u(3,:),'u3',0)
     call save1Ddata(Nx,res_num,t,x,rho0,'rho0',0)

!     call save1tDdata(Nx,res_num,t,x,u(1,:),'phi',0)


! --->  EVOLUTION

     rk3 = .true.

     if (rk3) then
        num_steps = 3
     end if

! Main iteration loop.

  do l=1,Nt

!    Time.

     t = t + dt

!    Time step informtaion to screen.

     if (mod(l,every_1D).eq.0) then
        write(*,"(A5,I6,A6,ES9.2,A3)") ' |   ',l,'    | ',t,'  |'
     end if

!    Recycle variables.

     u_p = u
     a_p = a

! Use RK3

     do k=1,num_steps

!       Calling the value of the RHS for each of the equations

        call rhs_evol( )

!       Applying the Runge-Kutta integration

        if (k.eq.1) then
              dt_temp = dt
              u = u_p + dt_temp*rhs_u
              a = a_p + dt_temp*rhs_a
           else if (k.eq.2) then
              dt_temp = 0.25*dt
              u = 0.75*u_p + 0.25*u + dt_temp*rhs_u
              a = 0.75*a_p + 0.25*a + dt_temp*rhs_a
           else
              dt_temp = 2.0D0*dt/3.0D0
              u = u_p/3.0D0 + 2.0D0*u/3.0D0 + dt_temp*rhs_u
              a = a_p/3.0D0 + 2.0D0*a/3.0D0 + dt_temp*rhs_a
        end if


! Applying boundary conditions

! Pure extrapolation DEBUG (check with [1])

        u(:,0)  = 3.0D0*u(:,1) - 3.0D0*u(:,2) + u(:,3)
        u(:,Nx) = 3.0D0*u(:,Nx-1) - 3.0D0*u(:,Nx-2) + u(:,Nx-3)

        u(2,0) = zero

      end do

! ---> Some diagnostics

     call diagnostics()

! --->   SAVE DATA TO FILE   <---

!    Save 0D data (only every every_0D time steps).

     if (mod(l,every_0D).eq.0) then
        call save0Ddata(Nx,res_num,t,dx,a,'a',1)
        call save0Ddata(Nx,res_num,t,dx,alpha,'alpha',1)
        call save0Ddata(Nx,res_num,t,dx,rho0,'rho0',1)
        call save0Ddata(Nx,res_num,t,dx,m,'m',1)
       call save0Ddata(Nx,res_num,t,dx,Ham,'Ham',1)
     end if

!    Save 1D data (only every every_1D time steps).

     if (mod(l,every_1D).eq.0) then

       call save1Ddata(Nx,res_num,t,x,a,'a',1)
       call save1Ddata(Nx,res_num,t,x,m,'m',1)
       call save1Ddata(Nx,res_num,t,x,vel,'vel',1)
       call save1Ddata(Nx,res_num,t,x,alpha,'alpha',1)
       call save1Ddata(Nx,res_num,t,x,p,'p',1)
       call save1Ddata(Nx,res_num,t,x,Ham,'Ham',1)

       call save1Ddata(Nx,res_num,t,x,u(1,:),'u1',1)
       call save1Ddata(Nx,res_num,t,x,u(2,:),'u2',1)
       call save1Ddata(Nx,res_num,t,x,u(3,:),'u3',1)
       call save1Ddata(Nx,res_num,t,x,rho0,'rho0',1)
print *, 'Star_radius=', Star_radius
print *, 'Star_radius_etiqueta=', Star_radius_etiqueta
print *, 'Star_radius=', x(Star_radius_etiqueta)
     end if

  end do

  print *,'-----------------------------'


! --->   END   <---

  end subroutine evolve

