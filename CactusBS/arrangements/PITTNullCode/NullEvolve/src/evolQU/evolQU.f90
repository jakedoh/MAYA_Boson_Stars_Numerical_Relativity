subroutine evolQU(Npts)
  implicit none
  integer, intent(in) :: Npts

  double precision, parameter :: Nrhs = 1.1;

! double precision, parameter :: qE = 1.9413402166757308;
! double precision, parameter :: uE = 1.8603493325402201;
  double precision :: qE, uE
  double precision, parameter :: xin = 0.45
  double precision, parameter :: xmax = 0.99
  double precision, parameter :: rwt = 20.0
  double precision, parameter :: xE = sqrt(3.0)/3.45;
  double precision, parameter :: rE = rwt*xE/(1-xE);
  double precision, parameter :: rinvE = 1.0/rE;

  integer :: mask(Npts)
  double precision :: xb(Npts), xbh(Npts), dx, dx_U, dx_UE, dxx_UE, RhsEph, dx_QE, dx_Uc, xc, rhsE, drhs
  double precision :: Q(Npts), U(Npts), Qexact(Npts), Uexact(Npts), rinvh(Npts), rinv(Npts), rhs(Npts)

  logical, parameter :: skipQ = .true.

  integer i

  ! grid setup

  dx = (xmax-xin)/dble(Npts-1)
  do i = 1, Npts
    xb(i) = xin + (i-1)*dx
  end do
  xbh = xb + 0.5*dx
  rinvh = (1-xbh)/rwt/xbh
  rinv  = (1-xb)/rwt/xb

  ! setup of the testbed
  rhs = rinvh**Nrhs;
  rhsE = rinvE**Nrhs;
  Qexact = rinv**Nrhs/dble(2-Nrhs)
  Uexact = rinvh**(1+Nrhs)/dble( (1+Nrhs)*(Nrhs-2) );
  qE = rinvE**Nrhs/dble(2-Nrhs)
  uE = rinvE**(1+Nrhs)/dble( (1+Nrhs)*(Nrhs-2) );

  ! integration algorithm

  mask = 0; where (xb>xE+1.25*dx) mask = 1

  if(skipQ) then
     Q = Qexact
  else
     Q = 0
     do i = 2, Npts
       if(mask(i).eq.0) then
          ! the RHS is known at xE and x(i+1/2)
          ! we interpolate onto the mid-point beween xE and x(i)
          drhs   = (rhs(i)-rhsE)/(xbh(i)-xE)
          RhsEph = rhsE + 0.5*(xb(i)-xE)*drhs
          Q(i)   = 0.5*RhsEph + ( qE -0.5*RhsEph ) * ((xb(i)-1)*xE/xb(i)/(xE-1))**2 
       end if
       Q(i) =    mask(i)  * (0.5*rhs(i) + ( Q(i-1) -0.5*rhs(i) ) * ((xb(i)-1)*xb(i-1)/xb(i)/(xb(i-1)-1))**2 )&
            + (1-mask(i)) * Q(i)
     end do
  end if

  U = Uexact
  do i = 2, Npts

    dx_U  = Q(i) / (rwt * xbh(i-1)*xbh(i))
    dx_UE = qE / (rwt*xE*xE)

    if(mask(i).eq.1) then
       U(i)  =  U(i-1) + dx * dx_U 
    else

       if(xE < xbh(i-1)) then
          ! we update the points i-1/2 based on dx_UE, xE and dx_U(i)
          ! then we update i+1/2 based on i-1/2 and dx_U(i)
   
          ! x(i-1/2) is outside the WT -- we update it
          xc = (xE + xbh(i-1))/2
          dx_Uc =  dx_UE * (xb(i)      -xc)/(xb(i)-xE) &
                +  dx_U          * (xE-xc)/(xE-xb(i))
          U(i-1) = uE    + dx_Uc * ( xbh(i-1)-xE )
          U(i)   = U(i-1) + dx_U * dx
       else
          ! we update the points i+1/2 based on dx_UE, xE and dx_U(i)
          if(abs(xE-xb(i)).gt.1.e-5*dx) then
             ! x(i-1/2) is outside the WT -- we update it
             xc = (xE + xbh(i))/2
             dx_Uc =  dx_UE * (xb(i)-xc)/(xb(i)-xE)&
                    + dx_U  * (xE-xc)/(xE-xb(i))
          else
             dx_Uc = dx_UE
          end if
          U(i) = uE + dx_Uc * ( xbh(i)-xE )
       end if
    end if


  end do

  if(.not.skipQ) then
     write (10000+Npts,"(A)") "# evolQU"
     write (10000+Npts,"(A)") "# [1] = x"
     write (10000+Npts,"(A)") "# [2] = Q"
     write (10000+Npts,"(A)") "# [3] = Qerr"
     write (10000+Npts,*) xE, qE, 0
     do i = 1, Npts
       if(mask(i).ne.0) then
         write (10000+Npts,*) xb(i), Q(i), mask(i)*(Q(i)-Qexact(i))
       end if
     end do
  end if
 
  write (30000+Npts,"(A)") "# evolQU"
  write (30000+Npts,"(A)") "# [1] = xh"
  write (30000+Npts,"(A)") "# [2] = U"
  write (30000+Npts,"(A)") "# [3] = Uerr"
  write (30000+Npts,*) xE, uE, 0
  do i = 1, Npts
    if(mask(i).ne.0) then
      write (30000+Npts,*) xbh(i), U(i), mask(i)*(U(i)-Uexact(i))
    end if
  end do

end subroutine evolQU
    
program quevol
  implicit none
! call evolQU(1001); call evolQU(2001); call evolQU(4001)
! call evolQU(1001); call evolQU(1201); call evolQU(1441)
! call evolQU(100+1); call evolQU(120+1); call evolQU(144+1)
  call evolQU(101); call evolQU(201); call evolQU(401)

end program
