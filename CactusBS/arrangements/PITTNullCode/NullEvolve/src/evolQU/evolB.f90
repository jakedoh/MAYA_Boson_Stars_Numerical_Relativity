subroutine evolB(Npts)
  implicit none
  integer, intent(in) :: Npts

  double precision, parameter :: xin = 0.45
  double precision, parameter :: xE = 0.5005

  integer :: mask(Npts)
  double precision :: xb(Npts), B(Npts), Bexact(Npts), BE, dx, rhs(Npts)

  integer i

  dx = (1-xin)/dble(Npts-1)
  do i = 1, Npts
    xb(i) = xin + (i-1)*dx
  end do

  rhs = cos(xb+0.5*dx)
  Bexact = sin(xb)
  BE = sin(xE)
  B = 0

  ! 'regular' evolution starts at x > xE + dx
  ! all interior points will be filled with the boundary stencil

  mask = 0; where (xb>xE+dx) mask = 1 ! evolution starts

  do i = 1, Npts
    if(mask(i).ne.0) then
      B(i) = B(i-1) + dx*rhs(i-1)
    else
      ! here i == B+1, B = i-1
      B(i) = BE + (xb(i)-xE) * cos(0.5*(xE+xb(i)))
    end if
  end do

  write (10000+Npts,*) xE, BE, 0
  do i = 1, Npts
    if(mask(i).ne.0) then
       write (10000+Npts,*) xb(i), B(i),         (B(i)-Bexact(i))
    end if
  end do
 
end subroutine evolB
    
program quevol
  implicit none

  call evolB(1001); call evolB(2001); call evolB(4001)
! call evolB(1001); call evolB(1201); call evolB(1441)
! call evolB(100+1); call evolB(120+1); call evolB(144+1)
! call evolB(101); call evolB(201); call evolB(401)

end program
