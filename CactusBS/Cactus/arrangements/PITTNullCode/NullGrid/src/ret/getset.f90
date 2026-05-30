subroutine getset(x,y)
  implicit none
  integer, intent(inout) :: x,y
  integer, save :: xo, yo
  logical, save :: was_called = .false.

  if (was_called) then
    x = xo
    y = yo
  else
    was_called = .true.
    xo = x
    yo = y
  endif
end subroutine getset
