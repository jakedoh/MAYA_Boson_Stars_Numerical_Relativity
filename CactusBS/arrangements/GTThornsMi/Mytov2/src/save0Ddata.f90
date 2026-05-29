
  subroutine save0Ddata(Nx,res_num,t,dspace,yval,base_name,first_index)

! --->   SAVE 1D DATA TO FILE   <---

! This routine saves 1D data to the files,
! that is, the whole spatial arrays at a
! given time.

  use arrays

  implicit none

  character(len=20) filestatus
  logical firstcall
  data firstcall / .true. /
  save firstcall

  character(len=*), intent(IN) :: base_name
  real(kind=8), dimension(Nx), intent(IN) :: yval

!  logical, save :: first = .true.
!  integer, save :: first_index
  character(len=256) :: filename

  integer i,Nx,first_index,res_num
  real(kind=8) t,dspace
  real(kind=8) min, max, nm1, nm2

  if (res_num.eq.1) then
    filename = base_name // '_1.tl'
  else if (res_num.eq.2) then
    filename = base_name // '_2.tl'
  else if (res_num.eq.3) then
    filename = base_name // '_3.tl'
  else if (res_num.eq.4) then
    filename = base_name // '_4.tl'
  else if (res_num.eq.5) then
    filename = base_name // '_5.tl'
  end if

  if (first_index.eq.0) then
     filestatus = 'replace'
  else
     filestatus = 'old'
  end if


! --->   Calculating scalars

  max = yval(1)
  min = yval(1)
  nm1 = 0.0D0
  nm2 = 0.0D0

  do i=2,Nx+1
     if (yval(i)>max) max = yval(i)
     if (yval(i)<min) min = yval(i)
  end do

  do i=2,Nx+1
     nm1 = nm1 + 0.5D0*(dabs(yval(i-1)) + dabs(yval(i)))*dspace
     nm2 = nm2 + 0.5D0*(yval(i-1)**2 + yval(i)**2)*dspace
     if (yval(i)>max) max = yval(i)
     if (yval(i)<min) min = yval(i)
  end do

  nm2 = dsqrt(nm2)

! --->   Saving data   <---

  if (filestatus=='replace') then
     open(1,file=filename,form='formatted',status=filestatus)
  write(1,*) '#Time           maximum         minimum         nm1         nm2'
  else
     open(1,file=filename,form='formatted',status=filestatus,position='append')
  end if
     write(1,"(5ES14.6)") t,max,min,nm1,nm2
  close(1)

! --->   END   <---

  end subroutine save0Ddata

