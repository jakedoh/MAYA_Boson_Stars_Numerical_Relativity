
  subroutine save1Ddata(Nx,res_num,t,xval,yval,base_name,first_index)

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
  real(kind=8), dimension(Nx), intent(IN) :: xval, yval

!  logical, save :: first = .true.
!  integer, save :: first_index
  character(len=256) :: filename

  integer i,Nx, first_index,res_num
  real(kind=8) t

  if (res_num.eq.1) then
    filename = base_name // '_1.xl'
  else if (res_num.eq.2) then
    filename = base_name // '_2.xl'
  else if (res_num.eq.3) then
    filename = base_name // '_3.xl'
  else if (res_num.eq.4) then
    filename = base_name // '_4.xl'
  else if (res_num.eq.5) then
    filename = base_name // '_5.xl'
  end if

  if (first_index.eq.0) then
     filestatus = 'replace'
  else
     filestatus = 'old'
  end if


! --->   Saving data   <---

  if (filestatus=='replace') then
     open(1,file=filename,form='formatted',status=filestatus)
  else
     open(1,file=filename,form='formatted',status=filestatus,position='append')
  end if
  write(1,*) ''
  write(1,*) '#Time = ',t
  do i=1,Nx+1,2**(res_num-1)
     write(1,"(2ES14.6)") xval(i),yval(i)
  end do
  write(1,*)
  close(1)

! --->   END   <---

  end subroutine save1Ddata

