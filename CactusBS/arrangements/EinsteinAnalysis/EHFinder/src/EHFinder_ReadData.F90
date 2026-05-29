! Read metric, lapse, shift and conformal factor from files.
! $Header$

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

!This routine reads in the metric (from ADMBase).
subroutine EHFinder_Read_Metric(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  character(len=512) :: in_files, in_vars
  character(len=64), dimension(6)  :: var_names
  character(len=10) :: iteration_string
  CCTK_INT :: i, nc, ntot, res

! Figure out which iteration number to read, based on the parameters
! last_iteration_number (the last iteration in the numerical evolution
! producing the metric) and saved_iteration_every (how often was the metric
! saved) and the current iteration and save it in a string variable.
  i = min ( last_iteration_number - saved_iteration_every * cctk_iteration, &
            last_iteration_number )
  write(iteration_string,'(i10)') i

! Trim the string variable.
  iteration_string = adjustl(iteration_string)
  nc = len_trim(iteration_string)

! Generate a string with the filenames. Note at present the requirement
! therefore is that all iterations of one variable is written in the
! same file.

  in_files = 'gxx gxy gxz gyy gyz gzz'
  if ( CCTK_EQUALS(file_type,'sep_time_files') ) then
    call create_filenames ( in_files, i )
  end if

! Fill in the string array, used to requesting the metric components at
! the specified iteration number.
  var_names(1) = 'admbase::gxx{cctk_iteration='//iteration_string(1:nc)//'}'
  var_names(2) = 'admbase::gxy{cctk_iteration='//iteration_string(1:nc)//'}'
  var_names(3) = 'admbase::gxz{cctk_iteration='//iteration_string(1:nc)//'}'
  var_names(4) = 'admbase::gyy{cctk_iteration='//iteration_string(1:nc)//'}'
  var_names(5) = 'admbase::gyz{cctk_iteration='//iteration_string(1:nc)//'}'
  var_names(6) = 'admbase::gzz{cctk_iteration='//iteration_string(1:nc)//'}'

! merge all the variable names into a single string.
  in_vars = ' '
  ntot = 0
  do i = 1, 6
    nc = len_trim(var_names(i))
    in_vars(ntot+1:ntot+1+nc+1) = var_names(i)(1:nc+1)
    ntot = ntot + nc + 1
  end do
   
! Call the routine that actualle does the file access. Note failures are
! just silently ignored.
  call IOUtil_RecoverVarsFromDatafiles ( res, cctkGH, in_files, in_vars )

end subroutine EHFinder_Read_Metric


! This routine reads in the lapse (from ADMBase).
subroutine EHFinder_Read_Lapse(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  character(len=128) :: in_files, in_vars
  character(len=10) :: iteration_string
  CCTK_INT :: i, nc, res

! Figure out which iteration number to read, based on the parameters
! last_iteration_number (the last iteration in the numerical evolution
! producing the metric) and saved_iteration_every (how often was the lapse
! saved) and the current iteration and save it in a string variable.
  i = min ( last_iteration_number - saved_iteration_every * cctk_iteration, &
            last_iteration_number )
  write(iteration_string,'(i10)') i 

! Trim the string variable.
  iteration_string = adjustl(iteration_string)
  nc = len_trim(iteration_string)

! Generate a string with the filename. Note at present the requirement
! therefore is that all iterations of one variable is written in the
! same file.
  in_files = 'alp'
  if ( CCTK_EQUALS(file_type,'sep_time_files') ) then
    call create_filenames ( in_files, i )
  end if

! Fill in the string used to requesting the lapse at the specified
! iteration number.
  in_vars = 'admbase::alp{cctk_iteration='//iteration_string(1:nc)//'}'

! Call the routine that actualle does the file access. Note failures are
! just silently ignored.
  call IOUtil_RecoverVarsFromDatafiles ( res, cctkGH, in_files, in_vars )

end subroutine EHFinder_Read_Lapse


! This routine reads in all the shift components (from ADMBAse).
subroutine EHFinder_Read_Shift(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  character(len=256) :: in_files, in_vars
  character(len=64), dimension(3)  :: var_names
  character(len=10) :: iteration_string
  CCTK_INT :: i, nc, ntot, res

! Figure out which iteration number to read, based on the parameters
! last_iteration_number (the last iteration in the numerical evolution
! producing the metric) and saved_iteration_every (how often was the metric
! saved) and the current iteration and save it in a string variable.
  i = min ( last_iteration_number - saved_iteration_every * cctk_iteration, &
            last_iteration_number )
  write(iteration_string,'(i10)') i 

! Trim the string variable.
  iteration_string = adjustl(iteration_string)
  nc = len_trim(iteration_string)

! Generate a string with the filenames. Note at present the requirement
! therefore is that all iterations of one variable is written in the
! same file.
  in_files = 'betax betay betaz'
  if ( CCTK_EQUALS(file_type,'sep_time_files') ) then
    call create_filenames ( in_files, i )
  end if

! Fill in the string array, used to requesting the shift components at
! the specified iteration number.
  var_names(1) = 'admbase::betax{cctk_iteration='//iteration_string(1:nc)//'}'
  var_names(2) = 'admbase::betay{cctk_iteration='//iteration_string(1:nc)//'}'
  var_names(3) = 'admbase::betaz{cctk_iteration='//iteration_string(1:nc)//'}'

! merge all the variable names into a single string.
  in_vars = ' '
  ntot = 0
  do i = 1, 3
    nc = len_trim(var_names(i))
    in_vars(ntot+1:ntot+1+nc+1) = var_names(i)(1:nc+1)
    ntot = ntot + nc + 1
  end do
 
! Call the routine that actualle does the file access. Note failures are
! just silently ignored.
  call IOUtil_RecoverVarsFromDatafiles ( res, cctkGH, in_files, in_vars )

end subroutine EHFinder_Read_Shift


! This routine reads in the conformal factor (from StaticConformal).
subroutine EHFinder_Read_Conformal(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  character(len=128) :: in_files, in_vars
  character(len=10) :: iteration_string
  CCTK_INT :: i, nc, res

! Figure out which iteration number to read, based on the parameters
! last_iteration_number (the last iteration in the numerical evolution
! producing the metric) and saved_iteration_every (how often was the metric
! saved) and the current iteration and save it in a string variable. Note
! that for this variable the default option is to only to read it in once at
! the first timestep.
  if ( read_conformal_factor_once .gt. 0 ) then
    write(iteration_string,'(i10)') 0
  else
    i = min ( last_iteration_number - saved_iteration_every * cctk_iteration, &
              last_iteration_number )
    write(iteration_string,'(i10)') i 
  end if

! Trim the string variable.
  iteration_string = adjustl(iteration_string)
  nc = len_trim(iteration_string)

! Generate a string with the filenames. Note at present the requirement
! therefore is that all iterations of one variable is written in the
! same file.
  in_files = 'psi'
  if ( CCTK_EQUALS(file_type,'sep_time_files') ) then
    call create_filenames ( in_files, i )
  end if

! Fill in the string array, used to requesting the conformal factor at
! the specified iteration number.
  in_vars = 'staticconformal::psi{cctk_iteration='//iteration_string(1:nc)//'}'

! Call the routine that actualle does the file access. Note failures are
! just silently ignored.
  call IOUtil_RecoverVarsFromDatafiles ( res, cctkGH, in_files, in_vars )

end subroutine EHFinder_Read_Conformal


! This routine reads in the excision mask (from SpaceMask).
subroutine EHFinder_Read_Mask(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  character(len=128) :: in_files, in_vars
  character(len=10) :: iteration_string
  CCTK_INT :: i, nc, res

! Figure out which iteration number to read, based on the parameters
! last_iteration_number (the last iteration in the numerical evolution
! producing the metric) and saved_iteration_every (how often was the metric
! saved) and the current iteration and save it in a string variable.
  i = min ( last_iteration_number - saved_iteration_every * cctk_iteration, &
            last_iteration_number )
  write(iteration_string,'(i10)') i 

! Trim the string variable.
  iteration_string = adjustl(iteration_string)
  nc = len_trim(iteration_string)

! Generate a string with the filename. Note at present the requirement
! therefore is that all iterations of one variable is written in the
! same file.
  in_files = 'emask'
  if ( CCTK_EQUALS(file_type,'sep_time_files') ) then
    call create_filenames ( in_files, i )
  end if

! Fill in the string array, used to requesting the old style mask at
! the specified iteration number.
  in_vars = 'spacemask::emask{cctk_iteration='//iteration_string(1:nc)//'}'

! Call the routine that actualle does the file access. Note failures are
! just silently ignored.
  call IOUtil_RecoverVarsFromDatafiles ( res, cctkGH, in_files, in_vars )

end subroutine EHFinder_Read_Mask

subroutine create_filenames ( in_files, i )

  implicit none

  character(len=*), intent(INOUT) :: in_files
  CCTK_INT, intent(IN) :: i

  character(len=len(in_files)) :: tmp1, tmp2, istr
  character(len=len_trim(in_files)) :: old_files
  character(len=7), parameter :: template1 = '_000000'
  character(len=7) :: template2
  logical :: last
  CCTK_INT :: j, nc, lnew_file, ltmp, ltmp2

  last = .false.
  old_files = in_files
  lnew_file = 0
  ltmp = len(old_files)
  ltmp2 = 0
  tmp2 = ''

  write(istr,'(i10)') i
  istr = adjustl(istr)
  nc = len_trim(istr)
  template2 = template1
  template2(8-nc:7) = istr(1:nc)

  j = scan (old_files(1:ltmp), ' ')
  name_loop: do
    if ( j == 0 ) last = .true.
!    print*,'last = ', last
    tmp1 = ''
    if ( .not. last  ) then
      tmp1 = old_files(1:j-1)
      old_files(1:ltmp-j) = old_files(j+1:ltmp)
    else
      tmp1 = old_files(1:ltmp)
      old_files = ''
    end if
!    print*,'tmp1 = ',tmp1(1:j-1)
!    print*,'old_files = ',old_files(1:ltmp)
    ltmp = ltmp - j
!    print*,'ltmp = ',ltmp
    tmp2 = adjustl(trim(tmp2)//' '//trim(tmp1)//template2)
!    print*,'tmp2 = ',trim(tmp2)
!    print* 
    if ( last ) exit name_loop
    j = scan (old_files(1:ltmp), ' ') 
  end do name_loop 

  in_files = tmp2
!  pause
end subroutine create_filenames 
