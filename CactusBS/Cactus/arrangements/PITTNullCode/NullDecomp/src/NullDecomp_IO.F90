! vim: syntax=fortran

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

module NullDecomp_IO
   implicit none

 contains

subroutine NullDecomp_WriteCoefFile(FieldName, TruncateFile, Lmin, Lmax, YlmCoef, time)
  character(*), intent(in) :: FieldName
  logical,      intent(in) :: TruncateFile 
  CCTK_INT,     intent(in) :: Lmin, Lmax
  CCTK_REAL,    intent(in) :: time
  CCTK_COMPLEX, dimension(Lmin:Lmax, -Lmax:Lmax) :: YlmCoef

  character(1000) DirName, FileName, lName, mName

  integer :: flen, nchar, l, m
  logical :: file_exists

  DECLARE_CCTK_PARAMETERS

  call CCTK_FortranString(nchar, out_dir, DirName)
  DirName = trim(DirName)

  do l = Lmin, Lmax
     do m = -l, l

       write(lName, '(i2.2)') l
       write(mName, '(i2.2)') abs(m)
       if(m.ge.0) then
         mName='p'//mName(1:2)
       else
         mName='m'//mName(1:2)
       end if

       flen = nchar+1+len(FieldName) +2+2+1+3+4

       FileName = DirName(1:nchar) //'/' //FieldName&
           //'.L'//lName(1:2)//'M'//mName(1:3) // '.asc'

       inquire(FILE=FileName(1:flen), EXIST=file_exists) 
       if (TruncateFile .or. .not. file_exists) then
         open(10, FILE=FileName(1:flen), STATUS='REPLACE')
       else 
         open(10, FILE=FileName(1:flen), STATUS='OLD', POSITION='APPEND')
       end if
       write(10, '(30g25.17)', advance='NO') time
       write(10, '(30g25.17, 30g25.17)', advance='NO') dble (YlmCoef(l, m)), dimag(YlmCoef(l, m))
       write(10,*)
       close(10)
            
     end do
  end do

end subroutine NullDecomp_WriteCoefFile


subroutine NullDecomp_WriteCoefRadiusFile(FieldName, TruncateFile, Lmin, Lmax, YlmCoef, time, radius)
  character(*), intent(in) :: FieldName
  logical,      intent(in) :: TruncateFile 
  CCTK_INT,     intent(in) :: Lmin, Lmax
  CCTK_REAL,    intent(in) :: time, radius
  CCTK_COMPLEX, dimension(Lmin:Lmax, -Lmax:Lmax) :: YlmCoef

  character(1000) DirName, FileName, lName, mName

  integer :: flen, nchar, l, m
  logical :: file_exists

  DECLARE_CCTK_PARAMETERS

  call CCTK_FortranString(nchar, out_dir, DirName)
  DirName = trim(DirName)

  do l = Lmin, Lmax
     do m = -l, l

       write(lName, '(i2.2)') l
       write(mName, '(i2.2)') abs(m)
       if(m.ge.0) then
         mName='p'//mName(1:2)
       else
         mName='m'//mName(1:2)
       end if

       flen = nchar+1+len(FieldName) +2+2+1+3+4

       FileName = DirName(1:nchar) //'/' //FieldName&
           //'.L'//lName(1:2)//'M'//mName(1:3) // '.asc'
 
       inquire(FILE=FileName(1:flen), EXIST=file_exists) 
       if (TruncateFile .or. .not. file_exists) then
         open(10, FILE=FileName(1:flen), STATUS='REPLACE')
       else 
         open(10, FILE=FileName(1:flen), STATUS='OLD', POSITION='APPEND')
       end if
       write(10, '(30g25.17)', advance='NO') time
       write(10, '(30g25.17)', advance='NO') radius
       write(10, '(30g25.17, 30g25.17)', advance='NO') dble (YlmCoef(l, m)), dimag(YlmCoef(l, m))
       write(10,*)
       close(10)
            
     end do
  end do

end subroutine NullDecomp_WriteCoefRadiusFile


end module
