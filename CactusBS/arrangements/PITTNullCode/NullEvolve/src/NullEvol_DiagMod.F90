! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"


module NullEvol_DiagMod
   implicit none

 contains

  subroutine NullEvol_DiagReArray(cctkGH, FieldName, TruncateFile, F3Dn, F3Ds, tmp, xgrid, time)
    use NullDecomp_Vars, only: Lmax
    use NullDecomp_SpinDecomp, only: SpinDecompCoefs
    use NullDecomp_IO
    use NullGrid_Vars, only: lsh, nx, zz
 
    implicit none

    CCTK_INT, parameter :: izero = 0
    integer, parameter :: ik = kind(izero)

    CCTK_POINTER, intent(in) :: cctkGH
    character(*), intent(in) :: FieldName
    logical,      intent(in) :: TruncateFile
    CCTK_REAL,    intent(in) :: xgrid(nx), time
    CCTK_REAL,    dimension(lsh(1), lsh(2), nx) :: F3Dn, F3Ds
    CCTK_COMPLEX, dimension(lsh(1), lsh(2), 2)  :: tmp
 
    CCTK_COMPLEX, dimension(0:Lmax, -Lmax:Lmax) :: coef

    ! bracket x_diag

    integer i
    double precision c1,c2,c3

    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_PARAMETERS

    i = nint((Diagnostics_Coord_x-xgrid(1)) / (xgrid(2)-xgrid(1)))
    if(i.eq.1) i=2
    if(i.eq.nx) i=nx-1

    c1 = (Diagnostics_Coord_x-xgrid(i))*(Diagnostics_Coord_x-xgrid(i+1))&
            /(xgrid(i-1)-xgrid(i))/(xgrid(i-1)-xgrid(i+1))
    c2 = (Diagnostics_Coord_x-xgrid(i-1))*(Diagnostics_Coord_x-xgrid(i+1))&
             /(xgrid(i)-xgrid(i-1))/(xgrid(i)-xgrid(i+1))
    c3 = (Diagnostics_Coord_x-xgrid(i-1))*(Diagnostics_Coord_x-xgrid(i))&
             /(xgrid(i+1)-xgrid(i-1))/(xgrid(i+1)-xgrid(i))

    tmp(:,:,1) = F3Dn(:,:,i-1)*c1 + F3Dn(:,:,i)*c2 + F3Dn(:,:,i+1)*c3
    tmp(:,:,2) = F3Ds(:,:,i-1)*c1 + F3Ds(:,:,i)*c2 + F3Ds(:,:,i+1)*c3
    call SpinDecompCoefs(cctkGH, lsh(1), lsh(2), 0_ik, zz, tmp, coef)

    ! Only output if this is processor 0
    if (CCTK_MyProc(cctkGH) .gt. 0) then
       return
    endif

   ! write on the output file
    call NullDecomp_WriteCoefFile(FieldName, TruncateFile, 0_ik, Lmax, coef, time)

  end subroutine 

  subroutine NullEvol_DiagImArray(cctkGH, FieldName, TruncateFile, SpinWeight, F3Dn, F3Ds, tmp, xgrid, time)
    use NullDecomp_Vars, only: Lmax
    use NullDecomp_SpinDecomp, only: SpinDecompCoefs
    use NullDecomp_IO
    use NullGrid_Vars, only: lsh, nx, zz
 
    implicit none

    CCTK_POINTER, intent(in) :: cctkGH
    character(*), intent(in) :: FieldName
    logical,      intent(in) :: TruncateFile
    CCTK_INT,     intent(in) :: SpinWeight
    CCTK_REAL,    intent(in) :: xgrid(nx), time
    CCTK_COMPLEX, dimension(lsh(1), lsh(2), nx) :: F3Dn, F3Ds
    CCTK_COMPLEX, dimension(lsh(1), lsh(2), 2)  :: tmp
 
    CCTK_COMPLEX, dimension(SpinWeight:Lmax, -Lmax:Lmax) :: coef

    ! bracket x_diag

    integer i
    double precision c1,c2,c3

    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_PARAMETERS

    i = nint((Diagnostics_Coord_x-xgrid(1)) / (xgrid(2)-xgrid(1)))
    if(i.eq.1) i=2
    if(i.eq.nx) i=nx-1

    c1 = (Diagnostics_Coord_x-xgrid(i))*(Diagnostics_Coord_x-xgrid(i+1))&
            /(xgrid(i-1)-xgrid(i))/(xgrid(i-1)-xgrid(i+1))
    c2 = (Diagnostics_Coord_x-xgrid(i-1))*(Diagnostics_Coord_x-xgrid(i+1))&
             /(xgrid(i)-xgrid(i-1))/(xgrid(i)-xgrid(i+1))
    c3 = (Diagnostics_Coord_x-xgrid(i-1))*(Diagnostics_Coord_x-xgrid(i))&
             /(xgrid(i+1)-xgrid(i-1))/(xgrid(i+1)-xgrid(i))

    tmp(:,:,1) = F3Dn(:,:,i-1)*c1 + F3Dn(:,:,i)*c2 + F3Dn(:,:,i+1)*c3
    tmp(:,:,2) = F3Ds(:,:,i-1)*c1 + F3Ds(:,:,i)*c2 + F3Ds(:,:,i+1)*c3
    call SpinDecompCoefs(cctkGH, lsh(1), lsh(2), SpinWeight, zz, tmp, coef)

    ! Only output if this is processor 0
    if (CCTK_MyProc(cctkGH) .gt. 0) then
       return
    endif

   ! write on the output file
    call NullDecomp_WriteCoefFile(FieldName, TruncateFile, SpinWeight, Lmax, coef, time)

  end subroutine 

end module NullEvol_DiagMod
