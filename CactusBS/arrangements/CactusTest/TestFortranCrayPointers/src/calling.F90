#include "cctk.h"
#include "cctk_Arguments.h"



module cctk_core
  implicit none

  ! This definition must stay in sync with the C definition in the flesh
  type, bind(C) :: cGH
     integer cctk_dim
     integer cctk_iteration

     CCTK_POINTER_TO_CONST cctk_gsh_ptr
     CCTK_POINTER_TO_CONST cctk_lsh_ptr
     CCTK_POINTER_TO_CONST cctk_lbnd_ptr
     CCTK_POINTER_TO_CONST cctk_ubnd_ptr

     CCTK_POINTER_TO_CONST cctk_ash_ptr
     integer cctk_alignment, cctk_alignment_offset
  end type cGH

contains

  pure function cctki_array_index(array, idx) result(value)
    implicit none
    CCTK_POINTER_TO_CONST, intent(in) :: array
    integer, intent(in) :: idx
    integer values(idx)
    pointer (array, values)
    integer value
    value = values(idx)
  end function cctki_array_index

  pure function cctki_array(array, dim) result(values)
    implicit none
    CCTK_POINTER_TO_CONST, intent(in) :: array
    integer, intent(in) :: dim
    integer array_values(dim)
    pointer (array, array_values)
    integer values(dim)
    values = array_values
  end function cctki_array

end module cctk_core



module cctk_core2
  use cctk
  use cctk_core
  implicit none

contains

  pure function cctki_ash(cctkGH, idx) result(value)
    implicit none
    CCTK_POINTER_TO_CONST, intent(in) :: cctkGH
    type(cGH) cctkGH_value
    pointer (cctkGH, cctkGH_value)
    integer, intent(in) :: idx
    integer value
    value = cctki_array_index(cctkGH_value % cctk_ash_ptr, idx)
  end function cctki_ash

  ! pure function cctki_ash2(cctkGH, idx) result(ash)
  !   implicit none
  !   CCTK_POINTER_TO_CONST, intent(in) :: cctkGH
  !   type(cGH) cctkGH_value
  !   pointer (cctkGH, cctkGH_value)
  !   integer, intent(in) :: idx
  !   integer ash(3)
  !   ash = cctki_array(cctkGH_value % cctk_ash_ptr, 3)
  ! end function cctki_ash2

  pure function cctki_group_ash(cctkGH, groupname, dim, idx) result(value)
    implicit none
    CCTK_POINTER_TO_CONST, intent(in) :: cctkGH
    character(*), intent(in) :: groupname
    integer, intent(in) :: dim
    integer, intent(in) :: idx
    integer value
    integer ash(dim)
    integer ierr
    call CCTK_GroupashGN(ierr, cctkGH, dim, ash, groupname)
    value = ash(idx)
  end function cctki_group_ash

end module cctk_core2



subroutine TestFortranCrayPointers_callee(cctkGH)
  use cctk
  use cctk_core
  use cctk_core2
  implicit none

  CCTK_POINTER_TO_CONST cctkGH
  type(cGH) cctkGH_value
  pointer (cctkGH, cctkGH_value)

  integer cctki_dim
  integer cctk_lsh
  cctk_lsh(cctki_dim) = &
       cctki_array_index(cctkGH_value % cctk_lsh_ptr, cctki_dim)

  integer, save :: cctki_varindex_gf = -1
  CCTK_POINTER_TO_CONST gf_ptr
  CCTK_REAL gf(cctki_ash(cctkGH, 1), cctki_ash(cctkGH, 2), cctki_ash(cctkGH, 3))
  pointer (gf_ptr, gf)

  integer, save :: cctki_varindex_ga = -1
  CCTK_POINTER_TO_CONST ga_ptr
  CCTK_REAL ga(cctki_group_ash(cctkGH, "TestFortranCrayPointers::ga", 2, 1), &
       cctki_group_ash(cctkGH, "TestFortranCrayPointers::ga", 2, 2))
  pointer (ga_ptr, ga)

  integer, save :: cctki_varindex_gs = -1
  CCTK_POINTER_TO_CONST gs_ptr
  CCTK_REAL gs
  pointer (gs_ptr, gs)

  integer ga_lsh(2)

  integer i, j, k
  integer ierr

  if (cctki_varindex_gf < 0) then
     call CCTK_VarIndex(cctki_varindex_gf, "TestFortranCrayPointers::gf")
  end if
  call CCTK_VarDataPtrI(gf_ptr, cctkGH, 0, cctki_varindex_gf)

  if (cctki_varindex_ga < 0) then
     call CCTK_VarIndex(cctki_varindex_ga, "TestFortranCrayPointers::ga")
  end if
  call CCTK_VarDataPtrI(ga_ptr, cctkGH, 0, cctki_varindex_ga)

  if (cctki_varindex_gs < 0) then
     call CCTK_VarIndex(cctki_varindex_gs, "TestFortranCrayPointers::gs")
  end if
  call CCTK_VarDataPtrI(gs_ptr, cctkGH, 0, cctki_varindex_gs)

  do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)
           gf(i,j,k) = 1
        end do
     end do
  end do

  call CCTK_GrouplshGN(ierr, cctkGH, 2, ga_lsh, "TestFortranCrayPointers::ga")
  if (ierr < 0) then
     call CCTK_ERROR("Cannot determine ga_lsh")
  end if

  do j = 1, ga_lsh(2)
     do i = 1, ga_lsh(1)
        ga(i,j) = 1
     end do
  end do

  gs = 1
end subroutine TestFortranCrayPointers_callee



subroutine TestFortranCrayPointers_caller(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  gf = 0
  ga = 0
  gs = 0
  call TestFortranCrayPointers_callee(cctkGH)
  if (any(gf /= 1) .or. any(ga /= 1) .or. gs /= 1) then
     call CCTK_ERROR("Grid variables are incorrect")
  end if
end subroutine TestFortranCrayPointers_caller
