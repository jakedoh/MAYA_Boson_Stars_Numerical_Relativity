! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

module NullInterp_Interp
  use cctk
  implicit none

contains

  !! WARNING WARNING WARNING
  !! allways pass tmp_cgfn and tmp_cgfs in that order
  !! these must be the distributed arrays of the same name
  subroutine NullInterp_cinterp(cctkGH, tmp_cgfn, tmp_cgfs, jn, js, spin) 
    use NullGrid_Vars
    use NullInterp_InterpUtil
    implicit none

    CCTK_POINTER, intent(in) :: cctkGH
    CCTK_INT,     intent(in) :: spin

    CCTK_COMPLEX, dimension(:,:), intent(inout) :: tmp_cgfn, tmp_cgfs
    CCTK_COMPLEX, dimension(:,:), intent(inout) :: jn, js

    CCTK_INT, parameter :: izero = 0
    integer, parameter :: ik = kind(izero)

    integer, dimension(2), save :: gf_indices

    logical, save :: NotCalled = .true.
    integer :: ierr, i, j
    CCTK_COMPLEX, save :: inactive_value = (0.0, 0.0)

    DECLARE_CCTK_PARAMETERS

    if (NotCalled ) then
       NotCalled  = .false.

       ! note the order of gf_indices and the call to cinterp is important
       ! we want to interpolate from north to south and south to north.
       call CCTK_VarIndex(gf_indices(1), "NullInterp::tmp_cgfs")
       call CCTK_VarIndex(gf_indices(2), "NullInterp::tmp_cgfn")
    endif

    tmp_cgfn = jn(1:lsh(1),1:lsh(2))
    tmp_cgfs = js(1:lsh(1),1:lsh(2))

    call CCTK_SyncGroup(ierr, cctkGH, "NullInterp::tmp_cgf")
    if (ierr .ne. 0 ) then
       call CCTK_WARN(0, "Error syncing group tmp_cgf")
    endif

    if (poison_test.ne.0) inactive_value = poison_value * dcmplx(1.0, 1.0)

    call NullInterp_Util_cinterp(cctkGH, lsh, tmp_cgfn, tmp_cgfs, int(gf_indices,ik), spin, guard)

    do j = 1, lsh(2)
       do i = 1, lsh(1)
          if (EG(i,j).ne.0) then
             jn(i,j) = tmp_cgfn(i,j)
             js(i,j) = tmp_cgfs(i,j)
          else
             jn(i,j) = inactive_value
             js(i,j) = inactive_value
          end if
       end do
    end do
 
  end subroutine NullInterp_cinterp

  !! WARNING WARNING
  !! always pass tmp_rgfn and tmp_rgfs in that order. These must be the sistributed
  !! arrays of the same name
  subroutine NullInterp_rinterp(cctkGH, tmp_rgfn, tmp_rgfs, bn, bs)
    use NullGrid_Vars
    use NullInterp_InterpUtil
    implicit none

    CCTK_POINTER, intent(in) :: cctkGH

    CCTK_REAL, dimension(:,:), intent(inout) :: tmp_rgfn, tmp_rgfs
    CCTK_REAL, dimension(:,:), intent(inout) :: bn, bs

    CCTK_INT, parameter :: izero = 0
    integer, parameter :: ik = kind(izero)

    integer, save, dimension(2) :: gf_indices
    logical, save :: NotCalled = .true.
    integer :: ierr, i, j
    CCTK_REAL, save :: inactive_value = 0.0

    DECLARE_CCTK_PARAMETERS

    if (NotCalled ) then
       NotCalled  = .false.
       ! note the order of gf_indices and the call to cinterp is important
       ! we want to interpolate from north to south and south to north.
       call CCTK_VarIndex(gf_indices(1), "NullInterp::tmp_rgfs")
       call CCTK_VarIndex(gf_indices(2), "NullInterp::tmp_rgfn")
    endif
    tmp_rgfn = bn(1:lsh(1),1:lsh(2))
    tmp_rgfs = bs(1:lsh(1),1:lsh(2))
    call CCTK_SyncGroup(ierr, cctkGH, "NullInterp::tmp_rgf")
    if (ierr .ne. 0 ) then
       call CCTK_WARN(0, "Error syncing group tmp_rgf")
    endif

    if (poison_test.ne.0) inactive_value = poison_value

    call NullInterp_Util_rinterp(cctkGH, lsh, tmp_rgfn, tmp_rgfs, int(gf_indices,ik), guard)

    do j = 1, lsh(2)
       do i = 1, lsh(1)
          if (EG(i,j).ne.0) then
             bn(i,j) = tmp_rgfn(i,j)
             bs(i,j) = tmp_rgfs(i,j)
          else
             bn(i,j) = inactive_value
             bs(i,j) = inactive_value
          end if
       end do
    end do
 
  end subroutine NullInterp_rinterp

  !! WARNING WARNING WARNING
  !! allways pass tmp_cgfn and tmp_cgfs in that order
  !! these must be the distributed arrays of the same name
  subroutine NullInterp_3cinterp(cctkGH,&
       tmp_cgfn1, tmp_cgfs1,&
       tmp_cgfn2, tmp_cgfs2,&
       tmp_cgfn3, tmp_cgfs3,&
       jn1, js1, jn2, js2, jn3, js3,&
       spin1, spin2, spin3) 
    use NullGrid_Vars
    use NullInterp_InterpUtil
    implicit none

    CCTK_INT, parameter :: izero = 0
    integer, parameter :: ik = kind(izero)

    CCTK_POINTER, intent(in) :: cctkGH
    CCTK_INT,     intent(in) :: spin1, spin2, spin3

    CCTK_COMPLEX, dimension(:,:), intent(inout) ::&
         tmp_cgfn1, tmp_cgfs1,& 
         tmp_cgfn2, tmp_cgfs2,&
         tmp_cgfn3, tmp_cgfs3
    CCTK_COMPLEX, dimension(:,:), intent(inout) ::&
         jn1, js1, jn2, js2, jn3, js3

    integer, save, dimension(6) :: gf_indices
    logical, save :: NotCalled = .true.
    integer :: ierr, i, j
    CCTK_COMPLEX, save :: inactive_value = (0.0, 0.0)

    DECLARE_CCTK_PARAMETERS

!   if (spin1.ne.1.or.spin2.ne.1.or.spin3.ne.1) then
!      call CCTK_WARN(0, "this routine only knows about spin-one quantities")
!   end if

    if (NotCalled ) then
       NotCalled  = .false.

       ! note the order of gf_indices and the call to cinterp is important
       ! we want to interpolate from north to south and south to north.
       call CCTK_VarIndex(gf_indices(1), "NullInterp::tmp_cgfs1")
       call CCTK_VarIndex(gf_indices(2), "NullInterp::tmp_cgfn1")
       call CCTK_VarIndex(gf_indices(3), "NullInterp::tmp_cgfs2")
       call CCTK_VarIndex(gf_indices(4), "NullInterp::tmp_cgfn2")
       call CCTK_VarIndex(gf_indices(5), "NullInterp::tmp_cgfs3")
       call CCTK_VarIndex(gf_indices(6), "NullInterp::tmp_cgfn3")
    endif

    tmp_cgfn1 = jn1(1:lsh(1),1:lsh(2))
    tmp_cgfs1 = js1(1:lsh(1),1:lsh(2))
    tmp_cgfn2 = jn2(1:lsh(1),1:lsh(2))
    tmp_cgfs2 = js2(1:lsh(1),1:lsh(2))
    tmp_cgfn3 = jn3(1:lsh(1),1:lsh(2))
    tmp_cgfs3 = js3(1:lsh(1),1:lsh(2))

    call CCTK_SyncGroup(ierr, cctkGH, "NullInterp::tmp_cgf3")
    if (ierr .ne. 0 ) then
       call CCTK_WARN(0, "Error syncing group tmp_cgf3")
    endif

    if (poison_test.ne.0) inactive_value = poison_value * dcmplx(1.0, 1.0)

    call NullInterp_Util_3cinterp(cctkGH,lsh,&
         tmp_cgfn1, tmp_cgfs1,&
         tmp_cgfn2, tmp_cgfs2,&
         tmp_cgfn3, tmp_cgfs3,&
         int(gf_indices,ik), spin1, spin2, spin3, guard)

    do j = 1, lsh(2)
       do i = 1, lsh(1)
          if (EG(i,j).ne.0) then

             jn1(i,j) = tmp_cgfn1(i,j)
             js1(i,j) = tmp_cgfs1(i,j)

             jn2(i,j) = tmp_cgfn2(i,j)
             js2(i,j) = tmp_cgfs2(i,j)

             jn3(i,j) = tmp_cgfn3(i,j)
             js3(i,j) = tmp_cgfs3(i,j)

          else

             jn1(i,j) = inactive_value
             js1(i,j) = inactive_value

             jn2(i,j) = inactive_value
             js2(i,j) = inactive_value

             jn3(i,j) = inactive_value
             js3(i,j) = inactive_value

          end if
       end do
    end do
 
  end subroutine NullInterp_3cinterp

end module NullInterp_Interp
