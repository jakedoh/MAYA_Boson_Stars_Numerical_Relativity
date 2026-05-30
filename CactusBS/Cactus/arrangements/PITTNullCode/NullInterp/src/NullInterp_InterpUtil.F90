! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

module NullInterp_InterpUtil
  use cctk
  implicit none

contains

  subroutine NullInterp_Util_cinterp(cctkGH, lsh_, Jn, Js, gf_indices, spin, stereo_mask)

    use NullGrid_Vars

    implicit none

    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_PARAMETERS

    !input variables
    CCTK_POINTER, intent(in) :: cctkGH
    CCTK_INT, dimension(2), intent(in) :: lsh_
    CCTK_COMPLEX, dimension(lsh_(1),lsh_(2)), intent(inout) :: Jn, Js
    CCTK_INT, dimension(lsh_(1),lsh_(2)), intent(in) :: stereo_mask
    CCTK_INT, dimension(2), intent(in) :: gf_indices
    CCTK_INT, intent(in) :: spin
    CCTK_INT :: i, j, k

    logical, save :: FirstTime = .true.
    CCTK_REAL, dimension(:), allocatable, save :: qcoord, pcoord
    CCTK_COMPLEX, dimension(:), allocatable, save :: Iarr1, Iarr2
    CCTK_INT, dimension(:), allocatable, save :: reindexi, reindexj
    CCTK_INT, save :: tot_nbp, nbp, mask_count

    ! interpolator variables
    integer :: ierror, N_in_arrays, N_out_arrays, N_interp_points, N_dims
    integer, save :: interp_handle, coord_system_handle
    integer, save :: param_table_handle = -1
    CCTK_POINTER, dimension(2) :: interp_coords
    integer :: interp_coords_type_code
    CCTK_INT, dimension(2) :: in_array_indices
    CCTK_INT, dimension(2) :: out_array_type_codes
    CCTK_POINTER, dimension(2) ::  out_arrays


    if (skip_interpolation.ne.0) then
       if (FirstTime) call CCTK_INFO("will skip interpolation");
       FirstTime = .false.
       do i = 1, lsh_(1)
          do j = 1, lsh_(2)
             if (stereo_mask(i,j).eq.1) then
                if (poison_test.ne.0) then
                   Jn(i,j) = poison_value
                   Js(i,j) = poison_value
                else
                   Jn(i,j) = 0
                   Js(i,j) = 0
                end if
             endif
          end do
       end do
       return;
    end if

    !    CCTK_REAL :: qs_min, qs_max,ps_min, ps_max
    if ( FirstTime ) then
       FirstTime = .false.

       mask_count = count(stereo_mask.eq.1)
       tot_nbp = mask_count
       allocate(qcoord(tot_nbp), pcoord(tot_nbp), reindexi(tot_nbp),&
               reindexj(tot_nbp), Iarr1(tot_nbp), Iarr2(tot_nbp))

       nbp = 0
       do i = 1, lsh_(1)
          do j = 1, lsh_(2)
             if (stereo_mask(i,j).eq.1) then
                nbp = nbp + 1
                qcoord(nbp) = dble(qs(i,j)/(qs(i,j)**2 + ps(i,j)**2))
                pcoord(nbp) = dble(-ps(i,j)/(qs(i,j)**2 + ps(i,j)**2))
                reindexi(nbp) = i
                reindexj(nbp) = j
             endif
          end do
       end do

       if ( nbp > tot_nbp ) then
          call CCTK_WARN(0,"Boundary point determination error")
       endif

    end if ! FirstTime

    if (param_table_handle < 0) then
       call  Util_TableCreateFromString(param_table_handle,&
            "order = " // char(ichar('0') + interpolation_order))

       if (param_table_handle < 0) then
          call  CCTK_WARN(-1, "can t create parameter table!")
       end if
    end if

    !interpolator call

    interp_handle =-1
    coord_system_handle = -1

! choice between Lagrange and Hermite interpolator    
    if (CCTK_EQUALS(interpolation_name, "Lagrange")) then
       call CCTK_InterpHandle (interp_handle, "generalized polynomial interpolation")
    else if (CCTK_EQUALS(interpolation_name, "Hermite")) then
       call CCTK_InterpHandle (interp_handle, "Hermite polynomial interpolation")
    end if

    if (interp_handle < 0) then
       call CCTK_WARN(0,"Interpolation operator not found")
    end if

    call CCTK_CoordSystemHandle (coord_system_handle, "stereo")
    if (coord_system_handle < 0) then
       call CCTK_WARN(0,"Coordinate system 'stereo' not registered")
    end if

    interp_coords(1)        = CCTK_PointerTo(qcoord)
    interp_coords(2)        = CCTK_PointerTo(pcoord)
    interp_coords_type_code = CCTK_VARIABLE_REAL
    in_array_indices        = gf_indices
    out_arrays(1)           = CCTK_PointerTo(Iarr1)
    out_arrays(2)           = CCTK_PointerTo(Iarr2)
    out_array_type_codes    = CCTK_VARIABLE_COMPLEX

    N_dims=2
    N_in_arrays=2
    N_out_arrays=2
    N_interp_points = nbp

    call CCTK_InterpGridArrays (ierror, cctkGH, N_dims,&
         interp_handle, param_table_handle, coord_system_handle,&
         N_interp_points, interp_coords_type_code,&
         interp_coords,&
         N_in_arrays, in_array_indices,&
         N_out_arrays, out_array_type_codes, out_arrays)

    if (ierror < 0) then
       call CCTK_WARN(1,"Interpolation error")
    endif

    do k = 1, nbp
       i = reindexi(k) 
       j = reindexj(k) 

       if (spin .ne. 0) then
          Jn(i,j) = Iarr1(k) * (-zz(i,j)/conjg(zz(i,j)))**spin
          Js(i,j) = Iarr2(k) * (-zz(i,j)/conjg(zz(i,j)))**spin
       else
          Jn(i,j) = Iarr1(k)
          Js(i,j) = Iarr2(k)
       endif

    end do

  end subroutine NullInterp_Util_cinterp

  subroutine NullInterp_Util_rinterp(cctkGH, lsh_, Jn, Js, gf_indices, stereo_mask)
 
    use NullGrid_Vars

    implicit none

    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_PARAMETERS

    CCTK_POINTER, intent(in) :: cctkGH
    CCTK_INT, dimension(2), intent(in) :: lsh_
    CCTK_REAL, dimension(lsh_(1),lsh_(2)), intent(inout) :: Jn, Js
    CCTK_INT, dimension(lsh_(1),lsh_(2)), intent(in) :: stereo_mask
    CCTK_INT, dimension(2), intent(in) :: gf_indices

    CCTK_INT :: i, j, k
    logical, save :: FirstTime = .true.
    CCTK_REAL, dimension(:), allocatable, save :: qcoord, pcoord
    CCTK_REAL, dimension(:), allocatable, save :: Iarr1, Iarr2
    CCTK_INT, dimension(:), allocatable, save :: reindexi, reindexj
    CCTK_INT, save :: tot_nbp, nbp, mask_count

    ! interpolator variables
    integer :: ierror, N_in_arrays, N_out_arrays, N_interp_points, N_dims
    integer, save :: interp_handle, coord_system_handle
    integer, save :: param_table_handle = -1
    CCTK_POINTER, dimension(2) :: interp_coords
    integer :: interp_coords_type_code
    CCTK_INT, dimension(2) :: in_array_indices
    CCTK_INT, dimension(2) :: out_array_type_codes
    CCTK_POINTER, dimension(2) ::  out_arrays

    if (skip_interpolation.ne.0) then
       if (FirstTime) call CCTK_INFO("will skip interpolation");
       FirstTime = .false.
       do i = 1, lsh_(1)
          do j = 1, lsh_(2)
             if (stereo_mask(i,j).eq.1) then
                if (poison_test.ne.0) then
                   Jn(i,j) = poison_value
                   Js(i,j) = poison_value
                else
                   Jn(i,j) = 0
                   Js(i,j) = 0
                end if
             endif
          end do
       end do
       return;
    end if

    if ( FirstTime ) then
       FirstTime = .false.

       mask_count = count(stereo_mask.eq.1)
       tot_nbp = mask_count 
!       tot_nbp = 2*lsh_(1)+2*lsh_(2) -4
       allocate(qcoord(tot_nbp), pcoord(tot_nbp), reindexi(tot_nbp), &
               reindexj(tot_nbp), Iarr1(tot_nbp), Iarr2(tot_nbp))

       nbp = 0
       do i = 1, lsh_(1) 
          do j = 1, lsh_(2) 
             if (stereo_mask(i,j).eq.1) then
                nbp = nbp + 1
                qcoord(nbp) = dble(qs(i,j)/(qs(i,j)**2 + ps(i,j)**2))
                pcoord(nbp) = dble(-ps(i,j)/(qs(i,j)**2 + ps(i,j)**2))
                reindexi(nbp) = i
                reindexj(nbp) = j
             endif
          end do
       end do

       if ( nbp > tot_nbp ) then
          call CCTK_WARN(0,"Boundary point determination error")
       end if

    end if ! FirstTime

    if (param_table_handle < 0) then
       call  Util_TableCreateFromString(param_table_handle,&
            "order = " // char(ichar('0') + interpolation_order))

       if (param_table_handle < 0) then
          call  CCTK_WARN(-1, "can t create parameter table!")
       end if
    end if

    !interpolator call

    interp_handle =-1
    coord_system_handle = -1

! choice between Lagrange and Hermite interpolator    
    if (CCTK_EQUALS(interpolation_name, "Lagrange")) then
       call CCTK_InterpHandle (interp_handle, "generalized polynomial interpolation")
    else if (CCTK_EQUALS(interpolation_name, "Hermite")) then
       call CCTK_InterpHandle (interp_handle, "Hermite polynomial interpolation")
    end if

    if (interp_handle < 0) then
       call CCTK_WARN(0,"Interpolation operator not found")
    end if

    call CCTK_CoordSystemHandle (coord_system_handle, "stereo")
    if (coord_system_handle < 0) then
       call CCTK_WARN(0,"Coordinate system 'stereo' not registered")
    end if

    interp_coords(1)        = CCTK_PointerTo(qcoord)
    interp_coords(2)        = CCTK_PointerTo(pcoord)
    interp_coords_type_code = CCTK_VARIABLE_REAL
    in_array_indices        = gf_indices
    out_arrays(1)           = CCTK_PointerTo(Iarr1)
    out_arrays(2)           = CCTK_PointerTo(Iarr2)
    out_array_type_codes    = CCTK_VARIABLE_REAL

    N_dims=2
    N_in_arrays=2
    N_out_arrays=2
    N_interp_points = nbp

    call CCTK_InterpGridArrays (ierror, cctkGH, N_dims,&
         interp_handle, param_table_handle, coord_system_handle,&
         N_interp_points, interp_coords_type_code,&
         interp_coords,&
         N_in_arrays, in_array_indices,&
         N_out_arrays, out_array_type_codes, out_arrays)

    if (ierror < 0) then
       call CCTK_WARN(1,"Interpolation error")
    endif

    do k = 1, nbp
       i = reindexi(k) 
       j = reindexj(k) 
       Jn(i,j) = Iarr1(k)
       Js(i,j) = Iarr2(k)
    end do

  end subroutine NullInterp_Util_rinterp

  subroutine NullInterp_Util_3cinterp(cctkGH, lsh_, J1n, J1s, J2n, J2s, J3n, J3s, gf_indices, spin1, spin2, spin3, stereo_mask)

    use NullGrid_Vars

    implicit none

    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_PARAMETERS

    !input variables
    CCTK_POINTER,               intent(in) :: cctkGH
    CCTK_INT,     dimension(2), intent(in) :: lsh_
    CCTK_INT,                   intent(in) :: spin1, spin2, spin3

    CCTK_COMPLEX, dimension(lsh_(1),lsh_(2)), intent(inout) :: J1n, J1s,&
         J2n, J2s, J3n, J3s
    CCTK_INT, dimension(lsh_(1),lsh_(2)), intent(in) :: stereo_mask
    CCTK_INT, dimension(6), intent(in) :: gf_indices
    CCTK_INT :: i, j, k

    logical, save :: FirstTime = .true.
    CCTK_REAL, dimension(:), allocatable, save :: qcoord, pcoord
    CCTK_COMPLEX, dimension(:), allocatable, save :: Iarr1, Iarr2,&
         Iarr3, Iarr4, Iarr5, Iarr6
    CCTK_INT, dimension(:), allocatable, save :: reindexi, reindexj
    CCTK_INT, save :: tot_nbp, nbp, mask_count

    !CCTK_INT, intent(in) :: spin

    ! interpolator variables
    integer :: ierror, N_in_arrays, N_out_arrays, N_interp_points, N_dims
    integer, save :: interp_handle, param_table_handle, coord_system_handle
    CCTK_POINTER, dimension(2) :: interp_coords
    integer :: interp_coords_type_code
    CCTK_INT, dimension(6) :: in_array_indices
    CCTK_INT, dimension(6) :: out_array_type_codes
    CCTK_POINTER, dimension(6) ::  out_arrays

    if (skip_interpolation.ne.0) then
       if (FirstTime) call CCTK_INFO("will skip interpolation");
       FirstTime = .false.
       do i = 1, lsh_(1)
          do j = 1, lsh_(2)
             if (stereo_mask(i,j).eq.1) then
                if (poison_test.ne.0) then
                   J1n(i,j) = poison_value
                   J1s(i,j) = poison_value
                else
                   J1n(i,j) = 0
                   J1s(i,j) = 0
                end if
             endif
          end do
       end do
       return;
    end if

    if ( FirstTime ) then
       FirstTime = .false.

       mask_count = count(stereo_mask.eq.1)
       tot_nbp = mask_count
       allocate(qcoord(tot_nbp), pcoord(tot_nbp),&
               reindexi(tot_nbp), reindexj(tot_nbp),&
               Iarr1(tot_nbp), Iarr2(tot_nbp),&
               Iarr3(tot_nbp), Iarr4(tot_nbp),&
               Iarr5(tot_nbp), Iarr6(tot_nbp))

       nbp = 0
       do i = 1, lsh_(1)
          do j = 1, lsh_(2)
             if (stereo_mask(i,j).eq.1) then
                nbp = nbp + 1
                qcoord(nbp) = qs(i,j)/(qs(i,j)**2 + ps(i,j)**2)
                pcoord(nbp) = -ps(i,j)/(qs(i,j)**2 + ps(i,j)**2)
                reindexi(nbp) = i
                reindexj(nbp) = j
             endif
          end do
       end do

       if ( nbp > tot_nbp ) then
          call CCTK_WARN(0,"Boundary point determination error")
       endif

       param_table_handle = -1
       call  Util_TableCreateFromString(param_table_handle,&
            "order = " // char(ichar('0') + interpolation_order))

       if (param_table_handle < 0) then
          call  CCTK_WARN(-1, "can t create parameter table!")
       end if

    end if ! FirstTime

    !interpolator call

!    interp_handle =-1
!    coord_system_handle = -1

! choice between Lagrange and Hermite interpolator    
!    if (CCTK_EQUALS(interpolation_name, "Lagrange")) then
!       call CCTK_InterpHandle (interp_handle, "generalized polynomial interpolation")
!    else if (CCTK_EQUALS(interpolation_name, "Hermite")) then
!       call CCTK_InterpHandle (interp_handle, "Hermite polynomial interpolation")
!    end if
!
!    if (interp_handle < 0) then
!       call CCTK_WARN(0,"Interpolation operator not found")
!    end if

!    call CCTK_CoordSystemHandle (coord_system_handle, "stereo")

    !interpolator call

    interp_handle =-1
    coord_system_handle = -1

! choice between Lagrange and Hermite interpolator    
    if (CCTK_EQUALS(interpolation_name, "Lagrange")) then
       call CCTK_InterpHandle (interp_handle, "generalized polynomial interpolation")
    else if (CCTK_EQUALS(interpolation_name, "Hermite")) then
       call CCTK_InterpHandle (interp_handle, "Hermite polynomial interpolation")
    end if

    if (interp_handle < 0) then
       call CCTK_WARN(0,"Interpolation operator not found")
    end if

    call CCTK_CoordSystemHandle (coord_system_handle, "stereo")
    if (coord_system_handle < 0) then
       call CCTK_WARN(0,"Coordinate system 'stereo' not registered")
    end if

    interp_coords(1)        = CCTK_PointerTo(qcoord)
    interp_coords(2)        = CCTK_PointerTo(pcoord)
    interp_coords_type_code = CCTK_VARIABLE_REAL
    in_array_indices        = gf_indices
    out_arrays(1)           = CCTK_PointerTo(Iarr1)
    out_arrays(2)           = CCTK_PointerTo(Iarr2)
    out_arrays(3)           = CCTK_PointerTo(Iarr3)
    out_arrays(4)           = CCTK_PointerTo(Iarr4)
    out_arrays(5)           = CCTK_PointerTo(Iarr5)
    out_arrays(6)           = CCTK_PointerTo(Iarr6)
    out_array_type_codes    = CCTK_VARIABLE_COMPLEX

    N_dims=2
    N_in_arrays=6
    N_out_arrays=6
    N_interp_points = nbp

    call CCTK_InterpGridArrays (ierror, cctkGH, N_dims,&
         interp_handle, param_table_handle, coord_system_handle,&
         N_interp_points, interp_coords_type_code,&
         interp_coords,&
         N_in_arrays, in_array_indices,&
         N_out_arrays, out_array_type_codes, out_arrays)

    if (ierror < 0) then
       call CCTK_WARN(1,"Interpolation error")
    endif

    do k = 1, nbp
       i = reindexi(k) 
       j = reindexj(k) 

       if (spin1 .ne. 0) then
          J1n(i,j) = Iarr1(k) * (-zz(i,j)/conjg(zz(i,j)))**spin1
          J1s(i,j) = Iarr2(k) * (-zz(i,j)/conjg(zz(i,j)))**spin1
       else
          J1n(i,j) = Iarr1(k)
          J1s(i,j) = Iarr2(k)
       endif

       if (spin2 .ne. 0) then
          J2n(i,j) = Iarr3(k) * (-zz(i,j)/conjg(zz(i,j)))**spin2
          J2s(i,j) = Iarr4(k) * (-zz(i,j)/conjg(zz(i,j)))**spin2
       else
          J2n(i,j) = Iarr3(k)
          J2s(i,j) = Iarr4(k)
       endif

       if (spin3 .ne. 0) then
          J3n(i,j) = Iarr5(k) * (-zz(i,j)/conjg(zz(i,j)))**spin3
          J3s(i,j) = Iarr6(k) * (-zz(i,j)/conjg(zz(i,j)))**spin3
       else
          J3n(i,j) = Iarr5(k)
          J3s(i,j) = Iarr6(k)
       endif


       ! J1(i,j) = Iarr1(k) * (-zz(i,j)/conjg(zz(i,j))) !**spin
       ! J2(i,j) = Iarr2(k) * (-zz(i,j)/conjg(zz(i,j))) !**spin
       ! J3(i,j) = Iarr3(k) * (-zz(i,j)/conjg(zz(i,j))) !**spin
       ! J4(i,j) = Iarr4(k) * (-zz(i,j)/conjg(zz(i,j))) !**spin
       ! J5(i,j) = Iarr5(k) * (-zz(i,j)/conjg(zz(i,j))) !**spin
       ! J6(i,j) = Iarr6(k) * (-zz(i,j)/conjg(zz(i,j))) !**spin

    end do

  end subroutine NullInterp_Util_3cinterp

end module NullInterp_InterpUtil
