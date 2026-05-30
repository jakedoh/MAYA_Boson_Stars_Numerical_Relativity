! Registration of symmetries for the necessary grid functions.
! $Header$

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine EHFinder_SetSym(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, dimension(3) :: sym

! All grid functions have even symmetries in all directions.
  sym = 1

! Set up symmetries for the level set function.
  call SetCartSymGN ( ierr, cctkGH, sym, 'ehfinder::f' )
  if ( ierr .gt. 0 ) then
    call CCTK_WARN(1,'Failed to register symmetry for the level set function')
  end if

! Set up symmetries for the temporary source function.
  call SetCartSymGN ( ierr, cctkGH, sym, 'ehfinder::sftmp' )
  if ( ierr .gt. 0 ) then
    call CCTK_WARN(1,'Failed to register symmetry for the temp source function')
  end if

! Set up symmetries for the mask function.
  call SetCartSymGN( ierr, cctkGH, sym, 'ehfinder::eh_mask' )
  if ( ierr .gt. 0 ) then
    call CCTK_WARN(1,'Failed to register symmetry for eh_mask')
  end if

! Set up symmetries for the surface index function.
  call SetCartSymGN(ierr,cctkGH,sym,'ehfinder::surface_index')
  if ( ierr .gt. 0 ) then
    call CCTK_WARN(1,'Failed to register symmetry for sc')
  end if

  return
end subroutine EHFinder_SetSym


subroutine EHFinder_ApplySymAll(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

! Select the variable groups for the None boundary condition. When the
! physical boundary condition (a noop in this case) is applied the
! symmetry boundary conditions will be applied automatically.

  ierr = Boundary_SelectGroupForBC ( cctkGH, &
       int(CCTK_ALL_FACES,ik), 1_ik, -1_ik, 'ehfinder::f', 'None' )
  if ( ierr /= 0 ) then
    call CCTK_WARN ( 0, 'Could not select f for boundary condition' )
  end if
  ierr = Boundary_SelectGroupForBC ( cctkGH, &
       int(CCTK_ALL_FACES,ik), 1_ik, -1_ik, 'ehfinder::eh_mask', 'None' )
  if ( ierr /= 0 ) then
    call CCTK_WARN ( 0, 'Could not select eh_mask for boundary condition' )
  end if

  return
end subroutine EHFinder_ApplySymAll


subroutine EHFinder_ApplySymF(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  ierr = Boundary_SelectGroupForBC ( cctkGH, &
       int(CCTK_ALL_FACES,ik), 1_ik, -1_ik, 'ehfinder::f', 'None' )
  if ( ierr /= 0 ) then
    call CCTK_WARN ( 0, 'Could not select f for boundary condition' )
  end if

  return
end subroutine EHFinder_ApplySymF


subroutine EHFinder_ApplySymFSFTMP(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  ierr = Boundary_SelectGroupForBC ( cctkGH, &
       int(CCTK_ALL_FACES,ik), 1_ik, -1_ik, 'ehfinder::f', 'None' )
  if ( ierr /= 0 ) then
    call CCTK_WARN ( 0, 'Could not select f for boundary condition' )
  end if

  ierr = Boundary_SelectGroupForBC ( cctkGH, &
       int(CCTK_ALL_FACES,ik), 1_ik, -1_ik, 'ehfinder::sftmp', 'None' )
  if ( ierr /= 0 ) then
    call CCTK_WARN ( 0, 'Could not select sftmp for boundary condition' )
  end if

  return
end subroutine EHFinder_ApplySymFSFTMP

subroutine EHFinder_ApplySymMask(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  ierr = Boundary_SelectGroupForBC ( cctkGH, &
       int(CCTK_ALL_FACES,ik), 1_ik, -1_ik, 'ehfinder::eh_mask', 'None' )

  if ( ierr /= 0 ) then
    call CCTK_WARN ( 0, 'Could not select eh_mask for boundary condition' )
  end if

  return
end subroutine EHFinder_ApplySymMask


subroutine EHFinder_ApplySymSC(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  ierr = Boundary_SelectGroupForBC ( cctkGH, &
       int(CCTK_ALL_FACES,ik), 1_ik, -1_ik, 'ehfinder::surface_index', 'None' )

  if ( ierr /= 0 ) then
    call CCTK_WARN ( 0, 'Could not select sc for boundary condition' )
  end if

  return
end subroutine EHFinder_ApplySymSC
