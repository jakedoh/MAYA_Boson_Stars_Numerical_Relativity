 /*@@
   @file      SourceData.F90
   @date      
   @author    Gabrielle Allen
   @desc 
              Elliptic initial data for wave equation

              Originally written in F77 fixed format. Converted to F90
              free format by Peter Diener.
   @enddesc 
   @version $Header$
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

#include "EllBase.h"

subroutine UniformCharge(CCTK_ARGUMENTS)

! Find static field for a uniformly charge sphere
!
! That is, solve Nabla^2 phi = - 4 pi rho
! where rho = Q/(4/3 * Pi * R^3)
! where Q is the total charge and R is the sphere radius

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL pi
  CCTK_REAL AbsTol(3), RelTol(3)
  CCTK_REAL charge_factor

  integer iphi,iMcoeff,iNcoeff
  integer ierr
  CCTK_INT length

  character*30 fsolver
  character*200 infoline

! Get variable indices for all grid functions

  call CCTK_VarIndex (iMcoeff,  "idscalarwaveelliptic::Mcoeff")
  if (iMcoeff .lt. 0) then
    call CCTK_WARN(0,"Grid variable index for Mcoeff not found")
  end if

  call CCTK_VarIndex (iNcoeff,  "idscalarwaveelliptic::Ncoeff")
  if (iNcoeff .lt. 0) then
    call CCTK_WARN(0,"Grid variable index for Ncoeff not found")
  end if

  call CCTK_VarIndex (iphi,"wavetoy::phi")
  if (iphi .lt. 0) then
    call CCTK_WARN(0,"Grid variable index for iphi not found")
  end if


! Set up all coefficients and initial guess for solution

  pi = 4.0d0*atan(1.0d0)
  charge_factor = 4.0d0*pi*charge*3.0d0/(4.0d0*pi*radius**3)

  phi = 0.0d0
  Mcoeff = 0.0d0
  
  where ( r <= radius )
    Ncoeff = charge_factor
  elsewhere
    Ncoeff = 0.0d0
  end where

! Set tolerance for stopping elliptic solve

  AbsTol(1)=1.0d-5
  AbsTol(2)=1.0d-5
  AbsTol(3)=1.0d-5

  RelTol(1)=-1
  RelTol(2)=-1
  RelTol(3)=-1

! Set any options needed for the elliptic solve

  call Ell_SetStrKey (ierr, "yes", "EllLinFlat::Bnd::Robin")
  call Ell_SetRealKey(ierr, 0.0d0, "EllLinFlat::Bnd::Robin::inf")
  call Ell_SetIntKey (ierr, 1,     "EllLinFlat::Bnd::Robin::falloff")

! Set parameters for specific solvers
  if (CCTK_EQUALS(solver,"sor")) then
    call Ell_SetIntKey(ierr, sor_maxit,"Ell::SORmaxit")
  end if
     
  call CCTK_FortranString(length,solver,fsolver)
  write(infoline,'("Going into elliptic solver ",A)') fsolver
  call CCTK_INFO(infoline)

! Call elliptic solver to fill out phi

  call Ell_LinFlatSolver( &
         ierr, &
         cctkGH, & 
         iphi,  &
         iMcoeff, iNcoeff, &
         AbsTol, RelTol, &
         fsolver)

  if (ierr .eq. ELL_SUCCESS) then
    call CCTK_INFO("Leaving elliptic solver: solve successful")
  else if (ierr .eq. ELL_NOCONVERGENCE) then
    call CCTK_INFO("Leaving elliptic solver: solver failed to converge")
  else if (ierr .eq. ELL_NOSOLVER) then
    call CCTK_INFO("Elliptic solver not found")
  else 
    write(infoline,'("Leaving elliptic solver: solve failed (Error ",I3,")")') ierr
    call CCTK_INFO(infoline)
  end if

! Set up last timestep ... assume (first order) time symmetry

  phi_p = phi

! Output exact solution if required

  if (output_tmp .eq. 1) then

    where ( r .ge. radius )
      temp = charge / r
    elsewhere
      temp = charge/(2.0d0*radius**3)*(3.0d0*radius**2-r**2)
    end where

    call CCTK_OutputVarAsByMethod(ierr,cctkGH, &
            "idscalarwaveelliptic::temp", &
            "IOASCII_1D","phi_exact")
  end if      

  return
end
