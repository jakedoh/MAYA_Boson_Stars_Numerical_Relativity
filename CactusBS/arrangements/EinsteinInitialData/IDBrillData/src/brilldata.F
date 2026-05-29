/*@@
  @file      brilldata.F
  @date
  @author    Carsten Gundlach (Cactus 4, Miguel Alcubierre)
  @desc
             Construct Brill wave initial data.
  @enddesc
  @version   $Header$
@@*/

#include "cctk.h" 
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#include "CactusElliptic/EllBase/src/EllBase.h"

      subroutine brilldata(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

      integer ipsi,iMcoeff,iNcoeff
      integer metric_index(6)
      integer ierr

      CCTK_REAL AbsTol(3),RelTol(3)

c     Get indices for metric.

      call CCTK_VarIndex(metric_index(1), "admbase::gxx")
      call CCTK_VarIndex(metric_index(2), "admbase::gxy")
      call CCTK_VarIndex(metric_index(3), "admbase::gxz")
      call CCTK_VarIndex(metric_index(4), "admbase::gyy")
      call CCTK_VarIndex(metric_index(5), "admbase::gyz")
      call CCTK_VarIndex(metric_index(6), "admbase::gzz")

c     Get indices for grid functions.

      call CCTK_VarIndex(ipsi,"idbrilldata::brillpsi")
      if (ipsi.lt.0) then
        call CCTK_WARN(0,"Grid variable index for ipsi not found")
      end if

      call CCTK_VarIndex(iMcoeff,"idbrilldata::brillMlinear")
      if (iMcoeff.lt.0) then
        call CCTK_WARN(0,"Grid variable index for brillMlinear not found")
      end if

      call CCTK_VarIndex(iNcoeff,"idbrilldata::brillNsource")
      if (iNcoeff.lt.0) then
        call CCTK_WARN(0,"Grid variable index for brillNsource not found")
      end if

c     Set up background metric and coefficients for linear solve.

      if (CCTK_EQUALS(initial_data,"brilldata2D")) then
         call setupbrilldata2D(CCTK_ARGUMENTS)
      else
         call setupbrilldata3D(CCTK_ARGUMENTS)
      end if

      if (output_coeffs .eq. 1) then
         call CCTK_OutputVarAs(ierr,cctkGH,"IDBrillData::brillMlinear","Debug_brillMlinear")
         call CCTK_OutputVarAs(ierr,cctkGH,"IDBrillData::brillNsource","Debug_brillNsource")
         call CCTK_OutputVarAs(ierr,cctkGH,"ADMBase::gxx","Debug_gxx")
         call CCTK_OutputVarAs(ierr,cctkGH,"ADMBase::gxy","Debug_gxy")
         call CCTK_OutputVarAs(ierr,cctkGH,"ADMBase::gxz","Debug_gxz")
         call CCTK_OutputVarAs(ierr,cctkGH,"ADMBase::gyy","Debug_gyy")
         call CCTK_OutputVarAs(ierr,cctkGH,"ADMBase::gyz","Debug_gyz")
         call CCTK_OutputVarAs(ierr,cctkGH,"ADMBase::gzz","Debug_gzz")
      end if

c     Tolerances for elliptic solve.

      AbsTol(1)= thresh
      AbsTol(2)= -1.0D0
      AbsTol(3)= -1.0D0
 
      RelTol(1)= -1.0D0
      RelTol(2)= -1.0D0
      RelTol(3)= -1.0D0

c     Boundaries.

      call Ell_SetStrKey(ierr,"yes","EllLinMetric::Bnd::Robin")
      call Ell_SetRealKey(ierr,1.0D0,"EllLinMetric::Bnd::Robin::inf")
      call Ell_SetIntKey(ierr,1,"EllLinMetric::Bnd::Robin::falloff")

c     Elliptic solver 

c     Just in case some solvers do not use the standard interface
      conformal_state = 0

      if (CCTK_EQUALS(solver,"sor")) then
         call Ell_SetIntKey(ierr,sor_maxit,"Ell::SORmaxit")
         call Ell_LinMetricSolver(ierr,cctkGH, 
     .   metric_index,ipsi,iMcoeff,iNcoeff,AbsTol,RelTol,"sor")
      else if (CCTK_EQUALS(solver,"petsc")) then
         call Ell_LinMetricSolver(ierr,cctkGH, 
     .   metric_index,ipsi,iMcoeff,iNcoeff,AbsTol,RelTol,"petsc")
      else if (CCTK_EQUALS(solver,"bam")) then
         call Ell_LinMetricSolver(ierr,cctkGH, 
     .   metric_index,ipsi,iMcoeff,iNcoeff,AbsTol,RelTol,"bam")
      end if

c     Check for errors.

      if (ierr.eq.ELL_SUCCESS) then
        call CCTK_INFO("Leaving elliptic solver: solve successful")
      else if (ierr.eq.ELL_NOCONVERGENCE) then
        call CCTK_INFO("Leaving elliptic solver: solver failed to converge")
      else if (ierr.eq.ELL_NOSOLVER) then
        call CCTK_INFO("Elliptic solver not found")
      else 
        call CCTK_INFO("Leaving elliptic solver: solve failed") 
      end if

c     Synchronization and symmetry boundaries.

      call CCTK_SyncGroup(ierr,cctkGH,"idbrilldata::brillconf")
      call CartSymGN(ierr,cctkGH,"idbrilldata::brillconf")

c     Construct final metric.

      call IDBrillData_Finish(CCTK_ARGUMENTS)

      return
      end
