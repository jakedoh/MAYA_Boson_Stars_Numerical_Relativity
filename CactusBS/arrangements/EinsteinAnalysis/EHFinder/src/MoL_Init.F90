! Registration of variables with MoL
! $Header$

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine EHFinder_MoLRegister(CCTK_ARGUMENTS)

  implicit none

  CCTK_INT :: ierr_cum, varindex, rhsindex

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  ierr_cum = 0

  call CCTK_GroupIndex ( varindex, 'ehfinder::f' )
  call CCTK_GroupIndex ( rhsindex, 'ehfinder::sf' )

  ierr_cum = ierr_cum + MolRegisterEvolvedGroup ( varindex, rhsindex )

  if ( evolve_generators .gt. 0 ) then

    if ( CCTK_EQUALS( generator_distribution, 'line' ) ) then

      call CCTK_GroupIndex (varindex, 'ehfinder::xg')
      call CCTK_GroupIndex(rhsindex, 'ehfinder::dxg')

      ierr_cum = ierr_cum + MoLRegisterEvolvedGroup(varindex, rhsindex)

      call CCTK_GroupIndex (varindex, 'ehfinder::yg')
      call CCTK_GroupIndex(rhsindex, 'ehfinder::dyg')

      ierr_cum = ierr_cum + MoLRegisterEvolvedGroup(varindex, rhsindex)

      call CCTK_GroupIndex (varindex, 'ehfinder::zg')
      call CCTK_GroupIndex(rhsindex, 'ehfinder::dzg')

      ierr_cum = ierr_cum + MoLRegisterEvolvedGroup(varindex, rhsindex)

    else if ( CCTK_EQUALS( generator_distribution, '2D array' ) ) then

      call CCTK_GroupIndex (varindex, 'ehfinder::xg2')
      call CCTK_GroupIndex(rhsindex, 'ehfinder::dxg2')

      ierr_cum = ierr_cum + MoLRegisterEvolvedGroup(varindex, rhsindex)

      call CCTK_GroupIndex (varindex, 'ehfinder::yg2')
      call CCTK_GroupIndex(rhsindex, 'ehfinder::dyg2')

      ierr_cum = ierr_cum + MoLRegisterEvolvedGroup(varindex, rhsindex)

      call CCTK_GroupIndex (varindex, 'ehfinder::zg2')
      call CCTK_GroupIndex(rhsindex, 'ehfinder::dzg2')

      ierr_cum = ierr_cum + MoLRegisterEvolvedGroup(varindex, rhsindex)

    end if

  end if
 
  if ( ierr_cum .gt. 0 ) then
    call CCTK_WARN(0,'Problems registering variables with MoL')
  end if
end subroutine EHFinder_MoLRegister


