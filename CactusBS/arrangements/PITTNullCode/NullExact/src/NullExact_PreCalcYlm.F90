! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

  ! set everything at the first few radial points

  subroutine NullExact_PreCalcYlm(CCTK_ARGUMENTS)

    use cctk
    use NullEvol_sYlm
    use NullGrid_Vars
    implicit none

    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    DECLARE_CCTK_FUNCTIONS

    CCTK_INT, parameter :: izero = 0
    integer, parameter :: ik = kind(izero)

    ! local variables, as needed:

    if (verbose.ne.0) then
      call CCTK_INFO("calculating Ylm's")
    endif

    call rsYlm(2_ik,l_in_Ylm,m_in_Ylm,null_lsh(1),null_lsh(2),zz,Ylm_2)
    call rsYlm(1_ik,l_in_Ylm,m_in_Ylm,null_lsh(1),null_lsh(2),zz,Ylm_1)
    call rsYlm(0_ik,l_in_Ylm,m_in_Ylm,null_lsh(1),null_lsh(2),zz,Ylm_0)

  end subroutine NullExact_PreCalcYlm
