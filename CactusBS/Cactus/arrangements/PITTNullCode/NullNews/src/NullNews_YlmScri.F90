! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"


subroutine NullNews_YlmScri(CCTK_ARGUMENTS)

!  use cctk
  use  NullDecomp_sYlm

  implicit none

  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

    CCTK_INT, parameter :: izero = 0
    integer, parameter :: ik = kind(izero)

    call sYlm(2_ik,l_YlmScri,m_YlmScri,null_lsh(1),null_lsh(2),zeta,YlmScri_2)
    call sYlm(1_ik,l_YlmScri,m_YlmScri,null_lsh(1),null_lsh(2),zeta,YlmScri_1)
    call sYlm(0_ik,l_YlmScri,m_YlmScri,null_lsh(1),null_lsh(2),zeta,YlmScri_0)

end subroutine NullNews_YlmScri



