! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"


subroutine NullSHRE_SchPreCalc(CCTK_ARGUMENTS)

  use cctk
  use NullGrid_Vars
  use NullSHRE_modSchYlm 
  implicit none

    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    DECLARE_CCTK_FUNCTIONS
    CCTK_INT :: n, l, m, k

    CCTK_COMPLEX, dimension(lsh(1), lsh(2), 2) :: SchY
    CCTK_COMPLEX, dimension(lsh(1), lsh(2), 2) :: Sch1Y

#define R cr
#define M mass 
 
    do l = 0, lmax
      do m = -l, l

         n = l*l + l + m + 1

         call SchYlm (l, m, lsh(1), lsh(2), qs, ps, pp, SchY)
         call Sch1Ylm (l, m, lsh(1), lsh(2), qs, ps, pp, Sch1Y)
               
              reTYN(:,:,n) = real(SchY(:,:,ip_n))
              imTYN(:,:,n) = imag(SchY(:,:,ip_n))
              reTYS(:,:,n) = real(SchY(:,:,ip_s))
              imTYS(:,:,n) = imag(SchY(:,:,ip_s))
 
              re1TYN(:,:,n) = real(Sch1Y(:,:,ip_n))
              im1TYN(:,:,n) = imag(Sch1Y(:,:,ip_n))
              re1TYS(:,:,n) = real(Sch1Y(:,:,ip_s))
              im1TYS(:,:,n) = imag(Sch1Y(:,:,ip_s))
 
      end do
    end do
 
 end subroutine NullSHRE_SchPreCalc

