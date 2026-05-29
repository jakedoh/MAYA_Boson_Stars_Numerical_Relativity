! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

  subroutine NullSHRE_PreCalcNSYlm (CCTK_ARGUMENTS)

    use cctk
    use NullEvol_sYlm
    use NullGrid_Vars
    implicit none

    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    DECLARE_CCTK_FUNCTIONS

    CCTK_INT, parameter :: izero = 0
    integer, parameter :: ik = kind(izero)

    CCTK_INT:: l, m, n 
    CCTK_COMPLEX, dimension(lsh(1), lsh(2), 2) :: NSYlm, NS1Ylm

      do l = 0, l_max
        do m = -l, l

           call sYlm(0_ik, l, m, lsh(1), lsh(2), zz, NSYlm)

           if (l.ge.1) then  
             call sYlm(1_ik, l, m, lsh(1), lsh(2), zz, NS1Ylm)
           else if (l.eq.0) then 
             NS1Ylm = 0.d0
           end if

           n = l*l + l + m + 1
             
           reYN(:,:,n) = dble(NSYlm(:,:,ip_n))
           imYN(:,:,n) = dimag(NSYlm(:,:,ip_n))
           reYS(:,:,n) = dble(NSYlm(:,:,ip_s))
           imYS(:,:,n) = dimag(NSYlm(:,:,ip_s))
             
           re1YN(:,:,n) = dble(NS1Ylm(:,:,ip_n))
           im1YN(:,:,n) = dimag(NS1Ylm(:,:,ip_n))
           re1YS(:,:,n) = dble(NS1Ylm(:,:,ip_s))
           im1YS(:,:,n) = dimag(NS1Ylm(:,:,ip_s))

       end do
    end do


  end subroutine NullSHRE_PreCalcNSYlm
