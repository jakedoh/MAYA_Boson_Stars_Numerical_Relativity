#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"


      
c     This alternative version avoids the CalcTmunu.inc 
c     mechanism altogether. We just make sure Tmunu and T2munu
c     are initialized as zero.

      subroutine TmunuBase_AltSetTmunu (CCTK_ARGUMENTS)
      implicit none
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
      
      integer i, j, k
      
      do k = 1, cctk_lsh(3)
         do j = 1, cctk_lsh(2)
            do i = 1, cctk_lsh(1)
               
               eTtt(i,j,k) = 0
               
               eTtx(i,j,k) = 0
               eTty(i,j,k) = 0
               eTtz(i,j,k) = 0
               
               eTxx(i,j,k) = 0
               eTxy(i,j,k) = 0
               eTxz(i,j,k) = 0
               eTyy(i,j,k) = 0
               eTyz(i,j,k) = 0
               eTzz(i,j,k) = 0
               
            end do
         end do
      end do
      
      if (stress_energy_2_storage .ne. 0) then

        eT2tt = eTtt
        
        eT2tx = eTtx
        eT2ty = eTty
        eT2tz = eTtz
        
        eT2xx = eTxx
        eT2xy = eTxy
        eT2xz = eTxz
        eT2yy = eTyy
        eT2yz = eTyz
        eT2zz = eTzz
      end if
      
      end
