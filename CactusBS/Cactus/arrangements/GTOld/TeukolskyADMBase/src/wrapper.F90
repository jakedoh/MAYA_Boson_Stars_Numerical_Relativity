#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine add_teukolsky( CCTK_ARGUMENTS )
      implicit none
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

    integer                  i, j, k
    CCTK_REAL                xx(3)
    CCTK_REAL                gg(3,3)        
    CCTK_REAL                kk(3,3)        
    CCTK_REAL                alph
    CCTK_REAL                beta(3)
    CCTK_REAL                psi
    CCTK_REAL                psi4
    CCTK_REAL                rr
    REAL                     random
    logical,save        ::   first = .true.
    character(len=100)       ::   buf

    do k = 1, cctk_lsh(3)
      do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)

           xx(1) = x(i,j,k)
           xx(2) = y(i,j,k)
           xx(3) = z(i,j,k)

           call calcTeukolskyWave( cctk_time, xx, gg, kk, alph, beta )

!            if (add_schwarzschild > 0) then
!               rr = (r(i,j,k)**4.d0 + epsilon**4.d0)**0.25.d0
!               psi = 1.d0 + mass / (2.d0 * rr)
!               psi4 = psi**4.d0
!               gg(1,1) = gg(1,1) + psi4
!               gg(2,2) = gg(2,2) + psi4
!               gg(3,3) = gg(3,3) + psi4
!            end if

           gxx(i,j,k) = addtomet*(gxx(i,j,k) - 1.d0) + gg(1,1)
           gxy(i,j,k) = addtomet*gxy(i,j,k) + gg(1,2)
           gxz(i,j,k) = addtomet*gxz(i,j,k) + gg(1,3)
           gyy(i,j,k) = addtomet*(gyy(i,j,k) - 1.d0) + gg(2,2)
           gyz(i,j,k) = addtomet*gyz(i,j,k) + gg(2,3)
           gzz(i,j,k) = addtomet*(gzz(i,j,k) - 1.d0) + gg(3,3)

           kxx(i,j,k) = addtomet*kxx(i,j,k) + kk(1,1)
           kxy(i,j,k) = addtomet*kxy(i,j,k) + kk(1,2)
           kxz(i,j,k) = addtomet*kxz(i,j,k) + kk(1,3)
           kyy(i,j,k) = addtomet*kyy(i,j,k) + kk(2,2)
           kyz(i,j,k) = addtomet*kyz(i,j,k) + kk(2,3)
           kzz(i,j,k) = addtomet*kzz(i,j,k) + kk(3,3)

           if(add_noise.ne.0) then
             call RANDOM_NUMBER(random);gxx(i,j,k) = 1.d0 + (gxx(i,j,k) - 1.d0) * (1d0 +  (random-0.5)*noise_scale * CCTK_DELTA_SPACE(1)**noise_power)
             if(first) then
               write(buf,*) 'Random start ',random
               call CCTK_INFO(buf)
               first = .false.
             end if
             call RANDOM_NUMBER(random);gxy(i,j,k) = gxy(i,j,k) * (1d0 +  (random-0.5)*noise_scale * CCTK_DELTA_SPACE(1)**noise_power)
             call RANDOM_NUMBER(random);gxz(i,j,k) = gxz(i,j,k) * (1d0 +  (random-0.5)*noise_scale * CCTK_DELTA_SPACE(1)**noise_power)
             call RANDOM_NUMBER(random);gyy(i,j,k) = 1.d0 + (gyy(i,j,k) - 1.d0) * (1d0 +  (random-0.5)*noise_scale * CCTK_DELTA_SPACE(1)**noise_power)
             call RANDOM_NUMBER(random);gyz(i,j,k) = gyz(i,j,k) * (1d0 +  (random-0.5)*noise_scale * CCTK_DELTA_SPACE(1)**noise_power)
             call RANDOM_NUMBER(random);gzz(i,j,k) = 1.d0 + (gzz(i,j,k) - 1.d0) * (1d0 +  (random-0.5)*noise_scale * CCTK_DELTA_SPACE(1)**noise_power)
             call RANDOM_NUMBER(random);kxx(i,j,k) = kxx(i,j,k) * (1d0 +  (random-0.5)*noise_scale * CCTK_DELTA_SPACE(1)**noise_power)
             call RANDOM_NUMBER(random);kxy(i,j,k) = kxy(i,j,k) * (1d0 +  (random-0.5)*noise_scale * CCTK_DELTA_SPACE(1)**noise_power)
             call RANDOM_NUMBER(random);kxz(i,j,k) = kxz(i,j,k) * (1d0 +  (random-0.5)*noise_scale * CCTK_DELTA_SPACE(1)**noise_power)
             call RANDOM_NUMBER(random);kyy(i,j,k) = kyy(i,j,k) * (1d0 +  (random-0.5)*noise_scale * CCTK_DELTA_SPACE(1)**noise_power)
             call RANDOM_NUMBER(random);kyz(i,j,k) = kyz(i,j,k) * (1d0 +  (random-0.5)*noise_scale * CCTK_DELTA_SPACE(1)**noise_power)
             call RANDOM_NUMBER(random);kzz(i,j,k) = kzz(i,j,k) * (1d0 +  (random-0.5)*noise_scale * CCTK_DELTA_SPACE(1)**noise_power)
           end if

!            alp(i,j,k) = alph
!            betax(i,j,k) = beta(1)
!            betay(i,j,k) = beta(2)
!            betaz(i,j,k) = beta(3)

        end do
      end do
    end do

end subroutine add_teukolsky
