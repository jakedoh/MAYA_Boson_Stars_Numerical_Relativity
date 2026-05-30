! Calculation of centered differences.
! $Header$

# include "physical_part.h"

idx = half / cctk_delta_space(1)
idy = half / cctk_delta_space(2)
idz = half / cctk_delta_space(3)

do l = 1, eh_number_level_sets
  do k = kzl, kzr
    do j = jyl, jyr
      do i = ixl, ixr
        if ( eh_mask(i,j,k,l) .ge. 0 ) then
          if ( ( .not. btest(eh_mask(i,j,k,l),0) ) .and. &
               ( .not. btest(eh_mask(i,j,k,l),1) ) ) then
            dfx(i,j,k,l) = idx * ( f(i+1,j,k,l) - f(i-1,j,k,l) )
          else if ( btest(eh_mask(i,j,k,l),0) ) then
            dfx(i,j,k,l) = idx * ( -three * f(i,j,k,l) + &
                                  four * f(i+1,j,k,l) - f(i+2,j,k,l) )
          else
            dfx(i,j,k,l) = idx * ( three * f(i,j,k,l) - &
                                 four * f(i-1,j,k,l) + f(i-2,j,k,l) )
          end if
          if ( ( .not. btest(eh_mask(i,j,k,l),2) ) .and. &
               ( .not. btest(eh_mask(i,j,k,l),3) ) ) then
            dfy(i,j,k,l) = idy * ( f(i,j+1,k,l) - f(i,j-1,k,l) )
          else if ( btest(eh_mask(i,j,k,l),2) ) then
            dfy(i,j,k,l) = idy * ( -three * f(i,j,k,l) + &
                                  four * f(i,j+1,k,l) - f(i,j+2,k,l) )
          else
            dfy(i,j,k,l) = idy * ( three * f(i,j,k,l) - &
                                  four * f(i,j-1,k,l) + f(i,j-2,k,l) )
          end if
          if ( ( .not. btest(eh_mask(i,j,k,l),4) ) .and. &
               ( .not. btest(eh_mask(i,j,k,l),5) ) ) then
            dfz(i,j,k,l) = idz * ( f(i,j,k+1,l) - f(i,j,k-1,l) )
          else if ( btest(eh_mask(i,j,k,l),4) ) then
            dfz(i,j,k,l) = idz * ( -three * f(i,j,k,l) + &
                                   four * f(i,j,k+1,l) - f(i,j,k+2,l) )
          else
            dfz(i,j,k,l) = idz * ( three * f(i,j,k,l) - &
                                  four * f(i,j,k-1,l) + f(i,j,k-2,l) )
          end if
        else
          dfx(i,j,k,l) = zero
          dfy(i,j,k,l) = zero
          dfz(i,j,k,l) = zero
        end if
      end do
    end do
  end do
end do
