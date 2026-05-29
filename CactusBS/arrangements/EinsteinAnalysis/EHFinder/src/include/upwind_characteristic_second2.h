! Calculation of characteristic upwinded differences except at boundaries.
! $Header$

! If this is an active cell we have to compute its derivative.
if ( eh_mask(i,j,k,l) .ge. 0 ) then

  ! If it is not a boundary cell in the x-direction then calculate as a
  ! starting point the centered derivative in the x-direction.
  if ( ( .not. btest(eh_mask(i,j,k,l),0) ) .and. &
       ( .not. btest(eh_mask(i,j,k,l),1) ) ) then
    dfx(i,j,k,l) = idx * ( f(i+1,j,k,l) - f(i-1,j,k,l) )

!   Otherwise we will have to calculate the appropriate one sided derivative.
  else if ( btest(eh_mask(i,j,k,l),0) ) then
    dfx(i,j,k,l) = idx * ( - three * f(i,j,k,l) + four * f(i+1,j,k,l) - f(i+2,j,k,l) )
  else
    dfx(i,j,k,l) = idx * ( three * f(i,j,k,l) - four * f(i-1,j,k,l) + f(i-2,j,k,l) )
  end if

  ! Do the same for the y direction.
  if ( ( .not. btest(eh_mask(i,j,k,l),2) ) .and. &
       ( .not. btest(eh_mask(i,j,k,l),3) ) ) then
    dfy(i,j,k,l) = idy * ( f(i,j+1,k,l) - f(i,j-1,k,l) )
  else if ( btest(eh_mask(i,j,k,l),2) ) then
    dfy(i,j,k,l) = idy * ( - three * f(i,j,k,l) + four * f(i,j+1,k,l) - f(i,j+2,k,l) )
  else
    dfy(i,j,k,l) = idy * ( three * f(i,j,k,l) - four * f(i,j-1,k,l) + f(i,j-2,k,l) )
  end if

  ! And for the z direction.
  if ( ( .not. btest(eh_mask(i,j,k,l),4) ) .and. &
       ( .not. btest(eh_mask(i,j,k,l),5) ) ) then
    dfz(i,j,k,l) = idz * ( f(i,j,k+1,l) - f(i,j,k-1,l) )
  else if ( btest(eh_mask(i,j,k,l),4) ) then
    dfz(i,j,k,l) = idz * ( - three * f(i,j,k,l) + four * f(i,j,k+1,l) - f(i,j,k+2,l) )
  else
    dfz(i,j,k,l) = idz * ( three * f(i,j,k,l) - four * f(i,j,k-1,l) + f(i,j,k-2,l) )
  end if

! Calculate the characteristic direction using these preliminary
! derivatives and multiply it with the timestep.

  dfup(1) = g3xx(i,j,k) * dfx(i,j,k,l) + &
            g3xy(i,j,k) * dfy(i,j,k,l) + &
            g3xz(i,j,k) * dfz(i,j,k,l)
  dfup(2) = g3xy(i,j,k) * dfx(i,j,k,l) + &
            g3yy(i,j,k) * dfy(i,j,k,l) + &
            g3yz(i,j,k) * dfz(i,j,k,l)
  dfup(3) = g3xz(i,j,k) * dfx(i,j,k,l) + &
            g3yz(i,j,k) * dfy(i,j,k,l) + &
            g3zz(i,j,k) * dfz(i,j,k,l)

  alp2 = alp(i,j,k)**2

  cfactor = one / sqrt ( alp2 * ( dfup(1) * dfx(i,j,k,l) + &
                                  dfup(2) * dfy(i,j,k,l) + &
                                  dfup(3) * dfz(i,j,k,l) ) )

  cdx(1) = ( -betax(i,j,k) + ssign * alp2 * dfup(1) * cfactor ) * &
                                                      cctk_delta_time
  cdx(2) = ( -betay(i,j,k) + ssign * alp2 * dfup(2) * cfactor ) * &
                                                      cctk_delta_time
  cdx(3) = ( -betaz(i,j,k) + ssign * alp2 * dfup(3) * cfactor ) * &
                                                      cctk_delta_time

  if ( ( .not. btest(eh_mask(i,j,k,l),0) ) .and. &
       ( .not. btest(eh_mask(i,j,k,l),1) ) ) then
    ! Look at the characteristic direction to figure out in which direction
    ! the stencil should be chosen.
    if ( cdx(1) .gt. zero ) then           ! Choose the stencil to the left.

      ! If the neighbour to the left is also an active non-boundary cell its
      ! safe to assume that i-2, i-1,i and i+1 have good values.
      ! Otherwise the neigbour must be a boundary cell, so we can only use a
      ! centered derivative, which we have already computed.
      if ( eh_mask(i-1,j,k,l) .eq. 0 ) then

        ! Choose the stencil, that gives the smalles second derivative.
        ! If it happens to be the centered, do not do anything since we have
        ! already calculated the centered derivative.
!        if ( f(i-2,j,k,l) - two * f(i-1,j,k,l) + f(i,j,k,l) .le. &
!             f(i-1,j,k,l) - two * f(i,j,k,l) + f(i+1,j,k,l) ) then
          dfx(i,j,k,l) = idx * ( three * f(i,j,k,l) - &
                               four * f(i-1,j,k,l) + f(i-2,j,k,l) )
!        end if
      end if

    else if ( cdx(1) .lt. zero ) then      ! Choose the stencil to the right.

      ! If the neighbour to the right is also an active non-boundary cell its
      ! safe to assume that i-1, i,i+1 and i+2 have good values.
      ! Otherwise the neigbour must be a boundary cell, so we can only use a
      ! centered derivative, which we have already computed.
      if ( eh_mask(i+1,j,k,l) .eq. 0 ) then

        ! Choose the stencil, that gives the smalles second derivative.
        ! If it happens to be the centered, do not do anything since we have
        ! already calculated the centered derivative.
!        if ( f(i-1,j,k,l) - two * f(i,j,k,l) + f(i+1,j,k,l) .ge. &
!             f(i,j,k,l) - two * f(i+1,j,k,l) + f(i+2,j,k,l) ) then
          dfx(i,j,k,l) = idx * ( -three * f(i,j,k,l) + &
                               four * f(i+1,j,k,l) - f(i+2,j,k,l) )
!        end if
      end if
    end if
  end if
  
  ! Do the same for the y direction.
  if ( ( .not. btest(eh_mask(i,j,k,l),2) ) .and. &
       ( .not. btest(eh_mask(i,j,k,l),3) ) ) then

    if ( cdx(2) .gt. zero ) then           ! Choose the stencil to the left.
      if ( eh_mask(i,j-1,k,l) .eq. 0 ) then
!        if ( f(i,j-2,k,l) - two * f(i,j-1,k,l) + f(i,j,k,l) .le. &
!             f(i,j-1,k,l) - two * f(i,j,k,l) + f(i,j+1,k,l) ) then
          dfy(i,j,k,l) = idy * ( three * f(i,j,k,l) - &
                               four * f(i,j-1,k,l) + f(i,j-2,k,l) )
!        end if
      end if
    else if ( cdx(2) .lt. zero ) then      ! Choose the stencil to the left.
      if ( eh_mask(i,j+1,k,l) .eq. 0 ) then
!        if ( f(i,j-1,k,l) - two * f(i,j,k,l) + f(i,j+1,k,l) .ge. &
!             f(i,j,k,l) - two * f(i,j+1,k,l) + f(i,j+2,k,l) ) then
          dfy(i,j,k,l) = idy * ( -three * f(i,j,k,l) + &
                               four * f(i,j+1,k,l) - f(i,j+2,k,l) )
!        end if
      end if
    end if
  end if

  ! And for the z direction.
  if ( ( .not. btest(eh_mask(i,j,k,l),4) ) .and. &
       ( .not. btest(eh_mask(i,j,k,l),5) ) ) then

    if ( cdx(3) .gt. zero ) then
      if ( eh_mask(i,j,k-1,l) .eq. 0 ) then
!        if ( f(i,j,k-2,l) - two * f(i,j,k-1,l) + f(i,j,k,l) .le. &
!             f(i,j,k-1,l) - two * f(i,j,k,l) + f(i,j,k+1,l) ) then
          dfz(i,j,k,l) = idz * ( three * f(i,j,k,l) - &
                               four * f(i,j,k-1,l) + f(i,j,k-2,l) )
!        end if
      end if
    else if ( cdx(3) .lt. zero ) then
      if ( eh_mask(i,j,k+1,l) .eq. 0 ) then
!        if ( f(i,j,k-1,l) - two * f(i,j,k,l) + f(i,j,k+1,l) .ge. &
!             f(i,j,k,l) - two * f(i,j,k+1,l) + f(i,j,k+2,l) ) then
          dfz(i,j,k,l) = idz * ( -three * f(i,j,k,l) + &
                               four * f(i,j,k+1,l) - f(i,j,k+2,l) )
!        end if
      end if
    end if
  end if
else
  dfx(i,j,k,l) = zero
  dfy(i,j,k,l) = zero
  dfz(i,j,k,l) = zero
end if
