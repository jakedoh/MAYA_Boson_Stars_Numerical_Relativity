! Calculation of shift upwinded differences except at boundaries.
! $Header$

! If this is an active cell...
if ( eh_mask(i,j,k,l) .ge. 0 ) then

  ! If the shift in the x-direction is <0.
  if ( betax(i,j,k) .lt. zero ) then

    ! If both of the two nearest cells to the right are active cells...
    if ( ( eh_mask(i+1,j,k,l) .ge. 0 ) .and. ( eh_mask(i+2,j,k,l) .ge. 0 ) ) then

      ! Calculate the right handed one sided derivative.
      dfx(i,j,k,l) = idx * ( -three * f(i,j,k,l) + &
                           four * f(i+1,j,k,l) - f(i+2,j,k,l) )

    ! Else if only the nearest neighbour to the right is active...
    else if ( eh_mask(i+1,j,k,l) .ge. 0 ) then

      ! Calculate the centered derivative.
      dfx(i,j,k,l) = idx * ( f(i+1,j,k,l) - f(i-1,j,k,l) )

    ! Else it must be a boundary cell...
    else

      ! So calculate the left handed one sided derivative.
      dfx(i,j,k,l) = idx * ( three * f(i,j,k,l) - &
                           four * f(i-1,j,k,l) + f(i-2,j,k,l) )
    end if

  ! Else if the shift is >= 0.
  else

    ! If both of the two nearest cells to the left are active cells...
    if ( ( eh_mask(i-1,j,k,l) .ge. 0 ) .and. ( eh_mask(i-2,j,k,l) .ge. 0 ) ) then

      ! Calculate the left handed one sided derivative.
      dfx(i,j,k,l) = idx * ( three * f(i,j,k,l) - &
                           four * f(i-1,j,k,l) + f(i-2,j,k,l) )

    ! Else if only the nearest neighbour to the left is active...
    else if ( eh_mask(i-1,j,k,l) .ge. 0 ) then

      ! Calculate the centered derivative.
      dfx(i,j,k,l) = idx * ( f(i+1,j,k,l) - f(i-1,j,k,l) )

    ! Else it must be a boundary cell...
    else

      ! So calculate the left handed one sided derivative.
      dfx(i,j,k,l) = idx * ( -three * f(i,j,k,l) + &
                           four * f(i+1,j,k,l) - f(i+2,j,k,l) )
    end if
  end if

  ! Ditto for the y-derivative.
  if ( betay(i,j,k) .lt. zero ) then
    if ( ( eh_mask(i,j+1,k,l) .ge. 0 ) .and. ( eh_mask(i,j+2,k,l) .ge. 0 ) ) then
      dfy(i,j,k,l) = idy * ( -three * f(i,j,k,l) + &
                           four * f(i,j+1,k,l) - f(i,j+2,k,l) )
    else if ( eh_mask(i,j+1,k,l) .ge. 0 ) then
      dfy(i,j,k,l) = idy * ( f(i,j+1,k,l) - f(i,j-1,k,l) )
    else
      dfy(i,j,k,l) = idy * ( three * f(i,j,k,l) - &
                           four * f(i,j-1,k,l) + f(i,j-2,k,l) )
    end if
  else
    if ( ( eh_mask(i,j-1,k,l) .ge. 0 ) .and. ( eh_mask(i,j-2,k,l) .ge. 0 ) ) then
      dfy(i,j,k,l) = idy * ( three * f(i,j,k,l) - &
                           four * f(i,j-1,k,l) + f(i,j-2,k,l) )
    else if ( eh_mask(i-1,j,k,l) .ge. 0 ) then
      dfy(i,j,k,l) = idy * ( f(i,j+1,k,l) - f(i,j-1,k,l) )
    else
      dfy(i,j,k,l) = idy * ( -three * f(i,j,k,l) + &
                           four * f(i,j+1,k,l) - f(i,j+2,k,l) )
    end if
  end if

  ! Ditto for the z-derivative.
  if ( betaz(i,j,k) .lt. zero ) then
    if ( ( eh_mask(i,j,k+1,l) .ge. 0 ) .and. ( eh_mask(i,j,k+2,l) .ge. 0 ) ) then
      dfz(i,j,k,l) = idz * ( -three * f(i,j,k,l) + &
                           four * f(i,j,k+1,l) - f(i,j,k+2,l) )
    else if ( eh_mask(i,j,k+1,l) .ge. 0 ) then
      dfz(i,j,k,l) = idz * ( f(i,j,k+1,l) - f(i,j,k-1,l) )
    else
      dfz(i,j,k,l) = idz * ( three * f(i,j,k,l) - &
                           four * f(i,j,k-1,l) + f(i,j,k-2,l) )
    end if
  else
	    if ( ( eh_mask(i,j,k-1,l) .ge. 0 ) .and. ( eh_mask(i,j,k-2,l) .ge. 0 ) ) then
      dfz(i,j,k,l) = idz * ( three * f(i,j,k,l) - &
                           four * f(i,j,k-1,l) + f(i,j,k-2,l) )
    else if ( eh_mask(i-1,j,k,l) .ge. 0 ) then
      dfz(i,j,k,l) = idz * ( f(i,j,k+1,l) - f(i,j,k-1,l) )
    else
      dfz(i,j,k,l) = idz * ( -three * f(i,j,k,l) + &
                           four * f(i,j,k+1,l) - f(i,j,k+2,l) )
    end if
  end if

! If the cell is not active set the derivatives to zero.
else
  dfx(i,j,k,l) = zero
  dfy(i,j,k,l) = zero
  dfz(i,j,k,l) = zero
end if
