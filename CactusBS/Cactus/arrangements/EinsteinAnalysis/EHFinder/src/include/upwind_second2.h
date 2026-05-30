! Calculate second order upwinded differences.
! The only difference from upwind_second.h is that here the loop over cells
! must be done in the including routine.
! $Header$

! If this is an active cell...
if ( eh_mask(i,j,k,l) .ge. 0 ) then

  ! If this is not a boundary cell...
  if ( ( .not. btest(eh_mask(i,j,k,l),0) ) .and. &
       ( .not. btest(eh_mask(i,j,k,l),1) ) ) then

    ! If the cell to the left is not a boundary cell...
    if ( ( .not. btest(eh_mask(i-1,j,k,l),0) ) .or. &
               ( eh_mask(i-1,j,k,l) .eq. -2 ) ) then

      ! Calculate left handed one sided derivative
      al = idx * ( three * f(i,j,k,l) - four * f(i-1,j,k,l) + f(i-2,j,k,l) )
    else
      al = zero
    end if

    ! If the cell to the right is not a boundary cell...
    if ( ( .not. btest(eh_mask(i+1,j,k,l),1) ) .or. &
               ( eh_mask(i+1,j,k,l) .eq. -2 ) ) then 

      ! Calculate right handed one sided derivative
      ar = idx * ( - three * f(i,j,k,l) + four * f(i+1,j,k,l) - f(i+2,j,k,l) )
    else
      ar = zero
    end if

  ! Else if the cell is a left boundary cell.
  else if ( btest(eh_mask(i,j,k,l),0) ) then

    ! calculate only right handed one sided difference.
    al = zero
    ar = idx * ( - three * f(i,j,k,l) + four * f(i+1,j,k,l) - f(i+2,j,k,l) )

  ! Else it must be a right boundary cell.
  else

    ! So calculate only the left handed one sided difference.
    al = idx * ( three * f(i,j,k,l) - four * f(i-1,j,k,l) + f(i-2,j,k,l) )
    ar = zero
  end if

  ! Do the same for the y-derivative.
  if ( ( .not. btest(eh_mask(i,j,k,l),2) ) .and. &
       ( .not. btest(eh_mask(i,j,k,l),3) ) ) then
    if ( ( .not. btest(eh_mask(i,j-1,k,l),2) ) .or. &
               ( eh_mask(i,j-1,k,l) .eq. -2 ) ) then
      bl = idy * ( three * f(i,j,k,l) - four * f(i,j-1,k,l) + f(i,j-2,k,l) )
    else
      bl = zero
    end if
    if ( ( .not. btest(eh_mask(i,j+1,k,l),3) ) .or. &
               ( eh_mask(i,j+1,k,l) .eq. -2 ) ) then 
      br = idy * ( - three * f(i,j,k,l) + four * f(i,j+1,k,l) - f(i,j+2,k,l) )
    else
      br = zero
    end if
  else if ( btest(eh_mask(i,j,k,l),2) ) then
    bl = zero
    br = idy * ( - three * f(i,j,k,l) + four * f(i,j+1,k,l) - f(i,j+2,k,l) )
  else
    bl = idy * ( three * f(i,j,k,l) - four * f(i,j-1,k,l) + f(i,j-2,k,l) )
    br = zero
  end if

  ! And the z-derivative
  if ( ( .not. btest(eh_mask(i,j,k,l),4) ) .and. &
       ( .not. btest(eh_mask(i,j,k,l),5) ) ) then
    if ( ( .not. btest(eh_mask(i,j,k-1,l),4) ) .or. &
               ( eh_mask(i,j,k-1,l) .eq. -2 ) ) then
      cl = idz * ( three * f(i,j,k,l) - four * f(i,j,k-1,l) + f(i,j,k-2,l) )
    else
      cl = zero
    end if
    if ( ( .not. btest(eh_mask(i,j,k+1,l),5) ) .or. &
               ( eh_mask(i,j,k+1,l) .eq. -2 ) ) then 
      cr = idz * ( - three * f(i,j,k,l) + four * f(i,j,k+1,l) - f(i,j,k+2,l) )
    else
      cr = zero
    end if
  else if ( btest(eh_mask(i,j,k,l),4) ) then
    cl = zero
    cr = idz * ( - three * f(i,j,k,l) + four * f(i,j,k+1,l) - f(i,j,k+2,l) )
  else
    cl = idz * ( three * f(i,j,k,l) - four * f(i,j,k-1,l) + f(i,j,k-2,l) )
    cr = zero
  end if

  ! Get the negative and positive parts of the one sided derivatives in
  ! the x-direction.
  alminus = half*(al-abs(al))
  alplus = half*(al+abs(al))
  arminus = half*(ar-abs(ar))
  arplus = half*(ar+abs(ar))

  ! Ditto for the y-direction.
  blminus = half*(bl-abs(bl))
  blplus = half*(bl+abs(bl))
  brminus = half*(br-abs(br))
  brplus = half*(br+abs(br))

  ! Ditto for the z-direction.
  clminus = half*(cl-abs(cl))
  clplus = half*(cl+abs(cl))
  crminus = half*(cr-abs(cr))
  crplus = half*(cr+abs(cr))

  ! If f>0 choose the correct one sided derivative depending on the
  ! magnitude of the negative and positive parts
  if ( f(i,j,k,l) .gt. 0 ) then
    if ( abs(alplus) .gt. abs(arminus) ) then
      dfx(i,j,k,l) = alplus
    else
      dfx(i,j,k,l) = arminus
    end if
    if ( abs(blplus) .gt. abs(brminus) ) then
      dfy(i,j,k,l) = blplus
    else
      dfy(i,j,k,l) = brminus
    end if
    if ( abs(clplus) .gt. abs(crminus) ) then
      dfz(i,j,k,l) = clplus
    else
      dfz(i,j,k,l) = crminus
    end if
  endif

  ! Ditto if f<0.
  if ( f(i,j,k,l) .lt. 0 ) then
    if ( abs(alminus) .gt. abs(arplus) ) then
      dfx(i,j,k,l) = alminus
    else
      dfx(i,j,k,l) = arplus
    end if
    if ( abs(blminus) .gt. abs(brplus) ) then
      dfy(i,j,k,l) = blminus
    else
      dfy(i,j,k,l) = brplus
    end if
    if ( abs(clminus) .gt. abs(crplus) ) then
      dfz(i,j,k,l) = clminus
    else
      dfz(i,j,k,l) = crplus
    end if
  end if

! If the cell is not active set the derivatives to zero.
else
  dfx(i,j,k,l) = zero
  dfy(i,j,k,l) = zero
  dfz(i,j,k,l) = zero
end if
