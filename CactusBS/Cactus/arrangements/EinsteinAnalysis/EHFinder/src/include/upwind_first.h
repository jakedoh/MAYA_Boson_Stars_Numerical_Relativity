! Calculate first order upwinded differences.
! $Header$

if ( eh_mask(i,j,k,l) .ge. 0 ) then
  if ( ( .not. btest(eh_mask(i,j,k,l),0) ) .and. &
       ( .not. btest(eh_mask(i,j,k,l),1) ) ) then
    al = idx * ( f(i,j,k,l) - f(i-1,j,k,l) )
    ar = idx * ( f(i+1,j,k,l) - f(i,j,k,l) )
  else if ( btest(eh_mask(i,j,k,l),0) ) then
    al = zero
    ar = idx * ( f(i+1,j,k,l) - f(i,j,k,l) )
  else
    al = idx * ( f(i,j,k,l) - f(i-1,j,k,l) )
    ar = zero
  end if
  if ( ( .not. btest(eh_mask(i,j,k,l),2) ) .and. &
       ( .not. btest(eh_mask(i,j,k,l),3) ) ) then
    bl = idy * ( f(i,j,k,l) - f(i,j-1,k,l) )
    br = idy * ( f(i,j+1,k,l) - f(i,j,k,l) )
  else if ( btest(eh_mask(i,j,k,l),2) ) then
    bl = zero
    br = idy * ( f(i,j+1,k,l) - f(i,j,k,l) )
  else
    bl = idy * ( f(i,j,k,l) - f(i,j-1,k,l) )
    br = zero
  end if
  if ( ( .not. btest(eh_mask(i,j,k,l),4) ) .and. &
       ( .not. btest(eh_mask(i,j,k,l),5) ) ) then
    cl = idz * ( f(i,j,k,l) - f(i,j,k-1,l) )
    cr = idz * ( f(i,j,k+1,l) - f(i,j,k,l) )
  else if ( btest(eh_mask(i,j,k,l),4) ) then
    cl = zero
    cr = idz * ( f(i,j,k+1,l) - f(i,j,k,l) )
  else
    cl = idz * ( f(i,j,k,l) - f(i,j,k-1,l) )
    cr = zero
  end if

  alminus = half*(al-abs(al))
  alplus = half*(al+abs(al))
  arminus = half*(ar-abs(ar))
  arplus = half*(ar+abs(ar))

  blminus = half*(bl-abs(bl))
  blplus = half*(bl+abs(bl))
  brminus = half*(br-abs(br))
  brplus = half*(br+abs(br))

  clminus = half*(cl-abs(cl))
  clplus = half*(cl+abs(cl))
  crminus = half*(cr-abs(cr))
  crplus = half*(cr+abs(cr))

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
else
  dfx(i,j,k,l) = zero
  dfy(i,j,k,l) = zero
  dfz(i,j,k,l) = zero
end if
