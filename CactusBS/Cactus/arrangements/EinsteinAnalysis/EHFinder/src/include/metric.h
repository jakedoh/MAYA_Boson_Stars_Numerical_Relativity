! Calculation of inverse three metric. At this point there is no distinction
! between the case of a physical or static conformal metric.
! $Header$

if ( eh_mask(i,j,k,l) .ge. 0 ) then
  tmp1 = gyy(i,j,k)*gzz(i,j,k) - gyz(i,j,k)**2
  tmp2 = gxz(i,j,k)*gyz(i,j,k) - gxy(i,j,k)*gzz(i,j,k)
  tmp3 = gxy(i,j,k)*gyz(i,j,k) - gxz(i,j,k)*gyy(i,j,k)

  idetg = one / ( gxx(i,j,k)*tmp1 + gxy(i,j,k)*tmp2 + gxz(i,j,k)*tmp3 )

  g3xx(i,j,k) = tmp1 * idetg
  g3xy(i,j,k) = tmp2 * idetg
  g3xz(i,j,k) = tmp3 * idetg
  g3yy(i,j,k) = ( gxx(i,j,k)*gzz(i,j,k) - gxz(i,j,k)**2 ) * idetg
  g3yz(i,j,k) = ( gxy(i,j,k)*gxz(i,j,k) - gxx(i,j,k)*gyz(i,j,k) ) * idetg
  g3zz(i,j,k) = ( gxx(i,j,k)*gyy(i,j,k) - gxy(i,j,k)**2 ) * idetg
end if
