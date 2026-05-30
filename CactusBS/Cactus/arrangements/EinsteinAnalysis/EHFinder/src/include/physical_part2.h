! Figure out the extent of the physical part of a chunk of a grid.
! This version to be used in EHFinder_Reparametrize_Check, where it does not
! make sense to go all the way to the physical boundary.
! $Header$

! nx, ny and nz contains the full size of the chunk.

nx = cctk_lsh(1)
ny = cctk_lsh(2)
nz = cctk_lsh(3)

! ngx, ngy, ngz contains the size of any ghostzones.

ngx = cctk_nghostzones(1)
ngy = cctk_nghostzones(2)
ngz = cctk_nghostzones(3)

! ixl and ixr is set to the minimum and maximum local cell index in the
! x-direction. The corresponding is done for the y- and z-directions.

symx = .false.; symy = .false.; symz = .false.;

ixl = 1
ixr = nx
jyl = 1
jyr = ny
kzl = 1
kzr = nz

! I check to see if we have any processor boundaries. If we have, adjust the
! minimum and maximum cell index appropriately with the size of the ghostzone.

if ( cctk_bbox(1) .eq. 0 ) ixl = ixl + ngx
if ( cctk_bbox(2) .eq. 0 ) ixr = ixr - ngx

if ( cctk_bbox(3) .eq. 0 ) jyl = jyl + ngy
if ( cctk_bbox(4) .eq. 0 ) jyr = jyr - ngy

if ( cctk_bbox(5) .eq. 0 ) kzl = kzl + ngz
if ( cctk_bbox(6) .eq. 0 ) kzr = kzr - ngz

! In octant I check to see if the chunk of grid contains any symmetry
! boundaries. This is done by checking if the lowest index of the local grid
! as seen on the global grid is equal to zero. If this is the case adjust the
! minimum and maximum cell index appropriately with the size of the ghostzone.

if ( CCTK_EQUALS ( domain, 'octant' ) ) then
  if (cctk_lbnd(1) .eq. 0 ) then
    ixl = ixl + ngx
    symx = .true.
  end if
  if (cctk_lbnd(2) .eq. 0 ) then
    jyl = jyl + ngy
    symy = .true.
  end if
  if (cctk_lbnd(3) .eq. 0 ) then
    kzl = kzl + ngz
    symz = .true.
  end if
end if

! In quadrant I do the same, except here I have to check in which direction
! the quadrant is directed. After that I check if the local grid is at the
! edge of the global grid without being a symmetry boundary.

if ( CCTK_EQUALS ( domain, 'quadrant' ) ) then
  if ( CCTK_EQUALS ( quadrant_direction, 'x' ) ) then
    if (cctk_lbnd(2) .eq. 0 ) then
      jyl = jyl + ngy
      symy = .true.
    end if
    if (cctk_lbnd(3) .eq. 0 ) then
      kzl = kzl + ngz
      symz = .true.
    end if
    if (cctk_lbnd(1) .eq. 0 ) ixl = ixl + 1
  end if
  if ( CCTK_EQUALS ( quadrant_direction, 'y' ) ) then
    if (cctk_lbnd(1) .eq. 0 ) then
      ixl = ixl + ngx
      symx = .true.
    end if
    if (cctk_lbnd(3) .eq. 0 ) then
      kzl = kzl + ngz
      symz = .true.
    end if
    if (cctk_lbnd(2) .eq. 0 ) jyl = jyl + 1
  end if
  if ( CCTK_EQUALS ( quadrant_direction,'z' ) ) then
    if (cctk_lbnd(1) .eq. 0 ) then
      ixl = ixl + ngx
      symx = .true.
    end if
    if (cctk_lbnd(2) .eq. 0 ) then
      jyl = jyl + ngy
      symy = .true.
    end if
    if (cctk_lbnd(3) .eq. 0 ) kzl = kzl + 1
  end if
end if

! Ditto for bitant. Here I have to check which plane is the symmetry-plane.

if ( CCTK_EQUALS ( domain, 'bitant' ) ) then
  if ( CCTK_EQUALS ( bitant_plane, 'xy' ) ) then
    if (cctk_lbnd(3) .eq. 0 ) then
      kzl = kzl + ngz
      symz = .true.
    end if
    if (cctk_lbnd(1) .eq. 0 ) ixl = ixl + 1
    if (cctk_lbnd(2) .eq. 0 ) jyl = jyl + 1
  endif
  if ( CCTK_EQUALS ( bitant_plane, 'xz' ) ) then
    if (cctk_lbnd(2) .eq. 0 ) then
      jyl = jyl + ngy
      symy = .true.
    end if
    if (cctk_lbnd(1) .eq. 0 ) ixl = ixl + 1
    if (cctk_lbnd(3) .eq. 0 ) kzl = kzl + 1
  endif
  if ( CCTK_EQUALS ( bitant_plane, 'yz' ) ) then
    if (cctk_lbnd(1) .eq. 0 ) then
      ixl = ixl + ngx
      symx = .true.
    end if
    if (cctk_lbnd(2) .eq. 0 ) jyl = jyl + 1
    if (cctk_lbnd(3) .eq. 0 ) kzl = kzl + 1
  endif
end if

if ( CCTK_EQUALS ( domain, 'full' ) ) then
  if (cctk_lbnd(1) .eq. 0 ) ixl = ixl + 1
  if (cctk_lbnd(2) .eq. 0 ) jyl = jyl + 1
  if (cctk_lbnd(3) .eq. 0 ) kzl = kzl + 1
end if

if (cctk_ubnd(1)+1 .eq. cctk_gsh(1) ) ixr = ixr - 1
if (cctk_ubnd(2)+1 .eq. cctk_gsh(2) ) jyr = jyr - 1
if (cctk_ubnd(3)+1 .eq. cctk_gsh(3) ) kzr = kzr - 1
