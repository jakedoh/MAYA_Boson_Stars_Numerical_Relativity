// -*-C++-*-

// This is an OpenCL kernel, without the surrounding function
// declaration. Such a declaration will be added by OpenCLRunTime,
// passing the Cactus cctkGH in a way such that (most) standard Cactus
// grid information (the cctk_* variables) can be accessed. That is,
// the code below will be placed into a function.



// Cactus parameters are evaluated when the OpenCL code is compiled
CCTK_REAL const kx      = 2*M_PI / wavelength;
CCTK_REAL const ky      = 2*M_PI / wavelength;
CCTK_REAL const kz      = 2*M_PI / wavelength;
CCTK_REAL const omega   = sqrt(pown(kx,2) + pown(ky,2) + pown(kz,2));
CCTK_REAL const cos_t   = cos(omega *  cctk_time                     );
CCTK_REAL const cos_t1  = cos(omega * (cctk_time -   cctk_delta_time));
CCTK_REAL const cos_t2  = cos(omega * (cctk_time - 2*cctk_delta_time));

// Grid points are index in the same way as for a CPU
// Using ptrdiff_t instead of int is more efficient on 64-bit
// architectures
ptrdiff_t const di = 1;
ptrdiff_t const dj =
  CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
ptrdiff_t const dk =
  CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);

// Coordinates are calculated in the same as as for a CPU
CCTK_REAL const x0 = CCTK_ORIGIN_SPACE(0);
CCTK_REAL const y0 = CCTK_ORIGIN_SPACE(1);
CCTK_REAL const z0 = CCTK_ORIGIN_SPACE(2);

CCTK_REAL const dx = CCTK_DELTA_SPACE(0);
CCTK_REAL const dy = CCTK_DELTA_SPACE(1);
CCTK_REAL const dz = CCTK_DELTA_SPACE(2);

ptrdiff_t const i0 = cctk_lbnd[0];
ptrdiff_t const j0 = cctk_lbnd[1];
ptrdiff_t const k0 = cctk_lbnd[2];


  
// Note: The kernel below is not vectorised (since it doesn't use
// CCTK_REAL_VEC). Therefore, vectorisation must be switched off in
// the paramter file (via OpenCLRunTime::vector_size_x = 1).



// This loop macro automatically parallelizes the code
// imin[] and imax[] are passed from the host
LC_LOOP3(init,
         i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
         cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
{
  // Calculate index of current point
  ptrdiff_t const ijk = di*i + dj*j + dk*k;
  
  // Calculate coordinate of local grid point
  CCTK_REAL const xval = x0 + dx * (i0 + i);
  CCTK_REAL const yval = y0 + dy * (j0 + j);
  CCTK_REAL const zval = z0 + dz * (k0 + k);
  
  // Set up a standing wave
  CCTK_REAL const cos_x = cos(kx*xval);
  CCTK_REAL const cos_y = cos(ky*yval);
  CCTK_REAL const cos_z = cos(kz*zval);
  
  CCTK_REAL const uval   = cos_x * cos_y * cos_z * cos_t;
  CCTK_REAL const upval  = cos_x * cos_y * cos_z * cos_t1;
  CCTK_REAL const uppval = cos_x * cos_y * cos_z * cos_t2;
  
  // Initialise grid functions
  u    [ijk] = uval;
  u_p  [ijk] = upval;
  u_p_p[ijk] = uppval;
  
} LC_ENDLOOP3(init);
