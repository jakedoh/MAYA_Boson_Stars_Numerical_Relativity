// -*-C++-*-

// This is an OpenCL kernel, without the surrounding function
// declaration. Such a declaration will be added by OpenCLRunTime,
// passing the Cactus cctkGH in a way such that (most) standard Cactus
// grid information (the cctk_* variables) can be accessed. That is,
// the code below will be placed into a function.



// Grid points are index in the same way as for a CPU
// Using ptrdiff_t instead of int is more efficient on 64-bit
// architectures
ptrdiff_t const di = 1;
ptrdiff_t const dj =
  CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
ptrdiff_t const dk =
  CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);


  
// Note: The kernel below is not vectorised (since it doesn't use
// CCTK_REAL_VEC). Therefore, vectorisation must be switched off in
// the paramter file (via OpenCLRunTime::vector_size_x = 1).



// This loop macro automatically parallelizes the code
LC_LOOP3(boundary,
         i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
         cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
{
  // Calculate index of current point
  ptrdiff_t const ijk = di*i + dj*j + dk*k;
  
  CCTK_REAL const uval = 0.0;
  
  u[ijk] = uval;
} LC_ENDLOOP3(boundary);
