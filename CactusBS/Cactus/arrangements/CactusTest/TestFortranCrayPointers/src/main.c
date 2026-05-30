#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include <cctk.h>
#include <cctk_Arguments.h>

CCTK_FCALL void
CCTK_FNAME(TestFortranCrayPointers_sub)(uintptr_t const* restrict pointers,
                                        int const* restrict n);

void TestFortranCrayPointers_main(CCTK_ARGUMENTS)
{
  /* Set the array size */
  int const n = 10;
  
  /* Allocate memory */
  CCTK_REAL* restrict const a = malloc(n * n * sizeof *a);
  CCTK_REAL* restrict const b = malloc(n * n * sizeof *b);
  CCTK_REAL* restrict const c = malloc(n * n * sizeof *c);
  
  /* Initialise the arrays */
  for (ptrdiff_t i=0; i<n; ++i) {
    for (ptrdiff_t j=0; j<n; ++j) {
      b[i+n*j] = 1.0;
      c[i+n*j] = 1.0;
    }
  }
  
  /* Put the arrays into a pointer array */
  uintptr_t const pointers[3] = { (uintptr_t)a, (uintptr_t)b, (uintptr_t)c, };
  
  /* Call the Fortran subroutine */
  CCTK_FNAME(TestFortranCrayPointers_sub)(pointers, &n);
  
  /* Print the result */
  for (ptrdiff_t i=0; i<n; ++i) {
    printf("%2d:", (int)i);
    for (ptrdiff_t j=0; j<n; ++j) {
      printf(" %6.3f", a[i+n*j]);
    }
    printf("\n");
  }
}
