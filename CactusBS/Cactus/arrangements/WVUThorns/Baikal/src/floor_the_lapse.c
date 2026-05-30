
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

#ifndef MAX
#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )

void Baikal_floor_the_lapse(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
#pragma omp parallel for
for (int i2 = 0; i2 < cctk_lsh[2]; i2++) {
    for (int i1 = 0; i1 < cctk_lsh[1]; i1++) {
        for (int i0 = 0; i0 < cctk_lsh[0]; i0++) {
            
            alphaGF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)] = MAX(alphaGF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)], lapse_floor);
            
        } // END LOOP: for (int i0 = 0; i0 < cctk_lsh[0]; i0++)
    } // END LOOP: for (int i1 = 0; i1 < cctk_lsh[1]; i1++)
} // END LOOP: for (int i2 = 0; i2 < cctk_lsh[2]; i2++)

}
#undef MAX
#endif
