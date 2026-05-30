
#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "SIMD/SIMD_intrinsics.h" // Contains needed definition of REAL_SIMD_ARRAY
#include "finite_difference_functions.h"

void Baikal_BSSN_constraints(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    const CCTK_REAL invdx0 = 1.0/CCTK_DELTA_SPACE(0);
    const CCTK_REAL invdx1 = 1.0/CCTK_DELTA_SPACE(1);
    const CCTK_REAL invdx2 = 1.0/CCTK_DELTA_SPACE(2);
    if(FD_order == 2) {
        #include "BSSN_constraints_enable_Tmunu_True_FD_order_2.h"
    }
    if(FD_order == 4) {
        #include "BSSN_constraints_enable_Tmunu_True_FD_order_4.h"
    }
}
