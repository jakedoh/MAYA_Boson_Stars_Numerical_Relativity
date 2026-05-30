
#include <math.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "SIMD/SIMD_intrinsics.h"
#include "finite_difference_functions.h"
extern void Baikal_BSSN_RHSs_enable_Tmunu_True_FD_order_2(CCTK_ARGUMENTS);
extern void Baikal_BSSN_RHSs_enable_Tmunu_True_FD_order_4(CCTK_ARGUMENTS);

void Baikal_driver_pt2_BSSN_RHSs(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;

    const CCTK_INT *FD_order = CCTK_ParameterGet("FD_order","Baikal",NULL);

    if(*FD_order == 2) {
        Baikal_BSSN_RHSs_enable_Tmunu_True_FD_order_2(CCTK_PASS_CTOC);
    }
    if(*FD_order == 4) {
        Baikal_BSSN_RHSs_enable_Tmunu_True_FD_order_4(CCTK_PASS_CTOC);
    }
} // END FUNCTION
