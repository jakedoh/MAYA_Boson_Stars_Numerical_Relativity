
#include <math.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "SIMD/SIMD_intrinsics.h"
#include "finite_difference_functions.h"
extern void BaikalVacuum_BSSN_Ricci_FD_order_4(CCTK_ARGUMENTS);
extern void BaikalVacuum_BSSN_Ricci_FD_order_6(CCTK_ARGUMENTS);
extern void BaikalVacuum_BSSN_Ricci_FD_order_8(CCTK_ARGUMENTS);

void BaikalVacuum_driver_pt1_BSSN_Ricci(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    const CCTK_INT *FD_order = CCTK_ParameterGet("FD_order","BaikalVacuum",NULL);
    if(*FD_order == 4) {
        BaikalVacuum_BSSN_Ricci_FD_order_4(CCTK_PASS_CTOC);
    }
    if(*FD_order == 6) {
        BaikalVacuum_BSSN_Ricci_FD_order_6(CCTK_PASS_CTOC);
    }
    if(*FD_order == 8) {
        BaikalVacuum_BSSN_Ricci_FD_order_8(CCTK_PASS_CTOC);
    }
} // END FUNCTION
