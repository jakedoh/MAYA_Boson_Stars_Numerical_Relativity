
#include <math.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "SIMD/SIMD_intrinsics.h"
#include "finite_difference_functions.h"

void Baikal_enforce_detgammabar_constraint(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    #include "enforcedetgammabar_constraint.h"
}
