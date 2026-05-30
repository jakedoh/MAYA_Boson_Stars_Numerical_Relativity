
#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void Baikal_NewRad(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  NewRad_Apply(cctkGH, aDD00GF, aDD00_rhsGF, 0.0, 1.0, 2.0);
  NewRad_Apply(cctkGH, aDD01GF, aDD01_rhsGF, 0.0, 1.0, 2.0);
  NewRad_Apply(cctkGH, aDD02GF, aDD02_rhsGF, 0.0, 1.0, 2.0);
  NewRad_Apply(cctkGH, aDD11GF, aDD11_rhsGF, 0.0, 1.0, 2.0);
  NewRad_Apply(cctkGH, aDD12GF, aDD12_rhsGF, 0.0, 1.0, 2.0);
  NewRad_Apply(cctkGH, aDD22GF, aDD22_rhsGF, 0.0, 1.0, 2.0);
  NewRad_Apply(cctkGH, alphaGF, alpha_rhsGF, 1.0, sqrt(2.0), 1.0);
  NewRad_Apply(cctkGH, betU0GF, betU0_rhsGF, 0.0, 1.0, 1.0);
  NewRad_Apply(cctkGH, betU1GF, betU1_rhsGF, 0.0, 1.0, 1.0);
  NewRad_Apply(cctkGH, betU2GF, betU2_rhsGF, 0.0, 1.0, 1.0);
  NewRad_Apply(cctkGH, cfGF, cf_rhsGF, 0.0, 1.0, 1.0);
  NewRad_Apply(cctkGH, hDD00GF, hDD00_rhsGF, 0.0, 1.0, 1.0);
  NewRad_Apply(cctkGH, hDD01GF, hDD01_rhsGF, 0.0, 1.0, 1.0);
  NewRad_Apply(cctkGH, hDD02GF, hDD02_rhsGF, 0.0, 1.0, 1.0);
  NewRad_Apply(cctkGH, hDD11GF, hDD11_rhsGF, 0.0, 1.0, 1.0);
  NewRad_Apply(cctkGH, hDD12GF, hDD12_rhsGF, 0.0, 1.0, 1.0);
  NewRad_Apply(cctkGH, hDD22GF, hDD22_rhsGF, 0.0, 1.0, 1.0);
  NewRad_Apply(cctkGH, lambdaU0GF, lambdaU0_rhsGF, 0.0, 1.0, 1.0);
  NewRad_Apply(cctkGH, lambdaU1GF, lambdaU1_rhsGF, 0.0, 1.0, 1.0);
  NewRad_Apply(cctkGH, lambdaU2GF, lambdaU2_rhsGF, 0.0, 1.0, 1.0);
  NewRad_Apply(cctkGH, trKGF, trK_rhsGF, 0.0, 1.0, 2.0);
  NewRad_Apply(cctkGH, vetU0GF, vetU0_rhsGF, 0.0, 1.0, 1.0);
  NewRad_Apply(cctkGH, vetU1GF, vetU1_rhsGF, 0.0, 1.0, 1.0);
  NewRad_Apply(cctkGH, vetU2GF, vetU2_rhsGF, 0.0, 1.0, 1.0);
}
