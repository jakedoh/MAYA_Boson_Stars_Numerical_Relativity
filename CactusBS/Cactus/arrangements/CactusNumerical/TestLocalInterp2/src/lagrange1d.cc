#include <cmath>
#include <cstdlib>

#include "LagrangeInterp.hh"

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#define MAXORDER  3

namespace {

CCTK_REAL * make_rand_poly(int order) {
  CCTK_REAL * coeff = new CCTK_REAL[order+1];
  for(int i = 0; i <= order; ++i) {
    CCTK_REAL ran = rand()/(CCTK_REAL)RAND_MAX;
    coeff[i] = 2*ran - 1;
  }
  return coeff;
}

void free_rand_poly(CCTK_REAL ** coeff) {
  delete[] *coeff;
}

CCTK_REAL eval_poly(CCTK_REAL const * const coeff,
    int order, CCTK_REAL const point) {
  CCTK_REAL out = 0;
  CCTK_REAL xp = 1.0;
  for(int i = 0; i <= order; ++i) {
    out += coeff[i] * xp;
    xp = xp * point;
  }
  return out;
}

}

void TestLocalInterp2_Lagrange1D(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL const delta = CCTK_DELTA_SPACE(0);
  CCTK_REAL const origin = CCTK_ORIGIN_SPACE(0);
  CCTK_INT const dim = cctk_lsh[0];
  CCTK_INT const stride = CCTK_GFINDEX3D(cctkGH, 1, 0, 0) -
                          CCTK_GFINDEX3D(cctkGH, 0, 0, 0);

  *lagrange_1d = true;
  for(int order = 0; order <= MAXORDER; ++order) {
    CCTK_REAL * pcoeff = make_rand_poly(order);

    for(int i = 0; i < cctk_lsh[0]; ++i)
    for(int j = 0; j < cctk_lsh[1]; ++j)
    for(int k = 0; k < cctk_lsh[2]; ++k) {
      int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
      poly[ijk] = eval_poly(pcoeff, order, x[ijk]);
    }

    CCTK_REAL point = M_PI/6.0*dim*delta + origin;
    CCTK_REAL val = eval_poly(pcoeff, order, point);

    switch(order) {
#define TYPECASE(O)                                                           \
    case O:                                                                   \
      {                                                                       \
        LagrangeInterp1D<O> const interp(origin, delta, dim, point);          \
        assert(!interp.out_of_bounds);                                        \
        CCTK_INT const offset = CCTK_GFINDEX3D(cctkGH, interp.point, 0, 0);   \
        CCTK_REAL const ival = interp.eval(&poly[offset], stride);            \
        if(std::abs(ival - val) > tolerance) {                                \
          CCTK_VWarn(warn_level, __LINE__, __FILE__, CCTK_THORNSTRING,        \
              "1D: Error ABOVE tolerance for order %d:\n  %.8e > %.8e", order,\
              std::abs(ival - val), tolerance);                               \
          *lagrange_1d = false;                                               \
        }                                                                     \
        else {                                                                \
          CCTK_VInfo(CCTK_THORNSTRING, "1D: Error below tolerance for order " \
              "%d:\n  %.8e <= %.8e", order, std::abs(ival - val), tolerance); \
        }                                                                     \
      }                                                                       \
      break

    TYPECASE(0);
    TYPECASE(1);
    TYPECASE(2);
    TYPECASE(3);
#undef TYPECASE
      default:
        assert(false);
    }

    free_rand_poly(&pcoeff);
  }
}
