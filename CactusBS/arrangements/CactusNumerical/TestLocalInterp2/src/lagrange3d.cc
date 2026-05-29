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

CCTK_REAL eval_tp_poly(CCTK_REAL const * const coeff[3],
    int order, CCTK_REAL const point[3]) {
  return
    eval_poly(coeff[0], order, point[0]) *
    eval_poly(coeff[1], order, point[1]) *
    eval_poly(coeff[2], order, point[2]);
}

}

void TestLocalInterp2_Lagrange3D(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL const delta[3] = {
    CCTK_DELTA_SPACE(0), CCTK_DELTA_SPACE(1), CCTK_DELTA_SPACE(2)
  };
  CCTK_REAL const origin[3] = {
    CCTK_ORIGIN_SPACE(0), CCTK_ORIGIN_SPACE(2), CCTK_ORIGIN_SPACE(2)
  };
  CCTK_INT const siz[3] = { cctk_ash[0], cctk_ash[1], cctk_ash[2] };

  *lagrange_3d = true;
  for(int order = 0; order <= MAXORDER; ++order) {
    CCTK_REAL * pcoeff[3] = {
      make_rand_poly(order),
      make_rand_poly(order),
      make_rand_poly(order)
    };

    for(int i = 0; i < cctk_lsh[0]; ++i)
    for(int j = 0; j < cctk_lsh[1]; ++j)
    for(int k = 0; k < cctk_lsh[2]; ++k) {
      int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
      CCTK_REAL const point[3] = {x[ijk], y[ijk], z[ijk]};
      poly[ijk] = eval_tp_poly(pcoeff, order, point);
    }

    CCTK_REAL point[3] = {
      M_PI/6.0*siz[0]*delta[0] + origin[0],
      0.5*siz[1]*delta[1] + origin[1],
      0.5*siz[2]*delta[2] + origin[2] + 0.5*delta[2]
    };
    CCTK_REAL const val = eval_tp_poly(pcoeff, order, point);

    switch(order) {
#define TYPECASE(O)                                                           \
    case O:                                                                   \
      {                                                                       \
        LagrangeInterpND<O,3> const interp(origin, delta, siz, point);        \
        assert(!interp.out_of_bounds);                                        \
        CCTK_REAL const ival = interp.eval(poly);                             \
        if(std::abs(ival - val) > tolerance) {                                \
          CCTK_VWarn(warn_level, __LINE__, __FILE__, CCTK_THORNSTRING,        \
              "3D: Error ABOVE tolerance for order %d:\n  %.8e > %.8e", order,\
              std::abs(ival - val), tolerance);                               \
          *lagrange_3d = false;                                               \
        }                                                                     \
        else {                                                                \
          CCTK_VInfo(CCTK_THORNSTRING, "3D: Error below tolerance for order " \
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

    free_rand_poly(&pcoeff[0]);
    free_rand_poly(&pcoeff[1]);
    free_rand_poly(&pcoeff[2]);
  }
}
