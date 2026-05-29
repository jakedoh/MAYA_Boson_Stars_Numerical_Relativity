#include <cmath>
#include <cstdlib>

#include "LagrangeInterp.hh"

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#define MAXORDER  3
#define SQ(X) ((X)*(X))

void TestLocalInterp2_Symmetric(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL const delta[3] = {
    CCTK_DELTA_SPACE(0), CCTK_DELTA_SPACE(1), CCTK_DELTA_SPACE(2)
  };
  CCTK_REAL const origin[3] = {
    CCTK_ORIGIN_SPACE(0), CCTK_ORIGIN_SPACE(2), CCTK_ORIGIN_SPACE(2)
  };
  CCTK_INT const siz[3] = { cctk_ash[0], cctk_ash[1], cctk_ash[2] };

  *symmetry = true;
  for(int order = 0; order <= MAXORDER; ++order) {
    for(int i = 0; i < cctk_lsh[0]; ++i)
    for(int j = 0; j < cctk_lsh[1]; ++j)
    for(int k = 0; k < cctk_lsh[2]; ++k) {
      int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
      func[ijk] = cos(x[ijk]) * sin(y[ijk]*y[ijk]) * SQ(z[ijk]);
    }

    CCTK_REAL point[2][3] = {
        {M_PI/10.0, 0.5, 0.125},
        {-M_PI/10.0, -0.5, -0.125}
    };

    switch(order) {
#define TYPECASE(O)                                                           \
    case O:                                                                   \
      {                                                                       \
        LagrangeInterpND<O,3> const interp_1(origin, delta, siz, point[0]);   \
        assert(!interp_1.out_of_bounds);                                      \
        LagrangeInterpND<O,3> const interp_2(origin, delta, siz, point[1]);   \
        assert(!interp_2.out_of_bounds);                                      \
        CCTK_REAL const val_1 = interp_1.eval(func);                          \
        CCTK_REAL const val_2 = interp_2.eval(func);                          \
        if(std::abs(val_1 - val_2) > 0) {                                     \
          CCTK_VWarn(warn_level, __LINE__, __FILE__, CCTK_THORNSTRING,        \
              "SYMM: Error ABOVE tolerance for order %d:\n  %.8e > %.8e",     \
              order, std::abs(val_1 - val_2), 0.0);                           \
          *symmetry = false;                                                  \
        }                                                                     \
        else {                                                                \
          CCTK_VInfo(CCTK_THORNSTRING, "SYMM: Error below tolerance for order"\
              " %d:\n  %.8e <= %.8e", order, std::abs(val_1 - val_2), 0.0);   \
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
  }
}
