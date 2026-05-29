#include <cmath>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "copy_mask.hh"

namespace NoExcision {

#ifdef HAVE_CARPET
  using namespace Carpet;
#endif

  /**
   * Modify the mask according to the CarpetReduce mask in order to be able to do the local reductions.
   */

  extern "C" 
    void CopyMask (CCTK_ARGUMENTS)
    {
      DECLARE_CCTK_ARGUMENTS;
      DECLARE_CCTK_PARAMETERS;

      CCTK_REAL * const weight =
        static_cast <CCTK_REAL *>
        (CCTK_VarDataPtr (cctkGH, 0, "CarpetReduce::weight"));

      if (not weight) {
        CCTK_WARN (CCTK_WARN_ABORT,
                   "CarpetReduce is not active, or CarpetReduce::mask does not have storage");
      }

      for (int k = 0; k < cctk_lsh[2]; ++ k) {
        for (int j = 0; j < cctk_lsh[1]; ++ j) {
          for (int i = 0; i < cctk_lsh[0]; ++ i) {
            int const ind = CCTK_GFINDEX3D (cctkGH, i, j, k);

            red_mask[ind] = weight[ind];
          }
        }
      }

    }


} // namespace NoExcision
