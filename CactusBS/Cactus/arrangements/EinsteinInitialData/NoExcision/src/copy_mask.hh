#include <cctk.h>
#include <cctk_Arguments.h>

namespace CarpetMask {

  extern "C" {
    void
    CopyMask (CCTK_ARGUMENTS);
  }

} // namespace CarpetMask
