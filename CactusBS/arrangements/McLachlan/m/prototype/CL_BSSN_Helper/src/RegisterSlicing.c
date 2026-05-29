#include <cctk.h>
#include <Slicing.h>

int
CL_BSSN_RegisterSlicing (void)
{
  Einstein_RegisterSlicing ("CL_BSSN");
  return 0;
}
