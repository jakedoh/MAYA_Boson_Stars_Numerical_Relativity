/*  File produced by Kranc */

#include "cctk.h"

extern "C" int ScalarTensor_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "ScalarTensor";
  CCTK_RegisterBanner(banner);
  return 0;
}
