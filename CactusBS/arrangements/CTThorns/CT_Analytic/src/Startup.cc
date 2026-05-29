/*  File produced by Kranc */

#include "cctk.h"

extern "C" int CT_Analytic_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "CT_Analytic";
  CCTK_RegisterBanner(banner);
  return 0;
}
