/*  File produced by Kranc */

#include "cctk.h"

extern "C" int ExternalEtaBeta_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "ExternalEtaBeta";
  CCTK_RegisterBanner(banner);
  return 0;
}
