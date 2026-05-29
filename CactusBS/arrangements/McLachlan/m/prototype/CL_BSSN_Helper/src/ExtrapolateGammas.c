#include <cctk.h>
#include <cctk_Arguments.h>

static void
extrap (cGH const * restrict cctkGH,
        CCTK_REAL * restrict var);

void
CL_BSSN_ExtrapolateGammas (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  
  extrap (cctkGH, dphi1);
  extrap (cctkGH, dphi2);
  extrap (cctkGH, dphi3);
  
  extrap (cctkGH, dgt111);
  extrap (cctkGH, dgt112);
  extrap (cctkGH, dgt113);
  extrap (cctkGH, dgt122);
  extrap (cctkGH, dgt123);
  extrap (cctkGH, dgt133);
  extrap (cctkGH, dgt211);
  extrap (cctkGH, dgt212);
  extrap (cctkGH, dgt213);
  extrap (cctkGH, dgt222);
  extrap (cctkGH, dgt223);
  extrap (cctkGH, dgt233);
  extrap (cctkGH, dgt311);
  extrap (cctkGH, dgt312);
  extrap (cctkGH, dgt313);
  extrap (cctkGH, dgt322);
  extrap (cctkGH, dgt323);
  extrap (cctkGH, dgt333);
  
  extrap (cctkGH, Xt1);
  extrap (cctkGH, Xt2);
  extrap (cctkGH, Xt3);
  
  extrap (cctkGH, dalpha1);
  extrap (cctkGH, dalpha2);
  extrap (cctkGH, dalpha3);
  
  extrap (cctkGH, dbeta11);
  extrap (cctkGH, dbeta21);
  extrap (cctkGH, dbeta31);
  extrap (cctkGH, dbeta12);
  extrap (cctkGH, dbeta22);
  extrap (cctkGH, dbeta32);
  extrap (cctkGH, dbeta13);
  extrap (cctkGH, dbeta23);
  extrap (cctkGH, dbeta33);
  
  extrap (cctkGH, B1);
  extrap (cctkGH, B2);
  extrap (cctkGH, B3);
}

static void
extrap (cGH const * restrict const cctkGH,
        CCTK_REAL * restrict const var)
{
  ExtrapolateGammas (cctkGH, var);
}
