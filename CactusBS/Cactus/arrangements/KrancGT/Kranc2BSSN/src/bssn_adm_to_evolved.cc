/*  File produced by Kranc */

#define KRANC_C

#include <algorithm>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Kranc.hh"
#include "Differencing.h"
#include "loopcontrol.h"

namespace Kranc2BSSN {


static void bssn_adm_to_evolved_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  /* Include user-supplied include files */
  /* Initialise finite differencing variables */
  const ptrdiff_t di CCTK_ATTRIBUTE_UNUSED = 1;
  const ptrdiff_t dj CCTK_ATTRIBUTE_UNUSED = 
    CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  const ptrdiff_t dk CCTK_ATTRIBUTE_UNUSED = 
    CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * di;
  const ptrdiff_t cdj CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dj;
  const ptrdiff_t cdk CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dk;
  const ptrdiff_t cctkLbnd1 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[0];
  const ptrdiff_t cctkLbnd2 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[1];
  const ptrdiff_t cctkLbnd3 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[2];
  const CCTK_REAL t CCTK_ATTRIBUTE_UNUSED = cctk_time;
  const CCTK_REAL cctkOriginSpace1 CCTK_ATTRIBUTE_UNUSED = 
    CCTK_ORIGIN_SPACE(0);
  const CCTK_REAL cctkOriginSpace2 CCTK_ATTRIBUTE_UNUSED = 
    CCTK_ORIGIN_SPACE(1);
  const CCTK_REAL cctkOriginSpace3 CCTK_ATTRIBUTE_UNUSED = 
    CCTK_ORIGIN_SPACE(2);
  const CCTK_REAL dt CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_TIME;
  const CCTK_REAL dx CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_SPACE(2);
  const CCTK_REAL dxi CCTK_ATTRIBUTE_UNUSED = pow(dx,-1);
  const CCTK_REAL dyi CCTK_ATTRIBUTE_UNUSED = pow(dy,-1);
  const CCTK_REAL dzi CCTK_ATTRIBUTE_UNUSED = pow(dz,-1);
  const CCTK_REAL khalf CCTK_ATTRIBUTE_UNUSED = 0.5;
  const CCTK_REAL kthird CCTK_ATTRIBUTE_UNUSED = 
    0.333333333333333333333333333333;
  const CCTK_REAL ktwothird CCTK_ATTRIBUTE_UNUSED = 
    0.666666666666666666666666666667;
  const CCTK_REAL kfourthird CCTK_ATTRIBUTE_UNUSED = 
    1.33333333333333333333333333333;
  const CCTK_REAL hdxi CCTK_ATTRIBUTE_UNUSED = 0.5*dxi;
  const CCTK_REAL hdyi CCTK_ATTRIBUTE_UNUSED = 0.5*dyi;
  const CCTK_REAL hdzi CCTK_ATTRIBUTE_UNUSED = 0.5*dzi;
  /* Initialize predefined quantities */
  const CCTK_REAL p1o12dx CCTK_ATTRIBUTE_UNUSED = 0.0833333333333333333333333333333*pow(dx,-1);
  const CCTK_REAL p1o12dy CCTK_ATTRIBUTE_UNUSED = 0.0833333333333333333333333333333*pow(dy,-1);
  const CCTK_REAL p1o12dz CCTK_ATTRIBUTE_UNUSED = 0.0833333333333333333333333333333*pow(dz,-1);
  const CCTK_REAL p1o144dxdy CCTK_ATTRIBUTE_UNUSED = 0.00694444444444444444444444444444*pow(dx,-1)*pow(dy,-1);
  const CCTK_REAL p1o144dxdz CCTK_ATTRIBUTE_UNUSED = 0.00694444444444444444444444444444*pow(dx,-1)*pow(dz,-1);
  const CCTK_REAL p1o144dydz CCTK_ATTRIBUTE_UNUSED = 0.00694444444444444444444444444444*pow(dy,-1)*pow(dz,-1);
  const CCTK_REAL p1o180dx2 CCTK_ATTRIBUTE_UNUSED = 0.00555555555555555555555555555556*pow(dx,-2);
  const CCTK_REAL p1o180dy2 CCTK_ATTRIBUTE_UNUSED = 0.00555555555555555555555555555556*pow(dy,-2);
  const CCTK_REAL p1o180dz2 CCTK_ATTRIBUTE_UNUSED = 0.00555555555555555555555555555556*pow(dz,-2);
  const CCTK_REAL p1o2dx CCTK_ATTRIBUTE_UNUSED = 0.5*pow(dx,-1);
  const CCTK_REAL p1o2dy CCTK_ATTRIBUTE_UNUSED = 0.5*pow(dy,-1);
  const CCTK_REAL p1o2dz CCTK_ATTRIBUTE_UNUSED = 0.5*pow(dz,-1);
  const CCTK_REAL p1o3600dxdy CCTK_ATTRIBUTE_UNUSED = 0.000277777777777777777777777777778*pow(dx,-1)*pow(dy,-1);
  const CCTK_REAL p1o3600dxdz CCTK_ATTRIBUTE_UNUSED = 0.000277777777777777777777777777778*pow(dx,-1)*pow(dz,-1);
  const CCTK_REAL p1o3600dydz CCTK_ATTRIBUTE_UNUSED = 0.000277777777777777777777777777778*pow(dy,-1)*pow(dz,-1);
  const CCTK_REAL p1o4dxdy CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dx,-1)*pow(dy,-1);
  const CCTK_REAL p1o4dxdz CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dx,-1)*pow(dz,-1);
  const CCTK_REAL p1o4dydz CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dy,-1)*pow(dz,-1);
  const CCTK_REAL p1o60dx CCTK_ATTRIBUTE_UNUSED = 0.0166666666666666666666666666667*pow(dx,-1);
  const CCTK_REAL p1o60dy CCTK_ATTRIBUTE_UNUSED = 0.0166666666666666666666666666667*pow(dy,-1);
  const CCTK_REAL p1o60dz CCTK_ATTRIBUTE_UNUSED = 0.0166666666666666666666666666667*pow(dz,-1);
  const CCTK_REAL p1odx CCTK_ATTRIBUTE_UNUSED = pow(dx,-1);
  const CCTK_REAL p1odx2 CCTK_ATTRIBUTE_UNUSED = pow(dx,-2);
  const CCTK_REAL p1odxdy CCTK_ATTRIBUTE_UNUSED = pow(dx,-1)*pow(dy,-1);
  const CCTK_REAL p1odxdz CCTK_ATTRIBUTE_UNUSED = pow(dx,-1)*pow(dz,-1);
  const CCTK_REAL p1ody CCTK_ATTRIBUTE_UNUSED = pow(dy,-1);
  const CCTK_REAL p1ody2 CCTK_ATTRIBUTE_UNUSED = pow(dy,-2);
  const CCTK_REAL p1odydz CCTK_ATTRIBUTE_UNUSED = pow(dy,-1)*pow(dz,-1);
  const CCTK_REAL p1odz CCTK_ATTRIBUTE_UNUSED = pow(dz,-1);
  const CCTK_REAL p1odz2 CCTK_ATTRIBUTE_UNUSED = pow(dz,-2);
  const CCTK_REAL pm1o12dx2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dx,-2);
  const CCTK_REAL pm1o12dy2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dy,-2);
  const CCTK_REAL pm1o12dz2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dz,-2);
  const CCTK_REAL pm1o2dx CCTK_ATTRIBUTE_UNUSED = -0.5*pow(dx,-1);
  const CCTK_REAL pm1o2dy CCTK_ATTRIBUTE_UNUSED = -0.5*pow(dy,-1);
  const CCTK_REAL pm1o2dz CCTK_ATTRIBUTE_UNUSED = -0.5*pow(dz,-1);
  const CCTK_REAL pm1o60dx CCTK_ATTRIBUTE_UNUSED = -0.0166666666666666666666666666667*pow(dx,-1);
  const CCTK_REAL pm1o60dy CCTK_ATTRIBUTE_UNUSED = -0.0166666666666666666666666666667*pow(dy,-1);
  const CCTK_REAL pm1o60dz CCTK_ATTRIBUTE_UNUSED = -0.0166666666666666666666666666667*pow(dz,-1);
  /* Assign local copies of arrays functions */
  
  
  /* Calculate temporaries and arrays functions */
  /* Copy local copies back to grid functions */
  /* Loop over the grid points */
  const int imin0=imin[0];
  const int imin1=imin[1];
  const int imin2=imin[2];
  const int imax0=imax[0];
  const int imax1=imax[1];
  const int imax2=imax[2];
  #pragma omp parallel
  CCTK_LOOP3(bssn_adm_to_evolved,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL alpL CCTK_ATTRIBUTE_UNUSED = alp[index];
    CCTK_REAL betaxL CCTK_ATTRIBUTE_UNUSED = betax[index];
    CCTK_REAL betayL CCTK_ATTRIBUTE_UNUSED = betay[index];
    CCTK_REAL betazL CCTK_ATTRIBUTE_UNUSED = betaz[index];
    CCTK_REAL gxxL CCTK_ATTRIBUTE_UNUSED = gxx[index];
    CCTK_REAL gxyL CCTK_ATTRIBUTE_UNUSED = gxy[index];
    CCTK_REAL gxzL CCTK_ATTRIBUTE_UNUSED = gxz[index];
    CCTK_REAL gyyL CCTK_ATTRIBUTE_UNUSED = gyy[index];
    CCTK_REAL gyzL CCTK_ATTRIBUTE_UNUSED = gyz[index];
    CCTK_REAL gzzL CCTK_ATTRIBUTE_UNUSED = gzz[index];
    CCTK_REAL KL CCTK_ATTRIBUTE_UNUSED = K[index];
    CCTK_REAL kxxL CCTK_ATTRIBUTE_UNUSED = kxx[index];
    CCTK_REAL kxyL CCTK_ATTRIBUTE_UNUSED = kxy[index];
    CCTK_REAL kxzL CCTK_ATTRIBUTE_UNUSED = kxz[index];
    CCTK_REAL kyyL CCTK_ATTRIBUTE_UNUSED = kyy[index];
    CCTK_REAL kyzL CCTK_ATTRIBUTE_UNUSED = kyz[index];
    CCTK_REAL kzzL CCTK_ATTRIBUTE_UNUSED = kzz[index];
    CCTK_REAL phiL CCTK_ATTRIBUTE_UNUSED = phi[index];
    
    /* Include user supplied include files */
    /* Precompute derivatives */
    
    switch (constraintsfdorder)
    {
      case 2:
      {
        break;
      }
      
      case 4:
      {
        break;
      }
      
      case 6:
      {
        break;
      }
      default:
        CCTK_BUILTIN_UNREACHABLE();
    }
    /* Calculate temporaries and grid functions */
    CCTK_REAL g11 CCTK_ATTRIBUTE_UNUSED = gxxL;
    
    CCTK_REAL g21 CCTK_ATTRIBUTE_UNUSED = gxyL;
    
    CCTK_REAL g22 CCTK_ATTRIBUTE_UNUSED = gyyL;
    
    CCTK_REAL g31 CCTK_ATTRIBUTE_UNUSED = gxzL;
    
    CCTK_REAL g32 CCTK_ATTRIBUTE_UNUSED = gyzL;
    
    CCTK_REAL g33 CCTK_ATTRIBUTE_UNUSED = gzzL;
    
    CCTK_REAL twelfth CCTK_ATTRIBUTE_UNUSED = 
      0.0833333333333333333333333333333;
    
    CCTK_REAL detg CCTK_ATTRIBUTE_UNUSED = 2*g21*g31*g32 - g33*pow(g21,2) 
      + g22*(g11*g33 - pow(g31,2)) - g11*pow(g32,2);
    
    CCTK_REAL detgmthirdroot CCTK_ATTRIBUTE_UNUSED = 
      pow(detg,-0.333333333333333333333333333333);
    
    CCTK_REAL invdetg CCTK_ATTRIBUTE_UNUSED = pow(detg,-1);
    
    CCTK_REAL gInv11 CCTK_ATTRIBUTE_UNUSED = invdetg*(g22*g33 - 
      pow(g32,2));
    
    CCTK_REAL gInv12 CCTK_ATTRIBUTE_UNUSED = (g31*g32 - g21*g33)*invdetg;
    
    CCTK_REAL gInv13 CCTK_ATTRIBUTE_UNUSED = (-(g22*g31) + 
      g21*g32)*invdetg;
    
    CCTK_REAL gInv21 CCTK_ATTRIBUTE_UNUSED = (g31*g32 - g21*g33)*invdetg;
    
    CCTK_REAL gInv22 CCTK_ATTRIBUTE_UNUSED = invdetg*(g11*g33 - 
      pow(g31,2));
    
    CCTK_REAL gInv23 CCTK_ATTRIBUTE_UNUSED = (g21*g31 - g11*g32)*invdetg;
    
    CCTK_REAL gInv31 CCTK_ATTRIBUTE_UNUSED = (-(g22*g31) + 
      g21*g32)*invdetg;
    
    CCTK_REAL gInv32 CCTK_ATTRIBUTE_UNUSED = (g21*g31 - g11*g32)*invdetg;
    
    CCTK_REAL gInv33 CCTK_ATTRIBUTE_UNUSED = invdetg*(g11*g22 - 
      pow(g21,2));
    
    CCTK_REAL k11 CCTK_ATTRIBUTE_UNUSED = kxxL;
    
    CCTK_REAL k21 CCTK_ATTRIBUTE_UNUSED = kxyL;
    
    CCTK_REAL k22 CCTK_ATTRIBUTE_UNUSED = kyyL;
    
    CCTK_REAL k31 CCTK_ATTRIBUTE_UNUSED = kxzL;
    
    CCTK_REAL k32 CCTK_ATTRIBUTE_UNUSED = kyzL;
    
    CCTK_REAL k33 CCTK_ATTRIBUTE_UNUSED = kzzL;
    
    CCTK_REAL h11L CCTK_ATTRIBUTE_UNUSED = detgmthirdroot*g11;
    
    CCTK_REAL h21L CCTK_ATTRIBUTE_UNUSED = detgmthirdroot*g21;
    
    CCTK_REAL h31L CCTK_ATTRIBUTE_UNUSED = detgmthirdroot*g31;
    
    CCTK_REAL h22L CCTK_ATTRIBUTE_UNUSED = detgmthirdroot*g22;
    
    CCTK_REAL h32L CCTK_ATTRIBUTE_UNUSED = detgmthirdroot*g32;
    
    CCTK_REAL h33L CCTK_ATTRIBUTE_UNUSED = detgmthirdroot*g33;
    
    phiL = twelfth*log(detg);
    
    KL = gInv11*k11 + (gInv12 + gInv21)*k21 + gInv22*k22 + (gInv13 + 
      gInv31)*k31 + (gInv23 + gInv32)*k32 + gInv33*k33;
    
    CCTK_REAL em4phi CCTK_ATTRIBUTE_UNUSED = exp(-4*phiL);
    
    CCTK_REAL A11L CCTK_ATTRIBUTE_UNUSED = 
      em4phi*(-0.333333333333333333333333333333*KL*g11 + k11);
    
    CCTK_REAL A21L CCTK_ATTRIBUTE_UNUSED = 
      em4phi*(-0.333333333333333333333333333333*KL*g21 + k21);
    
    CCTK_REAL A31L CCTK_ATTRIBUTE_UNUSED = 
      em4phi*(-0.333333333333333333333333333333*KL*g31 + k31);
    
    CCTK_REAL A22L CCTK_ATTRIBUTE_UNUSED = 
      em4phi*(-0.333333333333333333333333333333*KL*g22 + k22);
    
    CCTK_REAL A32L CCTK_ATTRIBUTE_UNUSED = 
      em4phi*(-0.333333333333333333333333333333*KL*g32 + k32);
    
    CCTK_REAL A33L CCTK_ATTRIBUTE_UNUSED = 
      em4phi*(-0.333333333333333333333333333333*KL*g33 + k33);
    
    CCTK_REAL alphaL CCTK_ATTRIBUTE_UNUSED = alpL;
    
    CCTK_REAL beta1L CCTK_ATTRIBUTE_UNUSED = betaxL;
    
    CCTK_REAL beta2L CCTK_ATTRIBUTE_UNUSED = betayL;
    
    CCTK_REAL beta3L CCTK_ATTRIBUTE_UNUSED = betazL;
    
    CCTK_REAL betat1L CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL betat2L CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL betat3L CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL Gam1L CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL Gam2L CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL Gam3L CCTK_ATTRIBUTE_UNUSED = 0;
    /* Copy local copies back to grid functions */
    A11[index] = A11L;
    A21[index] = A21L;
    A22[index] = A22L;
    A31[index] = A31L;
    A32[index] = A32L;
    A33[index] = A33L;
    alpha[index] = alphaL;
    beta1[index] = beta1L;
    beta2[index] = beta2L;
    beta3[index] = beta3L;
    betat1[index] = betat1L;
    betat2[index] = betat2L;
    betat3[index] = betat3L;
    Gam1[index] = Gam1L;
    Gam2[index] = Gam2L;
    Gam3[index] = Gam3L;
    h11[index] = h11L;
    h21[index] = h21L;
    h22[index] = h22L;
    h31[index] = h31L;
    h32[index] = h32L;
    h33[index] = h33L;
    K[index] = KL;
    phi[index] = phiL;
  }
  CCTK_ENDLOOP3(bssn_adm_to_evolved);
}
extern "C" void bssn_adm_to_evolved(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering bssn_adm_to_evolved_Body");
  }
  if (cctk_iteration % bssn_adm_to_evolved_calc_every != bssn_adm_to_evolved_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "admbase::curv",
    "admbase::lapse",
    "admbase::metric",
    "admbase::shift",
    "Kranc2BSSN::A_group",
    "Kranc2BSSN::alpha_group",
    "Kranc2BSSN::beta_group",
    "Kranc2BSSN::betat_group",
    "Kranc2BSSN::Gam_group",
    "Kranc2BSSN::h_group",
    "Kranc2BSSN::K_group",
    "Kranc2BSSN::phi_group"};
  AssertGroupStorage(cctkGH, "bssn_adm_to_evolved", 12, groups);
  
  switch (constraintsfdorder)
  {
    case 2:
    {
      break;
    }
    
    case 4:
    {
      break;
    }
    
    case 6:
    {
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  LoopOverEverything(cctkGH, bssn_adm_to_evolved_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving bssn_adm_to_evolved_Body");
  }
}

} // namespace Kranc2BSSN
