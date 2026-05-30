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

namespace EvolveScalarFields {

extern "C" void EvolveScalarFields_ScalarFieldTmunu_none_full4th_SelectBCs(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_EvolveScalarFields_ScalarFieldTmunu_none_full4th_SelectBCs
  DECLARE_CCTK_ARGUMENTS_CHECKED(EvolveScalarFields_ScalarFieldTmunu_none_full4th_SelectBCs);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % EvolveScalarFields_ScalarFieldTmunu_none_full4th_calc_every != EvolveScalarFields_ScalarFieldTmunu_none_full4th_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "EvolveScalarFields::SFTtt_group","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for EvolveScalarFields::SFTtt_group.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "TmunuBase::stress_energy_scalar","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for TmunuBase::stress_energy_scalar.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "TmunuBase::stress_energy_tensor","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for TmunuBase::stress_energy_tensor.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "TmunuBase::stress_energy_vector","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for TmunuBase::stress_energy_vector.");
  return;
}

static void EvolveScalarFields_ScalarFieldTmunu_none_full4th_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  CCTK_LOOP3(EvolveScalarFields_ScalarFieldTmunu_none_full4th,
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
    CCTK_REAL phiaL CCTK_ATTRIBUTE_UNUSED = phia[index];
    CCTK_REAL phibL CCTK_ATTRIBUTE_UNUSED = phib[index];
    CCTK_REAL piaL CCTK_ATTRIBUTE_UNUSED = pia[index];
    CCTK_REAL pibL CCTK_ATTRIBUTE_UNUSED = pib[index];
    CCTK_REAL SFTttL CCTK_ATTRIBUTE_UNUSED = SFTtt[index];
    
    CCTK_REAL eTttL, eTtxL, eTtyL, eTtzL, eTxxL, eTxyL, eTxzL, eTyyL, eTyzL, eTzzL CCTK_ATTRIBUTE_UNUSED ;
    
    if (assume_stress_energy_state>=0 ? assume_stress_energy_state : *stress_energy_state)
    {
      eTttL = eTtt[index];
      eTtxL = eTtx[index];
      eTtyL = eTty[index];
      eTtzL = eTtz[index];
      eTxxL = eTxx[index];
      eTxyL = eTxy[index];
      eTxzL = eTxz[index];
      eTyyL = eTyy[index];
      eTyzL = eTyz[index];
      eTzzL = eTzz[index];
    }
    else
    {
      eTttL = 0.;
      eTtxL = 0.;
      eTtyL = 0.;
      eTtzL = 0.;
      eTxxL = 0.;
      eTxyL = 0.;
      eTxzL = 0.;
      eTyyL = 0.;
      eTyzL = 0.;
      eTzzL = 0.;
    }
    /* Include user supplied include files */
    /* Precompute derivatives */
    const CCTK_REAL PDstandard4th1phia CCTK_ATTRIBUTE_UNUSED = PDstandard4th1(&phia[index]);
    const CCTK_REAL PDstandard4th2phia CCTK_ATTRIBUTE_UNUSED = PDstandard4th2(&phia[index]);
    const CCTK_REAL PDstandard4th3phia CCTK_ATTRIBUTE_UNUSED = PDstandard4th3(&phia[index]);
    const CCTK_REAL PDstandard4th1phib CCTK_ATTRIBUTE_UNUSED = PDstandard4th1(&phib[index]);
    const CCTK_REAL PDstandard4th2phib CCTK_ATTRIBUTE_UNUSED = PDstandard4th2(&phib[index]);
    const CCTK_REAL PDstandard4th3phib CCTK_ATTRIBUTE_UNUSED = PDstandard4th3(&phib[index]);
    /* Calculate temporaries and grid functions */
    ptrdiff_t dir1 CCTK_ATTRIBUTE_UNUSED = isgn(betaxL);
    
    ptrdiff_t dir2 CCTK_ATTRIBUTE_UNUSED = isgn(betayL);
    
    ptrdiff_t dir3 CCTK_ATTRIBUTE_UNUSED = isgn(betazL);
    
    CCTK_REAL gInv11 CCTK_ATTRIBUTE_UNUSED = (-(gyyL*gzzL) + 
      pow(gyzL,2))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*(-(gxxL*gzzL) + pow(gxzL,2)) + gxxL*pow(gyzL,2),-1);
    
    CCTK_REAL gInv12 CCTK_ATTRIBUTE_UNUSED = (-(gxzL*gyzL) + 
      gxyL*gzzL)*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + gyyL*pow(gxzL,2) 
      + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-1);
    
    CCTK_REAL gInv13 CCTK_ATTRIBUTE_UNUSED = (gxzL*gyyL - 
      gxyL*gyzL)*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*(-(gxxL*gzzL) + pow(gxzL,2)) + gxxL*pow(gyzL,2),-1);
    
    CCTK_REAL gInv22 CCTK_ATTRIBUTE_UNUSED = (-(gxxL*gzzL) + 
      pow(gxzL,2))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*(-(gxxL*gzzL) + pow(gxzL,2)) + gxxL*pow(gyzL,2),-1);
    
    CCTK_REAL gInv23 CCTK_ATTRIBUTE_UNUSED = (-(gxyL*gxzL) + 
      gxxL*gyzL)*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + gyyL*pow(gxzL,2) 
      + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-1);
    
    CCTK_REAL gInv33 CCTK_ATTRIBUTE_UNUSED = (-(gxxL*gyyL) + 
      pow(gxyL,2))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*(-(gxxL*gzzL) + pow(gxzL,2)) + gxxL*pow(gyzL,2),-1);
    
    CCTK_REAL dtphia CCTK_ATTRIBUTE_UNUSED = 
      betaxL*PDlopsided4th1(&phia[index]) + 
      betayL*PDlopsided4th2(&phia[index]) + 
      betazL*PDlopsided4th3(&phia[index]) - alpL*piaL;
    
    CCTK_REAL dtphib CCTK_ATTRIBUTE_UNUSED = 
      betaxL*PDlopsided4th1(&phib[index]) + 
      betayL*PDlopsided4th2(&phib[index]) + 
      betazL*PDlopsided4th3(&phib[index]) - alpL*pibL;
    
    CCTK_REAL betaSqr CCTK_ATTRIBUTE_UNUSED = 2*(betaxL*(betayL*gxyL + 
      betazL*gxzL) + betayL*betazL*gyzL) + gxxL*pow(betaxL,2) + 
      gyyL*pow(betayL,2) + gzzL*pow(betazL,2);
    
    CCTK_REAL dphiaSqr CCTK_ATTRIBUTE_UNUSED = 
      2*(gInv23*PDstandard4th2phia*PDstandard4th3phia + 
      PDstandard4th1phia*(gInv12*PDstandard4th2phia + 
      gInv13*PDstandard4th3phia)) - pow(piaL,2) + 
      gInv11*pow(PDstandard4th1phia,2) + gInv22*pow(PDstandard4th2phia,2) + 
      gInv33*pow(PDstandard4th3phia,2);
    
    CCTK_REAL dphibSqr CCTK_ATTRIBUTE_UNUSED = 
      2*(gInv23*PDstandard4th2phib*PDstandard4th3phib + 
      PDstandard4th1phib*(gInv12*PDstandard4th2phib + 
      gInv13*PDstandard4th3phib)) - pow(pibL,2) + 
      gInv11*pow(PDstandard4th1phib,2) + gInv22*pow(PDstandard4th2phib,2) + 
      gInv33*pow(PDstandard4th3phib,2);
    
    SFTttL = 0.5*(dphiaSqr + dphibSqr)*(-betaSqr + pow(alpL,2)) + 
      pow(dtphia,2) + pow(dtphib,2);
    
    CCTK_REAL SFTtj1 CCTK_ATTRIBUTE_UNUSED = -0.5*(betaxL*gxxL + 
      betayL*gxyL + betazL*gxzL)*(dphiaSqr + dphibSqr) + 
      dtphia*PDstandard4th1phia + dtphib*PDstandard4th1phib;
    
    CCTK_REAL SFTtj2 CCTK_ATTRIBUTE_UNUSED = -0.5*(betaxL*gxyL + 
      betayL*gyyL + betazL*gyzL)*(dphiaSqr + dphibSqr) + 
      dtphia*PDstandard4th2phia + dtphib*PDstandard4th2phib;
    
    CCTK_REAL SFTtj3 CCTK_ATTRIBUTE_UNUSED = -0.5*(betaxL*gxzL + 
      betayL*gyzL + betazL*gzzL)*(dphiaSqr + dphibSqr) + 
      dtphia*PDstandard4th3phia + dtphib*PDstandard4th3phib;
    
    CCTK_REAL SFTij11 CCTK_ATTRIBUTE_UNUSED = -0.5*gxxL*(dphiaSqr + 
      dphibSqr) + pow(PDstandard4th1phia,2) + pow(PDstandard4th1phib,2);
    
    CCTK_REAL SFTij21 CCTK_ATTRIBUTE_UNUSED = -0.5*gxyL*(dphiaSqr + 
      dphibSqr) + PDstandard4th1phia*PDstandard4th2phia + 
      PDstandard4th1phib*PDstandard4th2phib;
    
    CCTK_REAL SFTij31 CCTK_ATTRIBUTE_UNUSED = -0.5*gxzL*(dphiaSqr + 
      dphibSqr) + PDstandard4th1phia*PDstandard4th3phia + 
      PDstandard4th1phib*PDstandard4th3phib;
    
    CCTK_REAL SFTij22 CCTK_ATTRIBUTE_UNUSED = -0.5*gyyL*(dphiaSqr + 
      dphibSqr) + pow(PDstandard4th2phia,2) + pow(PDstandard4th2phib,2);
    
    CCTK_REAL SFTij32 CCTK_ATTRIBUTE_UNUSED = -0.5*gyzL*(dphiaSqr + 
      dphibSqr) + PDstandard4th2phia*PDstandard4th3phia + 
      PDstandard4th2phib*PDstandard4th3phib;
    
    CCTK_REAL SFTij33 CCTK_ATTRIBUTE_UNUSED = -0.5*gzzL*(dphiaSqr + 
      dphibSqr) + pow(PDstandard4th3phia,2) + pow(PDstandard4th3phib,2);
    
    eTttL = eTttL + SFTttL;
    
    eTtxL = eTtxL + SFTtj1;
    
    eTtyL = eTtyL + SFTtj2;
    
    eTtzL = eTtzL + SFTtj3;
    
    eTxxL = eTxxL + SFTij11;
    
    eTxyL = eTxyL + SFTij21;
    
    eTxzL = eTxzL + SFTij31;
    
    eTyyL = eTyyL + SFTij22;
    
    eTyzL = eTyzL + SFTij32;
    
    eTzzL = eTzzL + SFTij33;
    /* Copy local copies back to grid functions */
    eTtt[index] = eTttL;
    eTtx[index] = eTtxL;
    eTty[index] = eTtyL;
    eTtz[index] = eTtzL;
    eTxx[index] = eTxxL;
    eTxy[index] = eTxyL;
    eTxz[index] = eTxzL;
    eTyy[index] = eTyyL;
    eTyz[index] = eTyzL;
    eTzz[index] = eTzzL;
    SFTtt[index] = SFTttL;
  }
  CCTK_ENDLOOP3(EvolveScalarFields_ScalarFieldTmunu_none_full4th);
}
extern "C" void EvolveScalarFields_ScalarFieldTmunu_none_full4th(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_EvolveScalarFields_ScalarFieldTmunu_none_full4th
  DECLARE_CCTK_ARGUMENTS_CHECKED(EvolveScalarFields_ScalarFieldTmunu_none_full4th);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering EvolveScalarFields_ScalarFieldTmunu_none_full4th_Body");
  }
  if (cctk_iteration % EvolveScalarFields_ScalarFieldTmunu_none_full4th_calc_every != EvolveScalarFields_ScalarFieldTmunu_none_full4th_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "admbase::lapse",
    "admbase::metric",
    "admbase::shift",
    "EvolveScalarFields::phia_group",
    "EvolveScalarFields::phib_group",
    "EvolveScalarFields::pia_group",
    "EvolveScalarFields::pib_group",
    "EvolveScalarFields::SFTtt_group"};
  AssertGroupStorage(cctkGH, "EvolveScalarFields_ScalarFieldTmunu_none_full4th", 8, groups);
  
  EnsureStencilFits(cctkGH, "EvolveScalarFields_ScalarFieldTmunu_none_full4th", 3, 3, 3);
  
  LoopOverInterior(cctkGH, EvolveScalarFields_ScalarFieldTmunu_none_full4th_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving EvolveScalarFields_ScalarFieldTmunu_none_full4th_Body");
  }
}

} // namespace EvolveScalarFields
