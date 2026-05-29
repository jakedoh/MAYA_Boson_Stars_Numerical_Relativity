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

namespace ScalarTensor {

extern "C" void ScalarTensor_calc_rho_phiBubble_6th_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % ScalarTensor_calc_rho_phiBubble_6th_calc_every != ScalarTensor_calc_rho_phiBubble_6th_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ScalarTensor::rhoE_group","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ScalarTensor::rhoE_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ScalarTensor::rhoJ_group","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ScalarTensor::rhoJ_group.");
  return;
}

static void ScalarTensor_calc_rho_phiBubble_6th_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  CCTK_LOOP3(ScalarTensor_calc_rho_phiBubble_6th,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL betaxL CCTK_ATTRIBUTE_UNUSED = betax[index];
    CCTK_REAL betayL CCTK_ATTRIBUTE_UNUSED = betay[index];
    CCTK_REAL betazL CCTK_ATTRIBUTE_UNUSED = betaz[index];
    CCTK_REAL gLocal11L CCTK_ATTRIBUTE_UNUSED = gLocal11[index];
    CCTK_REAL gLocal21L CCTK_ATTRIBUTE_UNUSED = gLocal21[index];
    CCTK_REAL gLocal22L CCTK_ATTRIBUTE_UNUSED = gLocal22[index];
    CCTK_REAL gLocal31L CCTK_ATTRIBUTE_UNUSED = gLocal31[index];
    CCTK_REAL gLocal32L CCTK_ATTRIBUTE_UNUSED = gLocal32[index];
    CCTK_REAL gLocal33L CCTK_ATTRIBUTE_UNUSED = gLocal33[index];
    CCTK_REAL rhoEL CCTK_ATTRIBUTE_UNUSED = rhoE[index];
    CCTK_REAL STphiL CCTK_ATTRIBUTE_UNUSED = STphi[index];
    CCTK_REAL STpiL CCTK_ATTRIBUTE_UNUSED = STpi[index];
    
    /* Include user supplied include files */
    /* Precompute derivatives */
    const CCTK_REAL PDstandard6th1STphi CCTK_ATTRIBUTE_UNUSED = PDstandard6th1(&STphi[index]);
    const CCTK_REAL PDstandard6th2STphi CCTK_ATTRIBUTE_UNUSED = PDstandard6th2(&STphi[index]);
    const CCTK_REAL PDstandard6th3STphi CCTK_ATTRIBUTE_UNUSED = PDstandard6th3(&STphi[index]);
    /* Calculate temporaries and grid functions */
    ptrdiff_t dir1 CCTK_ATTRIBUTE_UNUSED = isgn(betaxL);
    
    ptrdiff_t dir2 CCTK_ATTRIBUTE_UNUSED = isgn(betayL);
    
    ptrdiff_t dir3 CCTK_ATTRIBUTE_UNUSED = isgn(betazL);
    
    CCTK_REAL detg CCTK_ATTRIBUTE_UNUSED = 2*gLocal21L*gLocal31L*gLocal32L 
      - gLocal33L*pow(gLocal21L,2) + gLocal22L*(gLocal11L*gLocal33L - 
      pow(gLocal31L,2)) - gLocal11L*pow(gLocal32L,2);
    
    CCTK_REAL invdetg CCTK_ATTRIBUTE_UNUSED = pow(detg,-1);
    
    CCTK_REAL gLocalInv11 CCTK_ATTRIBUTE_UNUSED = 
      invdetg*(gLocal22L*gLocal33L - pow(gLocal32L,2));
    
    CCTK_REAL gLocalInv12 CCTK_ATTRIBUTE_UNUSED = (gLocal31L*gLocal32L - 
      gLocal21L*gLocal33L)*invdetg;
    
    CCTK_REAL gLocalInv13 CCTK_ATTRIBUTE_UNUSED = (-(gLocal22L*gLocal31L) 
      + gLocal21L*gLocal32L)*invdetg;
    
    CCTK_REAL gLocalInv22 CCTK_ATTRIBUTE_UNUSED = 
      invdetg*(gLocal11L*gLocal33L - pow(gLocal31L,2));
    
    CCTK_REAL gLocalInv23 CCTK_ATTRIBUTE_UNUSED = (gLocal21L*gLocal31L - 
      gLocal11L*gLocal32L)*invdetg;
    
    CCTK_REAL gLocalInv33 CCTK_ATTRIBUTE_UNUSED = 
      invdetg*(gLocal11L*gLocal22L - pow(gLocal21L,2));
    
    rhoEL = 1.*(gLocalInv23*PDstandard6th2STphi*PDstandard6th3STphi + 
      PDstandard6th1STphi*(gLocalInv12*PDstandard6th2STphi + 
      gLocalInv13*PDstandard6th3STphi)) + 0.5*(pow(STpiL,2) + 
      gLocalInv11*pow(PDstandard6th1STphi,2) + 
      gLocalInv22*pow(PDstandard6th2STphi,2) + 
      gLocalInv33*pow(PDstandard6th3STphi,2)) + 0.125*Vlambda*(4*(-STphiL + 
      phiBG)*Veps*pow(phiBG,3) + pow(-pow(STphiL,2) + pow(phiBG,2),2));
    
    CCTK_REAL ConfFactorA CCTK_ATTRIBUTE_UNUSED = exp(-(phiBG*STalpha0) + 
      STphiL*(STalpha0 + 0.5*STphiL*STBeta) - 0.5*STBeta*pow(phiBG,2));
    
    CCTK_REAL rhoJL CCTK_ATTRIBUTE_UNUSED = rhoEL*pow(ConfFactorA,-4);
    /* Copy local copies back to grid functions */
    rhoE[index] = rhoEL;
    rhoJ[index] = rhoJL;
  }
  CCTK_ENDLOOP3(ScalarTensor_calc_rho_phiBubble_6th);
}
extern "C" void ScalarTensor_calc_rho_phiBubble_6th(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ScalarTensor_calc_rho_phiBubble_6th_Body");
  }
  if (cctk_iteration % ScalarTensor_calc_rho_phiBubble_6th_calc_every != ScalarTensor_calc_rho_phiBubble_6th_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "admbase::shift",
    "ScalarTensor::gLocal_group",
    "ScalarTensor::rhoE_group",
    "ScalarTensor::rhoJ_group",
    "ScalarTensor::STphi_group",
    "ScalarTensor::STpi_group"};
  AssertGroupStorage(cctkGH, "ScalarTensor_calc_rho_phiBubble_6th", 6, groups);
  
  EnsureStencilFits(cctkGH, "ScalarTensor_calc_rho_phiBubble_6th", 3, 3, 3);
  
  LoopOverInterior(cctkGH, ScalarTensor_calc_rho_phiBubble_6th_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ScalarTensor_calc_rho_phiBubble_6th_Body");
  }
}

} // namespace ScalarTensor
