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

namespace ML_ADMQuantities {

extern "C" void ML_ADMQuantities_calculate_SelectBCs(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_ML_ADMQuantities_calculate_SelectBCs
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_ADMQuantities_calculate_SelectBCs);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % ML_ADMQuantities_calculate_calc_every != ML_ADMQuantities_calculate_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADMQuantities::ML_Jadm","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_ADMQuantities::ML_Jadm.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADMQuantities::ML_Madm","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_ADMQuantities::ML_Madm.");
  return;
}

static void ML_ADMQuantities_calculate_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  const CCTK_REAL p1odx CCTK_ATTRIBUTE_UNUSED = pow(dx,-1);
  const CCTK_REAL p1ody CCTK_ATTRIBUTE_UNUSED = pow(dy,-1);
  const CCTK_REAL p1odz CCTK_ATTRIBUTE_UNUSED = pow(dz,-1);
  const CCTK_REAL pm1o12dx2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dx,-2);
  const CCTK_REAL pm1o12dy2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dy,-2);
  const CCTK_REAL pm1o12dz2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dz,-2);
  /* Jacobian variable pointers */
  const bool usejacobian1 = (!CCTK_IsFunctionAliased("MultiPatch_GetMap") || MultiPatch_GetMap(cctkGH) != jacobian_identity_map)
                        && strlen(jacobian_group) > 0;
  const bool usejacobian = assume_use_jacobian>=0 ? assume_use_jacobian : usejacobian1;
  if (usejacobian && (strlen(jacobian_derivative_group) == 0))
  {
    CCTK_WARN(CCTK_WARN_ALERT, "GenericFD::jacobian_group and GenericFD::jacobian_derivative_group must both be set to valid group names");
  }
  
  const CCTK_REAL* restrict jacobian_ptrs[9];
  if (usejacobian) GroupDataPointers(cctkGH, jacobian_group,
                                                9, jacobian_ptrs);
  
  const CCTK_REAL* restrict const J11 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[0] : 0;
  const CCTK_REAL* restrict const J12 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[1] : 0;
  const CCTK_REAL* restrict const J13 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[2] : 0;
  const CCTK_REAL* restrict const J21 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[3] : 0;
  const CCTK_REAL* restrict const J22 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[4] : 0;
  const CCTK_REAL* restrict const J23 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[5] : 0;
  const CCTK_REAL* restrict const J31 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[6] : 0;
  const CCTK_REAL* restrict const J32 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[7] : 0;
  const CCTK_REAL* restrict const J33 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[8] : 0;
  
  const CCTK_REAL* restrict jacobian_determinant_ptrs[1] CCTK_ATTRIBUTE_UNUSED;
  if (usejacobian && strlen(jacobian_determinant_group) > 0) GroupDataPointers(cctkGH, jacobian_determinant_group,
                                                1, jacobian_determinant_ptrs);
  
  const CCTK_REAL* restrict const detJ CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_determinant_ptrs[0] : 0;
  
  const CCTK_REAL* restrict jacobian_inverse_ptrs[9] CCTK_ATTRIBUTE_UNUSED;
  if (usejacobian && strlen(jacobian_inverse_group) > 0) GroupDataPointers(cctkGH, jacobian_inverse_group,
                                                9, jacobian_inverse_ptrs);
  
  const CCTK_REAL* restrict const iJ11 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[0] : 0;
  const CCTK_REAL* restrict const iJ12 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[1] : 0;
  const CCTK_REAL* restrict const iJ13 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[2] : 0;
  const CCTK_REAL* restrict const iJ21 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[3] : 0;
  const CCTK_REAL* restrict const iJ22 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[4] : 0;
  const CCTK_REAL* restrict const iJ23 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[5] : 0;
  const CCTK_REAL* restrict const iJ31 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[6] : 0;
  const CCTK_REAL* restrict const iJ32 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[7] : 0;
  const CCTK_REAL* restrict const iJ33 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[8] : 0;
  
  const CCTK_REAL* restrict jacobian_derivative_ptrs[18] CCTK_ATTRIBUTE_UNUSED;
  if (usejacobian) GroupDataPointers(cctkGH, jacobian_derivative_group,
                                      18, jacobian_derivative_ptrs);
  
  const CCTK_REAL* restrict const dJ111 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[0] : 0;
  const CCTK_REAL* restrict const dJ112 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[1] : 0;
  const CCTK_REAL* restrict const dJ113 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[2] : 0;
  const CCTK_REAL* restrict const dJ122 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[3] : 0;
  const CCTK_REAL* restrict const dJ123 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[4] : 0;
  const CCTK_REAL* restrict const dJ133 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[5] : 0;
  const CCTK_REAL* restrict const dJ211 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[6] : 0;
  const CCTK_REAL* restrict const dJ212 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[7] : 0;
  const CCTK_REAL* restrict const dJ213 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[8] : 0;
  const CCTK_REAL* restrict const dJ222 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[9] : 0;
  const CCTK_REAL* restrict const dJ223 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[10] : 0;
  const CCTK_REAL* restrict const dJ233 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[11] : 0;
  const CCTK_REAL* restrict const dJ311 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[12] : 0;
  const CCTK_REAL* restrict const dJ312 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[13] : 0;
  const CCTK_REAL* restrict const dJ313 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[14] : 0;
  const CCTK_REAL* restrict const dJ322 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[15] : 0;
  const CCTK_REAL* restrict const dJ323 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[16] : 0;
  const CCTK_REAL* restrict const dJ333 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[17] : 0;
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
  CCTK_LOOP3(ML_ADMQuantities_calculate,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL alphaL CCTK_ATTRIBUTE_UNUSED = alpha[index];
    CCTK_REAL At11L CCTK_ATTRIBUTE_UNUSED = At11[index];
    CCTK_REAL At12L CCTK_ATTRIBUTE_UNUSED = At12[index];
    CCTK_REAL At13L CCTK_ATTRIBUTE_UNUSED = At13[index];
    CCTK_REAL At22L CCTK_ATTRIBUTE_UNUSED = At22[index];
    CCTK_REAL At23L CCTK_ATTRIBUTE_UNUSED = At23[index];
    CCTK_REAL At33L CCTK_ATTRIBUTE_UNUSED = At33[index];
    CCTK_REAL beta1L CCTK_ATTRIBUTE_UNUSED = beta1[index];
    CCTK_REAL beta2L CCTK_ATTRIBUTE_UNUSED = beta2[index];
    CCTK_REAL beta3L CCTK_ATTRIBUTE_UNUSED = beta3[index];
    CCTK_REAL gt11L CCTK_ATTRIBUTE_UNUSED = gt11[index];
    CCTK_REAL gt12L CCTK_ATTRIBUTE_UNUSED = gt12[index];
    CCTK_REAL gt13L CCTK_ATTRIBUTE_UNUSED = gt13[index];
    CCTK_REAL gt22L CCTK_ATTRIBUTE_UNUSED = gt22[index];
    CCTK_REAL gt23L CCTK_ATTRIBUTE_UNUSED = gt23[index];
    CCTK_REAL gt33L CCTK_ATTRIBUTE_UNUSED = gt33[index];
    CCTK_REAL phiL CCTK_ATTRIBUTE_UNUSED = phi[index];
    CCTK_REAL trKL CCTK_ATTRIBUTE_UNUSED = trK[index];
    CCTK_REAL xL CCTK_ATTRIBUTE_UNUSED = x[index];
    CCTK_REAL Xt1L CCTK_ATTRIBUTE_UNUSED = Xt1[index];
    CCTK_REAL Xt2L CCTK_ATTRIBUTE_UNUSED = Xt2[index];
    CCTK_REAL Xt3L CCTK_ATTRIBUTE_UNUSED = Xt3[index];
    CCTK_REAL yL CCTK_ATTRIBUTE_UNUSED = y[index];
    CCTK_REAL zL CCTK_ATTRIBUTE_UNUSED = z[index];
    
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
    
    CCTK_REAL dJ111L, dJ112L, dJ113L, dJ122L, dJ123L, dJ133L, dJ211L, dJ212L, dJ213L, dJ222L, dJ223L, dJ233L, dJ311L, dJ312L, dJ313L, dJ322L, dJ323L, dJ333L, J11L, J12L, J13L, J21L, J22L, J23L, J31L, J32L, J33L CCTK_ATTRIBUTE_UNUSED ;
    
    if (usejacobian)
    {
      dJ111L = dJ111[index];
      dJ112L = dJ112[index];
      dJ113L = dJ113[index];
      dJ122L = dJ122[index];
      dJ123L = dJ123[index];
      dJ133L = dJ133[index];
      dJ211L = dJ211[index];
      dJ212L = dJ212[index];
      dJ213L = dJ213[index];
      dJ222L = dJ222[index];
      dJ223L = dJ223[index];
      dJ233L = dJ233[index];
      dJ311L = dJ311[index];
      dJ312L = dJ312[index];
      dJ313L = dJ313[index];
      dJ322L = dJ322[index];
      dJ323L = dJ323[index];
      dJ333L = dJ333[index];
      J11L = J11[index];
      J12L = J12[index];
      J13L = J13[index];
      J21L = J21[index];
      J22L = J22[index];
      J23L = J23[index];
      J31L = J31[index];
      J32L = J32[index];
      J33L = J33[index];
    }
    /* Include user supplied include files */
    /* Precompute derivatives */
    const CCTK_REAL PDstandardNth1gt11 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&gt11[index]);
    const CCTK_REAL PDstandardNth2gt11 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&gt11[index]);
    const CCTK_REAL PDstandardNth3gt11 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&gt11[index]);
    const CCTK_REAL PDstandardNth11gt11 CCTK_ATTRIBUTE_UNUSED = PDstandardNth11(&gt11[index]);
    const CCTK_REAL PDstandardNth22gt11 CCTK_ATTRIBUTE_UNUSED = PDstandardNth22(&gt11[index]);
    const CCTK_REAL PDstandardNth33gt11 CCTK_ATTRIBUTE_UNUSED = PDstandardNth33(&gt11[index]);
    const CCTK_REAL PDstandardNth12gt11 CCTK_ATTRIBUTE_UNUSED = PDstandardNth12(&gt11[index]);
    const CCTK_REAL PDstandardNth13gt11 CCTK_ATTRIBUTE_UNUSED = PDstandardNth13(&gt11[index]);
    const CCTK_REAL PDstandardNth23gt11 CCTK_ATTRIBUTE_UNUSED = PDstandardNth23(&gt11[index]);
    const CCTK_REAL PDstandardNth1gt12 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&gt12[index]);
    const CCTK_REAL PDstandardNth2gt12 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&gt12[index]);
    const CCTK_REAL PDstandardNth3gt12 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&gt12[index]);
    const CCTK_REAL PDstandardNth11gt12 CCTK_ATTRIBUTE_UNUSED = PDstandardNth11(&gt12[index]);
    const CCTK_REAL PDstandardNth22gt12 CCTK_ATTRIBUTE_UNUSED = PDstandardNth22(&gt12[index]);
    const CCTK_REAL PDstandardNth33gt12 CCTK_ATTRIBUTE_UNUSED = PDstandardNth33(&gt12[index]);
    const CCTK_REAL PDstandardNth12gt12 CCTK_ATTRIBUTE_UNUSED = PDstandardNth12(&gt12[index]);
    const CCTK_REAL PDstandardNth13gt12 CCTK_ATTRIBUTE_UNUSED = PDstandardNth13(&gt12[index]);
    const CCTK_REAL PDstandardNth23gt12 CCTK_ATTRIBUTE_UNUSED = PDstandardNth23(&gt12[index]);
    const CCTK_REAL PDstandardNth1gt13 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&gt13[index]);
    const CCTK_REAL PDstandardNth2gt13 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&gt13[index]);
    const CCTK_REAL PDstandardNth3gt13 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&gt13[index]);
    const CCTK_REAL PDstandardNth11gt13 CCTK_ATTRIBUTE_UNUSED = PDstandardNth11(&gt13[index]);
    const CCTK_REAL PDstandardNth22gt13 CCTK_ATTRIBUTE_UNUSED = PDstandardNth22(&gt13[index]);
    const CCTK_REAL PDstandardNth33gt13 CCTK_ATTRIBUTE_UNUSED = PDstandardNth33(&gt13[index]);
    const CCTK_REAL PDstandardNth12gt13 CCTK_ATTRIBUTE_UNUSED = PDstandardNth12(&gt13[index]);
    const CCTK_REAL PDstandardNth13gt13 CCTK_ATTRIBUTE_UNUSED = PDstandardNth13(&gt13[index]);
    const CCTK_REAL PDstandardNth23gt13 CCTK_ATTRIBUTE_UNUSED = PDstandardNth23(&gt13[index]);
    const CCTK_REAL PDstandardNth1gt22 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&gt22[index]);
    const CCTK_REAL PDstandardNth2gt22 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&gt22[index]);
    const CCTK_REAL PDstandardNth3gt22 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&gt22[index]);
    const CCTK_REAL PDstandardNth11gt22 CCTK_ATTRIBUTE_UNUSED = PDstandardNth11(&gt22[index]);
    const CCTK_REAL PDstandardNth22gt22 CCTK_ATTRIBUTE_UNUSED = PDstandardNth22(&gt22[index]);
    const CCTK_REAL PDstandardNth33gt22 CCTK_ATTRIBUTE_UNUSED = PDstandardNth33(&gt22[index]);
    const CCTK_REAL PDstandardNth12gt22 CCTK_ATTRIBUTE_UNUSED = PDstandardNth12(&gt22[index]);
    const CCTK_REAL PDstandardNth13gt22 CCTK_ATTRIBUTE_UNUSED = PDstandardNth13(&gt22[index]);
    const CCTK_REAL PDstandardNth23gt22 CCTK_ATTRIBUTE_UNUSED = PDstandardNth23(&gt22[index]);
    const CCTK_REAL PDstandardNth1gt23 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&gt23[index]);
    const CCTK_REAL PDstandardNth2gt23 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&gt23[index]);
    const CCTK_REAL PDstandardNth3gt23 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&gt23[index]);
    const CCTK_REAL PDstandardNth11gt23 CCTK_ATTRIBUTE_UNUSED = PDstandardNth11(&gt23[index]);
    const CCTK_REAL PDstandardNth22gt23 CCTK_ATTRIBUTE_UNUSED = PDstandardNth22(&gt23[index]);
    const CCTK_REAL PDstandardNth33gt23 CCTK_ATTRIBUTE_UNUSED = PDstandardNth33(&gt23[index]);
    const CCTK_REAL PDstandardNth12gt23 CCTK_ATTRIBUTE_UNUSED = PDstandardNth12(&gt23[index]);
    const CCTK_REAL PDstandardNth13gt23 CCTK_ATTRIBUTE_UNUSED = PDstandardNth13(&gt23[index]);
    const CCTK_REAL PDstandardNth23gt23 CCTK_ATTRIBUTE_UNUSED = PDstandardNth23(&gt23[index]);
    const CCTK_REAL PDstandardNth1gt33 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&gt33[index]);
    const CCTK_REAL PDstandardNth2gt33 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&gt33[index]);
    const CCTK_REAL PDstandardNth3gt33 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&gt33[index]);
    const CCTK_REAL PDstandardNth11gt33 CCTK_ATTRIBUTE_UNUSED = PDstandardNth11(&gt33[index]);
    const CCTK_REAL PDstandardNth22gt33 CCTK_ATTRIBUTE_UNUSED = PDstandardNth22(&gt33[index]);
    const CCTK_REAL PDstandardNth33gt33 CCTK_ATTRIBUTE_UNUSED = PDstandardNth33(&gt33[index]);
    const CCTK_REAL PDstandardNth12gt33 CCTK_ATTRIBUTE_UNUSED = PDstandardNth12(&gt33[index]);
    const CCTK_REAL PDstandardNth13gt33 CCTK_ATTRIBUTE_UNUSED = PDstandardNth13(&gt33[index]);
    const CCTK_REAL PDstandardNth23gt33 CCTK_ATTRIBUTE_UNUSED = PDstandardNth23(&gt33[index]);
    const CCTK_REAL PDstandardNth1trK CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&trK[index]);
    const CCTK_REAL PDstandardNth2trK CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&trK[index]);
    const CCTK_REAL PDstandardNth3trK CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&trK[index]);
    const CCTK_REAL PDstandardNth1Xt1 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&Xt1[index]);
    const CCTK_REAL PDstandardNth2Xt1 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&Xt1[index]);
    const CCTK_REAL PDstandardNth3Xt1 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&Xt1[index]);
    const CCTK_REAL PDstandardNth1Xt2 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&Xt2[index]);
    const CCTK_REAL PDstandardNth2Xt2 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&Xt2[index]);
    const CCTK_REAL PDstandardNth3Xt2 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&Xt2[index]);
    const CCTK_REAL PDstandardNth1Xt3 CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&Xt3[index]);
    const CCTK_REAL PDstandardNth2Xt3 CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&Xt3[index]);
    const CCTK_REAL PDstandardNth3Xt3 CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&Xt3[index]);
    /* Calculate temporaries and grid functions */
    CCTK_REAL JacPDstandardNth11gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth11gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth11gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth11gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth11gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth11gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth12gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth12gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth12gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth12gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth12gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth12gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth13gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth13gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth13gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth13gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth13gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth13gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1trK CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth21gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth21gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth21gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth21gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth21gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth21gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth22gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth22gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth22gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth22gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth22gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth22gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth23gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth23gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth23gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth23gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth23gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth23gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2trK CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth31gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth31gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth31gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth31gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth31gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth31gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth32gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth32gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth32gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth32gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth32gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth32gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth33gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth33gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth33gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth33gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth33gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth33gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3trK CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3Xt3 CCTK_ATTRIBUTE_UNUSED;
    
    if (usejacobian)
    {
      JacPDstandardNth1gt11 = J11L*PDstandardNth1gt11 + 
        J21L*PDstandardNth2gt11 + J31L*PDstandardNth3gt11;
      
      JacPDstandardNth1gt12 = J11L*PDstandardNth1gt12 + 
        J21L*PDstandardNth2gt12 + J31L*PDstandardNth3gt12;
      
      JacPDstandardNth1gt13 = J11L*PDstandardNth1gt13 + 
        J21L*PDstandardNth2gt13 + J31L*PDstandardNth3gt13;
      
      JacPDstandardNth1gt22 = J11L*PDstandardNth1gt22 + 
        J21L*PDstandardNth2gt22 + J31L*PDstandardNth3gt22;
      
      JacPDstandardNth1gt23 = J11L*PDstandardNth1gt23 + 
        J21L*PDstandardNth2gt23 + J31L*PDstandardNth3gt23;
      
      JacPDstandardNth1gt33 = J11L*PDstandardNth1gt33 + 
        J21L*PDstandardNth2gt33 + J31L*PDstandardNth3gt33;
      
      JacPDstandardNth1trK = J11L*PDstandardNth1trK + J21L*PDstandardNth2trK 
        + J31L*PDstandardNth3trK;
      
      JacPDstandardNth1Xt1 = J11L*PDstandardNth1Xt1 + J21L*PDstandardNth2Xt1 
        + J31L*PDstandardNth3Xt1;
      
      JacPDstandardNth1Xt2 = J11L*PDstandardNth1Xt2 + J21L*PDstandardNth2Xt2 
        + J31L*PDstandardNth3Xt2;
      
      JacPDstandardNth1Xt3 = J11L*PDstandardNth1Xt3 + J21L*PDstandardNth2Xt3 
        + J31L*PDstandardNth3Xt3;
      
      JacPDstandardNth2gt11 = J12L*PDstandardNth1gt11 + 
        J22L*PDstandardNth2gt11 + J32L*PDstandardNth3gt11;
      
      JacPDstandardNth2gt12 = J12L*PDstandardNth1gt12 + 
        J22L*PDstandardNth2gt12 + J32L*PDstandardNth3gt12;
      
      JacPDstandardNth2gt13 = J12L*PDstandardNth1gt13 + 
        J22L*PDstandardNth2gt13 + J32L*PDstandardNth3gt13;
      
      JacPDstandardNth2gt22 = J12L*PDstandardNth1gt22 + 
        J22L*PDstandardNth2gt22 + J32L*PDstandardNth3gt22;
      
      JacPDstandardNth2gt23 = J12L*PDstandardNth1gt23 + 
        J22L*PDstandardNth2gt23 + J32L*PDstandardNth3gt23;
      
      JacPDstandardNth2gt33 = J12L*PDstandardNth1gt33 + 
        J22L*PDstandardNth2gt33 + J32L*PDstandardNth3gt33;
      
      JacPDstandardNth2trK = J12L*PDstandardNth1trK + J22L*PDstandardNth2trK 
        + J32L*PDstandardNth3trK;
      
      JacPDstandardNth2Xt1 = J12L*PDstandardNth1Xt1 + J22L*PDstandardNth2Xt1 
        + J32L*PDstandardNth3Xt1;
      
      JacPDstandardNth2Xt2 = J12L*PDstandardNth1Xt2 + J22L*PDstandardNth2Xt2 
        + J32L*PDstandardNth3Xt2;
      
      JacPDstandardNth2Xt3 = J12L*PDstandardNth1Xt3 + J22L*PDstandardNth2Xt3 
        + J32L*PDstandardNth3Xt3;
      
      JacPDstandardNth3gt11 = J13L*PDstandardNth1gt11 + 
        J23L*PDstandardNth2gt11 + J33L*PDstandardNth3gt11;
      
      JacPDstandardNth3gt12 = J13L*PDstandardNth1gt12 + 
        J23L*PDstandardNth2gt12 + J33L*PDstandardNth3gt12;
      
      JacPDstandardNth3gt13 = J13L*PDstandardNth1gt13 + 
        J23L*PDstandardNth2gt13 + J33L*PDstandardNth3gt13;
      
      JacPDstandardNth3gt22 = J13L*PDstandardNth1gt22 + 
        J23L*PDstandardNth2gt22 + J33L*PDstandardNth3gt22;
      
      JacPDstandardNth3gt23 = J13L*PDstandardNth1gt23 + 
        J23L*PDstandardNth2gt23 + J33L*PDstandardNth3gt23;
      
      JacPDstandardNth3gt33 = J13L*PDstandardNth1gt33 + 
        J23L*PDstandardNth2gt33 + J33L*PDstandardNth3gt33;
      
      JacPDstandardNth3trK = J13L*PDstandardNth1trK + J23L*PDstandardNth2trK 
        + J33L*PDstandardNth3trK;
      
      JacPDstandardNth3Xt1 = J13L*PDstandardNth1Xt1 + J23L*PDstandardNth2Xt1 
        + J33L*PDstandardNth3Xt1;
      
      JacPDstandardNth3Xt2 = J13L*PDstandardNth1Xt2 + J23L*PDstandardNth2Xt2 
        + J33L*PDstandardNth3Xt2;
      
      JacPDstandardNth3Xt3 = J13L*PDstandardNth1Xt3 + J23L*PDstandardNth2Xt3 
        + J33L*PDstandardNth3Xt3;
      
      JacPDstandardNth11gt11 = dJ111L*PDstandardNth1gt11 + 
        2*(J11L*(J21L*PDstandardNth12gt11 + J31L*PDstandardNth13gt11) + 
        J21L*J31L*PDstandardNth23gt11) + dJ211L*PDstandardNth2gt11 + 
        dJ311L*PDstandardNth3gt11 + PDstandardNth11gt11*pow(J11L,2) + 
        PDstandardNth22gt11*pow(J21L,2) + PDstandardNth33gt11*pow(J31L,2);
      
      JacPDstandardNth11gt12 = dJ111L*PDstandardNth1gt12 + 
        2*(J11L*(J21L*PDstandardNth12gt12 + J31L*PDstandardNth13gt12) + 
        J21L*J31L*PDstandardNth23gt12) + dJ211L*PDstandardNth2gt12 + 
        dJ311L*PDstandardNth3gt12 + PDstandardNth11gt12*pow(J11L,2) + 
        PDstandardNth22gt12*pow(J21L,2) + PDstandardNth33gt12*pow(J31L,2);
      
      JacPDstandardNth11gt13 = dJ111L*PDstandardNth1gt13 + 
        2*(J11L*(J21L*PDstandardNth12gt13 + J31L*PDstandardNth13gt13) + 
        J21L*J31L*PDstandardNth23gt13) + dJ211L*PDstandardNth2gt13 + 
        dJ311L*PDstandardNth3gt13 + PDstandardNth11gt13*pow(J11L,2) + 
        PDstandardNth22gt13*pow(J21L,2) + PDstandardNth33gt13*pow(J31L,2);
      
      JacPDstandardNth11gt22 = dJ111L*PDstandardNth1gt22 + 
        2*(J11L*(J21L*PDstandardNth12gt22 + J31L*PDstandardNth13gt22) + 
        J21L*J31L*PDstandardNth23gt22) + dJ211L*PDstandardNth2gt22 + 
        dJ311L*PDstandardNth3gt22 + PDstandardNth11gt22*pow(J11L,2) + 
        PDstandardNth22gt22*pow(J21L,2) + PDstandardNth33gt22*pow(J31L,2);
      
      JacPDstandardNth11gt23 = dJ111L*PDstandardNth1gt23 + 
        2*(J11L*(J21L*PDstandardNth12gt23 + J31L*PDstandardNth13gt23) + 
        J21L*J31L*PDstandardNth23gt23) + dJ211L*PDstandardNth2gt23 + 
        dJ311L*PDstandardNth3gt23 + PDstandardNth11gt23*pow(J11L,2) + 
        PDstandardNth22gt23*pow(J21L,2) + PDstandardNth33gt23*pow(J31L,2);
      
      JacPDstandardNth11gt33 = dJ111L*PDstandardNth1gt33 + 
        2*(J11L*(J21L*PDstandardNth12gt33 + J31L*PDstandardNth13gt33) + 
        J21L*J31L*PDstandardNth23gt33) + dJ211L*PDstandardNth2gt33 + 
        dJ311L*PDstandardNth3gt33 + PDstandardNth11gt33*pow(J11L,2) + 
        PDstandardNth22gt33*pow(J21L,2) + PDstandardNth33gt33*pow(J31L,2);
      
      JacPDstandardNth22gt11 = dJ122L*PDstandardNth1gt11 + 
        2*(J12L*(J22L*PDstandardNth12gt11 + J32L*PDstandardNth13gt11) + 
        J22L*J32L*PDstandardNth23gt11) + dJ222L*PDstandardNth2gt11 + 
        dJ322L*PDstandardNth3gt11 + PDstandardNth11gt11*pow(J12L,2) + 
        PDstandardNth22gt11*pow(J22L,2) + PDstandardNth33gt11*pow(J32L,2);
      
      JacPDstandardNth22gt12 = dJ122L*PDstandardNth1gt12 + 
        2*(J12L*(J22L*PDstandardNth12gt12 + J32L*PDstandardNth13gt12) + 
        J22L*J32L*PDstandardNth23gt12) + dJ222L*PDstandardNth2gt12 + 
        dJ322L*PDstandardNth3gt12 + PDstandardNth11gt12*pow(J12L,2) + 
        PDstandardNth22gt12*pow(J22L,2) + PDstandardNth33gt12*pow(J32L,2);
      
      JacPDstandardNth22gt13 = dJ122L*PDstandardNth1gt13 + 
        2*(J12L*(J22L*PDstandardNth12gt13 + J32L*PDstandardNth13gt13) + 
        J22L*J32L*PDstandardNth23gt13) + dJ222L*PDstandardNth2gt13 + 
        dJ322L*PDstandardNth3gt13 + PDstandardNth11gt13*pow(J12L,2) + 
        PDstandardNth22gt13*pow(J22L,2) + PDstandardNth33gt13*pow(J32L,2);
      
      JacPDstandardNth22gt22 = dJ122L*PDstandardNth1gt22 + 
        2*(J12L*(J22L*PDstandardNth12gt22 + J32L*PDstandardNth13gt22) + 
        J22L*J32L*PDstandardNth23gt22) + dJ222L*PDstandardNth2gt22 + 
        dJ322L*PDstandardNth3gt22 + PDstandardNth11gt22*pow(J12L,2) + 
        PDstandardNth22gt22*pow(J22L,2) + PDstandardNth33gt22*pow(J32L,2);
      
      JacPDstandardNth22gt23 = dJ122L*PDstandardNth1gt23 + 
        2*(J12L*(J22L*PDstandardNth12gt23 + J32L*PDstandardNth13gt23) + 
        J22L*J32L*PDstandardNth23gt23) + dJ222L*PDstandardNth2gt23 + 
        dJ322L*PDstandardNth3gt23 + PDstandardNth11gt23*pow(J12L,2) + 
        PDstandardNth22gt23*pow(J22L,2) + PDstandardNth33gt23*pow(J32L,2);
      
      JacPDstandardNth22gt33 = dJ122L*PDstandardNth1gt33 + 
        2*(J12L*(J22L*PDstandardNth12gt33 + J32L*PDstandardNth13gt33) + 
        J22L*J32L*PDstandardNth23gt33) + dJ222L*PDstandardNth2gt33 + 
        dJ322L*PDstandardNth3gt33 + PDstandardNth11gt33*pow(J12L,2) + 
        PDstandardNth22gt33*pow(J22L,2) + PDstandardNth33gt33*pow(J32L,2);
      
      JacPDstandardNth33gt11 = dJ133L*PDstandardNth1gt11 + 
        2*(J13L*(J23L*PDstandardNth12gt11 + J33L*PDstandardNth13gt11) + 
        J23L*J33L*PDstandardNth23gt11) + dJ233L*PDstandardNth2gt11 + 
        dJ333L*PDstandardNth3gt11 + PDstandardNth11gt11*pow(J13L,2) + 
        PDstandardNth22gt11*pow(J23L,2) + PDstandardNth33gt11*pow(J33L,2);
      
      JacPDstandardNth33gt12 = dJ133L*PDstandardNth1gt12 + 
        2*(J13L*(J23L*PDstandardNth12gt12 + J33L*PDstandardNth13gt12) + 
        J23L*J33L*PDstandardNth23gt12) + dJ233L*PDstandardNth2gt12 + 
        dJ333L*PDstandardNth3gt12 + PDstandardNth11gt12*pow(J13L,2) + 
        PDstandardNth22gt12*pow(J23L,2) + PDstandardNth33gt12*pow(J33L,2);
      
      JacPDstandardNth33gt13 = dJ133L*PDstandardNth1gt13 + 
        2*(J13L*(J23L*PDstandardNth12gt13 + J33L*PDstandardNth13gt13) + 
        J23L*J33L*PDstandardNth23gt13) + dJ233L*PDstandardNth2gt13 + 
        dJ333L*PDstandardNth3gt13 + PDstandardNth11gt13*pow(J13L,2) + 
        PDstandardNth22gt13*pow(J23L,2) + PDstandardNth33gt13*pow(J33L,2);
      
      JacPDstandardNth33gt22 = dJ133L*PDstandardNth1gt22 + 
        2*(J13L*(J23L*PDstandardNth12gt22 + J33L*PDstandardNth13gt22) + 
        J23L*J33L*PDstandardNth23gt22) + dJ233L*PDstandardNth2gt22 + 
        dJ333L*PDstandardNth3gt22 + PDstandardNth11gt22*pow(J13L,2) + 
        PDstandardNth22gt22*pow(J23L,2) + PDstandardNth33gt22*pow(J33L,2);
      
      JacPDstandardNth33gt23 = dJ133L*PDstandardNth1gt23 + 
        2*(J13L*(J23L*PDstandardNth12gt23 + J33L*PDstandardNth13gt23) + 
        J23L*J33L*PDstandardNth23gt23) + dJ233L*PDstandardNth2gt23 + 
        dJ333L*PDstandardNth3gt23 + PDstandardNth11gt23*pow(J13L,2) + 
        PDstandardNth22gt23*pow(J23L,2) + PDstandardNth33gt23*pow(J33L,2);
      
      JacPDstandardNth33gt33 = dJ133L*PDstandardNth1gt33 + 
        2*(J13L*(J23L*PDstandardNth12gt33 + J33L*PDstandardNth13gt33) + 
        J23L*J33L*PDstandardNth23gt33) + dJ233L*PDstandardNth2gt33 + 
        dJ333L*PDstandardNth3gt33 + PDstandardNth11gt33*pow(J13L,2) + 
        PDstandardNth22gt33*pow(J23L,2) + PDstandardNth33gt33*pow(J33L,2);
      
      JacPDstandardNth12gt11 = J12L*(J11L*PDstandardNth11gt11 + 
        J21L*PDstandardNth12gt11 + J31L*PDstandardNth13gt11) + 
        J11L*(J22L*PDstandardNth12gt11 + J32L*PDstandardNth13gt11) + 
        dJ112L*PDstandardNth1gt11 + J22L*(J21L*PDstandardNth22gt11 + 
        J31L*PDstandardNth23gt11) + dJ212L*PDstandardNth2gt11 + 
        J32L*(J21L*PDstandardNth23gt11 + J31L*PDstandardNth33gt11) + 
        dJ312L*PDstandardNth3gt11;
      
      JacPDstandardNth12gt12 = J12L*(J11L*PDstandardNth11gt12 + 
        J21L*PDstandardNth12gt12 + J31L*PDstandardNth13gt12) + 
        J11L*(J22L*PDstandardNth12gt12 + J32L*PDstandardNth13gt12) + 
        dJ112L*PDstandardNth1gt12 + J22L*(J21L*PDstandardNth22gt12 + 
        J31L*PDstandardNth23gt12) + dJ212L*PDstandardNth2gt12 + 
        J32L*(J21L*PDstandardNth23gt12 + J31L*PDstandardNth33gt12) + 
        dJ312L*PDstandardNth3gt12;
      
      JacPDstandardNth12gt13 = J12L*(J11L*PDstandardNth11gt13 + 
        J21L*PDstandardNth12gt13 + J31L*PDstandardNth13gt13) + 
        J11L*(J22L*PDstandardNth12gt13 + J32L*PDstandardNth13gt13) + 
        dJ112L*PDstandardNth1gt13 + J22L*(J21L*PDstandardNth22gt13 + 
        J31L*PDstandardNth23gt13) + dJ212L*PDstandardNth2gt13 + 
        J32L*(J21L*PDstandardNth23gt13 + J31L*PDstandardNth33gt13) + 
        dJ312L*PDstandardNth3gt13;
      
      JacPDstandardNth12gt22 = J12L*(J11L*PDstandardNth11gt22 + 
        J21L*PDstandardNth12gt22 + J31L*PDstandardNth13gt22) + 
        J11L*(J22L*PDstandardNth12gt22 + J32L*PDstandardNth13gt22) + 
        dJ112L*PDstandardNth1gt22 + J22L*(J21L*PDstandardNth22gt22 + 
        J31L*PDstandardNth23gt22) + dJ212L*PDstandardNth2gt22 + 
        J32L*(J21L*PDstandardNth23gt22 + J31L*PDstandardNth33gt22) + 
        dJ312L*PDstandardNth3gt22;
      
      JacPDstandardNth12gt23 = J12L*(J11L*PDstandardNth11gt23 + 
        J21L*PDstandardNth12gt23 + J31L*PDstandardNth13gt23) + 
        J11L*(J22L*PDstandardNth12gt23 + J32L*PDstandardNth13gt23) + 
        dJ112L*PDstandardNth1gt23 + J22L*(J21L*PDstandardNth22gt23 + 
        J31L*PDstandardNth23gt23) + dJ212L*PDstandardNth2gt23 + 
        J32L*(J21L*PDstandardNth23gt23 + J31L*PDstandardNth33gt23) + 
        dJ312L*PDstandardNth3gt23;
      
      JacPDstandardNth12gt33 = J12L*(J11L*PDstandardNth11gt33 + 
        J21L*PDstandardNth12gt33 + J31L*PDstandardNth13gt33) + 
        J11L*(J22L*PDstandardNth12gt33 + J32L*PDstandardNth13gt33) + 
        dJ112L*PDstandardNth1gt33 + J22L*(J21L*PDstandardNth22gt33 + 
        J31L*PDstandardNth23gt33) + dJ212L*PDstandardNth2gt33 + 
        J32L*(J21L*PDstandardNth23gt33 + J31L*PDstandardNth33gt33) + 
        dJ312L*PDstandardNth3gt33;
      
      JacPDstandardNth13gt11 = J13L*(J11L*PDstandardNth11gt11 + 
        J21L*PDstandardNth12gt11 + J31L*PDstandardNth13gt11) + 
        J11L*(J23L*PDstandardNth12gt11 + J33L*PDstandardNth13gt11) + 
        dJ113L*PDstandardNth1gt11 + J23L*(J21L*PDstandardNth22gt11 + 
        J31L*PDstandardNth23gt11) + dJ213L*PDstandardNth2gt11 + 
        J33L*(J21L*PDstandardNth23gt11 + J31L*PDstandardNth33gt11) + 
        dJ313L*PDstandardNth3gt11;
      
      JacPDstandardNth13gt12 = J13L*(J11L*PDstandardNth11gt12 + 
        J21L*PDstandardNth12gt12 + J31L*PDstandardNth13gt12) + 
        J11L*(J23L*PDstandardNth12gt12 + J33L*PDstandardNth13gt12) + 
        dJ113L*PDstandardNth1gt12 + J23L*(J21L*PDstandardNth22gt12 + 
        J31L*PDstandardNth23gt12) + dJ213L*PDstandardNth2gt12 + 
        J33L*(J21L*PDstandardNth23gt12 + J31L*PDstandardNth33gt12) + 
        dJ313L*PDstandardNth3gt12;
      
      JacPDstandardNth13gt13 = J13L*(J11L*PDstandardNth11gt13 + 
        J21L*PDstandardNth12gt13 + J31L*PDstandardNth13gt13) + 
        J11L*(J23L*PDstandardNth12gt13 + J33L*PDstandardNth13gt13) + 
        dJ113L*PDstandardNth1gt13 + J23L*(J21L*PDstandardNth22gt13 + 
        J31L*PDstandardNth23gt13) + dJ213L*PDstandardNth2gt13 + 
        J33L*(J21L*PDstandardNth23gt13 + J31L*PDstandardNth33gt13) + 
        dJ313L*PDstandardNth3gt13;
      
      JacPDstandardNth13gt22 = J13L*(J11L*PDstandardNth11gt22 + 
        J21L*PDstandardNth12gt22 + J31L*PDstandardNth13gt22) + 
        J11L*(J23L*PDstandardNth12gt22 + J33L*PDstandardNth13gt22) + 
        dJ113L*PDstandardNth1gt22 + J23L*(J21L*PDstandardNth22gt22 + 
        J31L*PDstandardNth23gt22) + dJ213L*PDstandardNth2gt22 + 
        J33L*(J21L*PDstandardNth23gt22 + J31L*PDstandardNth33gt22) + 
        dJ313L*PDstandardNth3gt22;
      
      JacPDstandardNth13gt23 = J13L*(J11L*PDstandardNth11gt23 + 
        J21L*PDstandardNth12gt23 + J31L*PDstandardNth13gt23) + 
        J11L*(J23L*PDstandardNth12gt23 + J33L*PDstandardNth13gt23) + 
        dJ113L*PDstandardNth1gt23 + J23L*(J21L*PDstandardNth22gt23 + 
        J31L*PDstandardNth23gt23) + dJ213L*PDstandardNth2gt23 + 
        J33L*(J21L*PDstandardNth23gt23 + J31L*PDstandardNth33gt23) + 
        dJ313L*PDstandardNth3gt23;
      
      JacPDstandardNth13gt33 = J13L*(J11L*PDstandardNth11gt33 + 
        J21L*PDstandardNth12gt33 + J31L*PDstandardNth13gt33) + 
        J11L*(J23L*PDstandardNth12gt33 + J33L*PDstandardNth13gt33) + 
        dJ113L*PDstandardNth1gt33 + J23L*(J21L*PDstandardNth22gt33 + 
        J31L*PDstandardNth23gt33) + dJ213L*PDstandardNth2gt33 + 
        J33L*(J21L*PDstandardNth23gt33 + J31L*PDstandardNth33gt33) + 
        dJ313L*PDstandardNth3gt33;
      
      JacPDstandardNth21gt11 = J12L*(J11L*PDstandardNth11gt11 + 
        J21L*PDstandardNth12gt11 + J31L*PDstandardNth13gt11) + 
        J11L*(J22L*PDstandardNth12gt11 + J32L*PDstandardNth13gt11) + 
        dJ112L*PDstandardNth1gt11 + J22L*(J21L*PDstandardNth22gt11 + 
        J31L*PDstandardNth23gt11) + dJ212L*PDstandardNth2gt11 + 
        J32L*(J21L*PDstandardNth23gt11 + J31L*PDstandardNth33gt11) + 
        dJ312L*PDstandardNth3gt11;
      
      JacPDstandardNth21gt12 = J12L*(J11L*PDstandardNth11gt12 + 
        J21L*PDstandardNth12gt12 + J31L*PDstandardNth13gt12) + 
        J11L*(J22L*PDstandardNth12gt12 + J32L*PDstandardNth13gt12) + 
        dJ112L*PDstandardNth1gt12 + J22L*(J21L*PDstandardNth22gt12 + 
        J31L*PDstandardNth23gt12) + dJ212L*PDstandardNth2gt12 + 
        J32L*(J21L*PDstandardNth23gt12 + J31L*PDstandardNth33gt12) + 
        dJ312L*PDstandardNth3gt12;
      
      JacPDstandardNth21gt13 = J12L*(J11L*PDstandardNth11gt13 + 
        J21L*PDstandardNth12gt13 + J31L*PDstandardNth13gt13) + 
        J11L*(J22L*PDstandardNth12gt13 + J32L*PDstandardNth13gt13) + 
        dJ112L*PDstandardNth1gt13 + J22L*(J21L*PDstandardNth22gt13 + 
        J31L*PDstandardNth23gt13) + dJ212L*PDstandardNth2gt13 + 
        J32L*(J21L*PDstandardNth23gt13 + J31L*PDstandardNth33gt13) + 
        dJ312L*PDstandardNth3gt13;
      
      JacPDstandardNth21gt22 = J12L*(J11L*PDstandardNth11gt22 + 
        J21L*PDstandardNth12gt22 + J31L*PDstandardNth13gt22) + 
        J11L*(J22L*PDstandardNth12gt22 + J32L*PDstandardNth13gt22) + 
        dJ112L*PDstandardNth1gt22 + J22L*(J21L*PDstandardNth22gt22 + 
        J31L*PDstandardNth23gt22) + dJ212L*PDstandardNth2gt22 + 
        J32L*(J21L*PDstandardNth23gt22 + J31L*PDstandardNth33gt22) + 
        dJ312L*PDstandardNth3gt22;
      
      JacPDstandardNth21gt23 = J12L*(J11L*PDstandardNth11gt23 + 
        J21L*PDstandardNth12gt23 + J31L*PDstandardNth13gt23) + 
        J11L*(J22L*PDstandardNth12gt23 + J32L*PDstandardNth13gt23) + 
        dJ112L*PDstandardNth1gt23 + J22L*(J21L*PDstandardNth22gt23 + 
        J31L*PDstandardNth23gt23) + dJ212L*PDstandardNth2gt23 + 
        J32L*(J21L*PDstandardNth23gt23 + J31L*PDstandardNth33gt23) + 
        dJ312L*PDstandardNth3gt23;
      
      JacPDstandardNth21gt33 = J12L*(J11L*PDstandardNth11gt33 + 
        J21L*PDstandardNth12gt33 + J31L*PDstandardNth13gt33) + 
        J11L*(J22L*PDstandardNth12gt33 + J32L*PDstandardNth13gt33) + 
        dJ112L*PDstandardNth1gt33 + J22L*(J21L*PDstandardNth22gt33 + 
        J31L*PDstandardNth23gt33) + dJ212L*PDstandardNth2gt33 + 
        J32L*(J21L*PDstandardNth23gt33 + J31L*PDstandardNth33gt33) + 
        dJ312L*PDstandardNth3gt33;
      
      JacPDstandardNth23gt11 = J13L*(J12L*PDstandardNth11gt11 + 
        J22L*PDstandardNth12gt11 + J32L*PDstandardNth13gt11) + 
        J12L*(J23L*PDstandardNth12gt11 + J33L*PDstandardNth13gt11) + 
        dJ123L*PDstandardNth1gt11 + J23L*(J22L*PDstandardNth22gt11 + 
        J32L*PDstandardNth23gt11) + dJ223L*PDstandardNth2gt11 + 
        J33L*(J22L*PDstandardNth23gt11 + J32L*PDstandardNth33gt11) + 
        dJ323L*PDstandardNth3gt11;
      
      JacPDstandardNth23gt12 = J13L*(J12L*PDstandardNth11gt12 + 
        J22L*PDstandardNth12gt12 + J32L*PDstandardNth13gt12) + 
        J12L*(J23L*PDstandardNth12gt12 + J33L*PDstandardNth13gt12) + 
        dJ123L*PDstandardNth1gt12 + J23L*(J22L*PDstandardNth22gt12 + 
        J32L*PDstandardNth23gt12) + dJ223L*PDstandardNth2gt12 + 
        J33L*(J22L*PDstandardNth23gt12 + J32L*PDstandardNth33gt12) + 
        dJ323L*PDstandardNth3gt12;
      
      JacPDstandardNth23gt13 = J13L*(J12L*PDstandardNth11gt13 + 
        J22L*PDstandardNth12gt13 + J32L*PDstandardNth13gt13) + 
        J12L*(J23L*PDstandardNth12gt13 + J33L*PDstandardNth13gt13) + 
        dJ123L*PDstandardNth1gt13 + J23L*(J22L*PDstandardNth22gt13 + 
        J32L*PDstandardNth23gt13) + dJ223L*PDstandardNth2gt13 + 
        J33L*(J22L*PDstandardNth23gt13 + J32L*PDstandardNth33gt13) + 
        dJ323L*PDstandardNth3gt13;
      
      JacPDstandardNth23gt22 = J13L*(J12L*PDstandardNth11gt22 + 
        J22L*PDstandardNth12gt22 + J32L*PDstandardNth13gt22) + 
        J12L*(J23L*PDstandardNth12gt22 + J33L*PDstandardNth13gt22) + 
        dJ123L*PDstandardNth1gt22 + J23L*(J22L*PDstandardNth22gt22 + 
        J32L*PDstandardNth23gt22) + dJ223L*PDstandardNth2gt22 + 
        J33L*(J22L*PDstandardNth23gt22 + J32L*PDstandardNth33gt22) + 
        dJ323L*PDstandardNth3gt22;
      
      JacPDstandardNth23gt23 = J13L*(J12L*PDstandardNth11gt23 + 
        J22L*PDstandardNth12gt23 + J32L*PDstandardNth13gt23) + 
        J12L*(J23L*PDstandardNth12gt23 + J33L*PDstandardNth13gt23) + 
        dJ123L*PDstandardNth1gt23 + J23L*(J22L*PDstandardNth22gt23 + 
        J32L*PDstandardNth23gt23) + dJ223L*PDstandardNth2gt23 + 
        J33L*(J22L*PDstandardNth23gt23 + J32L*PDstandardNth33gt23) + 
        dJ323L*PDstandardNth3gt23;
      
      JacPDstandardNth23gt33 = J13L*(J12L*PDstandardNth11gt33 + 
        J22L*PDstandardNth12gt33 + J32L*PDstandardNth13gt33) + 
        J12L*(J23L*PDstandardNth12gt33 + J33L*PDstandardNth13gt33) + 
        dJ123L*PDstandardNth1gt33 + J23L*(J22L*PDstandardNth22gt33 + 
        J32L*PDstandardNth23gt33) + dJ223L*PDstandardNth2gt33 + 
        J33L*(J22L*PDstandardNth23gt33 + J32L*PDstandardNth33gt33) + 
        dJ323L*PDstandardNth3gt33;
      
      JacPDstandardNth31gt11 = J13L*(J11L*PDstandardNth11gt11 + 
        J21L*PDstandardNth12gt11 + J31L*PDstandardNth13gt11) + 
        J11L*(J23L*PDstandardNth12gt11 + J33L*PDstandardNth13gt11) + 
        dJ113L*PDstandardNth1gt11 + J23L*(J21L*PDstandardNth22gt11 + 
        J31L*PDstandardNth23gt11) + dJ213L*PDstandardNth2gt11 + 
        J33L*(J21L*PDstandardNth23gt11 + J31L*PDstandardNth33gt11) + 
        dJ313L*PDstandardNth3gt11;
      
      JacPDstandardNth31gt12 = J13L*(J11L*PDstandardNth11gt12 + 
        J21L*PDstandardNth12gt12 + J31L*PDstandardNth13gt12) + 
        J11L*(J23L*PDstandardNth12gt12 + J33L*PDstandardNth13gt12) + 
        dJ113L*PDstandardNth1gt12 + J23L*(J21L*PDstandardNth22gt12 + 
        J31L*PDstandardNth23gt12) + dJ213L*PDstandardNth2gt12 + 
        J33L*(J21L*PDstandardNth23gt12 + J31L*PDstandardNth33gt12) + 
        dJ313L*PDstandardNth3gt12;
      
      JacPDstandardNth31gt13 = J13L*(J11L*PDstandardNth11gt13 + 
        J21L*PDstandardNth12gt13 + J31L*PDstandardNth13gt13) + 
        J11L*(J23L*PDstandardNth12gt13 + J33L*PDstandardNth13gt13) + 
        dJ113L*PDstandardNth1gt13 + J23L*(J21L*PDstandardNth22gt13 + 
        J31L*PDstandardNth23gt13) + dJ213L*PDstandardNth2gt13 + 
        J33L*(J21L*PDstandardNth23gt13 + J31L*PDstandardNth33gt13) + 
        dJ313L*PDstandardNth3gt13;
      
      JacPDstandardNth31gt22 = J13L*(J11L*PDstandardNth11gt22 + 
        J21L*PDstandardNth12gt22 + J31L*PDstandardNth13gt22) + 
        J11L*(J23L*PDstandardNth12gt22 + J33L*PDstandardNth13gt22) + 
        dJ113L*PDstandardNth1gt22 + J23L*(J21L*PDstandardNth22gt22 + 
        J31L*PDstandardNth23gt22) + dJ213L*PDstandardNth2gt22 + 
        J33L*(J21L*PDstandardNth23gt22 + J31L*PDstandardNth33gt22) + 
        dJ313L*PDstandardNth3gt22;
      
      JacPDstandardNth31gt23 = J13L*(J11L*PDstandardNth11gt23 + 
        J21L*PDstandardNth12gt23 + J31L*PDstandardNth13gt23) + 
        J11L*(J23L*PDstandardNth12gt23 + J33L*PDstandardNth13gt23) + 
        dJ113L*PDstandardNth1gt23 + J23L*(J21L*PDstandardNth22gt23 + 
        J31L*PDstandardNth23gt23) + dJ213L*PDstandardNth2gt23 + 
        J33L*(J21L*PDstandardNth23gt23 + J31L*PDstandardNth33gt23) + 
        dJ313L*PDstandardNth3gt23;
      
      JacPDstandardNth31gt33 = J13L*(J11L*PDstandardNth11gt33 + 
        J21L*PDstandardNth12gt33 + J31L*PDstandardNth13gt33) + 
        J11L*(J23L*PDstandardNth12gt33 + J33L*PDstandardNth13gt33) + 
        dJ113L*PDstandardNth1gt33 + J23L*(J21L*PDstandardNth22gt33 + 
        J31L*PDstandardNth23gt33) + dJ213L*PDstandardNth2gt33 + 
        J33L*(J21L*PDstandardNth23gt33 + J31L*PDstandardNth33gt33) + 
        dJ313L*PDstandardNth3gt33;
      
      JacPDstandardNth32gt11 = J13L*(J12L*PDstandardNth11gt11 + 
        J22L*PDstandardNth12gt11 + J32L*PDstandardNth13gt11) + 
        J12L*(J23L*PDstandardNth12gt11 + J33L*PDstandardNth13gt11) + 
        dJ123L*PDstandardNth1gt11 + J23L*(J22L*PDstandardNth22gt11 + 
        J32L*PDstandardNth23gt11) + dJ223L*PDstandardNth2gt11 + 
        J33L*(J22L*PDstandardNth23gt11 + J32L*PDstandardNth33gt11) + 
        dJ323L*PDstandardNth3gt11;
      
      JacPDstandardNth32gt12 = J13L*(J12L*PDstandardNth11gt12 + 
        J22L*PDstandardNth12gt12 + J32L*PDstandardNth13gt12) + 
        J12L*(J23L*PDstandardNth12gt12 + J33L*PDstandardNth13gt12) + 
        dJ123L*PDstandardNth1gt12 + J23L*(J22L*PDstandardNth22gt12 + 
        J32L*PDstandardNth23gt12) + dJ223L*PDstandardNth2gt12 + 
        J33L*(J22L*PDstandardNth23gt12 + J32L*PDstandardNth33gt12) + 
        dJ323L*PDstandardNth3gt12;
      
      JacPDstandardNth32gt13 = J13L*(J12L*PDstandardNth11gt13 + 
        J22L*PDstandardNth12gt13 + J32L*PDstandardNth13gt13) + 
        J12L*(J23L*PDstandardNth12gt13 + J33L*PDstandardNth13gt13) + 
        dJ123L*PDstandardNth1gt13 + J23L*(J22L*PDstandardNth22gt13 + 
        J32L*PDstandardNth23gt13) + dJ223L*PDstandardNth2gt13 + 
        J33L*(J22L*PDstandardNth23gt13 + J32L*PDstandardNth33gt13) + 
        dJ323L*PDstandardNth3gt13;
      
      JacPDstandardNth32gt22 = J13L*(J12L*PDstandardNth11gt22 + 
        J22L*PDstandardNth12gt22 + J32L*PDstandardNth13gt22) + 
        J12L*(J23L*PDstandardNth12gt22 + J33L*PDstandardNth13gt22) + 
        dJ123L*PDstandardNth1gt22 + J23L*(J22L*PDstandardNth22gt22 + 
        J32L*PDstandardNth23gt22) + dJ223L*PDstandardNth2gt22 + 
        J33L*(J22L*PDstandardNth23gt22 + J32L*PDstandardNth33gt22) + 
        dJ323L*PDstandardNth3gt22;
      
      JacPDstandardNth32gt23 = J13L*(J12L*PDstandardNth11gt23 + 
        J22L*PDstandardNth12gt23 + J32L*PDstandardNth13gt23) + 
        J12L*(J23L*PDstandardNth12gt23 + J33L*PDstandardNth13gt23) + 
        dJ123L*PDstandardNth1gt23 + J23L*(J22L*PDstandardNth22gt23 + 
        J32L*PDstandardNth23gt23) + dJ223L*PDstandardNth2gt23 + 
        J33L*(J22L*PDstandardNth23gt23 + J32L*PDstandardNth33gt23) + 
        dJ323L*PDstandardNth3gt23;
      
      JacPDstandardNth32gt33 = J13L*(J12L*PDstandardNth11gt33 + 
        J22L*PDstandardNth12gt33 + J32L*PDstandardNth13gt33) + 
        J12L*(J23L*PDstandardNth12gt33 + J33L*PDstandardNth13gt33) + 
        dJ123L*PDstandardNth1gt33 + J23L*(J22L*PDstandardNth22gt33 + 
        J32L*PDstandardNth23gt33) + dJ223L*PDstandardNth2gt33 + 
        J33L*(J22L*PDstandardNth23gt33 + J32L*PDstandardNth33gt33) + 
        dJ323L*PDstandardNth3gt33;
    }
    else
    {
      JacPDstandardNth1gt11 = PDstandardNth1gt11;
      
      JacPDstandardNth1gt12 = PDstandardNth1gt12;
      
      JacPDstandardNth1gt13 = PDstandardNth1gt13;
      
      JacPDstandardNth1gt22 = PDstandardNth1gt22;
      
      JacPDstandardNth1gt23 = PDstandardNth1gt23;
      
      JacPDstandardNth1gt33 = PDstandardNth1gt33;
      
      JacPDstandardNth1trK = PDstandardNth1trK;
      
      JacPDstandardNth1Xt1 = PDstandardNth1Xt1;
      
      JacPDstandardNth1Xt2 = PDstandardNth1Xt2;
      
      JacPDstandardNth1Xt3 = PDstandardNth1Xt3;
      
      JacPDstandardNth2gt11 = PDstandardNth2gt11;
      
      JacPDstandardNth2gt12 = PDstandardNth2gt12;
      
      JacPDstandardNth2gt13 = PDstandardNth2gt13;
      
      JacPDstandardNth2gt22 = PDstandardNth2gt22;
      
      JacPDstandardNth2gt23 = PDstandardNth2gt23;
      
      JacPDstandardNth2gt33 = PDstandardNth2gt33;
      
      JacPDstandardNth2trK = PDstandardNth2trK;
      
      JacPDstandardNth2Xt1 = PDstandardNth2Xt1;
      
      JacPDstandardNth2Xt2 = PDstandardNth2Xt2;
      
      JacPDstandardNth2Xt3 = PDstandardNth2Xt3;
      
      JacPDstandardNth3gt11 = PDstandardNth3gt11;
      
      JacPDstandardNth3gt12 = PDstandardNth3gt12;
      
      JacPDstandardNth3gt13 = PDstandardNth3gt13;
      
      JacPDstandardNth3gt22 = PDstandardNth3gt22;
      
      JacPDstandardNth3gt23 = PDstandardNth3gt23;
      
      JacPDstandardNth3gt33 = PDstandardNth3gt33;
      
      JacPDstandardNth3trK = PDstandardNth3trK;
      
      JacPDstandardNth3Xt1 = PDstandardNth3Xt1;
      
      JacPDstandardNth3Xt2 = PDstandardNth3Xt2;
      
      JacPDstandardNth3Xt3 = PDstandardNth3Xt3;
      
      JacPDstandardNth11gt11 = PDstandardNth11gt11;
      
      JacPDstandardNth11gt12 = PDstandardNth11gt12;
      
      JacPDstandardNth11gt13 = PDstandardNth11gt13;
      
      JacPDstandardNth11gt22 = PDstandardNth11gt22;
      
      JacPDstandardNth11gt23 = PDstandardNth11gt23;
      
      JacPDstandardNth11gt33 = PDstandardNth11gt33;
      
      JacPDstandardNth22gt11 = PDstandardNth22gt11;
      
      JacPDstandardNth22gt12 = PDstandardNth22gt12;
      
      JacPDstandardNth22gt13 = PDstandardNth22gt13;
      
      JacPDstandardNth22gt22 = PDstandardNth22gt22;
      
      JacPDstandardNth22gt23 = PDstandardNth22gt23;
      
      JacPDstandardNth22gt33 = PDstandardNth22gt33;
      
      JacPDstandardNth33gt11 = PDstandardNth33gt11;
      
      JacPDstandardNth33gt12 = PDstandardNth33gt12;
      
      JacPDstandardNth33gt13 = PDstandardNth33gt13;
      
      JacPDstandardNth33gt22 = PDstandardNth33gt22;
      
      JacPDstandardNth33gt23 = PDstandardNth33gt23;
      
      JacPDstandardNth33gt33 = PDstandardNth33gt33;
      
      JacPDstandardNth12gt11 = PDstandardNth12gt11;
      
      JacPDstandardNth12gt12 = PDstandardNth12gt12;
      
      JacPDstandardNth12gt13 = PDstandardNth12gt13;
      
      JacPDstandardNth12gt22 = PDstandardNth12gt22;
      
      JacPDstandardNth12gt23 = PDstandardNth12gt23;
      
      JacPDstandardNth12gt33 = PDstandardNth12gt33;
      
      JacPDstandardNth13gt11 = PDstandardNth13gt11;
      
      JacPDstandardNth13gt12 = PDstandardNth13gt12;
      
      JacPDstandardNth13gt13 = PDstandardNth13gt13;
      
      JacPDstandardNth13gt22 = PDstandardNth13gt22;
      
      JacPDstandardNth13gt23 = PDstandardNth13gt23;
      
      JacPDstandardNth13gt33 = PDstandardNth13gt33;
      
      JacPDstandardNth21gt11 = PDstandardNth12gt11;
      
      JacPDstandardNth21gt12 = PDstandardNth12gt12;
      
      JacPDstandardNth21gt13 = PDstandardNth12gt13;
      
      JacPDstandardNth21gt22 = PDstandardNth12gt22;
      
      JacPDstandardNth21gt23 = PDstandardNth12gt23;
      
      JacPDstandardNth21gt33 = PDstandardNth12gt33;
      
      JacPDstandardNth23gt11 = PDstandardNth23gt11;
      
      JacPDstandardNth23gt12 = PDstandardNth23gt12;
      
      JacPDstandardNth23gt13 = PDstandardNth23gt13;
      
      JacPDstandardNth23gt22 = PDstandardNth23gt22;
      
      JacPDstandardNth23gt23 = PDstandardNth23gt23;
      
      JacPDstandardNth23gt33 = PDstandardNth23gt33;
      
      JacPDstandardNth31gt11 = PDstandardNth13gt11;
      
      JacPDstandardNth31gt12 = PDstandardNth13gt12;
      
      JacPDstandardNth31gt13 = PDstandardNth13gt13;
      
      JacPDstandardNth31gt22 = PDstandardNth13gt22;
      
      JacPDstandardNth31gt23 = PDstandardNth13gt23;
      
      JacPDstandardNth31gt33 = PDstandardNth13gt33;
      
      JacPDstandardNth32gt11 = PDstandardNth23gt11;
      
      JacPDstandardNth32gt12 = PDstandardNth23gt12;
      
      JacPDstandardNth32gt13 = PDstandardNth23gt13;
      
      JacPDstandardNth32gt22 = PDstandardNth23gt22;
      
      JacPDstandardNth32gt23 = PDstandardNth23gt23;
      
      JacPDstandardNth32gt33 = PDstandardNth23gt33;
    }
    
    CCTK_REAL detgt CCTK_ATTRIBUTE_UNUSED = 1;
    
    CCTK_REAL gtu11 CCTK_ATTRIBUTE_UNUSED = (gt22L*gt33L - 
      pow(gt23L,2))*pow(detgt,-1);
    
    CCTK_REAL gtu21 CCTK_ATTRIBUTE_UNUSED = (gt13L*gt23L - 
      gt12L*gt33L)*pow(detgt,-1);
    
    CCTK_REAL gtu31 CCTK_ATTRIBUTE_UNUSED = (-(gt13L*gt22L) + 
      gt12L*gt23L)*pow(detgt,-1);
    
    CCTK_REAL gtu22 CCTK_ATTRIBUTE_UNUSED = (gt11L*gt33L - 
      pow(gt13L,2))*pow(detgt,-1);
    
    CCTK_REAL gtu32 CCTK_ATTRIBUTE_UNUSED = (gt12L*gt13L - 
      gt11L*gt23L)*pow(detgt,-1);
    
    CCTK_REAL gtu33 CCTK_ATTRIBUTE_UNUSED = (gt11L*gt22L - 
      pow(gt12L,2))*pow(detgt,-1);
    
    CCTK_REAL dgtu111 CCTK_ATTRIBUTE_UNUSED = 
      -2*gtu11*(gtu21*JacPDstandardNth1gt12 + gtu31*JacPDstandardNth1gt13) - 
      2*gtu21*gtu31*JacPDstandardNth1gt23 - 
      JacPDstandardNth1gt11*pow(gtu11,2) - JacPDstandardNth1gt22*pow(gtu21,2) 
      - JacPDstandardNth1gt33*pow(gtu31,2);
    
    CCTK_REAL dgtu211 CCTK_ATTRIBUTE_UNUSED = 
      -(gtu21*(gtu11*JacPDstandardNth1gt11 + gtu21*JacPDstandardNth1gt12 + 
      gtu31*JacPDstandardNth1gt13)) - gtu22*(gtu11*JacPDstandardNth1gt12 + 
      gtu21*JacPDstandardNth1gt22 + gtu31*JacPDstandardNth1gt23) - 
      gtu32*(gtu11*JacPDstandardNth1gt13 + gtu21*JacPDstandardNth1gt23 + 
      gtu31*JacPDstandardNth1gt33);
    
    CCTK_REAL dgtu311 CCTK_ATTRIBUTE_UNUSED = 
      -(gtu11*(gtu31*JacPDstandardNth1gt11 + gtu32*JacPDstandardNth1gt12 + 
      gtu33*JacPDstandardNth1gt13)) - gtu21*(gtu31*JacPDstandardNth1gt12 + 
      gtu32*JacPDstandardNth1gt22 + gtu33*JacPDstandardNth1gt23) - 
      gtu31*(gtu31*JacPDstandardNth1gt13 + gtu32*JacPDstandardNth1gt23 + 
      gtu33*JacPDstandardNth1gt33);
    
    CCTK_REAL dgtu221 CCTK_ATTRIBUTE_UNUSED = 
      -2*gtu21*(gtu22*JacPDstandardNth1gt12 + gtu32*JacPDstandardNth1gt13) - 
      2*gtu22*gtu32*JacPDstandardNth1gt23 - 
      JacPDstandardNth1gt11*pow(gtu21,2) - JacPDstandardNth1gt22*pow(gtu22,2) 
      - JacPDstandardNth1gt33*pow(gtu32,2);
    
    CCTK_REAL dgtu321 CCTK_ATTRIBUTE_UNUSED = 
      -(gtu21*(gtu31*JacPDstandardNth1gt11 + gtu32*JacPDstandardNth1gt12 + 
      gtu33*JacPDstandardNth1gt13)) - gtu22*(gtu31*JacPDstandardNth1gt12 + 
      gtu32*JacPDstandardNth1gt22 + gtu33*JacPDstandardNth1gt23) - 
      gtu32*(gtu31*JacPDstandardNth1gt13 + gtu32*JacPDstandardNth1gt23 + 
      gtu33*JacPDstandardNth1gt33);
    
    CCTK_REAL dgtu331 CCTK_ATTRIBUTE_UNUSED = 
      -2*gtu31*(gtu32*JacPDstandardNth1gt12 + gtu33*JacPDstandardNth1gt13) - 
      2*gtu32*gtu33*JacPDstandardNth1gt23 - 
      JacPDstandardNth1gt11*pow(gtu31,2) - JacPDstandardNth1gt22*pow(gtu32,2) 
      - JacPDstandardNth1gt33*pow(gtu33,2);
    
    CCTK_REAL dgtu112 CCTK_ATTRIBUTE_UNUSED = 
      -2*gtu11*(gtu21*JacPDstandardNth2gt12 + gtu31*JacPDstandardNth2gt13) - 
      2*gtu21*gtu31*JacPDstandardNth2gt23 - 
      JacPDstandardNth2gt11*pow(gtu11,2) - JacPDstandardNth2gt22*pow(gtu21,2) 
      - JacPDstandardNth2gt33*pow(gtu31,2);
    
    CCTK_REAL dgtu212 CCTK_ATTRIBUTE_UNUSED = 
      -(gtu21*(gtu11*JacPDstandardNth2gt11 + gtu21*JacPDstandardNth2gt12 + 
      gtu31*JacPDstandardNth2gt13)) - gtu22*(gtu11*JacPDstandardNth2gt12 + 
      gtu21*JacPDstandardNth2gt22 + gtu31*JacPDstandardNth2gt23) - 
      gtu32*(gtu11*JacPDstandardNth2gt13 + gtu21*JacPDstandardNth2gt23 + 
      gtu31*JacPDstandardNth2gt33);
    
    CCTK_REAL dgtu312 CCTK_ATTRIBUTE_UNUSED = 
      -(gtu11*(gtu31*JacPDstandardNth2gt11 + gtu32*JacPDstandardNth2gt12 + 
      gtu33*JacPDstandardNth2gt13)) - gtu21*(gtu31*JacPDstandardNth2gt12 + 
      gtu32*JacPDstandardNth2gt22 + gtu33*JacPDstandardNth2gt23) - 
      gtu31*(gtu31*JacPDstandardNth2gt13 + gtu32*JacPDstandardNth2gt23 + 
      gtu33*JacPDstandardNth2gt33);
    
    CCTK_REAL dgtu222 CCTK_ATTRIBUTE_UNUSED = 
      -2*gtu21*(gtu22*JacPDstandardNth2gt12 + gtu32*JacPDstandardNth2gt13) - 
      2*gtu22*gtu32*JacPDstandardNth2gt23 - 
      JacPDstandardNth2gt11*pow(gtu21,2) - JacPDstandardNth2gt22*pow(gtu22,2) 
      - JacPDstandardNth2gt33*pow(gtu32,2);
    
    CCTK_REAL dgtu322 CCTK_ATTRIBUTE_UNUSED = 
      -(gtu21*(gtu31*JacPDstandardNth2gt11 + gtu32*JacPDstandardNth2gt12 + 
      gtu33*JacPDstandardNth2gt13)) - gtu22*(gtu31*JacPDstandardNth2gt12 + 
      gtu32*JacPDstandardNth2gt22 + gtu33*JacPDstandardNth2gt23) - 
      gtu32*(gtu31*JacPDstandardNth2gt13 + gtu32*JacPDstandardNth2gt23 + 
      gtu33*JacPDstandardNth2gt33);
    
    CCTK_REAL dgtu332 CCTK_ATTRIBUTE_UNUSED = 
      -2*gtu31*(gtu32*JacPDstandardNth2gt12 + gtu33*JacPDstandardNth2gt13) - 
      2*gtu32*gtu33*JacPDstandardNth2gt23 - 
      JacPDstandardNth2gt11*pow(gtu31,2) - JacPDstandardNth2gt22*pow(gtu32,2) 
      - JacPDstandardNth2gt33*pow(gtu33,2);
    
    CCTK_REAL dgtu113 CCTK_ATTRIBUTE_UNUSED = 
      -2*gtu11*(gtu21*JacPDstandardNth3gt12 + gtu31*JacPDstandardNth3gt13) - 
      2*gtu21*gtu31*JacPDstandardNth3gt23 - 
      JacPDstandardNth3gt11*pow(gtu11,2) - JacPDstandardNth3gt22*pow(gtu21,2) 
      - JacPDstandardNth3gt33*pow(gtu31,2);
    
    CCTK_REAL dgtu213 CCTK_ATTRIBUTE_UNUSED = 
      -(gtu21*(gtu11*JacPDstandardNth3gt11 + gtu21*JacPDstandardNth3gt12 + 
      gtu31*JacPDstandardNth3gt13)) - gtu22*(gtu11*JacPDstandardNth3gt12 + 
      gtu21*JacPDstandardNth3gt22 + gtu31*JacPDstandardNth3gt23) - 
      gtu32*(gtu11*JacPDstandardNth3gt13 + gtu21*JacPDstandardNth3gt23 + 
      gtu31*JacPDstandardNth3gt33);
    
    CCTK_REAL dgtu313 CCTK_ATTRIBUTE_UNUSED = 
      -(gtu11*(gtu31*JacPDstandardNth3gt11 + gtu32*JacPDstandardNth3gt12 + 
      gtu33*JacPDstandardNth3gt13)) - gtu21*(gtu31*JacPDstandardNth3gt12 + 
      gtu32*JacPDstandardNth3gt22 + gtu33*JacPDstandardNth3gt23) - 
      gtu31*(gtu31*JacPDstandardNth3gt13 + gtu32*JacPDstandardNth3gt23 + 
      gtu33*JacPDstandardNth3gt33);
    
    CCTK_REAL dgtu223 CCTK_ATTRIBUTE_UNUSED = 
      -2*gtu21*(gtu22*JacPDstandardNth3gt12 + gtu32*JacPDstandardNth3gt13) - 
      2*gtu22*gtu32*JacPDstandardNth3gt23 - 
      JacPDstandardNth3gt11*pow(gtu21,2) - JacPDstandardNth3gt22*pow(gtu22,2) 
      - JacPDstandardNth3gt33*pow(gtu32,2);
    
    CCTK_REAL dgtu323 CCTK_ATTRIBUTE_UNUSED = 
      -(gtu21*(gtu31*JacPDstandardNth3gt11 + gtu32*JacPDstandardNth3gt12 + 
      gtu33*JacPDstandardNth3gt13)) - gtu22*(gtu31*JacPDstandardNth3gt12 + 
      gtu32*JacPDstandardNth3gt22 + gtu33*JacPDstandardNth3gt23) - 
      gtu32*(gtu31*JacPDstandardNth3gt13 + gtu32*JacPDstandardNth3gt23 + 
      gtu33*JacPDstandardNth3gt33);
    
    CCTK_REAL dgtu333 CCTK_ATTRIBUTE_UNUSED = 
      -2*gtu31*(gtu32*JacPDstandardNth3gt12 + gtu33*JacPDstandardNth3gt13) - 
      2*gtu32*gtu33*JacPDstandardNth3gt23 - 
      JacPDstandardNth3gt11*pow(gtu31,2) - JacPDstandardNth3gt22*pow(gtu32,2) 
      - JacPDstandardNth3gt33*pow(gtu33,2);
    
    CCTK_REAL Gtl111 CCTK_ATTRIBUTE_UNUSED = 0.5*JacPDstandardNth1gt11;
    
    CCTK_REAL Gtl112 CCTK_ATTRIBUTE_UNUSED = 0.5*JacPDstandardNth2gt11;
    
    CCTK_REAL Gtl113 CCTK_ATTRIBUTE_UNUSED = 0.5*JacPDstandardNth3gt11;
    
    CCTK_REAL Gtl122 CCTK_ATTRIBUTE_UNUSED = -0.5*JacPDstandardNth1gt22 + 
      JacPDstandardNth2gt12;
    
    CCTK_REAL Gtl123 CCTK_ATTRIBUTE_UNUSED = 0.5*(-JacPDstandardNth1gt23 + 
      JacPDstandardNth2gt13 + JacPDstandardNth3gt12);
    
    CCTK_REAL Gtl133 CCTK_ATTRIBUTE_UNUSED = -0.5*JacPDstandardNth1gt33 + 
      JacPDstandardNth3gt13;
    
    CCTK_REAL Gtl211 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth1gt12 - 
      0.5*JacPDstandardNth2gt11;
    
    CCTK_REAL Gtl212 CCTK_ATTRIBUTE_UNUSED = 0.5*JacPDstandardNth1gt22;
    
    CCTK_REAL Gtl213 CCTK_ATTRIBUTE_UNUSED = 0.5*(JacPDstandardNth1gt23 - 
      JacPDstandardNth2gt13 + JacPDstandardNth3gt12);
    
    CCTK_REAL Gtl222 CCTK_ATTRIBUTE_UNUSED = 0.5*JacPDstandardNth2gt22;
    
    CCTK_REAL Gtl223 CCTK_ATTRIBUTE_UNUSED = 0.5*JacPDstandardNth3gt22;
    
    CCTK_REAL Gtl233 CCTK_ATTRIBUTE_UNUSED = -0.5*JacPDstandardNth2gt33 + 
      JacPDstandardNth3gt23;
    
    CCTK_REAL Gtl311 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth1gt13 - 
      0.5*JacPDstandardNth3gt11;
    
    CCTK_REAL Gtl312 CCTK_ATTRIBUTE_UNUSED = 0.5*(JacPDstandardNth1gt23 + 
      JacPDstandardNth2gt13 - JacPDstandardNth3gt12);
    
    CCTK_REAL Gtl313 CCTK_ATTRIBUTE_UNUSED = 0.5*JacPDstandardNth1gt33;
    
    CCTK_REAL Gtl322 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth2gt23 - 
      0.5*JacPDstandardNth3gt22;
    
    CCTK_REAL Gtl323 CCTK_ATTRIBUTE_UNUSED = 0.5*JacPDstandardNth2gt33;
    
    CCTK_REAL Gtl333 CCTK_ATTRIBUTE_UNUSED = 0.5*JacPDstandardNth3gt33;
    
    CCTK_REAL Gtlu111 CCTK_ATTRIBUTE_UNUSED = Gtl111*gtu11 + Gtl112*gtu21 
      + Gtl113*gtu31;
    
    CCTK_REAL Gtlu112 CCTK_ATTRIBUTE_UNUSED = Gtl111*gtu21 + Gtl112*gtu22 
      + Gtl113*gtu32;
    
    CCTK_REAL Gtlu113 CCTK_ATTRIBUTE_UNUSED = Gtl111*gtu31 + Gtl112*gtu32 
      + Gtl113*gtu33;
    
    CCTK_REAL Gtlu121 CCTK_ATTRIBUTE_UNUSED = Gtl112*gtu11 + Gtl122*gtu21 
      + Gtl123*gtu31;
    
    CCTK_REAL Gtlu122 CCTK_ATTRIBUTE_UNUSED = Gtl112*gtu21 + Gtl122*gtu22 
      + Gtl123*gtu32;
    
    CCTK_REAL Gtlu123 CCTK_ATTRIBUTE_UNUSED = Gtl112*gtu31 + Gtl122*gtu32 
      + Gtl123*gtu33;
    
    CCTK_REAL Gtlu131 CCTK_ATTRIBUTE_UNUSED = Gtl113*gtu11 + Gtl123*gtu21 
      + Gtl133*gtu31;
    
    CCTK_REAL Gtlu132 CCTK_ATTRIBUTE_UNUSED = Gtl113*gtu21 + Gtl123*gtu22 
      + Gtl133*gtu32;
    
    CCTK_REAL Gtlu133 CCTK_ATTRIBUTE_UNUSED = Gtl113*gtu31 + Gtl123*gtu32 
      + Gtl133*gtu33;
    
    CCTK_REAL Gtlu211 CCTK_ATTRIBUTE_UNUSED = Gtl211*gtu11 + Gtl212*gtu21 
      + Gtl213*gtu31;
    
    CCTK_REAL Gtlu212 CCTK_ATTRIBUTE_UNUSED = Gtl211*gtu21 + Gtl212*gtu22 
      + Gtl213*gtu32;
    
    CCTK_REAL Gtlu213 CCTK_ATTRIBUTE_UNUSED = Gtl211*gtu31 + Gtl212*gtu32 
      + Gtl213*gtu33;
    
    CCTK_REAL Gtlu221 CCTK_ATTRIBUTE_UNUSED = Gtl212*gtu11 + Gtl222*gtu21 
      + Gtl223*gtu31;
    
    CCTK_REAL Gtlu222 CCTK_ATTRIBUTE_UNUSED = Gtl212*gtu21 + Gtl222*gtu22 
      + Gtl223*gtu32;
    
    CCTK_REAL Gtlu223 CCTK_ATTRIBUTE_UNUSED = Gtl212*gtu31 + Gtl222*gtu32 
      + Gtl223*gtu33;
    
    CCTK_REAL Gtlu231 CCTK_ATTRIBUTE_UNUSED = Gtl213*gtu11 + Gtl223*gtu21 
      + Gtl233*gtu31;
    
    CCTK_REAL Gtlu232 CCTK_ATTRIBUTE_UNUSED = Gtl213*gtu21 + Gtl223*gtu22 
      + Gtl233*gtu32;
    
    CCTK_REAL Gtlu233 CCTK_ATTRIBUTE_UNUSED = Gtl213*gtu31 + Gtl223*gtu32 
      + Gtl233*gtu33;
    
    CCTK_REAL Gtlu311 CCTK_ATTRIBUTE_UNUSED = Gtl311*gtu11 + Gtl312*gtu21 
      + Gtl313*gtu31;
    
    CCTK_REAL Gtlu312 CCTK_ATTRIBUTE_UNUSED = Gtl311*gtu21 + Gtl312*gtu22 
      + Gtl313*gtu32;
    
    CCTK_REAL Gtlu313 CCTK_ATTRIBUTE_UNUSED = Gtl311*gtu31 + Gtl312*gtu32 
      + Gtl313*gtu33;
    
    CCTK_REAL Gtlu321 CCTK_ATTRIBUTE_UNUSED = Gtl312*gtu11 + Gtl322*gtu21 
      + Gtl323*gtu31;
    
    CCTK_REAL Gtlu322 CCTK_ATTRIBUTE_UNUSED = Gtl312*gtu21 + Gtl322*gtu22 
      + Gtl323*gtu32;
    
    CCTK_REAL Gtlu323 CCTK_ATTRIBUTE_UNUSED = Gtl312*gtu31 + Gtl322*gtu32 
      + Gtl323*gtu33;
    
    CCTK_REAL Gtlu331 CCTK_ATTRIBUTE_UNUSED = Gtl313*gtu11 + Gtl323*gtu21 
      + Gtl333*gtu31;
    
    CCTK_REAL Gtlu332 CCTK_ATTRIBUTE_UNUSED = Gtl313*gtu21 + Gtl323*gtu22 
      + Gtl333*gtu32;
    
    CCTK_REAL Gtlu333 CCTK_ATTRIBUTE_UNUSED = Gtl313*gtu31 + Gtl323*gtu32 
      + Gtl333*gtu33;
    
    CCTK_REAL Gt111 CCTK_ATTRIBUTE_UNUSED = Gtl111*gtu11 + Gtl211*gtu21 + 
      Gtl311*gtu31;
    
    CCTK_REAL Gt211 CCTK_ATTRIBUTE_UNUSED = Gtl111*gtu21 + Gtl211*gtu22 + 
      Gtl311*gtu32;
    
    CCTK_REAL Gt311 CCTK_ATTRIBUTE_UNUSED = Gtl111*gtu31 + Gtl211*gtu32 + 
      Gtl311*gtu33;
    
    CCTK_REAL Gt112 CCTK_ATTRIBUTE_UNUSED = Gtl112*gtu11 + Gtl212*gtu21 + 
      Gtl312*gtu31;
    
    CCTK_REAL Gt212 CCTK_ATTRIBUTE_UNUSED = Gtl112*gtu21 + Gtl212*gtu22 + 
      Gtl312*gtu32;
    
    CCTK_REAL Gt312 CCTK_ATTRIBUTE_UNUSED = Gtl112*gtu31 + Gtl212*gtu32 + 
      Gtl312*gtu33;
    
    CCTK_REAL Gt113 CCTK_ATTRIBUTE_UNUSED = Gtl113*gtu11 + Gtl213*gtu21 + 
      Gtl313*gtu31;
    
    CCTK_REAL Gt213 CCTK_ATTRIBUTE_UNUSED = Gtl113*gtu21 + Gtl213*gtu22 + 
      Gtl313*gtu32;
    
    CCTK_REAL Gt313 CCTK_ATTRIBUTE_UNUSED = Gtl113*gtu31 + Gtl213*gtu32 + 
      Gtl313*gtu33;
    
    CCTK_REAL Gt122 CCTK_ATTRIBUTE_UNUSED = Gtl122*gtu11 + Gtl222*gtu21 + 
      Gtl322*gtu31;
    
    CCTK_REAL Gt222 CCTK_ATTRIBUTE_UNUSED = Gtl122*gtu21 + Gtl222*gtu22 + 
      Gtl322*gtu32;
    
    CCTK_REAL Gt322 CCTK_ATTRIBUTE_UNUSED = Gtl122*gtu31 + Gtl222*gtu32 + 
      Gtl322*gtu33;
    
    CCTK_REAL Gt123 CCTK_ATTRIBUTE_UNUSED = Gtl123*gtu11 + Gtl223*gtu21 + 
      Gtl323*gtu31;
    
    CCTK_REAL Gt223 CCTK_ATTRIBUTE_UNUSED = Gtl123*gtu21 + Gtl223*gtu22 + 
      Gtl323*gtu32;
    
    CCTK_REAL Gt323 CCTK_ATTRIBUTE_UNUSED = Gtl123*gtu31 + Gtl223*gtu32 + 
      Gtl323*gtu33;
    
    CCTK_REAL Gt133 CCTK_ATTRIBUTE_UNUSED = Gtl133*gtu11 + Gtl233*gtu21 + 
      Gtl333*gtu31;
    
    CCTK_REAL Gt233 CCTK_ATTRIBUTE_UNUSED = Gtl133*gtu21 + Gtl233*gtu22 + 
      Gtl333*gtu32;
    
    CCTK_REAL Gt333 CCTK_ATTRIBUTE_UNUSED = Gtl133*gtu31 + Gtl233*gtu32 + 
      Gtl333*gtu33;
    
    CCTK_REAL Xtn1 CCTK_ATTRIBUTE_UNUSED = Gt111*gtu11 + Gt122*gtu22 + 
      2*(Gt112*gtu21 + Gt113*gtu31 + Gt123*gtu32) + Gt133*gtu33;
    
    CCTK_REAL Xtn2 CCTK_ATTRIBUTE_UNUSED = Gt211*gtu11 + Gt222*gtu22 + 
      2*(Gt212*gtu21 + Gt213*gtu31 + Gt223*gtu32) + Gt233*gtu33;
    
    CCTK_REAL Xtn3 CCTK_ATTRIBUTE_UNUSED = Gt311*gtu11 + Gt322*gtu22 + 
      2*(Gt312*gtu21 + Gt313*gtu31 + Gt323*gtu32) + Gt333*gtu33;
    
    CCTK_REAL Rt11 CCTK_ATTRIBUTE_UNUSED = 3*(Gt111*Gtlu111 + 
      Gt112*Gtlu112 + Gt113*Gtlu113) + 2*(Gt211*Gtlu121 + Gt212*Gtlu122 + 
      Gt213*Gtlu123 + Gt311*Gtlu131 + Gt312*Gtlu132 + Gt313*Gtlu133) + 
      Gt211*Gtlu211 + Gt212*Gtlu212 + Gt213*Gtlu213 + Gt311*Gtlu311 + 
      Gt312*Gtlu312 + Gt313*Gtlu313 + gt11L*JacPDstandardNth1Xt1 + 
      gt12L*JacPDstandardNth1Xt2 + gt13L*JacPDstandardNth1Xt3 + 
      0.5*(-(gtu11*JacPDstandardNth11gt11) - gtu21*(JacPDstandardNth12gt11 + 
      JacPDstandardNth21gt11) - gtu22*JacPDstandardNth22gt11 + 
      gtu31*(-JacPDstandardNth13gt11 - JacPDstandardNth31gt11) + 
      gtu32*(-JacPDstandardNth23gt11 - JacPDstandardNth32gt11) - 
      gtu33*JacPDstandardNth33gt11) + Gtl111*Xtn1 + Gtl112*Xtn2 + 
      Gtl113*Xtn3;
    
    CCTK_REAL Rt12 CCTK_ATTRIBUTE_UNUSED = 0.5*(4*(Gt211*Gtlu221 + 
      Gt212*Gtlu222 + Gt213*Gtlu223) + 2*(Gt122*Gtlu112 + Gt123*Gtlu113 + 
      Gt111*Gtlu121 + Gt212*Gtlu121 + Gt222*Gtlu122 + Gt113*Gtlu123 + 
      Gt223*Gtlu123 + Gt312*Gtlu131 + Gt322*Gtlu132 + Gt323*Gtlu133 + 
      Gt111*Gtlu211 + Gt112*(Gtlu111 + Gtlu122 + Gtlu212) + Gt113*Gtlu213 + 
      Gt311*Gtlu231 + Gt312*Gtlu232 + Gt313*Gtlu233 + Gt311*Gtlu321 + 
      Gt312*Gtlu322 + Gt313*Gtlu323) - gtu11*JacPDstandardNth11gt12 + 
      gt12L*JacPDstandardNth1Xt1 + gt22L*JacPDstandardNth1Xt2 + 
      gt23L*JacPDstandardNth1Xt3 + gtu21*(-JacPDstandardNth12gt12 - 
      JacPDstandardNth21gt12) - gtu22*JacPDstandardNth22gt12 + 
      gt11L*JacPDstandardNth2Xt1 + gt12L*JacPDstandardNth2Xt2 + 
      gt13L*JacPDstandardNth2Xt3 + gtu31*(-JacPDstandardNth13gt12 - 
      JacPDstandardNth31gt12) + gtu32*(-JacPDstandardNth23gt12 - 
      JacPDstandardNth32gt12) - gtu33*JacPDstandardNth33gt12 + Gtl112*Xtn1 + 
      Gtl211*Xtn1 + Gtl122*Xtn2 + Gtl212*Xtn2 + Gtl123*Xtn3 + Gtl213*Xtn3);
    
    CCTK_REAL Rt13 CCTK_ATTRIBUTE_UNUSED = 0.5*(2*(Gt123*Gtlu112 + 
      Gt133*Gtlu113 + Gt213*Gtlu121 + Gt223*Gtlu122 + Gt233*Gtlu123 + 
      Gt111*Gtlu131 + Gt313*Gtlu131 + Gt112*Gtlu132 + Gt323*Gtlu132 + 
      Gt333*Gtlu133 + Gt211*Gtlu231 + Gt212*Gtlu232 + Gt213*Gtlu233 + 
      Gt111*Gtlu311 + Gt112*Gtlu312 + Gt113*(Gtlu111 + Gtlu133 + Gtlu313) + 
      Gt211*Gtlu321 + Gt212*Gtlu322 + Gt213*Gtlu323) + 4*(Gt311*Gtlu331 + 
      Gt312*Gtlu332 + Gt313*Gtlu333) - gtu11*JacPDstandardNth11gt13 + 
      gt13L*JacPDstandardNth1Xt1 + gt23L*JacPDstandardNth1Xt2 + 
      gt33L*JacPDstandardNth1Xt3 + gtu21*(-JacPDstandardNth12gt13 - 
      JacPDstandardNth21gt13) - gtu22*JacPDstandardNth22gt13 + 
      gtu31*(-JacPDstandardNth13gt13 - JacPDstandardNth31gt13) + 
      gtu32*(-JacPDstandardNth23gt13 - JacPDstandardNth32gt13) - 
      gtu33*JacPDstandardNth33gt13 + gt11L*JacPDstandardNth3Xt1 + 
      gt12L*JacPDstandardNth3Xt2 + gt13L*JacPDstandardNth3Xt3 + Gtl113*Xtn1 + 
      Gtl311*Xtn1 + Gtl123*Xtn2 + Gtl312*Xtn2 + Gtl133*Xtn3 + Gtl313*Xtn3);
    
    CCTK_REAL Rt22 CCTK_ATTRIBUTE_UNUSED = Gt112*(Gtlu121 + 2*Gtlu211) + 
      Gt122*(Gtlu122 + 2*Gtlu212) + Gt123*(Gtlu123 + 2*Gtlu213) + 
      3*(Gt212*Gtlu221 + Gt222*Gtlu222 + Gt223*Gtlu223) + 2*(Gt312*Gtlu231 + 
      Gt322*Gtlu232 + Gt323*Gtlu233) + Gt312*Gtlu321 + Gt322*Gtlu322 + 
      Gt323*Gtlu323 + gt12L*JacPDstandardNth2Xt1 + gt22L*JacPDstandardNth2Xt2 
      + gt23L*JacPDstandardNth2Xt3 + 0.5*(-(gtu11*JacPDstandardNth11gt22) - 
      gtu21*(JacPDstandardNth12gt22 + JacPDstandardNth21gt22) - 
      gtu22*JacPDstandardNth22gt22 + gtu31*(-JacPDstandardNth13gt22 - 
      JacPDstandardNth31gt22) + gtu32*(-JacPDstandardNth23gt22 - 
      JacPDstandardNth32gt22) - gtu33*JacPDstandardNth33gt22) + Gtl212*Xtn1 + 
      Gtl222*Xtn2 + Gtl223*Xtn3;
    
    CCTK_REAL Rt23 CCTK_ATTRIBUTE_UNUSED = 0.5*(2*(Gt123*Gtlu133 + 
      Gt113*Gtlu211 + Gt123*Gtlu212 + Gt133*Gtlu213 + Gt213*Gtlu221 + 
      Gt223*Gtlu222 + Gt233*Gtlu223 + Gt212*Gtlu231 + Gt313*Gtlu231 + 
      Gt222*Gtlu232 + Gt323*Gtlu232 + Gt223*Gtlu233 + Gt333*Gtlu233 + 
      Gt112*(Gtlu131 + Gtlu311) + Gt122*(Gtlu132 + Gtlu312) + Gt123*Gtlu313 + 
      Gt212*Gtlu321 + Gt222*Gtlu322 + Gt223*Gtlu323) + 4*(Gt312*Gtlu331 + 
      Gt322*Gtlu332 + Gt323*Gtlu333) - gtu11*JacPDstandardNth11gt23 + 
      gtu21*(-JacPDstandardNth12gt23 - JacPDstandardNth21gt23) - 
      gtu22*JacPDstandardNth22gt23 + gt13L*JacPDstandardNth2Xt1 + 
      gt23L*JacPDstandardNth2Xt2 + gt33L*JacPDstandardNth2Xt3 + 
      gtu31*(-JacPDstandardNth13gt23 - JacPDstandardNth31gt23) + 
      gtu32*(-JacPDstandardNth23gt23 - JacPDstandardNth32gt23) - 
      gtu33*JacPDstandardNth33gt23 + gt12L*JacPDstandardNth3Xt1 + 
      gt22L*JacPDstandardNth3Xt2 + gt23L*JacPDstandardNth3Xt3 + Gtl213*Xtn1 + 
      Gtl312*Xtn1 + Gtl223*Xtn2 + Gtl322*Xtn2 + Gtl233*Xtn3 + Gtl323*Xtn3);
    
    CCTK_REAL Rt33 CCTK_ATTRIBUTE_UNUSED = Gt113*(Gtlu131 + 2*Gtlu311) + 
      Gt123*(Gtlu132 + 2*Gtlu312) + Gt133*(Gtlu133 + 2*Gtlu313) + 
      Gt213*(Gtlu231 + 2*Gtlu321) + Gt223*(Gtlu232 + 2*Gtlu322) + 
      Gt233*(Gtlu233 + 2*Gtlu323) + 3*(Gt313*Gtlu331 + Gt323*Gtlu332 + 
      Gt333*Gtlu333) + 0.5*(-(gtu11*JacPDstandardNth11gt33) - 
      gtu21*(JacPDstandardNth12gt33 + JacPDstandardNth21gt33) - 
      gtu22*JacPDstandardNth22gt33 + gtu31*(-JacPDstandardNth13gt33 - 
      JacPDstandardNth31gt33) + gtu32*(-JacPDstandardNth23gt33 - 
      JacPDstandardNth32gt33) - gtu33*JacPDstandardNth33gt33) + 
      gt13L*JacPDstandardNth3Xt1 + gt23L*JacPDstandardNth3Xt2 + 
      gt33L*JacPDstandardNth3Xt3 + Gtl313*Xtn1 + Gtl323*Xtn2 + Gtl333*Xtn3;
    
    CCTK_REAL trRt CCTK_ATTRIBUTE_UNUSED = gtu11*Rt11 + gtu22*Rt22 + 
      2*(gtu21*Rt12 + gtu31*Rt13 + gtu32*Rt23) + gtu33*Rt33;
    
    CCTK_REAL ephi CCTK_ATTRIBUTE_UNUSED = 
      IfThen(conformalMethod,pow(phiL,-0.5),exp(phiL));
    
    CCTK_REAL Atm11 CCTK_ATTRIBUTE_UNUSED = At11L*gtu11 + At12L*gtu21 + 
      At13L*gtu31;
    
    CCTK_REAL Atm21 CCTK_ATTRIBUTE_UNUSED = At11L*gtu21 + At12L*gtu22 + 
      At13L*gtu32;
    
    CCTK_REAL Atm31 CCTK_ATTRIBUTE_UNUSED = At11L*gtu31 + At12L*gtu32 + 
      At13L*gtu33;
    
    CCTK_REAL Atm12 CCTK_ATTRIBUTE_UNUSED = At12L*gtu11 + At22L*gtu21 + 
      At23L*gtu31;
    
    CCTK_REAL Atm22 CCTK_ATTRIBUTE_UNUSED = At12L*gtu21 + At22L*gtu22 + 
      At23L*gtu32;
    
    CCTK_REAL Atm32 CCTK_ATTRIBUTE_UNUSED = At12L*gtu31 + At22L*gtu32 + 
      At23L*gtu33;
    
    CCTK_REAL Atm13 CCTK_ATTRIBUTE_UNUSED = At13L*gtu11 + At23L*gtu21 + 
      At33L*gtu31;
    
    CCTK_REAL Atm23 CCTK_ATTRIBUTE_UNUSED = At13L*gtu21 + At23L*gtu22 + 
      At33L*gtu32;
    
    CCTK_REAL Atm33 CCTK_ATTRIBUTE_UNUSED = At13L*gtu31 + At23L*gtu32 + 
      At33L*gtu33;
    
    CCTK_REAL rho CCTK_ATTRIBUTE_UNUSED = pow(alphaL,-2)*(eTttL - 
      2*(beta2L*eTtyL + beta3L*eTtzL) + 2*(beta1L*(-eTtxL + beta2L*eTxyL + 
      beta3L*eTxzL) + beta2L*beta3L*eTyzL) + eTxxL*pow(beta1L,2) + 
      eTyyL*pow(beta2L,2) + eTzzL*pow(beta3L,2));
    
    CCTK_REAL S1 CCTK_ATTRIBUTE_UNUSED = (-eTtxL + beta1L*eTxxL + 
      beta2L*eTxyL + beta3L*eTxzL)*pow(alphaL,-1);
    
    CCTK_REAL S2 CCTK_ATTRIBUTE_UNUSED = (-eTtyL + beta1L*eTxyL + 
      beta2L*eTyyL + beta3L*eTyzL)*pow(alphaL,-1);
    
    CCTK_REAL S3 CCTK_ATTRIBUTE_UNUSED = (-eTtzL + beta1L*eTxzL + 
      beta2L*eTyzL + beta3L*eTzzL)*pow(alphaL,-1);
    
    CCTK_REAL MadmL CCTK_ATTRIBUTE_UNUSED = 
      -0.0198943678864869169711104704216*(Gtlu111*(Gt111*gtu11 + Gt112*gtu21 
      + Gt113*gtu31) + Gtlu112*(Gt112*gtu11 + Gt122*gtu21 + Gt123*gtu31) + 
      Gtlu113*(Gt113*gtu11 + Gt123*gtu21 + Gt133*gtu31) + 
      Gtlu121*(Gt211*gtu11 + Gt212*gtu21 + Gt213*gtu31) + 
      Gtlu122*(Gt212*gtu11 + Gt222*gtu21 + Gt223*gtu31) + 
      Gtlu123*(Gt213*gtu11 + Gt223*gtu21 + Gt233*gtu31) + 
      Gtlu131*(Gt311*gtu11 + Gt312*gtu21 + Gt313*gtu31) + 
      Gtlu132*(Gt312*gtu11 + Gt322*gtu21 + Gt323*gtu31) + 
      Gtlu133*(Gt313*gtu11 + Gt323*gtu21 + Gt333*gtu31) + 
      Gtlu211*(Gt111*gtu21 + Gt112*gtu22 + Gt113*gtu32) + 
      Gtlu212*(Gt112*gtu21 + Gt122*gtu22 + Gt123*gtu32) + 
      Gtlu213*(Gt113*gtu21 + Gt123*gtu22 + Gt133*gtu32) + 
      Gtlu221*(Gt211*gtu21 + Gt212*gtu22 + Gt213*gtu32) + 
      Gtlu222*(Gt212*gtu21 + Gt222*gtu22 + Gt223*gtu32) + 
      Gtlu223*(Gt213*gtu21 + Gt223*gtu22 + Gt233*gtu32) + 
      Gtlu231*(Gt311*gtu21 + Gt312*gtu22 + Gt313*gtu32) + 
      Gtlu232*(Gt312*gtu21 + Gt322*gtu22 + Gt323*gtu32) + 
      Gtlu233*(Gt313*gtu21 + Gt323*gtu22 + Gt333*gtu32) + 
      Gtlu311*(Gt111*gtu31 + Gt112*gtu32 + Gt113*gtu33) + 
      Gtlu312*(Gt112*gtu31 + Gt122*gtu32 + Gt123*gtu33) + 
      Gtlu313*(Gt113*gtu31 + Gt123*gtu32 + Gt133*gtu33) + 
      Gtlu321*(Gt211*gtu31 + Gt212*gtu32 + Gt213*gtu33) + 
      Gtlu322*(Gt212*gtu31 + Gt222*gtu32 + Gt223*gtu33) + 
      Gtlu323*(Gt213*gtu31 + Gt223*gtu32 + Gt233*gtu33) + 
      Gtlu331*(Gt311*gtu31 + Gt312*gtu32 + Gt313*gtu33) + 
      Gtlu332*(Gt312*gtu31 + Gt322*gtu32 + Gt323*gtu33) + 
      Gtlu333*(Gt313*gtu31 + Gt323*gtu32 + Gt333*gtu33) - trRt + ephi*trRt - 
      (2*Atm12*Atm21 + 2*Atm13*Atm31 + 2*Atm23*Atm32 + 16*Pi*rho - 
      0.666666666666666666666666666667*pow(trKL,2) + pow(Atm11,2) + 
      pow(Atm22,2) + pow(Atm33,2))*pow(ephi,5));
    
    CCTK_REAL Jadm1L CCTK_ATTRIBUTE_UNUSED = 
      0.0198943678864869169711104704216*(2*Atm23 - 2*Atm32 + 
      zL*(At11L*dgtu112 + At22L*dgtu222 + 2*(At12L*dgtu212 + At13L*dgtu312 + 
      At23L*dgtu322) + At33L*dgtu332 - 
      1.33333333333333333333333333333*JacPDstandardNth2trK - 16*Pi*S2) + 
      yL*(-(At11L*dgtu113) - 2*At12L*dgtu213 - At22L*dgtu223 - 
      2*At13L*dgtu313 - 2*At23L*dgtu323 - At33L*dgtu333 + 
      1.33333333333333333333333333333*JacPDstandardNth3trK + 
      16*Pi*S3))*pow(ephi,6);
    
    CCTK_REAL Jadm2L CCTK_ATTRIBUTE_UNUSED = 
      0.0198943678864869169711104704216*(-2*Atm13 + 2*Atm31 + 
      zL*(-(At11L*dgtu111) - 2*At12L*dgtu211 - At22L*dgtu221 - 
      2*At13L*dgtu311 - 2*At23L*dgtu321 - At33L*dgtu331 + 
      1.33333333333333333333333333333*JacPDstandardNth1trK + 16*Pi*S1) + 
      xL*(At11L*dgtu113 + At22L*dgtu223 + 2*(At12L*dgtu213 + At13L*dgtu313 + 
      At23L*dgtu323) + At33L*dgtu333 - 
      1.33333333333333333333333333333*JacPDstandardNth3trK - 
      16*Pi*S3))*pow(ephi,6);
    
    CCTK_REAL Jadm3L CCTK_ATTRIBUTE_UNUSED = 
      0.0198943678864869169711104704216*(2*Atm12 - 2*Atm21 + 
      yL*(At11L*dgtu111 + At22L*dgtu221 + 2*(At12L*dgtu211 + At13L*dgtu311 + 
      At23L*dgtu321) + At33L*dgtu331 - 
      1.33333333333333333333333333333*JacPDstandardNth1trK - 16*Pi*S1) + 
      xL*(-(At11L*dgtu112) - 2*At12L*dgtu212 - At22L*dgtu222 - 
      2*At13L*dgtu312 - 2*At23L*dgtu322 - At33L*dgtu332 + 
      1.33333333333333333333333333333*JacPDstandardNth2trK + 
      16*Pi*S2))*pow(ephi,6);
    /* Copy local copies back to grid functions */
    Jadm1[index] = Jadm1L;
    Jadm2[index] = Jadm2L;
    Jadm3[index] = Jadm3L;
    Madm[index] = MadmL;
  }
  CCTK_ENDLOOP3(ML_ADMQuantities_calculate);
}
extern "C" void ML_ADMQuantities_calculate(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_ML_ADMQuantities_calculate
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_ADMQuantities_calculate);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_ADMQuantities_calculate_Body");
  }
  if (cctk_iteration % ML_ADMQuantities_calculate_calc_every != ML_ADMQuantities_calculate_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "grid::coordinates",
    "ML_BSSN::ML_curv",
    "ML_BSSN::ML_Gamma",
    "ML_BSSN::ML_lapse",
    "ML_BSSN::ML_log_confac",
    "ML_BSSN::ML_metric",
    "ML_BSSN::ML_shift",
    "ML_BSSN::ML_trace_curv",
    "ML_ADMQuantities::ML_Jadm",
    "ML_ADMQuantities::ML_Madm"};
  AssertGroupStorage(cctkGH, "ML_ADMQuantities_calculate", 10, groups);
  
  EnsureStencilFits(cctkGH, "ML_ADMQuantities_calculate", 2, 2, 2);
  
  LoopOverInterior(cctkGH, ML_ADMQuantities_calculate_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_ADMQuantities_calculate_Body");
  }
}

} // namespace ML_ADMQuantities
