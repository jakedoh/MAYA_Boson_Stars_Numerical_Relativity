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

namespace ML_ADMConstraints {

extern "C" void ML_ADMConstraints_evaluate_SelectBCs(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_ML_ADMConstraints_evaluate_SelectBCs
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_ADMConstraints_evaluate_SelectBCs);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % ML_ADMConstraints_evaluate_calc_every != ML_ADMConstraints_evaluate_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADMConstraints::ML_Ham","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_ADMConstraints::ML_Ham.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_ADMConstraints::ML_mom","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_ADMConstraints::ML_mom.");
  return;
}

static void ML_ADMConstraints_evaluate_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  CCTK_LOOP3(ML_ADMConstraints_evaluate,
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
    CCTK_REAL kxxL CCTK_ATTRIBUTE_UNUSED = kxx[index];
    CCTK_REAL kxyL CCTK_ATTRIBUTE_UNUSED = kxy[index];
    CCTK_REAL kxzL CCTK_ATTRIBUTE_UNUSED = kxz[index];
    CCTK_REAL kyyL CCTK_ATTRIBUTE_UNUSED = kyy[index];
    CCTK_REAL kyzL CCTK_ATTRIBUTE_UNUSED = kyz[index];
    CCTK_REAL kzzL CCTK_ATTRIBUTE_UNUSED = kzz[index];
    
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
    const CCTK_REAL PDstandardNth1gxx CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&gxx[index]);
    const CCTK_REAL PDstandardNth2gxx CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&gxx[index]);
    const CCTK_REAL PDstandardNth3gxx CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&gxx[index]);
    const CCTK_REAL PDstandardNth11gxx CCTK_ATTRIBUTE_UNUSED = PDstandardNth11(&gxx[index]);
    const CCTK_REAL PDstandardNth22gxx CCTK_ATTRIBUTE_UNUSED = PDstandardNth22(&gxx[index]);
    const CCTK_REAL PDstandardNth33gxx CCTK_ATTRIBUTE_UNUSED = PDstandardNth33(&gxx[index]);
    const CCTK_REAL PDstandardNth12gxx CCTK_ATTRIBUTE_UNUSED = PDstandardNth12(&gxx[index]);
    const CCTK_REAL PDstandardNth13gxx CCTK_ATTRIBUTE_UNUSED = PDstandardNth13(&gxx[index]);
    const CCTK_REAL PDstandardNth23gxx CCTK_ATTRIBUTE_UNUSED = PDstandardNth23(&gxx[index]);
    const CCTK_REAL PDstandardNth1gxy CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&gxy[index]);
    const CCTK_REAL PDstandardNth2gxy CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&gxy[index]);
    const CCTK_REAL PDstandardNth3gxy CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&gxy[index]);
    const CCTK_REAL PDstandardNth11gxy CCTK_ATTRIBUTE_UNUSED = PDstandardNth11(&gxy[index]);
    const CCTK_REAL PDstandardNth22gxy CCTK_ATTRIBUTE_UNUSED = PDstandardNth22(&gxy[index]);
    const CCTK_REAL PDstandardNth33gxy CCTK_ATTRIBUTE_UNUSED = PDstandardNth33(&gxy[index]);
    const CCTK_REAL PDstandardNth12gxy CCTK_ATTRIBUTE_UNUSED = PDstandardNth12(&gxy[index]);
    const CCTK_REAL PDstandardNth13gxy CCTK_ATTRIBUTE_UNUSED = PDstandardNth13(&gxy[index]);
    const CCTK_REAL PDstandardNth23gxy CCTK_ATTRIBUTE_UNUSED = PDstandardNth23(&gxy[index]);
    const CCTK_REAL PDstandardNth1gxz CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&gxz[index]);
    const CCTK_REAL PDstandardNth2gxz CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&gxz[index]);
    const CCTK_REAL PDstandardNth3gxz CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&gxz[index]);
    const CCTK_REAL PDstandardNth11gxz CCTK_ATTRIBUTE_UNUSED = PDstandardNth11(&gxz[index]);
    const CCTK_REAL PDstandardNth22gxz CCTK_ATTRIBUTE_UNUSED = PDstandardNth22(&gxz[index]);
    const CCTK_REAL PDstandardNth33gxz CCTK_ATTRIBUTE_UNUSED = PDstandardNth33(&gxz[index]);
    const CCTK_REAL PDstandardNth12gxz CCTK_ATTRIBUTE_UNUSED = PDstandardNth12(&gxz[index]);
    const CCTK_REAL PDstandardNth13gxz CCTK_ATTRIBUTE_UNUSED = PDstandardNth13(&gxz[index]);
    const CCTK_REAL PDstandardNth23gxz CCTK_ATTRIBUTE_UNUSED = PDstandardNth23(&gxz[index]);
    const CCTK_REAL PDstandardNth1gyy CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&gyy[index]);
    const CCTK_REAL PDstandardNth2gyy CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&gyy[index]);
    const CCTK_REAL PDstandardNth3gyy CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&gyy[index]);
    const CCTK_REAL PDstandardNth11gyy CCTK_ATTRIBUTE_UNUSED = PDstandardNth11(&gyy[index]);
    const CCTK_REAL PDstandardNth22gyy CCTK_ATTRIBUTE_UNUSED = PDstandardNth22(&gyy[index]);
    const CCTK_REAL PDstandardNth33gyy CCTK_ATTRIBUTE_UNUSED = PDstandardNth33(&gyy[index]);
    const CCTK_REAL PDstandardNth12gyy CCTK_ATTRIBUTE_UNUSED = PDstandardNth12(&gyy[index]);
    const CCTK_REAL PDstandardNth13gyy CCTK_ATTRIBUTE_UNUSED = PDstandardNth13(&gyy[index]);
    const CCTK_REAL PDstandardNth23gyy CCTK_ATTRIBUTE_UNUSED = PDstandardNth23(&gyy[index]);
    const CCTK_REAL PDstandardNth1gyz CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&gyz[index]);
    const CCTK_REAL PDstandardNth2gyz CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&gyz[index]);
    const CCTK_REAL PDstandardNth3gyz CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&gyz[index]);
    const CCTK_REAL PDstandardNth11gyz CCTK_ATTRIBUTE_UNUSED = PDstandardNth11(&gyz[index]);
    const CCTK_REAL PDstandardNth22gyz CCTK_ATTRIBUTE_UNUSED = PDstandardNth22(&gyz[index]);
    const CCTK_REAL PDstandardNth33gyz CCTK_ATTRIBUTE_UNUSED = PDstandardNth33(&gyz[index]);
    const CCTK_REAL PDstandardNth12gyz CCTK_ATTRIBUTE_UNUSED = PDstandardNth12(&gyz[index]);
    const CCTK_REAL PDstandardNth13gyz CCTK_ATTRIBUTE_UNUSED = PDstandardNth13(&gyz[index]);
    const CCTK_REAL PDstandardNth23gyz CCTK_ATTRIBUTE_UNUSED = PDstandardNth23(&gyz[index]);
    const CCTK_REAL PDstandardNth1gzz CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&gzz[index]);
    const CCTK_REAL PDstandardNth2gzz CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&gzz[index]);
    const CCTK_REAL PDstandardNth3gzz CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&gzz[index]);
    const CCTK_REAL PDstandardNth11gzz CCTK_ATTRIBUTE_UNUSED = PDstandardNth11(&gzz[index]);
    const CCTK_REAL PDstandardNth22gzz CCTK_ATTRIBUTE_UNUSED = PDstandardNth22(&gzz[index]);
    const CCTK_REAL PDstandardNth33gzz CCTK_ATTRIBUTE_UNUSED = PDstandardNth33(&gzz[index]);
    const CCTK_REAL PDstandardNth12gzz CCTK_ATTRIBUTE_UNUSED = PDstandardNth12(&gzz[index]);
    const CCTK_REAL PDstandardNth13gzz CCTK_ATTRIBUTE_UNUSED = PDstandardNth13(&gzz[index]);
    const CCTK_REAL PDstandardNth23gzz CCTK_ATTRIBUTE_UNUSED = PDstandardNth23(&gzz[index]);
    const CCTK_REAL PDstandardNth1kxx CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&kxx[index]);
    const CCTK_REAL PDstandardNth2kxx CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&kxx[index]);
    const CCTK_REAL PDstandardNth3kxx CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&kxx[index]);
    const CCTK_REAL PDstandardNth1kxy CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&kxy[index]);
    const CCTK_REAL PDstandardNth2kxy CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&kxy[index]);
    const CCTK_REAL PDstandardNth3kxy CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&kxy[index]);
    const CCTK_REAL PDstandardNth1kxz CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&kxz[index]);
    const CCTK_REAL PDstandardNth2kxz CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&kxz[index]);
    const CCTK_REAL PDstandardNth3kxz CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&kxz[index]);
    const CCTK_REAL PDstandardNth1kyy CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&kyy[index]);
    const CCTK_REAL PDstandardNth2kyy CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&kyy[index]);
    const CCTK_REAL PDstandardNth3kyy CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&kyy[index]);
    const CCTK_REAL PDstandardNth1kyz CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&kyz[index]);
    const CCTK_REAL PDstandardNth2kyz CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&kyz[index]);
    const CCTK_REAL PDstandardNth3kyz CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&kyz[index]);
    const CCTK_REAL PDstandardNth1kzz CCTK_ATTRIBUTE_UNUSED = PDstandardNth1(&kzz[index]);
    const CCTK_REAL PDstandardNth2kzz CCTK_ATTRIBUTE_UNUSED = PDstandardNth2(&kzz[index]);
    const CCTK_REAL PDstandardNth3kzz CCTK_ATTRIBUTE_UNUSED = PDstandardNth3(&kzz[index]);
    /* Calculate temporaries and grid functions */
    CCTK_REAL JacPDstandardNth11gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth11gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth11gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth12gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth12gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth12gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth12gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth12gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth12gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth13gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth13gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth13gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth13gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth13gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth13gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1kxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1kxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1kyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1kyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1kzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth21gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth21gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth21gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth21gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth21gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth21gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth22gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth22gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth22gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth23gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth23gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth23gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth23gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth23gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth23gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2kxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2kxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2kxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2kyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2kzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth31gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth31gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth31gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth31gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth31gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth31gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth32gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth32gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth32gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth32gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth32gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth32gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth33gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth33gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth33gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3kxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3kxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3kxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3kyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3kyz CCTK_ATTRIBUTE_UNUSED;
    
    if (usejacobian)
    {
      JacPDstandardNth1gxx = J11L*PDstandardNth1gxx + J21L*PDstandardNth2gxx 
        + J31L*PDstandardNth3gxx;
      
      JacPDstandardNth1gxy = J11L*PDstandardNth1gxy + J21L*PDstandardNth2gxy 
        + J31L*PDstandardNth3gxy;
      
      JacPDstandardNth1gxz = J11L*PDstandardNth1gxz + J21L*PDstandardNth2gxz 
        + J31L*PDstandardNth3gxz;
      
      JacPDstandardNth1gyy = J11L*PDstandardNth1gyy + J21L*PDstandardNth2gyy 
        + J31L*PDstandardNth3gyy;
      
      JacPDstandardNth1gyz = J11L*PDstandardNth1gyz + J21L*PDstandardNth2gyz 
        + J31L*PDstandardNth3gyz;
      
      JacPDstandardNth1gzz = J11L*PDstandardNth1gzz + J21L*PDstandardNth2gzz 
        + J31L*PDstandardNth3gzz;
      
      JacPDstandardNth1kxy = J11L*PDstandardNth1kxy + J21L*PDstandardNth2kxy 
        + J31L*PDstandardNth3kxy;
      
      JacPDstandardNth1kxz = J11L*PDstandardNth1kxz + J21L*PDstandardNth2kxz 
        + J31L*PDstandardNth3kxz;
      
      JacPDstandardNth1kyy = J11L*PDstandardNth1kyy + J21L*PDstandardNth2kyy 
        + J31L*PDstandardNth3kyy;
      
      JacPDstandardNth1kyz = J11L*PDstandardNth1kyz + J21L*PDstandardNth2kyz 
        + J31L*PDstandardNth3kyz;
      
      JacPDstandardNth1kzz = J11L*PDstandardNth1kzz + J21L*PDstandardNth2kzz 
        + J31L*PDstandardNth3kzz;
      
      JacPDstandardNth2gxx = J12L*PDstandardNth1gxx + J22L*PDstandardNth2gxx 
        + J32L*PDstandardNth3gxx;
      
      JacPDstandardNth2gxy = J12L*PDstandardNth1gxy + J22L*PDstandardNth2gxy 
        + J32L*PDstandardNth3gxy;
      
      JacPDstandardNth2gxz = J12L*PDstandardNth1gxz + J22L*PDstandardNth2gxz 
        + J32L*PDstandardNth3gxz;
      
      JacPDstandardNth2gyy = J12L*PDstandardNth1gyy + J22L*PDstandardNth2gyy 
        + J32L*PDstandardNth3gyy;
      
      JacPDstandardNth2gyz = J12L*PDstandardNth1gyz + J22L*PDstandardNth2gyz 
        + J32L*PDstandardNth3gyz;
      
      JacPDstandardNth2gzz = J12L*PDstandardNth1gzz + J22L*PDstandardNth2gzz 
        + J32L*PDstandardNth3gzz;
      
      JacPDstandardNth2kxx = J12L*PDstandardNth1kxx + J22L*PDstandardNth2kxx 
        + J32L*PDstandardNth3kxx;
      
      JacPDstandardNth2kxy = J12L*PDstandardNth1kxy + J22L*PDstandardNth2kxy 
        + J32L*PDstandardNth3kxy;
      
      JacPDstandardNth2kxz = J12L*PDstandardNth1kxz + J22L*PDstandardNth2kxz 
        + J32L*PDstandardNth3kxz;
      
      JacPDstandardNth2kyz = J12L*PDstandardNth1kyz + J22L*PDstandardNth2kyz 
        + J32L*PDstandardNth3kyz;
      
      JacPDstandardNth2kzz = J12L*PDstandardNth1kzz + J22L*PDstandardNth2kzz 
        + J32L*PDstandardNth3kzz;
      
      JacPDstandardNth3gxx = J13L*PDstandardNth1gxx + J23L*PDstandardNth2gxx 
        + J33L*PDstandardNth3gxx;
      
      JacPDstandardNth3gxy = J13L*PDstandardNth1gxy + J23L*PDstandardNth2gxy 
        + J33L*PDstandardNth3gxy;
      
      JacPDstandardNth3gxz = J13L*PDstandardNth1gxz + J23L*PDstandardNth2gxz 
        + J33L*PDstandardNth3gxz;
      
      JacPDstandardNth3gyy = J13L*PDstandardNth1gyy + J23L*PDstandardNth2gyy 
        + J33L*PDstandardNth3gyy;
      
      JacPDstandardNth3gyz = J13L*PDstandardNth1gyz + J23L*PDstandardNth2gyz 
        + J33L*PDstandardNth3gyz;
      
      JacPDstandardNth3gzz = J13L*PDstandardNth1gzz + J23L*PDstandardNth2gzz 
        + J33L*PDstandardNth3gzz;
      
      JacPDstandardNth3kxx = J13L*PDstandardNth1kxx + J23L*PDstandardNth2kxx 
        + J33L*PDstandardNth3kxx;
      
      JacPDstandardNth3kxy = J13L*PDstandardNth1kxy + J23L*PDstandardNth2kxy 
        + J33L*PDstandardNth3kxy;
      
      JacPDstandardNth3kxz = J13L*PDstandardNth1kxz + J23L*PDstandardNth2kxz 
        + J33L*PDstandardNth3kxz;
      
      JacPDstandardNth3kyy = J13L*PDstandardNth1kyy + J23L*PDstandardNth2kyy 
        + J33L*PDstandardNth3kyy;
      
      JacPDstandardNth3kyz = J13L*PDstandardNth1kyz + J23L*PDstandardNth2kyz 
        + J33L*PDstandardNth3kyz;
      
      JacPDstandardNth11gyy = dJ111L*PDstandardNth1gyy + 
        2*(J11L*(J21L*PDstandardNth12gyy + J31L*PDstandardNth13gyy) + 
        J21L*J31L*PDstandardNth23gyy) + dJ211L*PDstandardNth2gyy + 
        dJ311L*PDstandardNth3gyy + PDstandardNth11gyy*pow(J11L,2) + 
        PDstandardNth22gyy*pow(J21L,2) + PDstandardNth33gyy*pow(J31L,2);
      
      JacPDstandardNth11gyz = dJ111L*PDstandardNth1gyz + 
        2*(J11L*(J21L*PDstandardNth12gyz + J31L*PDstandardNth13gyz) + 
        J21L*J31L*PDstandardNth23gyz) + dJ211L*PDstandardNth2gyz + 
        dJ311L*PDstandardNth3gyz + PDstandardNth11gyz*pow(J11L,2) + 
        PDstandardNth22gyz*pow(J21L,2) + PDstandardNth33gyz*pow(J31L,2);
      
      JacPDstandardNth11gzz = dJ111L*PDstandardNth1gzz + 
        2*(J11L*(J21L*PDstandardNth12gzz + J31L*PDstandardNth13gzz) + 
        J21L*J31L*PDstandardNth23gzz) + dJ211L*PDstandardNth2gzz + 
        dJ311L*PDstandardNth3gzz + PDstandardNth11gzz*pow(J11L,2) + 
        PDstandardNth22gzz*pow(J21L,2) + PDstandardNth33gzz*pow(J31L,2);
      
      JacPDstandardNth22gxx = dJ122L*PDstandardNth1gxx + 
        2*(J12L*(J22L*PDstandardNth12gxx + J32L*PDstandardNth13gxx) + 
        J22L*J32L*PDstandardNth23gxx) + dJ222L*PDstandardNth2gxx + 
        dJ322L*PDstandardNth3gxx + PDstandardNth11gxx*pow(J12L,2) + 
        PDstandardNth22gxx*pow(J22L,2) + PDstandardNth33gxx*pow(J32L,2);
      
      JacPDstandardNth22gxz = dJ122L*PDstandardNth1gxz + 
        2*(J12L*(J22L*PDstandardNth12gxz + J32L*PDstandardNth13gxz) + 
        J22L*J32L*PDstandardNth23gxz) + dJ222L*PDstandardNth2gxz + 
        dJ322L*PDstandardNth3gxz + PDstandardNth11gxz*pow(J12L,2) + 
        PDstandardNth22gxz*pow(J22L,2) + PDstandardNth33gxz*pow(J32L,2);
      
      JacPDstandardNth22gzz = dJ122L*PDstandardNth1gzz + 
        2*(J12L*(J22L*PDstandardNth12gzz + J32L*PDstandardNth13gzz) + 
        J22L*J32L*PDstandardNth23gzz) + dJ222L*PDstandardNth2gzz + 
        dJ322L*PDstandardNth3gzz + PDstandardNth11gzz*pow(J12L,2) + 
        PDstandardNth22gzz*pow(J22L,2) + PDstandardNth33gzz*pow(J32L,2);
      
      JacPDstandardNth33gxx = dJ133L*PDstandardNth1gxx + 
        2*(J13L*(J23L*PDstandardNth12gxx + J33L*PDstandardNth13gxx) + 
        J23L*J33L*PDstandardNth23gxx) + dJ233L*PDstandardNth2gxx + 
        dJ333L*PDstandardNth3gxx + PDstandardNth11gxx*pow(J13L,2) + 
        PDstandardNth22gxx*pow(J23L,2) + PDstandardNth33gxx*pow(J33L,2);
      
      JacPDstandardNth33gxy = dJ133L*PDstandardNth1gxy + 
        2*(J13L*(J23L*PDstandardNth12gxy + J33L*PDstandardNth13gxy) + 
        J23L*J33L*PDstandardNth23gxy) + dJ233L*PDstandardNth2gxy + 
        dJ333L*PDstandardNth3gxy + PDstandardNth11gxy*pow(J13L,2) + 
        PDstandardNth22gxy*pow(J23L,2) + PDstandardNth33gxy*pow(J33L,2);
      
      JacPDstandardNth33gyy = dJ133L*PDstandardNth1gyy + 
        2*(J13L*(J23L*PDstandardNth12gyy + J33L*PDstandardNth13gyy) + 
        J23L*J33L*PDstandardNth23gyy) + dJ233L*PDstandardNth2gyy + 
        dJ333L*PDstandardNth3gyy + PDstandardNth11gyy*pow(J13L,2) + 
        PDstandardNth22gyy*pow(J23L,2) + PDstandardNth33gyy*pow(J33L,2);
      
      JacPDstandardNth12gxx = J12L*(J11L*PDstandardNth11gxx + 
        J21L*PDstandardNth12gxx + J31L*PDstandardNth13gxx) + 
        J11L*(J22L*PDstandardNth12gxx + J32L*PDstandardNth13gxx) + 
        dJ112L*PDstandardNth1gxx + J22L*(J21L*PDstandardNth22gxx + 
        J31L*PDstandardNth23gxx) + dJ212L*PDstandardNth2gxx + 
        J32L*(J21L*PDstandardNth23gxx + J31L*PDstandardNth33gxx) + 
        dJ312L*PDstandardNth3gxx;
      
      JacPDstandardNth12gxy = J12L*(J11L*PDstandardNth11gxy + 
        J21L*PDstandardNth12gxy + J31L*PDstandardNth13gxy) + 
        J11L*(J22L*PDstandardNth12gxy + J32L*PDstandardNth13gxy) + 
        dJ112L*PDstandardNth1gxy + J22L*(J21L*PDstandardNth22gxy + 
        J31L*PDstandardNth23gxy) + dJ212L*PDstandardNth2gxy + 
        J32L*(J21L*PDstandardNth23gxy + J31L*PDstandardNth33gxy) + 
        dJ312L*PDstandardNth3gxy;
      
      JacPDstandardNth12gxz = J12L*(J11L*PDstandardNth11gxz + 
        J21L*PDstandardNth12gxz + J31L*PDstandardNth13gxz) + 
        J11L*(J22L*PDstandardNth12gxz + J32L*PDstandardNth13gxz) + 
        dJ112L*PDstandardNth1gxz + J22L*(J21L*PDstandardNth22gxz + 
        J31L*PDstandardNth23gxz) + dJ212L*PDstandardNth2gxz + 
        J32L*(J21L*PDstandardNth23gxz + J31L*PDstandardNth33gxz) + 
        dJ312L*PDstandardNth3gxz;
      
      JacPDstandardNth12gyy = J12L*(J11L*PDstandardNth11gyy + 
        J21L*PDstandardNth12gyy + J31L*PDstandardNth13gyy) + 
        J11L*(J22L*PDstandardNth12gyy + J32L*PDstandardNth13gyy) + 
        dJ112L*PDstandardNth1gyy + J22L*(J21L*PDstandardNth22gyy + 
        J31L*PDstandardNth23gyy) + dJ212L*PDstandardNth2gyy + 
        J32L*(J21L*PDstandardNth23gyy + J31L*PDstandardNth33gyy) + 
        dJ312L*PDstandardNth3gyy;
      
      JacPDstandardNth12gyz = J12L*(J11L*PDstandardNth11gyz + 
        J21L*PDstandardNth12gyz + J31L*PDstandardNth13gyz) + 
        J11L*(J22L*PDstandardNth12gyz + J32L*PDstandardNth13gyz) + 
        dJ112L*PDstandardNth1gyz + J22L*(J21L*PDstandardNth22gyz + 
        J31L*PDstandardNth23gyz) + dJ212L*PDstandardNth2gyz + 
        J32L*(J21L*PDstandardNth23gyz + J31L*PDstandardNth33gyz) + 
        dJ312L*PDstandardNth3gyz;
      
      JacPDstandardNth12gzz = J12L*(J11L*PDstandardNth11gzz + 
        J21L*PDstandardNth12gzz + J31L*PDstandardNth13gzz) + 
        J11L*(J22L*PDstandardNth12gzz + J32L*PDstandardNth13gzz) + 
        dJ112L*PDstandardNth1gzz + J22L*(J21L*PDstandardNth22gzz + 
        J31L*PDstandardNth23gzz) + dJ212L*PDstandardNth2gzz + 
        J32L*(J21L*PDstandardNth23gzz + J31L*PDstandardNth33gzz) + 
        dJ312L*PDstandardNth3gzz;
      
      JacPDstandardNth13gxx = J13L*(J11L*PDstandardNth11gxx + 
        J21L*PDstandardNth12gxx + J31L*PDstandardNth13gxx) + 
        J11L*(J23L*PDstandardNth12gxx + J33L*PDstandardNth13gxx) + 
        dJ113L*PDstandardNth1gxx + J23L*(J21L*PDstandardNth22gxx + 
        J31L*PDstandardNth23gxx) + dJ213L*PDstandardNth2gxx + 
        J33L*(J21L*PDstandardNth23gxx + J31L*PDstandardNth33gxx) + 
        dJ313L*PDstandardNth3gxx;
      
      JacPDstandardNth13gxy = J13L*(J11L*PDstandardNth11gxy + 
        J21L*PDstandardNth12gxy + J31L*PDstandardNth13gxy) + 
        J11L*(J23L*PDstandardNth12gxy + J33L*PDstandardNth13gxy) + 
        dJ113L*PDstandardNth1gxy + J23L*(J21L*PDstandardNth22gxy + 
        J31L*PDstandardNth23gxy) + dJ213L*PDstandardNth2gxy + 
        J33L*(J21L*PDstandardNth23gxy + J31L*PDstandardNth33gxy) + 
        dJ313L*PDstandardNth3gxy;
      
      JacPDstandardNth13gxz = J13L*(J11L*PDstandardNth11gxz + 
        J21L*PDstandardNth12gxz + J31L*PDstandardNth13gxz) + 
        J11L*(J23L*PDstandardNth12gxz + J33L*PDstandardNth13gxz) + 
        dJ113L*PDstandardNth1gxz + J23L*(J21L*PDstandardNth22gxz + 
        J31L*PDstandardNth23gxz) + dJ213L*PDstandardNth2gxz + 
        J33L*(J21L*PDstandardNth23gxz + J31L*PDstandardNth33gxz) + 
        dJ313L*PDstandardNth3gxz;
      
      JacPDstandardNth13gyy = J13L*(J11L*PDstandardNth11gyy + 
        J21L*PDstandardNth12gyy + J31L*PDstandardNth13gyy) + 
        J11L*(J23L*PDstandardNth12gyy + J33L*PDstandardNth13gyy) + 
        dJ113L*PDstandardNth1gyy + J23L*(J21L*PDstandardNth22gyy + 
        J31L*PDstandardNth23gyy) + dJ213L*PDstandardNth2gyy + 
        J33L*(J21L*PDstandardNth23gyy + J31L*PDstandardNth33gyy) + 
        dJ313L*PDstandardNth3gyy;
      
      JacPDstandardNth13gyz = J13L*(J11L*PDstandardNth11gyz + 
        J21L*PDstandardNth12gyz + J31L*PDstandardNth13gyz) + 
        J11L*(J23L*PDstandardNth12gyz + J33L*PDstandardNth13gyz) + 
        dJ113L*PDstandardNth1gyz + J23L*(J21L*PDstandardNth22gyz + 
        J31L*PDstandardNth23gyz) + dJ213L*PDstandardNth2gyz + 
        J33L*(J21L*PDstandardNth23gyz + J31L*PDstandardNth33gyz) + 
        dJ313L*PDstandardNth3gyz;
      
      JacPDstandardNth13gzz = J13L*(J11L*PDstandardNth11gzz + 
        J21L*PDstandardNth12gzz + J31L*PDstandardNth13gzz) + 
        J11L*(J23L*PDstandardNth12gzz + J33L*PDstandardNth13gzz) + 
        dJ113L*PDstandardNth1gzz + J23L*(J21L*PDstandardNth22gzz + 
        J31L*PDstandardNth23gzz) + dJ213L*PDstandardNth2gzz + 
        J33L*(J21L*PDstandardNth23gzz + J31L*PDstandardNth33gzz) + 
        dJ313L*PDstandardNth3gzz;
      
      JacPDstandardNth21gxx = J12L*(J11L*PDstandardNth11gxx + 
        J21L*PDstandardNth12gxx + J31L*PDstandardNth13gxx) + 
        J11L*(J22L*PDstandardNth12gxx + J32L*PDstandardNth13gxx) + 
        dJ112L*PDstandardNth1gxx + J22L*(J21L*PDstandardNth22gxx + 
        J31L*PDstandardNth23gxx) + dJ212L*PDstandardNth2gxx + 
        J32L*(J21L*PDstandardNth23gxx + J31L*PDstandardNth33gxx) + 
        dJ312L*PDstandardNth3gxx;
      
      JacPDstandardNth21gxy = J12L*(J11L*PDstandardNth11gxy + 
        J21L*PDstandardNth12gxy + J31L*PDstandardNth13gxy) + 
        J11L*(J22L*PDstandardNth12gxy + J32L*PDstandardNth13gxy) + 
        dJ112L*PDstandardNth1gxy + J22L*(J21L*PDstandardNth22gxy + 
        J31L*PDstandardNth23gxy) + dJ212L*PDstandardNth2gxy + 
        J32L*(J21L*PDstandardNth23gxy + J31L*PDstandardNth33gxy) + 
        dJ312L*PDstandardNth3gxy;
      
      JacPDstandardNth21gxz = J12L*(J11L*PDstandardNth11gxz + 
        J21L*PDstandardNth12gxz + J31L*PDstandardNth13gxz) + 
        J11L*(J22L*PDstandardNth12gxz + J32L*PDstandardNth13gxz) + 
        dJ112L*PDstandardNth1gxz + J22L*(J21L*PDstandardNth22gxz + 
        J31L*PDstandardNth23gxz) + dJ212L*PDstandardNth2gxz + 
        J32L*(J21L*PDstandardNth23gxz + J31L*PDstandardNth33gxz) + 
        dJ312L*PDstandardNth3gxz;
      
      JacPDstandardNth21gyy = J12L*(J11L*PDstandardNth11gyy + 
        J21L*PDstandardNth12gyy + J31L*PDstandardNth13gyy) + 
        J11L*(J22L*PDstandardNth12gyy + J32L*PDstandardNth13gyy) + 
        dJ112L*PDstandardNth1gyy + J22L*(J21L*PDstandardNth22gyy + 
        J31L*PDstandardNth23gyy) + dJ212L*PDstandardNth2gyy + 
        J32L*(J21L*PDstandardNth23gyy + J31L*PDstandardNth33gyy) + 
        dJ312L*PDstandardNth3gyy;
      
      JacPDstandardNth21gyz = J12L*(J11L*PDstandardNth11gyz + 
        J21L*PDstandardNth12gyz + J31L*PDstandardNth13gyz) + 
        J11L*(J22L*PDstandardNth12gyz + J32L*PDstandardNth13gyz) + 
        dJ112L*PDstandardNth1gyz + J22L*(J21L*PDstandardNth22gyz + 
        J31L*PDstandardNth23gyz) + dJ212L*PDstandardNth2gyz + 
        J32L*(J21L*PDstandardNth23gyz + J31L*PDstandardNth33gyz) + 
        dJ312L*PDstandardNth3gyz;
      
      JacPDstandardNth21gzz = J12L*(J11L*PDstandardNth11gzz + 
        J21L*PDstandardNth12gzz + J31L*PDstandardNth13gzz) + 
        J11L*(J22L*PDstandardNth12gzz + J32L*PDstandardNth13gzz) + 
        dJ112L*PDstandardNth1gzz + J22L*(J21L*PDstandardNth22gzz + 
        J31L*PDstandardNth23gzz) + dJ212L*PDstandardNth2gzz + 
        J32L*(J21L*PDstandardNth23gzz + J31L*PDstandardNth33gzz) + 
        dJ312L*PDstandardNth3gzz;
      
      JacPDstandardNth23gxx = J13L*(J12L*PDstandardNth11gxx + 
        J22L*PDstandardNth12gxx + J32L*PDstandardNth13gxx) + 
        J12L*(J23L*PDstandardNth12gxx + J33L*PDstandardNth13gxx) + 
        dJ123L*PDstandardNth1gxx + J23L*(J22L*PDstandardNth22gxx + 
        J32L*PDstandardNth23gxx) + dJ223L*PDstandardNth2gxx + 
        J33L*(J22L*PDstandardNth23gxx + J32L*PDstandardNth33gxx) + 
        dJ323L*PDstandardNth3gxx;
      
      JacPDstandardNth23gxy = J13L*(J12L*PDstandardNth11gxy + 
        J22L*PDstandardNth12gxy + J32L*PDstandardNth13gxy) + 
        J12L*(J23L*PDstandardNth12gxy + J33L*PDstandardNth13gxy) + 
        dJ123L*PDstandardNth1gxy + J23L*(J22L*PDstandardNth22gxy + 
        J32L*PDstandardNth23gxy) + dJ223L*PDstandardNth2gxy + 
        J33L*(J22L*PDstandardNth23gxy + J32L*PDstandardNth33gxy) + 
        dJ323L*PDstandardNth3gxy;
      
      JacPDstandardNth23gxz = J13L*(J12L*PDstandardNth11gxz + 
        J22L*PDstandardNth12gxz + J32L*PDstandardNth13gxz) + 
        J12L*(J23L*PDstandardNth12gxz + J33L*PDstandardNth13gxz) + 
        dJ123L*PDstandardNth1gxz + J23L*(J22L*PDstandardNth22gxz + 
        J32L*PDstandardNth23gxz) + dJ223L*PDstandardNth2gxz + 
        J33L*(J22L*PDstandardNth23gxz + J32L*PDstandardNth33gxz) + 
        dJ323L*PDstandardNth3gxz;
      
      JacPDstandardNth23gyy = J13L*(J12L*PDstandardNth11gyy + 
        J22L*PDstandardNth12gyy + J32L*PDstandardNth13gyy) + 
        J12L*(J23L*PDstandardNth12gyy + J33L*PDstandardNth13gyy) + 
        dJ123L*PDstandardNth1gyy + J23L*(J22L*PDstandardNth22gyy + 
        J32L*PDstandardNth23gyy) + dJ223L*PDstandardNth2gyy + 
        J33L*(J22L*PDstandardNth23gyy + J32L*PDstandardNth33gyy) + 
        dJ323L*PDstandardNth3gyy;
      
      JacPDstandardNth23gyz = J13L*(J12L*PDstandardNth11gyz + 
        J22L*PDstandardNth12gyz + J32L*PDstandardNth13gyz) + 
        J12L*(J23L*PDstandardNth12gyz + J33L*PDstandardNth13gyz) + 
        dJ123L*PDstandardNth1gyz + J23L*(J22L*PDstandardNth22gyz + 
        J32L*PDstandardNth23gyz) + dJ223L*PDstandardNth2gyz + 
        J33L*(J22L*PDstandardNth23gyz + J32L*PDstandardNth33gyz) + 
        dJ323L*PDstandardNth3gyz;
      
      JacPDstandardNth23gzz = J13L*(J12L*PDstandardNth11gzz + 
        J22L*PDstandardNth12gzz + J32L*PDstandardNth13gzz) + 
        J12L*(J23L*PDstandardNth12gzz + J33L*PDstandardNth13gzz) + 
        dJ123L*PDstandardNth1gzz + J23L*(J22L*PDstandardNth22gzz + 
        J32L*PDstandardNth23gzz) + dJ223L*PDstandardNth2gzz + 
        J33L*(J22L*PDstandardNth23gzz + J32L*PDstandardNth33gzz) + 
        dJ323L*PDstandardNth3gzz;
      
      JacPDstandardNth31gxx = J13L*(J11L*PDstandardNth11gxx + 
        J21L*PDstandardNth12gxx + J31L*PDstandardNth13gxx) + 
        J11L*(J23L*PDstandardNth12gxx + J33L*PDstandardNth13gxx) + 
        dJ113L*PDstandardNth1gxx + J23L*(J21L*PDstandardNth22gxx + 
        J31L*PDstandardNth23gxx) + dJ213L*PDstandardNth2gxx + 
        J33L*(J21L*PDstandardNth23gxx + J31L*PDstandardNth33gxx) + 
        dJ313L*PDstandardNth3gxx;
      
      JacPDstandardNth31gxy = J13L*(J11L*PDstandardNth11gxy + 
        J21L*PDstandardNth12gxy + J31L*PDstandardNth13gxy) + 
        J11L*(J23L*PDstandardNth12gxy + J33L*PDstandardNth13gxy) + 
        dJ113L*PDstandardNth1gxy + J23L*(J21L*PDstandardNth22gxy + 
        J31L*PDstandardNth23gxy) + dJ213L*PDstandardNth2gxy + 
        J33L*(J21L*PDstandardNth23gxy + J31L*PDstandardNth33gxy) + 
        dJ313L*PDstandardNth3gxy;
      
      JacPDstandardNth31gxz = J13L*(J11L*PDstandardNth11gxz + 
        J21L*PDstandardNth12gxz + J31L*PDstandardNth13gxz) + 
        J11L*(J23L*PDstandardNth12gxz + J33L*PDstandardNth13gxz) + 
        dJ113L*PDstandardNth1gxz + J23L*(J21L*PDstandardNth22gxz + 
        J31L*PDstandardNth23gxz) + dJ213L*PDstandardNth2gxz + 
        J33L*(J21L*PDstandardNth23gxz + J31L*PDstandardNth33gxz) + 
        dJ313L*PDstandardNth3gxz;
      
      JacPDstandardNth31gyy = J13L*(J11L*PDstandardNth11gyy + 
        J21L*PDstandardNth12gyy + J31L*PDstandardNth13gyy) + 
        J11L*(J23L*PDstandardNth12gyy + J33L*PDstandardNth13gyy) + 
        dJ113L*PDstandardNth1gyy + J23L*(J21L*PDstandardNth22gyy + 
        J31L*PDstandardNth23gyy) + dJ213L*PDstandardNth2gyy + 
        J33L*(J21L*PDstandardNth23gyy + J31L*PDstandardNth33gyy) + 
        dJ313L*PDstandardNth3gyy;
      
      JacPDstandardNth31gyz = J13L*(J11L*PDstandardNth11gyz + 
        J21L*PDstandardNth12gyz + J31L*PDstandardNth13gyz) + 
        J11L*(J23L*PDstandardNth12gyz + J33L*PDstandardNth13gyz) + 
        dJ113L*PDstandardNth1gyz + J23L*(J21L*PDstandardNth22gyz + 
        J31L*PDstandardNth23gyz) + dJ213L*PDstandardNth2gyz + 
        J33L*(J21L*PDstandardNth23gyz + J31L*PDstandardNth33gyz) + 
        dJ313L*PDstandardNth3gyz;
      
      JacPDstandardNth31gzz = J13L*(J11L*PDstandardNth11gzz + 
        J21L*PDstandardNth12gzz + J31L*PDstandardNth13gzz) + 
        J11L*(J23L*PDstandardNth12gzz + J33L*PDstandardNth13gzz) + 
        dJ113L*PDstandardNth1gzz + J23L*(J21L*PDstandardNth22gzz + 
        J31L*PDstandardNth23gzz) + dJ213L*PDstandardNth2gzz + 
        J33L*(J21L*PDstandardNth23gzz + J31L*PDstandardNth33gzz) + 
        dJ313L*PDstandardNth3gzz;
      
      JacPDstandardNth32gxx = J13L*(J12L*PDstandardNth11gxx + 
        J22L*PDstandardNth12gxx + J32L*PDstandardNth13gxx) + 
        J12L*(J23L*PDstandardNth12gxx + J33L*PDstandardNth13gxx) + 
        dJ123L*PDstandardNth1gxx + J23L*(J22L*PDstandardNth22gxx + 
        J32L*PDstandardNth23gxx) + dJ223L*PDstandardNth2gxx + 
        J33L*(J22L*PDstandardNth23gxx + J32L*PDstandardNth33gxx) + 
        dJ323L*PDstandardNth3gxx;
      
      JacPDstandardNth32gxy = J13L*(J12L*PDstandardNth11gxy + 
        J22L*PDstandardNth12gxy + J32L*PDstandardNth13gxy) + 
        J12L*(J23L*PDstandardNth12gxy + J33L*PDstandardNth13gxy) + 
        dJ123L*PDstandardNth1gxy + J23L*(J22L*PDstandardNth22gxy + 
        J32L*PDstandardNth23gxy) + dJ223L*PDstandardNth2gxy + 
        J33L*(J22L*PDstandardNth23gxy + J32L*PDstandardNth33gxy) + 
        dJ323L*PDstandardNth3gxy;
      
      JacPDstandardNth32gxz = J13L*(J12L*PDstandardNth11gxz + 
        J22L*PDstandardNth12gxz + J32L*PDstandardNth13gxz) + 
        J12L*(J23L*PDstandardNth12gxz + J33L*PDstandardNth13gxz) + 
        dJ123L*PDstandardNth1gxz + J23L*(J22L*PDstandardNth22gxz + 
        J32L*PDstandardNth23gxz) + dJ223L*PDstandardNth2gxz + 
        J33L*(J22L*PDstandardNth23gxz + J32L*PDstandardNth33gxz) + 
        dJ323L*PDstandardNth3gxz;
      
      JacPDstandardNth32gyy = J13L*(J12L*PDstandardNth11gyy + 
        J22L*PDstandardNth12gyy + J32L*PDstandardNth13gyy) + 
        J12L*(J23L*PDstandardNth12gyy + J33L*PDstandardNth13gyy) + 
        dJ123L*PDstandardNth1gyy + J23L*(J22L*PDstandardNth22gyy + 
        J32L*PDstandardNth23gyy) + dJ223L*PDstandardNth2gyy + 
        J33L*(J22L*PDstandardNth23gyy + J32L*PDstandardNth33gyy) + 
        dJ323L*PDstandardNth3gyy;
      
      JacPDstandardNth32gyz = J13L*(J12L*PDstandardNth11gyz + 
        J22L*PDstandardNth12gyz + J32L*PDstandardNth13gyz) + 
        J12L*(J23L*PDstandardNth12gyz + J33L*PDstandardNth13gyz) + 
        dJ123L*PDstandardNth1gyz + J23L*(J22L*PDstandardNth22gyz + 
        J32L*PDstandardNth23gyz) + dJ223L*PDstandardNth2gyz + 
        J33L*(J22L*PDstandardNth23gyz + J32L*PDstandardNth33gyz) + 
        dJ323L*PDstandardNth3gyz;
      
      JacPDstandardNth32gzz = J13L*(J12L*PDstandardNth11gzz + 
        J22L*PDstandardNth12gzz + J32L*PDstandardNth13gzz) + 
        J12L*(J23L*PDstandardNth12gzz + J33L*PDstandardNth13gzz) + 
        dJ123L*PDstandardNth1gzz + J23L*(J22L*PDstandardNth22gzz + 
        J32L*PDstandardNth23gzz) + dJ223L*PDstandardNth2gzz + 
        J33L*(J22L*PDstandardNth23gzz + J32L*PDstandardNth33gzz) + 
        dJ323L*PDstandardNth3gzz;
    }
    else
    {
      JacPDstandardNth1gxx = PDstandardNth1gxx;
      
      JacPDstandardNth1gxy = PDstandardNth1gxy;
      
      JacPDstandardNth1gxz = PDstandardNth1gxz;
      
      JacPDstandardNth1gyy = PDstandardNth1gyy;
      
      JacPDstandardNth1gyz = PDstandardNth1gyz;
      
      JacPDstandardNth1gzz = PDstandardNth1gzz;
      
      JacPDstandardNth1kxy = PDstandardNth1kxy;
      
      JacPDstandardNth1kxz = PDstandardNth1kxz;
      
      JacPDstandardNth1kyy = PDstandardNth1kyy;
      
      JacPDstandardNth1kyz = PDstandardNth1kyz;
      
      JacPDstandardNth1kzz = PDstandardNth1kzz;
      
      JacPDstandardNth2gxx = PDstandardNth2gxx;
      
      JacPDstandardNth2gxy = PDstandardNth2gxy;
      
      JacPDstandardNth2gxz = PDstandardNth2gxz;
      
      JacPDstandardNth2gyy = PDstandardNth2gyy;
      
      JacPDstandardNth2gyz = PDstandardNth2gyz;
      
      JacPDstandardNth2gzz = PDstandardNth2gzz;
      
      JacPDstandardNth2kxx = PDstandardNth2kxx;
      
      JacPDstandardNth2kxy = PDstandardNth2kxy;
      
      JacPDstandardNth2kxz = PDstandardNth2kxz;
      
      JacPDstandardNth2kyz = PDstandardNth2kyz;
      
      JacPDstandardNth2kzz = PDstandardNth2kzz;
      
      JacPDstandardNth3gxx = PDstandardNth3gxx;
      
      JacPDstandardNth3gxy = PDstandardNth3gxy;
      
      JacPDstandardNth3gxz = PDstandardNth3gxz;
      
      JacPDstandardNth3gyy = PDstandardNth3gyy;
      
      JacPDstandardNth3gyz = PDstandardNth3gyz;
      
      JacPDstandardNth3gzz = PDstandardNth3gzz;
      
      JacPDstandardNth3kxx = PDstandardNth3kxx;
      
      JacPDstandardNth3kxy = PDstandardNth3kxy;
      
      JacPDstandardNth3kxz = PDstandardNth3kxz;
      
      JacPDstandardNth3kyy = PDstandardNth3kyy;
      
      JacPDstandardNth3kyz = PDstandardNth3kyz;
      
      JacPDstandardNth11gyy = PDstandardNth11gyy;
      
      JacPDstandardNth11gyz = PDstandardNth11gyz;
      
      JacPDstandardNth11gzz = PDstandardNth11gzz;
      
      JacPDstandardNth22gxx = PDstandardNth22gxx;
      
      JacPDstandardNth22gxz = PDstandardNth22gxz;
      
      JacPDstandardNth22gzz = PDstandardNth22gzz;
      
      JacPDstandardNth33gxx = PDstandardNth33gxx;
      
      JacPDstandardNth33gxy = PDstandardNth33gxy;
      
      JacPDstandardNth33gyy = PDstandardNth33gyy;
      
      JacPDstandardNth12gxx = PDstandardNth12gxx;
      
      JacPDstandardNth12gxy = PDstandardNth12gxy;
      
      JacPDstandardNth12gxz = PDstandardNth12gxz;
      
      JacPDstandardNth12gyy = PDstandardNth12gyy;
      
      JacPDstandardNth12gyz = PDstandardNth12gyz;
      
      JacPDstandardNth12gzz = PDstandardNth12gzz;
      
      JacPDstandardNth13gxx = PDstandardNth13gxx;
      
      JacPDstandardNth13gxy = PDstandardNth13gxy;
      
      JacPDstandardNth13gxz = PDstandardNth13gxz;
      
      JacPDstandardNth13gyy = PDstandardNth13gyy;
      
      JacPDstandardNth13gyz = PDstandardNth13gyz;
      
      JacPDstandardNth13gzz = PDstandardNth13gzz;
      
      JacPDstandardNth21gxx = PDstandardNth12gxx;
      
      JacPDstandardNth21gxy = PDstandardNth12gxy;
      
      JacPDstandardNth21gxz = PDstandardNth12gxz;
      
      JacPDstandardNth21gyy = PDstandardNth12gyy;
      
      JacPDstandardNth21gyz = PDstandardNth12gyz;
      
      JacPDstandardNth21gzz = PDstandardNth12gzz;
      
      JacPDstandardNth23gxx = PDstandardNth23gxx;
      
      JacPDstandardNth23gxy = PDstandardNth23gxy;
      
      JacPDstandardNth23gxz = PDstandardNth23gxz;
      
      JacPDstandardNth23gyy = PDstandardNth23gyy;
      
      JacPDstandardNth23gyz = PDstandardNth23gyz;
      
      JacPDstandardNth23gzz = PDstandardNth23gzz;
      
      JacPDstandardNth31gxx = PDstandardNth13gxx;
      
      JacPDstandardNth31gxy = PDstandardNth13gxy;
      
      JacPDstandardNth31gxz = PDstandardNth13gxz;
      
      JacPDstandardNth31gyy = PDstandardNth13gyy;
      
      JacPDstandardNth31gyz = PDstandardNth13gyz;
      
      JacPDstandardNth31gzz = PDstandardNth13gzz;
      
      JacPDstandardNth32gxx = PDstandardNth23gxx;
      
      JacPDstandardNth32gxy = PDstandardNth23gxy;
      
      JacPDstandardNth32gxz = PDstandardNth23gxz;
      
      JacPDstandardNth32gyy = PDstandardNth23gyy;
      
      JacPDstandardNth32gyz = PDstandardNth23gyz;
      
      JacPDstandardNth32gzz = PDstandardNth23gzz;
    }
    
    CCTK_REAL detg CCTK_ATTRIBUTE_UNUSED = 2*gxyL*gxzL*gyzL - 
      gzzL*pow(gxyL,2) + gyyL*(gxxL*gzzL - pow(gxzL,2)) - gxxL*pow(gyzL,2);
    
    CCTK_REAL gu11 CCTK_ATTRIBUTE_UNUSED = (gyyL*gzzL - 
      pow(gyzL,2))*pow(detg,-1);
    
    CCTK_REAL gu21 CCTK_ATTRIBUTE_UNUSED = (gxzL*gyzL - 
      gxyL*gzzL)*pow(detg,-1);
    
    CCTK_REAL gu31 CCTK_ATTRIBUTE_UNUSED = (-(gxzL*gyyL) + 
      gxyL*gyzL)*pow(detg,-1);
    
    CCTK_REAL gu22 CCTK_ATTRIBUTE_UNUSED = (gxxL*gzzL - 
      pow(gxzL,2))*pow(detg,-1);
    
    CCTK_REAL gu32 CCTK_ATTRIBUTE_UNUSED = (gxyL*gxzL - 
      gxxL*gyzL)*pow(detg,-1);
    
    CCTK_REAL gu33 CCTK_ATTRIBUTE_UNUSED = (gxxL*gyyL - 
      pow(gxyL,2))*pow(detg,-1);
    
    CCTK_REAL G111 CCTK_ATTRIBUTE_UNUSED = 0.5*(gu11*JacPDstandardNth1gxx 
      + 2*(gu21*JacPDstandardNth1gxy + gu31*JacPDstandardNth1gxz) - 
      gu21*JacPDstandardNth2gxx - gu31*JacPDstandardNth3gxx);
    
    CCTK_REAL G211 CCTK_ATTRIBUTE_UNUSED = 0.5*(gu21*JacPDstandardNth1gxx 
      + 2*(gu22*JacPDstandardNth1gxy + gu32*JacPDstandardNth1gxz) - 
      gu22*JacPDstandardNth2gxx - gu32*JacPDstandardNth3gxx);
    
    CCTK_REAL G311 CCTK_ATTRIBUTE_UNUSED = 0.5*(gu31*JacPDstandardNth1gxx 
      + 2*(gu32*JacPDstandardNth1gxy + gu33*JacPDstandardNth1gxz) - 
      gu32*JacPDstandardNth2gxx - gu33*JacPDstandardNth3gxx);
    
    CCTK_REAL G112 CCTK_ATTRIBUTE_UNUSED = 0.5*(gu21*JacPDstandardNth1gyy 
      + gu11*JacPDstandardNth2gxx + gu31*(JacPDstandardNth1gyz + 
      JacPDstandardNth2gxz - JacPDstandardNth3gxy));
    
    CCTK_REAL G212 CCTK_ATTRIBUTE_UNUSED = 0.5*(gu22*JacPDstandardNth1gyy 
      + gu21*JacPDstandardNth2gxx + gu32*(JacPDstandardNth1gyz + 
      JacPDstandardNth2gxz - JacPDstandardNth3gxy));
    
    CCTK_REAL G312 CCTK_ATTRIBUTE_UNUSED = 0.5*(gu32*JacPDstandardNth1gyy 
      + gu31*JacPDstandardNth2gxx + gu33*(JacPDstandardNth1gyz + 
      JacPDstandardNth2gxz - JacPDstandardNth3gxy));
    
    CCTK_REAL G113 CCTK_ATTRIBUTE_UNUSED = 0.5*(gu31*JacPDstandardNth1gzz 
      + gu11*JacPDstandardNth3gxx + gu21*(JacPDstandardNth1gyz - 
      JacPDstandardNth2gxz + JacPDstandardNth3gxy));
    
    CCTK_REAL G213 CCTK_ATTRIBUTE_UNUSED = 0.5*(gu32*JacPDstandardNth1gzz 
      + gu21*JacPDstandardNth3gxx + gu22*(JacPDstandardNth1gyz - 
      JacPDstandardNth2gxz + JacPDstandardNth3gxy));
    
    CCTK_REAL G313 CCTK_ATTRIBUTE_UNUSED = 0.5*(gu33*JacPDstandardNth1gzz 
      + gu31*JacPDstandardNth3gxx + gu32*(JacPDstandardNth1gyz - 
      JacPDstandardNth2gxz + JacPDstandardNth3gxy));
    
    CCTK_REAL G122 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gu11*(-JacPDstandardNth1gyy + 2*JacPDstandardNth2gxy) + 
      gu21*JacPDstandardNth2gyy + gu31*(2*JacPDstandardNth2gyz - 
      JacPDstandardNth3gyy));
    
    CCTK_REAL G222 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gu21*(-JacPDstandardNth1gyy + 2*JacPDstandardNth2gxy) + 
      gu22*JacPDstandardNth2gyy + gu32*(2*JacPDstandardNth2gyz - 
      JacPDstandardNth3gyy));
    
    CCTK_REAL G322 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gu31*(-JacPDstandardNth1gyy + 2*JacPDstandardNth2gxy) + 
      gu32*JacPDstandardNth2gyy + gu33*(2*JacPDstandardNth2gyz - 
      JacPDstandardNth3gyy));
    
    CCTK_REAL G123 CCTK_ATTRIBUTE_UNUSED = 0.5*(gu31*JacPDstandardNth2gzz 
      + gu11*(-JacPDstandardNth1gyz + JacPDstandardNth2gxz + 
      JacPDstandardNth3gxy) + gu21*JacPDstandardNth3gyy);
    
    CCTK_REAL G223 CCTK_ATTRIBUTE_UNUSED = 0.5*(gu32*JacPDstandardNth2gzz 
      + gu21*(-JacPDstandardNth1gyz + JacPDstandardNth2gxz + 
      JacPDstandardNth3gxy) + gu22*JacPDstandardNth3gyy);
    
    CCTK_REAL G323 CCTK_ATTRIBUTE_UNUSED = 0.5*(gu33*JacPDstandardNth2gzz 
      + gu31*(-JacPDstandardNth1gyz + JacPDstandardNth2gxz + 
      JacPDstandardNth3gxy) + gu32*JacPDstandardNth3gyy);
    
    CCTK_REAL G133 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gu11*(-JacPDstandardNth1gzz + 2*JacPDstandardNth3gxz) + 
      gu21*(-JacPDstandardNth2gzz + 2*JacPDstandardNth3gyz) + 
      gu31*JacPDstandardNth3gzz);
    
    CCTK_REAL G233 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gu21*(-JacPDstandardNth1gzz + 2*JacPDstandardNth3gxz) + 
      gu22*(-JacPDstandardNth2gzz + 2*JacPDstandardNth3gyz) + 
      gu32*JacPDstandardNth3gzz);
    
    CCTK_REAL G333 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gu31*(-JacPDstandardNth1gzz + 2*JacPDstandardNth3gxz) + 
      gu32*(-JacPDstandardNth2gzz + 2*JacPDstandardNth3gyz) + 
      gu33*JacPDstandardNth3gzz);
    
    CCTK_REAL R11 CCTK_ATTRIBUTE_UNUSED = 0.5*(-4*(gxzL*G123*G311*gu32 + 
      gyzL*G223*G311*gu32 + gyzL*G211*G323*gu32 + gzzL*G311*G323*gu32) + 
      gxyL*(-4*G123*G211*gu32 + 4*G113*G212*gu32 - 2*G133*G211*gu33) + 
      4*(gxzL*G113*G312*gu32 + gyzL*G213*G312*gu32 + gyzL*G212*G313*gu32 + 
      gzzL*G312*G313*gu32 + G112*(gxyL*G212*gu22 + gxzL*G312*gu22 + 
      gxxL*G113*gu32 + gxyL*G213*gu32 + gxzL*G313*gu32) + 
      gxyL*G113*G213*gu33) - 2*(G111*(gxxL*G122*gu22 + gxyL*G222*gu22 + 
      gxzL*G322*gu22 + 2*gxxL*G123*gu32 + 2*gxyL*G223*gu32 + 2*gxzL*G323*gu32 
      + gxxL*G133*gu33 + gxyL*G233*gu33 + gxzL*G333*gu33) + 
      gu32*JacPDstandardNth11gyz) + gu21*(-JacPDstandardNth12gxx + 
      JacPDstandardNth21gxx) + gu31*(-JacPDstandardNth13gxx + 
      JacPDstandardNth31gxx) + 2*gu32*JacPDstandardNth31gxy + 
      gu32*(2*JacPDstandardNth21gxz - JacPDstandardNth23gxx - 
      JacPDstandardNth32gxx) + gyyL*(4*G212*G213*gu32 - 4*G211*G223*gu32 + 
      2*gu33*pow(G213,2)) + gu22*(-2*gxyL*G122*G211 - 2*gyyL*G211*G222 - 
      2*gxzL*G122*G311 - 2*gyzL*G222*G311 + 4*gyzL*G212*G312 - 
      2*gyzL*G211*G322 - 2*gzzL*G311*G322 - JacPDstandardNth11gyy + 
      2*JacPDstandardNth21gxy - JacPDstandardNth22gxx + 2*gxxL*pow(G112,2) + 
      2*gyyL*pow(G212,2) + 2*gzzL*pow(G312,2)) + gu33*(-2*gyyL*G211*G233 - 
      2*gxzL*G133*G311 - 2*gyzL*G233*G311 + 4*gxzL*G113*G313 + 
      4*gyzL*G213*G313 - 2*gyzL*G211*G333 - 2*gzzL*G311*G333 - 
      JacPDstandardNth11gzz + 2*JacPDstandardNth31gxz - JacPDstandardNth33gxx 
      + 2*gxxL*pow(G113,2) + 2*gzzL*pow(G313,2)));
    
    CCTK_REAL R12 CCTK_ATTRIBUTE_UNUSED = 0.5*(2*((gxzL*G123 + 
      gyzL*G223)*G311*gu31 + G111*((gxxL*G122 + gxyL*G222 + gxzL*G322)*gu21 + 
      (gxxL*G123 + gxyL*G223 + gxzL*G323)*gu31) + (G122*(gxxL*G113 + 
      gxyL*G213 + gxzL*G313) + G222*(gxyL*G113 + gyzL*G313) + (gxzL*G113 + 
      gyzL*G213)*G322)*gu32 + (G123*(gxxL*G113 + gxyL*G213 + gxzL*G313) + 
      G223*(gxyL*G113 + gyzL*G313) + gyzL*G213*G323)*gu33 + G323*((gyzL*G211 
      + gzzL*G311)*gu31 + gxzL*G113*gu33)) + gxyL*((2*G123*G211 - 
      2*G113*G212)*gu31 - 2*G212*(G123*gu32 + G133*gu33)) + 
      gyyL*((-2*G212*G213 + 2*G211*G223)*gu31 + 2*G213*(G222*gu32 + 
      G223*gu33)) - 2*((gyzL*G212 + gzzL*G312)*G313*gu31 + (gyyL*G212 + 
      gyzL*G312)*(G223*gu32 + G233*gu33) + G112*(2*(gxyL*G212 + 
      gxzL*G312)*gu21 + (gxxL*G113 + gxyL*G213 + gxzL*G313)*gu31 + (gxxL*G123 
      + gxyL*G223 + gxzL*G323)*gu32 + (gxxL*G133 + gxyL*G233 + 
      gxzL*G333)*gu33) + G312*((gxzL*G113 + gyzL*G213)*gu31 + gxzL*(G123*gu32 
      + G133*gu33)) + gu21*JacPDstandardNth12gxy) + 
      gu22*(-JacPDstandardNth12gyy + JacPDstandardNth21gyy) + 
      gu32*(gzzL*(2*G313*G322 - 2*G312*G323) - 2*(gyzL*G212*G323 + 
      JacPDstandardNth12gyz) + JacPDstandardNth21gyz + JacPDstandardNth22gxz 
      - JacPDstandardNth23gxy + JacPDstandardNth31gyy) + 
      gu31*(JacPDstandardNth11gyz - JacPDstandardNth12gxz - 
      JacPDstandardNth13gxy + JacPDstandardNth32gxx) + 
      gu33*(-2*gyzL*G212*G333 + gzzL*(2*G313*G323 - 2*G312*G333) - 
      JacPDstandardNth12gzz + JacPDstandardNth31gyz + JacPDstandardNth32gxz - 
      JacPDstandardNth33gxy) + gu21*(-4*gyzL*G212*G312 + 2*(G211*(gxyL*G122 + 
      gyyL*G222 + gyzL*G322) + G311*(gxzL*G122 + gyzL*G222 + gzzL*G322)) + 
      JacPDstandardNth11gyy + JacPDstandardNth22gxx - 2*(gxxL*pow(G112,2) + 
      gyyL*pow(G212,2) + gzzL*pow(G312,2))));
    
    CCTK_REAL R13 CCTK_ATTRIBUTE_UNUSED = 0.5*(gxyL*(2*G123*G211*gu21 - 
      2*G113*G212*gu21 + 2*G123*G212*gu22 + 2*G133*G211*gu31 + 
      2*G133*G212*gu32) + 2*(gxzL*G123*G311*gu21 + gyzL*G223*G311*gu21 + 
      gyzL*G211*G323*gu21 + gzzL*G311*G323*gu21 + gyyL*G212*G223*gu22 + 
      gxzL*G123*G312*gu22 + gyzL*G223*G312*gu22 + gyyL*G211*G233*gu31 + 
      G111*(gxxL*G123*gu21 + gxyL*G223*gu21 + gxzL*G323*gu21 + gxxL*G133*gu31 
      + gxyL*G233*gu31 + gxzL*G333*gu31) + gyyL*G212*G233*gu32 + 
      gxzL*G133*G312*gu32 + gyzL*G233*G312*gu32 + G112*((-(gxxL*G113) - 
      gxyL*G213 - gxzL*G313)*gu21 + gxxL*G123*gu22 + gxyL*G223*gu22 + 
      gxzL*G323*gu22 + gxxL*G133*gu32 + gxyL*G233*gu32 + gxzL*G333*gu32)) + 
      gu22*(-2*gzzL*G313*G322 + 2*gyzL*G212*G323 + 2*gzzL*G312*G323 - 
      JacPDstandardNth13gyy + JacPDstandardNth21gyz - JacPDstandardNth22gxz) 
      + gu21*(JacPDstandardNth11gyz - JacPDstandardNth12gxz - 
      JacPDstandardNth13gxy + JacPDstandardNth23gxx) + 
      gu22*JacPDstandardNth23gxy + gu32*JacPDstandardNth31gyz + 
      gu33*(-JacPDstandardNth13gzz + JacPDstandardNth31gzz) + 
      gu32*(-2*gzzL*G313*G323 + 2*gyzL*G212*G333 + 2*gzzL*G312*G333 - 
      2*JacPDstandardNth13gyz + JacPDstandardNth21gzz - 
      JacPDstandardNth32gxz) + gu32*JacPDstandardNth33gxy - 
      2*(gxzL*G113*G312*gu21 + gyzL*G213*G312*gu21 + gyzL*G212*G313*gu21 + 
      gzzL*G312*G313*gu21 + gxxL*G113*G122*gu22 + gxyL*G122*G213*gu22 + 
      gxyL*G113*G222*gu22 + gxzL*G122*G313*gu22 + gyzL*G222*G313*gu22 + 
      gxzL*G113*G322*gu22 + gyzL*G213*G322*gu22 + gxxL*G113*G123*gu32 + 
      gxyL*G123*G213*gu32 + gxyL*G113*G223*gu32 + gxzL*G123*G313*gu32 + 
      gyzL*G223*G313*gu32 + gxzL*G113*G323*gu32 + gyzL*G213*G323*gu32 + 
      gu31*JacPDstandardNth13gxz + gxxL*gu31*pow(G113,2)) + 
      gyyL*(-2*G212*G213*gu21 + 2*G211*G223*gu21 - 2*G213*G222*gu22 - 
      2*G213*G223*gu32 - 2*gu31*pow(G213,2)) + gu31*(-4*gxyL*G113*G213 + 
      2*gxzL*G133*G311 + 2*gyzL*G233*G311 - 4*gxzL*G113*G313 - 
      4*gyzL*G213*G313 + 2*gyzL*G211*G333 + 2*gzzL*G311*G333 + 
      JacPDstandardNth11gzz + JacPDstandardNth33gxx - 2*gzzL*pow(G313,2)));
    
    CCTK_REAL R22 CCTK_ATTRIBUTE_UNUSED = 0.5*(4*(G112*((gxyL*G212 + 
      gxzL*G312)*gu11 + (gxxL*G123 + gxyL*G223 + gxzL*G323)*gu31) + 
      gxyL*G123*G223*gu33) + gu21*(JacPDstandardNth12gyy - 
      JacPDstandardNth21gyy) - 2*((gzzL*G311*G322 + G111*(gxxL*G122 + 
      gxyL*G222 + gxzL*G322))*gu11 + (gxxL*G122*G133 + (gyzL*G222 + 
      gzzL*G322)*G333)*gu33 + gu31*JacPDstandardNth22gxz) + 
      gu31*(gxyL*(4*G123*G212 - 4*G113*G222) - 4*(gxxL*G113*G122 + 
      gxyL*G122*G213 + gyyL*G213*G222 + gxzL*G122*G313 + gyzL*G222*G313 + 
      gxzL*G113*G322 + (gyzL*G213 + gzzL*G313)*G322) + 4*(gyyL*G212*G223 + 
      gxzL*G123*G312 + gyzL*G223*G312 + gyzL*G212*G323 + gzzL*G312*G323) - 
      JacPDstandardNth13gyy - JacPDstandardNth31gyy + 
      2*(JacPDstandardNth12gyz + JacPDstandardNth32gxy)) + 
      gu32*(-JacPDstandardNth23gyy + JacPDstandardNth32gyy) + 
      gu11*(4*gyzL*G212*G312 - 2*(gxyL*G122*G211 + gyyL*G211*G222 + 
      gxzL*G122*G311 + gyzL*G222*G311 + gyzL*G211*G322) - 
      JacPDstandardNth11gyy - JacPDstandardNth22gxx + 
      2*(JacPDstandardNth12gxy + gxxL*pow(G112,2) + gyyL*pow(G212,2) + 
      gzzL*pow(G312,2))) + gu33*(4*(gxzL*G123*G323 + gyzL*G223*G323) - 
      2*(gxyL*G133*G222 + gxyL*G122*G233 + gyyL*G222*G233 + gxzL*G133*G322 + 
      gyzL*G233*G322 + gxzL*G122*G333) - JacPDstandardNth22gzz - 
      JacPDstandardNth33gyy + 2*(JacPDstandardNth32gyz + gxxL*pow(G123,2) + 
      gyyL*pow(G223,2) + gzzL*pow(G323,2))));
    
    CCTK_REAL R23 CCTK_ATTRIBUTE_UNUSED = 0.5*(gxyL*((-2*G123*G211 + 
      2*G113*G212)*gu11 + G212*(-2*G123*gu21 + 2*G133*gu31) + (-4*G123*G223 + 
      2*(G133*G222 + G122*G233))*gu32) + 2*(G313*((gyzL*G212 + 
      gzzL*G312)*gu11 + gxzL*G122*gu21) + G233*(gyyL*G212 + gyzL*G312)*gu31 + 
      G312*((gxzL*G113 + gyzL*G213)*gu11 + gxzL*G133*gu31) + G112*((gxxL*G113 
      + gxyL*G213 + gxzL*G313)*gu11 + (-(gxxL*G123) - gxyL*G223 - 
      gxzL*G323)*gu21 + gxxL*G133*gu31 + gxyL*G233*gu31 + gxzL*G333*gu31) + 
      (gyzL*G222 + gzzL*G322)*G333*gu32 + G122*((gxxL*G113 + gxyL*G213)*gu21 
      + gxxL*G133*gu32) + G322*((gxzL*G113 + gyzL*G213)*gu21 + 
      gxzL*G133*gu32) + G222*((gxyL*G113 + gyzL*G313)*gu21 + gyyL*G233*gu32)) 
      + gu11*(-JacPDstandardNth11gyz + JacPDstandardNth12gxz + 
      JacPDstandardNth13gxy - JacPDstandardNth23gxx) + 
      gu21*(-2*gyzL*G212*G323 + gzzL*(2*G313*G322 - 2*G312*G323) + 
      JacPDstandardNth13gyy - JacPDstandardNth21gyz + JacPDstandardNth22gxz - 
      JacPDstandardNth23gxy) + gu33*(-JacPDstandardNth23gzz + 
      JacPDstandardNth32gzz) + gu31*(2*gyzL*G212*G333 + gzzL*(-2*G313*G323 + 
      2*G312*G333) + JacPDstandardNth12gzz - JacPDstandardNth31gyz + 
      JacPDstandardNth32gxz + JacPDstandardNth33gxy) - 2*(((gxzL*G123 + 
      gyzL*G223)*G311 + G111*(gxxL*G123 + gxyL*G223 + gxzL*G323))*gu11 + 
      (gxzL*G123*G312 + G223*(gyyL*G212 + gyzL*G312))*gu21 + G323*((gyzL*G211 
      + gzzL*G311)*gu11 + gxzL*G113*gu31) + gu31*(G123*(gxxL*G113 + gxyL*G213 
      + gxzL*G313) + G223*(gxyL*G113 + gyzL*G313) + gyzL*G213*G323 + 
      JacPDstandardNth23gxz) + gxxL*gu32*pow(G123,2)) + gyyL*((2*G212*G213 - 
      2*G211*G223)*gu11 + G213*(2*G222*gu21 - 2*G223*gu31) - 
      2*gu32*pow(G223,2)) + gu32*(gyzL*(2*G233*G322 - 4*G223*G323) + 
      gxzL*(-4*G123*G323 + 2*G122*G333) + JacPDstandardNth22gzz + 
      JacPDstandardNth33gyy - 2*(JacPDstandardNth23gyz + gzzL*pow(G323,2))));
    
    CCTK_REAL R33 CCTK_ATTRIBUTE_UNUSED = 0.5*(4*((G213*(gxyL*G123 + 
      gyyL*G223) + (gxzL*G123 + gyzL*G223)*G313 + (gyzL*G213 + 
      gzzL*G313)*G323)*gu21 + G113*((gxyL*G213 + gxzL*G313)*gu11 + (gxxL*G123 
      + gxyL*G223 + gxzL*G323)*gu21) + gxyL*G123*G223*gu22) + 
      gu21*(-4*(G133*(gxxL*G112 + gxyL*G212) + (gxyL*G112 + gyyL*G212)*G233 + 
      (gxzL*G133 + gyzL*G233)*G312 + (gxzL*G112 + gyzL*G212 + 
      gzzL*G312)*G333) - JacPDstandardNth12gzz + 2*JacPDstandardNth13gyz - 
      JacPDstandardNth21gzz + 2*JacPDstandardNth23gxz) + 
      gu31*(JacPDstandardNth13gzz - JacPDstandardNth31gzz) + 
      gu32*(JacPDstandardNth23gzz - JacPDstandardNth32gzz) - 
      2*((gzzL*G311*G333 + G111*(gxxL*G133 + gxyL*G233 + gxzL*G333))*gu11 + 
      (gxxL*G122*G133 + (gyzL*G222 + gzzL*G322)*G333)*gu22 + 
      gu21*JacPDstandardNth33gxy) + gu11*(4*gyzL*G213*G313 - 
      2*(gxyL*G133*G211 + gyyL*G211*G233 + gxzL*G133*G311 + gyzL*G233*G311 + 
      gyzL*G211*G333) - JacPDstandardNth11gzz - JacPDstandardNth33gxx + 
      2*(JacPDstandardNth13gxz + gxxL*pow(G113,2) + gyyL*pow(G213,2) + 
      gzzL*pow(G313,2))) + gu22*(4*(gxzL*G123*G323 + gyzL*G223*G323) - 
      2*(gxyL*G133*G222 + gxyL*G122*G233 + gyyL*G222*G233 + gxzL*G133*G322 + 
      gyzL*G233*G322 + gxzL*G122*G333) - JacPDstandardNth22gzz - 
      JacPDstandardNth33gyy + 2*(JacPDstandardNth23gyz + gxxL*pow(G123,2) + 
      gyyL*pow(G223,2) + gzzL*pow(G323,2))));
    
    CCTK_REAL trR CCTK_ATTRIBUTE_UNUSED = gu11*R11 + gu22*R22 + 
      2*(gu21*R12 + gu31*R13 + gu32*R23) + gu33*R33;
    
    CCTK_REAL Km11 CCTK_ATTRIBUTE_UNUSED = kxxL*gu11 + kxyL*gu21 + 
      kxzL*gu31;
    
    CCTK_REAL Km21 CCTK_ATTRIBUTE_UNUSED = kxxL*gu21 + kxyL*gu22 + 
      kxzL*gu32;
    
    CCTK_REAL Km31 CCTK_ATTRIBUTE_UNUSED = kxxL*gu31 + kxyL*gu32 + 
      kxzL*gu33;
    
    CCTK_REAL Km12 CCTK_ATTRIBUTE_UNUSED = kxyL*gu11 + kyyL*gu21 + 
      kyzL*gu31;
    
    CCTK_REAL Km22 CCTK_ATTRIBUTE_UNUSED = kxyL*gu21 + kyyL*gu22 + 
      kyzL*gu32;
    
    CCTK_REAL Km32 CCTK_ATTRIBUTE_UNUSED = kxyL*gu31 + kyyL*gu32 + 
      kyzL*gu33;
    
    CCTK_REAL Km13 CCTK_ATTRIBUTE_UNUSED = kxzL*gu11 + kyzL*gu21 + 
      kzzL*gu31;
    
    CCTK_REAL Km23 CCTK_ATTRIBUTE_UNUSED = kxzL*gu21 + kyzL*gu22 + 
      kzzL*gu32;
    
    CCTK_REAL Km33 CCTK_ATTRIBUTE_UNUSED = kxzL*gu31 + kyzL*gu32 + 
      kzzL*gu33;
    
    CCTK_REAL trK CCTK_ATTRIBUTE_UNUSED = Km11 + Km22 + Km33;
    
    CCTK_REAL rho CCTK_ATTRIBUTE_UNUSED = pow(alpL,-2)*(eTttL - 
      2*(betayL*eTtyL + betazL*eTtzL) + 2*(betaxL*(-eTtxL + betayL*eTxyL + 
      betazL*eTxzL) + betayL*betazL*eTyzL) + eTxxL*pow(betaxL,2) + 
      eTyyL*pow(betayL,2) + eTzzL*pow(betazL,2));
    
    CCTK_REAL S1 CCTK_ATTRIBUTE_UNUSED = (-eTtxL + betaxL*eTxxL + 
      betayL*eTxyL + betazL*eTxzL)*pow(alpL,-1);
    
    CCTK_REAL S2 CCTK_ATTRIBUTE_UNUSED = (-eTtyL + betaxL*eTxyL + 
      betayL*eTyyL + betazL*eTyzL)*pow(alpL,-1);
    
    CCTK_REAL S3 CCTK_ATTRIBUTE_UNUSED = (-eTtzL + betaxL*eTxzL + 
      betayL*eTyzL + betazL*eTzzL)*pow(alpL,-1);
    
    CCTK_REAL HL CCTK_ATTRIBUTE_UNUSED = -2*Km12*Km21 - 2*Km13*Km31 - 
      2*Km23*Km32 - 16*Pi*rho + trR - pow(Km11,2) - pow(Km22,2) - pow(Km33,2) 
      + pow(trK,2);
    
    CCTK_REAL M1L CCTK_ATTRIBUTE_UNUSED = gu21*(-(kxxL*G112) + kyyL*G211 + 
      kxyL*(G111 - G212) + kyzL*G311 - kxzL*G312 - JacPDstandardNth1kxy + 
      JacPDstandardNth2kxx) + gu22*(-(kxxL*G122) + kyyL*G212 + kxyL*(G112 - 
      G222) + kyzL*G312 - kxzL*G322 - JacPDstandardNth1kyy + 
      JacPDstandardNth2kxy) + gu31*(-(kxxL*G113) + kyzL*G211 - kxyL*G213 + 
      kzzL*G311 + kxzL*(G111 - G313) - JacPDstandardNth1kxz + 
      JacPDstandardNth3kxx) + gu32*(kyyL*G213 + kxyL*(G113 - 2*G223) + 
      kzzL*G312 + kyzL*(G212 + G313) + kxzL*(G112 - 2*G323) - 2*(kxxL*G123 + 
      JacPDstandardNth1kyz) + JacPDstandardNth2kxz + JacPDstandardNth3kxy) + 
      gu33*(-(kxxL*G133) + kyzL*G213 - kxyL*G233 + kzzL*G313 + kxzL*(G113 - 
      G333) - JacPDstandardNth1kzz + JacPDstandardNth3kxz) - 8*Pi*S1;
    
    CCTK_REAL M2L CCTK_ATTRIBUTE_UNUSED = gu11*(kxxL*G112 - kyyL*G211 + 
      kxyL*(-G111 + G212) - kyzL*G311 + kxzL*G312 + JacPDstandardNth1kxy - 
      JacPDstandardNth2kxx) + gu21*(kxxL*G122 - kyyL*G212 + kxyL*(-G112 + 
      G222) - kyzL*G312 + kxzL*G322 + JacPDstandardNth1kyy - 
      JacPDstandardNth2kxy) + gu31*(kxxL*G123 + kxyL*G223 + kzzL*G312 + 
      kyzL*(G212 - 2*G313) + kxzL*(G112 + G323) + JacPDstandardNth1kyz - 
      2*(kxyL*G113 + kyyL*G213 + JacPDstandardNth2kxz) + 
      JacPDstandardNth3kxy) + gu32*(kxzL*G122 - kxyL*G123 - kyyL*G223 + 
      kzzL*G322 + kyzL*(G222 - G323) - JacPDstandardNth2kyz + 
      JacPDstandardNth3kyy) + gu33*(kxzL*G123 - kxyL*G133 - kyyL*G233 + 
      kzzL*G323 + kyzL*(G223 - G333) - JacPDstandardNth2kzz + 
      JacPDstandardNth3kyz) - 8*Pi*S2;
    
    CCTK_REAL M3L CCTK_ATTRIBUTE_UNUSED = gu11*(kxxL*G113 - kyzL*G211 + 
      kxyL*G213 - kzzL*G311 + kxzL*(-G111 + G313) + JacPDstandardNth1kxz - 
      JacPDstandardNth3kxx) + gu21*(kxxL*G123 + kyyL*G213 + kxyL*(G113 + 
      G223) + kyzL*G313 + kxzL*G323 + JacPDstandardNth1kyz + 
      JacPDstandardNth2kxz - 2*(kxzL*G112 + kyzL*G212 + kzzL*G312 + 
      JacPDstandardNth3kxy)) + gu31*(kxxL*G133 - kyzL*G213 + kxyL*G233 - 
      kzzL*G313 + kxzL*(-G113 + G333) + JacPDstandardNth1kzz - 
      JacPDstandardNth3kxz) + gu22*(-(kxzL*G122) + kxyL*G123 + kyyL*G223 - 
      kzzL*G322 + kyzL*(-G222 + G323) + JacPDstandardNth2kyz - 
      JacPDstandardNth3kyy) + gu32*(-(kxzL*G123) + kxyL*G133 + kyyL*G233 - 
      kzzL*G323 + kyzL*(-G223 + G333) + JacPDstandardNth2kzz - 
      JacPDstandardNth3kyz) - 8*Pi*S3;
    /* Copy local copies back to grid functions */
    H[index] = HL;
    M1[index] = M1L;
    M2[index] = M2L;
    M3[index] = M3L;
  }
  CCTK_ENDLOOP3(ML_ADMConstraints_evaluate);
}
extern "C" void ML_ADMConstraints_evaluate(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_ML_ADMConstraints_evaluate
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_ADMConstraints_evaluate);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_ADMConstraints_evaluate_Body");
  }
  if (cctk_iteration % ML_ADMConstraints_evaluate_calc_every != ML_ADMConstraints_evaluate_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "ADMBase::curv",
    "ADMBase::lapse",
    "ADMBase::metric",
    "ADMBase::shift",
    "ML_ADMConstraints::ML_Ham",
    "ML_ADMConstraints::ML_mom"};
  AssertGroupStorage(cctkGH, "ML_ADMConstraints_evaluate", 6, groups);
  
  EnsureStencilFits(cctkGH, "ML_ADMConstraints_evaluate", 2, 2, 2);
  
  LoopOverInterior(cctkGH, ML_ADMConstraints_evaluate_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_ADMConstraints_evaluate_Body");
  }
}

} // namespace ML_ADMConstraints
