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

namespace Kranc2BSSNChiMatter {

extern "C" void bssnchimatter_calc_constraints_matter_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % bssnchimatter_calc_constraints_matter_calc_every != bssnchimatter_calc_constraints_matter_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "Kranc2BSSNChiMatter::bssnCcons_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for Kranc2BSSNChiMatter::bssnCcons_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "Kranc2BSSNChiMatter::bssnmom_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for Kranc2BSSNChiMatter::bssnmom_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "Kranc2BSSNChiMatter::scalarconstraints","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for Kranc2BSSNChiMatter::scalarconstraints.");
  return;
}

static void bssnchimatter_calc_constraints_matter_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  CCTK_LOOP3(bssnchimatter_calc_constraints_matter,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL A11L CCTK_ATTRIBUTE_UNUSED = A11[index];
    CCTK_REAL A21L CCTK_ATTRIBUTE_UNUSED = A21[index];
    CCTK_REAL A22L CCTK_ATTRIBUTE_UNUSED = A22[index];
    CCTK_REAL A31L CCTK_ATTRIBUTE_UNUSED = A31[index];
    CCTK_REAL A32L CCTK_ATTRIBUTE_UNUSED = A32[index];
    CCTK_REAL A33L CCTK_ATTRIBUTE_UNUSED = A33[index];
    CCTK_REAL alphaL CCTK_ATTRIBUTE_UNUSED = alpha[index];
    CCTK_REAL beta1L CCTK_ATTRIBUTE_UNUSED = beta1[index];
    CCTK_REAL beta2L CCTK_ATTRIBUTE_UNUSED = beta2[index];
    CCTK_REAL beta3L CCTK_ATTRIBUTE_UNUSED = beta3[index];
    CCTK_REAL bssnhamL CCTK_ATTRIBUTE_UNUSED = bssnham[index];
    CCTK_REAL bssnmom1L CCTK_ATTRIBUTE_UNUSED = bssnmom1[index];
    CCTK_REAL bssnmom2L CCTK_ATTRIBUTE_UNUSED = bssnmom2[index];
    CCTK_REAL bssnmom3L CCTK_ATTRIBUTE_UNUSED = bssnmom3[index];
    CCTK_REAL chiL CCTK_ATTRIBUTE_UNUSED = chi[index];
    CCTK_REAL Gam1L CCTK_ATTRIBUTE_UNUSED = Gam1[index];
    CCTK_REAL Gam2L CCTK_ATTRIBUTE_UNUSED = Gam2[index];
    CCTK_REAL Gam3L CCTK_ATTRIBUTE_UNUSED = Gam3[index];
    CCTK_REAL h11L CCTK_ATTRIBUTE_UNUSED = h11[index];
    CCTK_REAL h21L CCTK_ATTRIBUTE_UNUSED = h21[index];
    CCTK_REAL h22L CCTK_ATTRIBUTE_UNUSED = h22[index];
    CCTK_REAL h31L CCTK_ATTRIBUTE_UNUSED = h31[index];
    CCTK_REAL h32L CCTK_ATTRIBUTE_UNUSED = h32[index];
    CCTK_REAL h33L CCTK_ATTRIBUTE_UNUSED = h33[index];
    CCTK_REAL KL CCTK_ATTRIBUTE_UNUSED = K[index];
    
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
    CCTK_REAL PDstandardNth1A11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2A11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3A11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1A21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2A21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3A21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1A22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2A22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3A22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1A31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2A31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3A31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1A32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2A32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3A32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1A33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2A33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3A33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1chi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2chi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3chi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth11chi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth22chi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth33chi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth12chi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth13chi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth23chi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1Gam1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2Gam1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3Gam1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1Gam2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2Gam2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3Gam2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1Gam3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2Gam3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3Gam3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1h11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2h11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3h11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth11h11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth22h11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth33h11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth12h11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth13h11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth23h11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1h21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2h21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3h21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth11h21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth22h21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth33h21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth12h21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth13h21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth23h21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1h22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2h22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3h22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth11h22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth22h22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth33h22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth12h22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth13h22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth23h22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1h31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2h31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3h31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth11h31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth22h31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth33h31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth12h31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth13h31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth23h31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1h32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2h32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3h32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth11h32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth22h32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth33h32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth12h32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth13h32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth23h32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1h33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2h33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3h33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth11h33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth22h33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth33h33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth12h33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth13h33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth23h33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1K CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2K CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3K CCTK_ATTRIBUTE_UNUSED;
    
    switch (constraintsfdorder)
    {
      case 2:
      {
        PDstandardNth1A11 = PDstandardNthconstraintsfdorder21(&A11[index]);
        PDstandardNth2A11 = PDstandardNthconstraintsfdorder22(&A11[index]);
        PDstandardNth3A11 = PDstandardNthconstraintsfdorder23(&A11[index]);
        PDstandardNth1A21 = PDstandardNthconstraintsfdorder21(&A21[index]);
        PDstandardNth2A21 = PDstandardNthconstraintsfdorder22(&A21[index]);
        PDstandardNth3A21 = PDstandardNthconstraintsfdorder23(&A21[index]);
        PDstandardNth1A22 = PDstandardNthconstraintsfdorder21(&A22[index]);
        PDstandardNth2A22 = PDstandardNthconstraintsfdorder22(&A22[index]);
        PDstandardNth3A22 = PDstandardNthconstraintsfdorder23(&A22[index]);
        PDstandardNth1A31 = PDstandardNthconstraintsfdorder21(&A31[index]);
        PDstandardNth2A31 = PDstandardNthconstraintsfdorder22(&A31[index]);
        PDstandardNth3A31 = PDstandardNthconstraintsfdorder23(&A31[index]);
        PDstandardNth1A32 = PDstandardNthconstraintsfdorder21(&A32[index]);
        PDstandardNth2A32 = PDstandardNthconstraintsfdorder22(&A32[index]);
        PDstandardNth3A32 = PDstandardNthconstraintsfdorder23(&A32[index]);
        PDstandardNth1A33 = PDstandardNthconstraintsfdorder21(&A33[index]);
        PDstandardNth2A33 = PDstandardNthconstraintsfdorder22(&A33[index]);
        PDstandardNth3A33 = PDstandardNthconstraintsfdorder23(&A33[index]);
        PDstandardNth1chi = PDstandardNthconstraintsfdorder21(&chi[index]);
        PDstandardNth2chi = PDstandardNthconstraintsfdorder22(&chi[index]);
        PDstandardNth3chi = PDstandardNthconstraintsfdorder23(&chi[index]);
        PDstandardNth11chi = PDstandardNthconstraintsfdorder211(&chi[index]);
        PDstandardNth22chi = PDstandardNthconstraintsfdorder222(&chi[index]);
        PDstandardNth33chi = PDstandardNthconstraintsfdorder233(&chi[index]);
        PDstandardNth12chi = PDstandardNthconstraintsfdorder212(&chi[index]);
        PDstandardNth13chi = PDstandardNthconstraintsfdorder213(&chi[index]);
        PDstandardNth23chi = PDstandardNthconstraintsfdorder223(&chi[index]);
        PDstandardNth1Gam1 = PDstandardNthconstraintsfdorder21(&Gam1[index]);
        PDstandardNth2Gam1 = PDstandardNthconstraintsfdorder22(&Gam1[index]);
        PDstandardNth3Gam1 = PDstandardNthconstraintsfdorder23(&Gam1[index]);
        PDstandardNth1Gam2 = PDstandardNthconstraintsfdorder21(&Gam2[index]);
        PDstandardNth2Gam2 = PDstandardNthconstraintsfdorder22(&Gam2[index]);
        PDstandardNth3Gam2 = PDstandardNthconstraintsfdorder23(&Gam2[index]);
        PDstandardNth1Gam3 = PDstandardNthconstraintsfdorder21(&Gam3[index]);
        PDstandardNth2Gam3 = PDstandardNthconstraintsfdorder22(&Gam3[index]);
        PDstandardNth3Gam3 = PDstandardNthconstraintsfdorder23(&Gam3[index]);
        PDstandardNth1h11 = PDstandardNthconstraintsfdorder21(&h11[index]);
        PDstandardNth2h11 = PDstandardNthconstraintsfdorder22(&h11[index]);
        PDstandardNth3h11 = PDstandardNthconstraintsfdorder23(&h11[index]);
        PDstandardNth11h11 = PDstandardNthconstraintsfdorder211(&h11[index]);
        PDstandardNth22h11 = PDstandardNthconstraintsfdorder222(&h11[index]);
        PDstandardNth33h11 = PDstandardNthconstraintsfdorder233(&h11[index]);
        PDstandardNth12h11 = PDstandardNthconstraintsfdorder212(&h11[index]);
        PDstandardNth13h11 = PDstandardNthconstraintsfdorder213(&h11[index]);
        PDstandardNth23h11 = PDstandardNthconstraintsfdorder223(&h11[index]);
        PDstandardNth1h21 = PDstandardNthconstraintsfdorder21(&h21[index]);
        PDstandardNth2h21 = PDstandardNthconstraintsfdorder22(&h21[index]);
        PDstandardNth3h21 = PDstandardNthconstraintsfdorder23(&h21[index]);
        PDstandardNth11h21 = PDstandardNthconstraintsfdorder211(&h21[index]);
        PDstandardNth22h21 = PDstandardNthconstraintsfdorder222(&h21[index]);
        PDstandardNth33h21 = PDstandardNthconstraintsfdorder233(&h21[index]);
        PDstandardNth12h21 = PDstandardNthconstraintsfdorder212(&h21[index]);
        PDstandardNth13h21 = PDstandardNthconstraintsfdorder213(&h21[index]);
        PDstandardNth23h21 = PDstandardNthconstraintsfdorder223(&h21[index]);
        PDstandardNth1h22 = PDstandardNthconstraintsfdorder21(&h22[index]);
        PDstandardNth2h22 = PDstandardNthconstraintsfdorder22(&h22[index]);
        PDstandardNth3h22 = PDstandardNthconstraintsfdorder23(&h22[index]);
        PDstandardNth11h22 = PDstandardNthconstraintsfdorder211(&h22[index]);
        PDstandardNth22h22 = PDstandardNthconstraintsfdorder222(&h22[index]);
        PDstandardNth33h22 = PDstandardNthconstraintsfdorder233(&h22[index]);
        PDstandardNth12h22 = PDstandardNthconstraintsfdorder212(&h22[index]);
        PDstandardNth13h22 = PDstandardNthconstraintsfdorder213(&h22[index]);
        PDstandardNth23h22 = PDstandardNthconstraintsfdorder223(&h22[index]);
        PDstandardNth1h31 = PDstandardNthconstraintsfdorder21(&h31[index]);
        PDstandardNth2h31 = PDstandardNthconstraintsfdorder22(&h31[index]);
        PDstandardNth3h31 = PDstandardNthconstraintsfdorder23(&h31[index]);
        PDstandardNth11h31 = PDstandardNthconstraintsfdorder211(&h31[index]);
        PDstandardNth22h31 = PDstandardNthconstraintsfdorder222(&h31[index]);
        PDstandardNth33h31 = PDstandardNthconstraintsfdorder233(&h31[index]);
        PDstandardNth12h31 = PDstandardNthconstraintsfdorder212(&h31[index]);
        PDstandardNth13h31 = PDstandardNthconstraintsfdorder213(&h31[index]);
        PDstandardNth23h31 = PDstandardNthconstraintsfdorder223(&h31[index]);
        PDstandardNth1h32 = PDstandardNthconstraintsfdorder21(&h32[index]);
        PDstandardNth2h32 = PDstandardNthconstraintsfdorder22(&h32[index]);
        PDstandardNth3h32 = PDstandardNthconstraintsfdorder23(&h32[index]);
        PDstandardNth11h32 = PDstandardNthconstraintsfdorder211(&h32[index]);
        PDstandardNth22h32 = PDstandardNthconstraintsfdorder222(&h32[index]);
        PDstandardNth33h32 = PDstandardNthconstraintsfdorder233(&h32[index]);
        PDstandardNth12h32 = PDstandardNthconstraintsfdorder212(&h32[index]);
        PDstandardNth13h32 = PDstandardNthconstraintsfdorder213(&h32[index]);
        PDstandardNth23h32 = PDstandardNthconstraintsfdorder223(&h32[index]);
        PDstandardNth1h33 = PDstandardNthconstraintsfdorder21(&h33[index]);
        PDstandardNth2h33 = PDstandardNthconstraintsfdorder22(&h33[index]);
        PDstandardNth3h33 = PDstandardNthconstraintsfdorder23(&h33[index]);
        PDstandardNth11h33 = PDstandardNthconstraintsfdorder211(&h33[index]);
        PDstandardNth22h33 = PDstandardNthconstraintsfdorder222(&h33[index]);
        PDstandardNth33h33 = PDstandardNthconstraintsfdorder233(&h33[index]);
        PDstandardNth12h33 = PDstandardNthconstraintsfdorder212(&h33[index]);
        PDstandardNth13h33 = PDstandardNthconstraintsfdorder213(&h33[index]);
        PDstandardNth23h33 = PDstandardNthconstraintsfdorder223(&h33[index]);
        PDstandardNth1K = PDstandardNthconstraintsfdorder21(&K[index]);
        PDstandardNth2K = PDstandardNthconstraintsfdorder22(&K[index]);
        PDstandardNth3K = PDstandardNthconstraintsfdorder23(&K[index]);
        break;
      }
      
      case 4:
      {
        PDstandardNth1A11 = PDstandardNthconstraintsfdorder41(&A11[index]);
        PDstandardNth2A11 = PDstandardNthconstraintsfdorder42(&A11[index]);
        PDstandardNth3A11 = PDstandardNthconstraintsfdorder43(&A11[index]);
        PDstandardNth1A21 = PDstandardNthconstraintsfdorder41(&A21[index]);
        PDstandardNth2A21 = PDstandardNthconstraintsfdorder42(&A21[index]);
        PDstandardNth3A21 = PDstandardNthconstraintsfdorder43(&A21[index]);
        PDstandardNth1A22 = PDstandardNthconstraintsfdorder41(&A22[index]);
        PDstandardNth2A22 = PDstandardNthconstraintsfdorder42(&A22[index]);
        PDstandardNth3A22 = PDstandardNthconstraintsfdorder43(&A22[index]);
        PDstandardNth1A31 = PDstandardNthconstraintsfdorder41(&A31[index]);
        PDstandardNth2A31 = PDstandardNthconstraintsfdorder42(&A31[index]);
        PDstandardNth3A31 = PDstandardNthconstraintsfdorder43(&A31[index]);
        PDstandardNth1A32 = PDstandardNthconstraintsfdorder41(&A32[index]);
        PDstandardNth2A32 = PDstandardNthconstraintsfdorder42(&A32[index]);
        PDstandardNth3A32 = PDstandardNthconstraintsfdorder43(&A32[index]);
        PDstandardNth1A33 = PDstandardNthconstraintsfdorder41(&A33[index]);
        PDstandardNth2A33 = PDstandardNthconstraintsfdorder42(&A33[index]);
        PDstandardNth3A33 = PDstandardNthconstraintsfdorder43(&A33[index]);
        PDstandardNth1chi = PDstandardNthconstraintsfdorder41(&chi[index]);
        PDstandardNth2chi = PDstandardNthconstraintsfdorder42(&chi[index]);
        PDstandardNth3chi = PDstandardNthconstraintsfdorder43(&chi[index]);
        PDstandardNth11chi = PDstandardNthconstraintsfdorder411(&chi[index]);
        PDstandardNth22chi = PDstandardNthconstraintsfdorder422(&chi[index]);
        PDstandardNth33chi = PDstandardNthconstraintsfdorder433(&chi[index]);
        PDstandardNth12chi = PDstandardNthconstraintsfdorder412(&chi[index]);
        PDstandardNth13chi = PDstandardNthconstraintsfdorder413(&chi[index]);
        PDstandardNth23chi = PDstandardNthconstraintsfdorder423(&chi[index]);
        PDstandardNth1Gam1 = PDstandardNthconstraintsfdorder41(&Gam1[index]);
        PDstandardNth2Gam1 = PDstandardNthconstraintsfdorder42(&Gam1[index]);
        PDstandardNth3Gam1 = PDstandardNthconstraintsfdorder43(&Gam1[index]);
        PDstandardNth1Gam2 = PDstandardNthconstraintsfdorder41(&Gam2[index]);
        PDstandardNth2Gam2 = PDstandardNthconstraintsfdorder42(&Gam2[index]);
        PDstandardNth3Gam2 = PDstandardNthconstraintsfdorder43(&Gam2[index]);
        PDstandardNth1Gam3 = PDstandardNthconstraintsfdorder41(&Gam3[index]);
        PDstandardNth2Gam3 = PDstandardNthconstraintsfdorder42(&Gam3[index]);
        PDstandardNth3Gam3 = PDstandardNthconstraintsfdorder43(&Gam3[index]);
        PDstandardNth1h11 = PDstandardNthconstraintsfdorder41(&h11[index]);
        PDstandardNth2h11 = PDstandardNthconstraintsfdorder42(&h11[index]);
        PDstandardNth3h11 = PDstandardNthconstraintsfdorder43(&h11[index]);
        PDstandardNth11h11 = PDstandardNthconstraintsfdorder411(&h11[index]);
        PDstandardNth22h11 = PDstandardNthconstraintsfdorder422(&h11[index]);
        PDstandardNth33h11 = PDstandardNthconstraintsfdorder433(&h11[index]);
        PDstandardNth12h11 = PDstandardNthconstraintsfdorder412(&h11[index]);
        PDstandardNth13h11 = PDstandardNthconstraintsfdorder413(&h11[index]);
        PDstandardNth23h11 = PDstandardNthconstraintsfdorder423(&h11[index]);
        PDstandardNth1h21 = PDstandardNthconstraintsfdorder41(&h21[index]);
        PDstandardNth2h21 = PDstandardNthconstraintsfdorder42(&h21[index]);
        PDstandardNth3h21 = PDstandardNthconstraintsfdorder43(&h21[index]);
        PDstandardNth11h21 = PDstandardNthconstraintsfdorder411(&h21[index]);
        PDstandardNth22h21 = PDstandardNthconstraintsfdorder422(&h21[index]);
        PDstandardNth33h21 = PDstandardNthconstraintsfdorder433(&h21[index]);
        PDstandardNth12h21 = PDstandardNthconstraintsfdorder412(&h21[index]);
        PDstandardNth13h21 = PDstandardNthconstraintsfdorder413(&h21[index]);
        PDstandardNth23h21 = PDstandardNthconstraintsfdorder423(&h21[index]);
        PDstandardNth1h22 = PDstandardNthconstraintsfdorder41(&h22[index]);
        PDstandardNth2h22 = PDstandardNthconstraintsfdorder42(&h22[index]);
        PDstandardNth3h22 = PDstandardNthconstraintsfdorder43(&h22[index]);
        PDstandardNth11h22 = PDstandardNthconstraintsfdorder411(&h22[index]);
        PDstandardNth22h22 = PDstandardNthconstraintsfdorder422(&h22[index]);
        PDstandardNth33h22 = PDstandardNthconstraintsfdorder433(&h22[index]);
        PDstandardNth12h22 = PDstandardNthconstraintsfdorder412(&h22[index]);
        PDstandardNth13h22 = PDstandardNthconstraintsfdorder413(&h22[index]);
        PDstandardNth23h22 = PDstandardNthconstraintsfdorder423(&h22[index]);
        PDstandardNth1h31 = PDstandardNthconstraintsfdorder41(&h31[index]);
        PDstandardNth2h31 = PDstandardNthconstraintsfdorder42(&h31[index]);
        PDstandardNth3h31 = PDstandardNthconstraintsfdorder43(&h31[index]);
        PDstandardNth11h31 = PDstandardNthconstraintsfdorder411(&h31[index]);
        PDstandardNth22h31 = PDstandardNthconstraintsfdorder422(&h31[index]);
        PDstandardNth33h31 = PDstandardNthconstraintsfdorder433(&h31[index]);
        PDstandardNth12h31 = PDstandardNthconstraintsfdorder412(&h31[index]);
        PDstandardNth13h31 = PDstandardNthconstraintsfdorder413(&h31[index]);
        PDstandardNth23h31 = PDstandardNthconstraintsfdorder423(&h31[index]);
        PDstandardNth1h32 = PDstandardNthconstraintsfdorder41(&h32[index]);
        PDstandardNth2h32 = PDstandardNthconstraintsfdorder42(&h32[index]);
        PDstandardNth3h32 = PDstandardNthconstraintsfdorder43(&h32[index]);
        PDstandardNth11h32 = PDstandardNthconstraintsfdorder411(&h32[index]);
        PDstandardNth22h32 = PDstandardNthconstraintsfdorder422(&h32[index]);
        PDstandardNth33h32 = PDstandardNthconstraintsfdorder433(&h32[index]);
        PDstandardNth12h32 = PDstandardNthconstraintsfdorder412(&h32[index]);
        PDstandardNth13h32 = PDstandardNthconstraintsfdorder413(&h32[index]);
        PDstandardNth23h32 = PDstandardNthconstraintsfdorder423(&h32[index]);
        PDstandardNth1h33 = PDstandardNthconstraintsfdorder41(&h33[index]);
        PDstandardNth2h33 = PDstandardNthconstraintsfdorder42(&h33[index]);
        PDstandardNth3h33 = PDstandardNthconstraintsfdorder43(&h33[index]);
        PDstandardNth11h33 = PDstandardNthconstraintsfdorder411(&h33[index]);
        PDstandardNth22h33 = PDstandardNthconstraintsfdorder422(&h33[index]);
        PDstandardNth33h33 = PDstandardNthconstraintsfdorder433(&h33[index]);
        PDstandardNth12h33 = PDstandardNthconstraintsfdorder412(&h33[index]);
        PDstandardNth13h33 = PDstandardNthconstraintsfdorder413(&h33[index]);
        PDstandardNth23h33 = PDstandardNthconstraintsfdorder423(&h33[index]);
        PDstandardNth1K = PDstandardNthconstraintsfdorder41(&K[index]);
        PDstandardNth2K = PDstandardNthconstraintsfdorder42(&K[index]);
        PDstandardNth3K = PDstandardNthconstraintsfdorder43(&K[index]);
        break;
      }
      
      case 6:
      {
        PDstandardNth1A11 = PDstandardNthconstraintsfdorder61(&A11[index]);
        PDstandardNth2A11 = PDstandardNthconstraintsfdorder62(&A11[index]);
        PDstandardNth3A11 = PDstandardNthconstraintsfdorder63(&A11[index]);
        PDstandardNth1A21 = PDstandardNthconstraintsfdorder61(&A21[index]);
        PDstandardNth2A21 = PDstandardNthconstraintsfdorder62(&A21[index]);
        PDstandardNth3A21 = PDstandardNthconstraintsfdorder63(&A21[index]);
        PDstandardNth1A22 = PDstandardNthconstraintsfdorder61(&A22[index]);
        PDstandardNth2A22 = PDstandardNthconstraintsfdorder62(&A22[index]);
        PDstandardNth3A22 = PDstandardNthconstraintsfdorder63(&A22[index]);
        PDstandardNth1A31 = PDstandardNthconstraintsfdorder61(&A31[index]);
        PDstandardNth2A31 = PDstandardNthconstraintsfdorder62(&A31[index]);
        PDstandardNth3A31 = PDstandardNthconstraintsfdorder63(&A31[index]);
        PDstandardNth1A32 = PDstandardNthconstraintsfdorder61(&A32[index]);
        PDstandardNth2A32 = PDstandardNthconstraintsfdorder62(&A32[index]);
        PDstandardNth3A32 = PDstandardNthconstraintsfdorder63(&A32[index]);
        PDstandardNth1A33 = PDstandardNthconstraintsfdorder61(&A33[index]);
        PDstandardNth2A33 = PDstandardNthconstraintsfdorder62(&A33[index]);
        PDstandardNth3A33 = PDstandardNthconstraintsfdorder63(&A33[index]);
        PDstandardNth1chi = PDstandardNthconstraintsfdorder61(&chi[index]);
        PDstandardNth2chi = PDstandardNthconstraintsfdorder62(&chi[index]);
        PDstandardNth3chi = PDstandardNthconstraintsfdorder63(&chi[index]);
        PDstandardNth11chi = PDstandardNthconstraintsfdorder611(&chi[index]);
        PDstandardNth22chi = PDstandardNthconstraintsfdorder622(&chi[index]);
        PDstandardNth33chi = PDstandardNthconstraintsfdorder633(&chi[index]);
        PDstandardNth12chi = PDstandardNthconstraintsfdorder612(&chi[index]);
        PDstandardNth13chi = PDstandardNthconstraintsfdorder613(&chi[index]);
        PDstandardNth23chi = PDstandardNthconstraintsfdorder623(&chi[index]);
        PDstandardNth1Gam1 = PDstandardNthconstraintsfdorder61(&Gam1[index]);
        PDstandardNth2Gam1 = PDstandardNthconstraintsfdorder62(&Gam1[index]);
        PDstandardNth3Gam1 = PDstandardNthconstraintsfdorder63(&Gam1[index]);
        PDstandardNth1Gam2 = PDstandardNthconstraintsfdorder61(&Gam2[index]);
        PDstandardNth2Gam2 = PDstandardNthconstraintsfdorder62(&Gam2[index]);
        PDstandardNth3Gam2 = PDstandardNthconstraintsfdorder63(&Gam2[index]);
        PDstandardNth1Gam3 = PDstandardNthconstraintsfdorder61(&Gam3[index]);
        PDstandardNth2Gam3 = PDstandardNthconstraintsfdorder62(&Gam3[index]);
        PDstandardNth3Gam3 = PDstandardNthconstraintsfdorder63(&Gam3[index]);
        PDstandardNth1h11 = PDstandardNthconstraintsfdorder61(&h11[index]);
        PDstandardNth2h11 = PDstandardNthconstraintsfdorder62(&h11[index]);
        PDstandardNth3h11 = PDstandardNthconstraintsfdorder63(&h11[index]);
        PDstandardNth11h11 = PDstandardNthconstraintsfdorder611(&h11[index]);
        PDstandardNth22h11 = PDstandardNthconstraintsfdorder622(&h11[index]);
        PDstandardNth33h11 = PDstandardNthconstraintsfdorder633(&h11[index]);
        PDstandardNth12h11 = PDstandardNthconstraintsfdorder612(&h11[index]);
        PDstandardNth13h11 = PDstandardNthconstraintsfdorder613(&h11[index]);
        PDstandardNth23h11 = PDstandardNthconstraintsfdorder623(&h11[index]);
        PDstandardNth1h21 = PDstandardNthconstraintsfdorder61(&h21[index]);
        PDstandardNth2h21 = PDstandardNthconstraintsfdorder62(&h21[index]);
        PDstandardNth3h21 = PDstandardNthconstraintsfdorder63(&h21[index]);
        PDstandardNth11h21 = PDstandardNthconstraintsfdorder611(&h21[index]);
        PDstandardNth22h21 = PDstandardNthconstraintsfdorder622(&h21[index]);
        PDstandardNth33h21 = PDstandardNthconstraintsfdorder633(&h21[index]);
        PDstandardNth12h21 = PDstandardNthconstraintsfdorder612(&h21[index]);
        PDstandardNth13h21 = PDstandardNthconstraintsfdorder613(&h21[index]);
        PDstandardNth23h21 = PDstandardNthconstraintsfdorder623(&h21[index]);
        PDstandardNth1h22 = PDstandardNthconstraintsfdorder61(&h22[index]);
        PDstandardNth2h22 = PDstandardNthconstraintsfdorder62(&h22[index]);
        PDstandardNth3h22 = PDstandardNthconstraintsfdorder63(&h22[index]);
        PDstandardNth11h22 = PDstandardNthconstraintsfdorder611(&h22[index]);
        PDstandardNth22h22 = PDstandardNthconstraintsfdorder622(&h22[index]);
        PDstandardNth33h22 = PDstandardNthconstraintsfdorder633(&h22[index]);
        PDstandardNth12h22 = PDstandardNthconstraintsfdorder612(&h22[index]);
        PDstandardNth13h22 = PDstandardNthconstraintsfdorder613(&h22[index]);
        PDstandardNth23h22 = PDstandardNthconstraintsfdorder623(&h22[index]);
        PDstandardNth1h31 = PDstandardNthconstraintsfdorder61(&h31[index]);
        PDstandardNth2h31 = PDstandardNthconstraintsfdorder62(&h31[index]);
        PDstandardNth3h31 = PDstandardNthconstraintsfdorder63(&h31[index]);
        PDstandardNth11h31 = PDstandardNthconstraintsfdorder611(&h31[index]);
        PDstandardNth22h31 = PDstandardNthconstraintsfdorder622(&h31[index]);
        PDstandardNth33h31 = PDstandardNthconstraintsfdorder633(&h31[index]);
        PDstandardNth12h31 = PDstandardNthconstraintsfdorder612(&h31[index]);
        PDstandardNth13h31 = PDstandardNthconstraintsfdorder613(&h31[index]);
        PDstandardNth23h31 = PDstandardNthconstraintsfdorder623(&h31[index]);
        PDstandardNth1h32 = PDstandardNthconstraintsfdorder61(&h32[index]);
        PDstandardNth2h32 = PDstandardNthconstraintsfdorder62(&h32[index]);
        PDstandardNth3h32 = PDstandardNthconstraintsfdorder63(&h32[index]);
        PDstandardNth11h32 = PDstandardNthconstraintsfdorder611(&h32[index]);
        PDstandardNth22h32 = PDstandardNthconstraintsfdorder622(&h32[index]);
        PDstandardNth33h32 = PDstandardNthconstraintsfdorder633(&h32[index]);
        PDstandardNth12h32 = PDstandardNthconstraintsfdorder612(&h32[index]);
        PDstandardNth13h32 = PDstandardNthconstraintsfdorder613(&h32[index]);
        PDstandardNth23h32 = PDstandardNthconstraintsfdorder623(&h32[index]);
        PDstandardNth1h33 = PDstandardNthconstraintsfdorder61(&h33[index]);
        PDstandardNth2h33 = PDstandardNthconstraintsfdorder62(&h33[index]);
        PDstandardNth3h33 = PDstandardNthconstraintsfdorder63(&h33[index]);
        PDstandardNth11h33 = PDstandardNthconstraintsfdorder611(&h33[index]);
        PDstandardNth22h33 = PDstandardNthconstraintsfdorder622(&h33[index]);
        PDstandardNth33h33 = PDstandardNthconstraintsfdorder633(&h33[index]);
        PDstandardNth12h33 = PDstandardNthconstraintsfdorder612(&h33[index]);
        PDstandardNth13h33 = PDstandardNthconstraintsfdorder613(&h33[index]);
        PDstandardNth23h33 = PDstandardNthconstraintsfdorder623(&h33[index]);
        PDstandardNth1K = PDstandardNthconstraintsfdorder61(&K[index]);
        PDstandardNth2K = PDstandardNthconstraintsfdorder62(&K[index]);
        PDstandardNth3K = PDstandardNthconstraintsfdorder63(&K[index]);
        break;
      }
      default:
        CCTK_BUILTIN_UNREACHABLE();
    }
    /* Calculate temporaries and grid functions */
    CCTK_REAL Sij11 CCTK_ATTRIBUTE_UNUSED = eTxxL;
    
    CCTK_REAL Sij21 CCTK_ATTRIBUTE_UNUSED = eTxyL;
    
    CCTK_REAL Sij31 CCTK_ATTRIBUTE_UNUSED = eTxzL;
    
    CCTK_REAL Sij22 CCTK_ATTRIBUTE_UNUSED = eTyyL;
    
    CCTK_REAL Sij32 CCTK_ATTRIBUTE_UNUSED = eTyzL;
    
    CCTK_REAL Sij33 CCTK_ATTRIBUTE_UNUSED = eTzzL;
    
    CCTK_REAL Tj1 CCTK_ATTRIBUTE_UNUSED = eTtxL;
    
    CCTK_REAL Tj2 CCTK_ATTRIBUTE_UNUSED = eTtyL;
    
    CCTK_REAL Tj3 CCTK_ATTRIBUTE_UNUSED = eTtzL;
    
    CCTK_REAL Si1 CCTK_ATTRIBUTE_UNUSED = (beta1L*Sij11 + beta2L*Sij21 + 
      beta3L*Sij31 - Tj1)*pow(alphaL,-1);
    
    CCTK_REAL Si2 CCTK_ATTRIBUTE_UNUSED = (beta1L*Sij21 + beta2L*Sij22 + 
      beta3L*Sij32 - Tj2)*pow(alphaL,-1);
    
    CCTK_REAL Si3 CCTK_ATTRIBUTE_UNUSED = (beta1L*Sij31 + beta2L*Sij32 + 
      beta3L*Sij33 - Tj3)*pow(alphaL,-1);
    
    CCTK_REAL rho CCTK_ATTRIBUTE_UNUSED = pow(alphaL,-2)*(eTttL + 
      2*(beta2L*beta3L*Sij32 + beta1L*(beta2L*Sij21 + beta3L*Sij31 - Tj1)) - 
      2*(beta2L*Tj2 + beta3L*Tj3) + Sij11*pow(beta1L,2) + Sij22*pow(beta2L,2) 
      + Sij33*pow(beta3L,2));
    
    CCTK_REAL pi CCTK_ATTRIBUTE_UNUSED = 3.1415926535897932385;
    
    CCTK_REAL deth CCTK_ATTRIBUTE_UNUSED = 2*h21L*h31L*h32L - 
      h33L*pow(h21L,2) + h22L*(h11L*h33L - pow(h31L,2)) - h11L*pow(h32L,2);
    
    CCTK_REAL invdeth CCTK_ATTRIBUTE_UNUSED = pow(deth,-1);
    
    ptrdiff_t dir1 CCTK_ATTRIBUTE_UNUSED = isgn(beta1L);
    
    ptrdiff_t dir2 CCTK_ATTRIBUTE_UNUSED = isgn(beta2L);
    
    ptrdiff_t dir3 CCTK_ATTRIBUTE_UNUSED = isgn(beta3L);
    
    CCTK_REAL hInv11 CCTK_ATTRIBUTE_UNUSED = invdeth*(h22L*h33L - 
      pow(h32L,2));
    
    CCTK_REAL hInv12 CCTK_ATTRIBUTE_UNUSED = (h31L*h32L - 
      h21L*h33L)*invdeth;
    
    CCTK_REAL hInv13 CCTK_ATTRIBUTE_UNUSED = (-(h22L*h31L) + 
      h21L*h32L)*invdeth;
    
    CCTK_REAL hInv22 CCTK_ATTRIBUTE_UNUSED = invdeth*(h11L*h33L - 
      pow(h31L,2));
    
    CCTK_REAL hInv23 CCTK_ATTRIBUTE_UNUSED = (h21L*h31L - 
      h11L*h32L)*invdeth;
    
    CCTK_REAL hInv33 CCTK_ATTRIBUTE_UNUSED = invdeth*(h11L*h22L - 
      pow(h21L,2));
    
    CCTK_REAL gamma111 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv11*PDstandardNth1h11 + 2*(hInv12*PDstandardNth1h21 + 
      hInv13*PDstandardNth1h31) - hInv12*PDstandardNth2h11 - 
      hInv13*PDstandardNth3h11);
    
    CCTK_REAL gamma211 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv12*PDstandardNth1h11 + 2*(hInv22*PDstandardNth1h21 + 
      hInv23*PDstandardNth1h31) - hInv22*PDstandardNth2h11 - 
      hInv23*PDstandardNth3h11);
    
    CCTK_REAL gamma311 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv13*PDstandardNth1h11 + 2*(hInv23*PDstandardNth1h21 + 
      hInv33*PDstandardNth1h31) - hInv23*PDstandardNth2h11 - 
      hInv33*PDstandardNth3h11);
    
    CCTK_REAL gamma121 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv12*PDstandardNth1h22 + hInv11*PDstandardNth2h11 + 
      hInv13*(PDstandardNth1h32 + PDstandardNth2h31 - PDstandardNth3h21));
    
    CCTK_REAL gamma221 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv22*PDstandardNth1h22 + hInv12*PDstandardNth2h11 + 
      hInv23*(PDstandardNth1h32 + PDstandardNth2h31 - PDstandardNth3h21));
    
    CCTK_REAL gamma321 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv23*PDstandardNth1h22 + hInv13*PDstandardNth2h11 + 
      hInv33*(PDstandardNth1h32 + PDstandardNth2h31 - PDstandardNth3h21));
    
    CCTK_REAL gamma131 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv13*PDstandardNth1h33 + hInv11*PDstandardNth3h11 + 
      hInv12*(PDstandardNth1h32 - PDstandardNth2h31 + PDstandardNth3h21));
    
    CCTK_REAL gamma231 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv23*PDstandardNth1h33 + hInv12*PDstandardNth3h11 + 
      hInv22*(PDstandardNth1h32 - PDstandardNth2h31 + PDstandardNth3h21));
    
    CCTK_REAL gamma331 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv33*PDstandardNth1h33 + hInv13*PDstandardNth3h11 + 
      hInv23*(PDstandardNth1h32 - PDstandardNth2h31 + PDstandardNth3h21));
    
    CCTK_REAL gamma122 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv11*(-PDstandardNth1h22 + 2*PDstandardNth2h21) + 
      hInv12*PDstandardNth2h22 + hInv13*(2*PDstandardNth2h32 - 
      PDstandardNth3h22));
    
    CCTK_REAL gamma222 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv12*(-PDstandardNth1h22 + 2*PDstandardNth2h21) + 
      hInv22*PDstandardNth2h22 + hInv23*(2*PDstandardNth2h32 - 
      PDstandardNth3h22));
    
    CCTK_REAL gamma322 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv13*(-PDstandardNth1h22 + 2*PDstandardNth2h21) + 
      hInv23*PDstandardNth2h22 + hInv33*(2*PDstandardNth2h32 - 
      PDstandardNth3h22));
    
    CCTK_REAL gamma132 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv13*PDstandardNth2h33 + hInv11*(-PDstandardNth1h32 + 
      PDstandardNth2h31 + PDstandardNth3h21) + hInv12*PDstandardNth3h22);
    
    CCTK_REAL gamma232 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv23*PDstandardNth2h33 + hInv12*(-PDstandardNth1h32 + 
      PDstandardNth2h31 + PDstandardNth3h21) + hInv22*PDstandardNth3h22);
    
    CCTK_REAL gamma332 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv33*PDstandardNth2h33 + hInv13*(-PDstandardNth1h32 + 
      PDstandardNth2h31 + PDstandardNth3h21) + hInv23*PDstandardNth3h22);
    
    CCTK_REAL gamma133 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv11*(-PDstandardNth1h33 + 2*PDstandardNth3h31) + 
      hInv12*(-PDstandardNth2h33 + 2*PDstandardNth3h32) + 
      hInv13*PDstandardNth3h33);
    
    CCTK_REAL gamma233 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv12*(-PDstandardNth1h33 + 2*PDstandardNth3h31) + 
      hInv22*(-PDstandardNth2h33 + 2*PDstandardNth3h32) + 
      hInv23*PDstandardNth3h33);
    
    CCTK_REAL gamma333 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv13*(-PDstandardNth1h33 + 2*PDstandardNth3h31) + 
      hInv23*(-PDstandardNth2h33 + 2*PDstandardNth3h32) + 
      hInv33*PDstandardNth3h33);
    
    CCTK_REAL AInv11 CCTK_ATTRIBUTE_UNUSED = 2*(A32L*hInv12*hInv13 + 
      hInv11*(A21L*hInv12 + A31L*hInv13)) + A11L*pow(hInv11,2) + 
      A22L*pow(hInv12,2) + A33L*pow(hInv13,2);
    
    CCTK_REAL AInv12 CCTK_ATTRIBUTE_UNUSED = hInv12*(A11L*hInv11 + 
      A21L*hInv12 + A31L*hInv13) + (A21L*hInv11 + A22L*hInv12 + 
      A32L*hInv13)*hInv22 + (A31L*hInv11 + A32L*hInv12 + A33L*hInv13)*hInv23;
    
    CCTK_REAL AInv13 CCTK_ATTRIBUTE_UNUSED = hInv13*(A11L*hInv11 + 
      A21L*hInv12 + A31L*hInv13) + (A21L*hInv11 + A22L*hInv12 + 
      A32L*hInv13)*hInv23 + (A31L*hInv11 + A32L*hInv12 + A33L*hInv13)*hInv33;
    
    CCTK_REAL AInv22 CCTK_ATTRIBUTE_UNUSED = 2*(A32L*hInv22*hInv23 + 
      hInv12*(A21L*hInv22 + A31L*hInv23)) + A11L*pow(hInv12,2) + 
      A22L*pow(hInv22,2) + A33L*pow(hInv23,2);
    
    CCTK_REAL AInv23 CCTK_ATTRIBUTE_UNUSED = hInv13*(A11L*hInv12 + 
      A21L*hInv22 + A31L*hInv23) + hInv23*(A21L*hInv12 + A22L*hInv22 + 
      A32L*hInv23) + (A31L*hInv12 + A32L*hInv22 + A33L*hInv23)*hInv33;
    
    CCTK_REAL AInv33 CCTK_ATTRIBUTE_UNUSED = 2*(A32L*hInv23*hInv33 + 
      hInv13*(A21L*hInv23 + A31L*hInv33)) + A11L*pow(hInv13,2) + 
      A22L*pow(hInv23,2) + A33L*pow(hInv33,2);
    
    CCTK_REAL em4phi CCTK_ATTRIBUTE_UNUSED = chiL;
    
    CCTK_REAL chiCutoff CCTK_ATTRIBUTE_UNUSED = fmax(chiL,chiEps);
    
    CCTK_REAL DDphi11 CCTK_ATTRIBUTE_UNUSED = 
      0.25*pow(chiCutoff,-2)*(chiCutoff*(-PDstandardNth11chi + 
      gamma111*PDstandardNth1chi + gamma211*PDstandardNth2chi + 
      gamma311*PDstandardNth3chi) + pow(PDstandardNth1chi,2));
    
    CCTK_REAL DDphi21 CCTK_ATTRIBUTE_UNUSED = 
      0.25*(PDstandardNth1chi*(chiCutoff*gamma121 + PDstandardNth2chi) + 
      chiCutoff*(-PDstandardNth12chi + gamma221*PDstandardNth2chi + 
      gamma321*PDstandardNth3chi))*pow(chiCutoff,-2);
    
    CCTK_REAL DDphi31 CCTK_ATTRIBUTE_UNUSED = 
      0.25*(PDstandardNth1chi*(chiCutoff*gamma131 + PDstandardNth3chi) + 
      chiCutoff*(-PDstandardNth13chi + gamma231*PDstandardNth2chi + 
      gamma331*PDstandardNth3chi))*pow(chiCutoff,-2);
    
    CCTK_REAL DDphi12 CCTK_ATTRIBUTE_UNUSED = 
      0.25*(PDstandardNth1chi*(chiCutoff*gamma121 + PDstandardNth2chi) + 
      chiCutoff*(-PDstandardNth12chi + gamma221*PDstandardNth2chi + 
      gamma321*PDstandardNth3chi))*pow(chiCutoff,-2);
    
    CCTK_REAL DDphi22 CCTK_ATTRIBUTE_UNUSED = 
      0.25*pow(chiCutoff,-2)*(chiCutoff*(gamma122*PDstandardNth1chi - 
      PDstandardNth22chi + gamma222*PDstandardNth2chi + 
      gamma322*PDstandardNth3chi) + pow(PDstandardNth2chi,2));
    
    CCTK_REAL DDphi32 CCTK_ATTRIBUTE_UNUSED = 
      0.25*(PDstandardNth2chi*(chiCutoff*gamma232 + PDstandardNth3chi) + 
      chiCutoff*(gamma132*PDstandardNth1chi - PDstandardNth23chi + 
      gamma332*PDstandardNth3chi))*pow(chiCutoff,-2);
    
    CCTK_REAL DDphi13 CCTK_ATTRIBUTE_UNUSED = 
      0.25*(PDstandardNth1chi*(chiCutoff*gamma131 + PDstandardNth3chi) + 
      chiCutoff*(-PDstandardNth13chi + gamma231*PDstandardNth2chi + 
      gamma331*PDstandardNth3chi))*pow(chiCutoff,-2);
    
    CCTK_REAL DDphi23 CCTK_ATTRIBUTE_UNUSED = 
      0.25*(PDstandardNth2chi*(chiCutoff*gamma232 + PDstandardNth3chi) + 
      chiCutoff*(gamma132*PDstandardNth1chi - PDstandardNth23chi + 
      gamma332*PDstandardNth3chi))*pow(chiCutoff,-2);
    
    CCTK_REAL DDphi33 CCTK_ATTRIBUTE_UNUSED = 
      0.25*pow(chiCutoff,-2)*(chiCutoff*(gamma133*PDstandardNth1chi + 
      gamma233*PDstandardNth2chi - PDstandardNth33chi + 
      gamma333*PDstandardNth3chi) + pow(PDstandardNth3chi,2));
    
    CCTK_REAL Dphi1 CCTK_ATTRIBUTE_UNUSED = 
      -0.25*PDstandardNth1chi*pow(chiCutoff,-1);
    
    CCTK_REAL Dphi2 CCTK_ATTRIBUTE_UNUSED = 
      -0.25*PDstandardNth2chi*pow(chiCutoff,-1);
    
    CCTK_REAL Dphi3 CCTK_ATTRIBUTE_UNUSED = 
      -0.25*PDstandardNth3chi*pow(chiCutoff,-1);
    
    CCTK_REAL R11 CCTK_ATTRIBUTE_UNUSED = gamma221*(4*h21L*gamma111*hInv12 
      + 2*h11L*gamma121*hInv12 + 2*h21L*gamma231*hInv13 + 
      2*h11L*gamma122*hInv22 + 4*h21L*gamma131*hInv23) + 
      gamma321*(4*h31L*gamma111*hInv12 + 2*h11L*gamma131*hInv12 + 
      2*h31L*gamma231*hInv13 + 4*h31L*gamma121*hInv22 + 
      4*h31L*gamma131*hInv23) + gamma231*(4*h21L*gamma111*hInv13 + 
      2*h11L*gamma121*hInv13 + 4*h21L*gamma121*hInv23 + 
      4*h21L*gamma131*hInv33) + gamma331*(4*h31L*gamma121*hInv23 + 
      2*h11L*gamma132*hInv23 + 4*h31L*gamma131*hInv33) - 
      hInv12*PDstandardNth12h11 + h21L*(Gam1L*gamma211 + Gam2L*gamma221 + 
      Gam3L*gamma231 + 4*gamma111*gamma211*hInv11 + 
      4*gamma121*gamma211*hInv12 + 4*gamma131*gamma211*hInv13 + 
      4*gamma121*gamma221*hInv22 + 2*gamma222*gamma231*hInv23 + 
      2*gamma231*gamma232*hInv33 + PDstandardNth1Gam2) + h31L*(Gam1L*gamma311 
      + Gam2L*gamma321 + Gam3L*gamma331 + 4*gamma111*gamma311*hInv11 + 
      4*gamma121*gamma311*hInv12 + 4*gamma131*gamma311*hInv13 + 
      2*gamma221*gamma322*hInv22 + 2*gamma231*gamma322*hInv23 + 
      2*gamma231*gamma332*hInv33 + PDstandardNth1Gam3) - 
      hInv23*PDstandardNth23h11 + h11L*(Gam1L*gamma111 + Gam2L*gamma121 + 
      Gam3L*gamma131 + 6*gamma111*gamma121*hInv12 + 
      6*gamma111*gamma131*hInv13 + 6*gamma121*gamma131*hInv23 + 
      PDstandardNth1Gam1 + 3*hInv11*pow(gamma111,2) + 
      3*hInv22*pow(gamma121,2) + 3*hInv33*pow(gamma131,2)) + 
      2*(h32L*gamma211*gamma311*hInv11 + h21L*gamma231*gamma311*hInv11 + 
      h11L*gamma122*gamma211*hInv12 + h22L*gamma211*gamma221*hInv12 + 
      h21L*gamma211*gamma222*hInv12 + h11L*gamma132*gamma311*hInv12 + 
      h32L*gamma221*gamma311*hInv12 + h21L*gamma232*gamma311*hInv12 + 
      h32L*gamma211*gamma321*hInv12 + h31L*gamma221*gamma321*hInv12 + 
      h21L*gamma231*gamma321*hInv12 + h33L*gamma311*gamma321*hInv12 + 
      h31L*gamma211*gamma322*hInv12 + h31L*gamma321*gamma331*hInv12 + 
      h31L*gamma311*gamma332*hInv12 + h11L*gamma132*gamma211*hInv13 + 
      h22L*gamma211*gamma231*hInv13 + h21L*gamma211*gamma232*hInv13 + 
      h11L*gamma133*gamma311*hInv13 + h32L*gamma231*gamma311*hInv13 + 
      h11L*gamma131*gamma331*hInv13 + h32L*gamma211*gamma331*hInv13 + 
      h21L*gamma231*gamma331*hInv13 + h33L*gamma311*gamma331*hInv13 + 
      h21L*gamma221*gamma222*hInv22 + h32L*gamma221*gamma321*hInv22 + 
      h21L*gamma232*gamma321*hInv22 + h11L*gamma132*gamma221*hInv23 + 
      h11L*gamma122*gamma231*hInv23 + h22L*gamma221*gamma231*hInv23 + 
      h21L*gamma221*gamma232*hInv23 + h11L*gamma133*gamma321*hInv23 + 
      h32L*gamma231*gamma321*hInv23 + h21L*gamma233*gamma321*hInv23 + 
      h32L*gamma221*gamma331*hInv23 + h21L*gamma232*gamma331*hInv23 + 
      h33L*gamma321*gamma331*hInv23 + h31L*gamma221*gamma332*hInv23 + 
      h31L*gamma331*gamma332*hInv23 + h31L*gamma321*gamma333*hInv23 + 
      h11L*gamma132*gamma231*hInv33 + h32L*gamma231*gamma331*hInv33 + 
      h21L*gamma233*gamma331*hInv33 + h21L*hInv12*pow(gamma221,2)) + 
      hInv11*(2*h11L*gamma121*gamma211 + 2*h21L*gamma211*gamma221 + 
      2*h11L*gamma131*gamma311 + 2*h31L*gamma211*gamma321 + 
      2*h31L*gamma311*gamma331 - 0.5*PDstandardNth11h11 + 
      h22L*pow(gamma211,2) + h33L*pow(gamma311,2)) + 
      hInv22*(2*h11L*gamma132*gamma321 + 2*h31L*gamma321*gamma332 - 
      0.5*PDstandardNth22h11 + h22L*pow(gamma221,2) + h33L*pow(gamma321,2)) + 
      hInv13*(2*h21L*gamma233*gamma311 + 4*h31L*gamma111*gamma331 + 
      2*h31L*gamma211*gamma332 + 2*h31L*gamma311*gamma333 - 
      PDstandardNth13h11 + 2*h31L*pow(gamma331,2)) + 
      hInv33*(2*h11L*gamma133*gamma331 + 2*h31L*gamma331*gamma333 - 
      0.5*PDstandardNth33h11 + h22L*pow(gamma231,2) + h33L*pow(gamma331,2));
    
    CCTK_REAL R21 CCTK_ATTRIBUTE_UNUSED = 0.5*(Gam1L*(h21L*gamma111 + 
      h11L*gamma121 + h22L*gamma211 + h21L*gamma221 + h32L*gamma311 + 
      h31L*gamma321) + Gam2L*(h21L*gamma121 + h11L*gamma122 + h22L*gamma221 + 
      h21L*gamma222 + h32L*gamma321 + h31L*gamma322) + Gam3L*(h21L*gamma131 + 
      h11L*gamma132 + h22L*gamma231 + h21L*gamma232 + h32L*gamma331 + 
      h31L*gamma332) - hInv11*PDstandardNth11h21 - 
      2*hInv12*PDstandardNth12h21 - 2*hInv13*PDstandardNth13h21 + 
      h21L*PDstandardNth1Gam1 + h22L*PDstandardNth1Gam2 + 
      h32L*PDstandardNth1Gam3 - hInv22*PDstandardNth22h21 - 
      2*hInv23*PDstandardNth23h21 + h11L*PDstandardNth2Gam1 + 
      h21L*PDstandardNth2Gam2 + h31L*PDstandardNth2Gam3 - 
      hInv33*PDstandardNth33h21 + 2*((h21L*gamma111*gamma131 + 
      h21L*gamma132*gamma211 + h21L*gamma131*gamma221 + 
      h11L*gamma132*gamma221 + h22L*gamma111*gamma231 + 
      h22L*gamma221*gamma231 + h22L*gamma211*gamma232 + 
      h21L*gamma221*gamma232 + h21L*gamma133*gamma311 + 
      h22L*gamma233*gamma311 + h31L*gamma131*gamma321 + 
      h11L*gamma133*gamma321 + h32L*gamma231*gamma321 + 
      h21L*gamma233*gamma321 + h32L*gamma111*gamma331 + 
      h32L*gamma221*gamma331 + h33L*gamma321*gamma331 + 
      2*gamma121*(h11L*gamma131 + h21L*gamma231 + h31L*gamma331) + 
      h32L*gamma211*gamma332 + h31L*gamma221*gamma332 + 
      h32L*gamma311*gamma333 + h31L*gamma321*gamma333)*hInv13 + 
      (h21L*gamma121*gamma131 + h21L*gamma132*gamma221 + 
      h21L*gamma131*gamma222 + h11L*gamma132*gamma222 + 
      h22L*gamma121*gamma231 + h22L*gamma222*gamma231 + 
      h22L*gamma221*gamma232 + h21L*gamma222*gamma232 + 
      h21L*gamma133*gamma321 + h22L*gamma233*gamma321 + 
      h31L*gamma131*gamma322 + h11L*gamma133*gamma322 + 
      h32L*gamma231*gamma322 + h21L*gamma233*gamma322 + 
      h32L*gamma121*gamma331 + h32L*gamma222*gamma331 + 
      h33L*gamma322*gamma331 + 2*gamma122*(h11L*gamma131 + h21L*gamma231 + 
      h31L*gamma331) + h32L*gamma221*gamma332 + h31L*gamma222*gamma332 + 
      h32L*gamma321*gamma333 + h31L*gamma322*gamma333)*hInv23 + 
      hInv11*(2*h22L*gamma211*gamma221 + h21L*gamma131*gamma311 + 
      h32L*gamma221*gamma311 + h22L*gamma231*gamma311 + 
      gamma121*(3*h21L*gamma211 + h11L*gamma221 + 2*h31L*gamma311) + 
      h11L*gamma131*gamma321 + 2*h32L*gamma211*gamma321 + 
      h31L*gamma221*gamma321 + h21L*gamma231*gamma321 + 
      h33L*gamma311*gamma321 + gamma111*(2*h11L*gamma121 + h22L*gamma211 + 
      h21L*gamma221 + h32L*gamma311 + h31L*gamma321) + h32L*gamma311*gamma331 
      + h31L*gamma321*gamma331 + h21L*pow(gamma111,2) + h21L*pow(gamma221,2)) 
      + hInv12*(2*h21L*gamma122*gamma211 + h22L*gamma211*gamma222 + 
      h21L*gamma221*gamma222 + 2*h31L*gamma122*gamma311 + 
      h32L*gamma222*gamma311 + gamma121*(h22L*gamma211 + h21L*gamma221 + 
      h11L*gamma222 + h32L*gamma311) + h21L*gamma131*gamma321 + 
      h32L*gamma221*gamma321 + h31L*gamma222*gamma321 + 
      h22L*gamma231*gamma321 + h11L*gamma131*gamma322 + 
      h32L*gamma211*gamma322 + h21L*gamma231*gamma322 + 
      h33L*gamma311*gamma322 + gamma111*(h21L*gamma121 + 2*h11L*gamma122 + 
      h21L*gamma222 + h31L*gamma322) + h32L*gamma321*gamma331 + 
      h31L*gamma322*gamma331 + h22L*pow(gamma221,2)) + 
      hInv22*(2*h22L*gamma221*gamma222 + h21L*gamma132*gamma321 + 
      h32L*gamma222*gamma321 + h22L*gamma232*gamma321 + 
      gamma122*(3*h21L*gamma221 + h11L*gamma222 + 2*h31L*gamma321) + 
      h11L*gamma132*gamma322 + 2*h32L*gamma221*gamma322 + 
      h31L*gamma222*gamma322 + h21L*gamma232*gamma322 + 
      h33L*gamma321*gamma322 + gamma121*(2*h11L*gamma122 + h22L*gamma221 + 
      h21L*gamma222 + h32L*gamma321 + h31L*gamma322) + h32L*gamma321*gamma332 
      + h31L*gamma322*gamma332 + h21L*pow(gamma121,2) + h21L*pow(gamma222,2)) 
      + hInv33*(2*h22L*gamma231*gamma232 + h21L*gamma133*gamma331 + 
      h32L*gamma232*gamma331 + h22L*gamma233*gamma331 + 
      gamma132*(3*h21L*gamma231 + h11L*gamma232 + 2*h31L*gamma331) + 
      h11L*gamma133*gamma332 + 2*h32L*gamma231*gamma332 + 
      h31L*gamma232*gamma332 + h21L*gamma233*gamma332 + 
      h33L*gamma331*gamma332 + gamma131*(2*h11L*gamma132 + h22L*gamma231 + 
      h21L*gamma232 + h32L*gamma331 + h31L*gamma332) + h32L*gamma331*gamma333 
      + h31L*gamma332*gamma333 + h21L*pow(gamma131,2) + h21L*pow(gamma232,2)) 
      + hInv12*(h21L*gamma122*gamma211 + h22L*gamma111*gamma221 + 
      h11L*gamma122*gamma221 + h22L*gamma211*gamma222 + 
      h21L*gamma221*gamma222 + h21L*gamma132*gamma311 + 
      h22L*gamma232*gamma311 + h32L*gamma111*gamma321 + 
      h11L*gamma132*gamma321 + 2*h32L*gamma221*gamma321 + 
      h21L*gamma232*gamma321 + gamma121*(h21L*gamma111 + 3*h21L*gamma221 + 
      3*h31L*gamma321) + h32L*gamma211*gamma322 + h31L*gamma221*gamma322 + 
      h32L*gamma311*gamma332 + h31L*gamma321*gamma332 + 
      2*h11L*pow(gamma121,2) + h22L*pow(gamma221,2) + h33L*pow(gamma321,2)) + 
      hInv13*(h22L*gamma131*gamma211 + 2*h21L*gamma132*gamma211 + 
      h21L*gamma121*gamma231 + h22L*gamma221*gamma231 + 
      h11L*gamma121*gamma232 + h22L*gamma211*gamma232 + 
      h21L*gamma221*gamma232 + h32L*gamma131*gamma311 + 
      2*h31L*gamma132*gamma311 + h32L*gamma232*gamma311 + 
      h32L*gamma231*gamma321 + h31L*gamma232*gamma321 + 
      h21L*gamma131*gamma331 + h22L*gamma231*gamma331 + 
      h11L*gamma131*gamma332 + h32L*gamma211*gamma332 + 
      h21L*gamma231*gamma332 + h33L*gamma311*gamma332 + 
      h31L*gamma331*gamma332 + gamma111*(h21L*gamma131 + 2*h11L*gamma132 + 
      h21L*gamma232 + h31L*gamma332) + h32L*pow(gamma331,2)) + 
      hInv23*(h22L*gamma131*gamma221 + 2*h21L*gamma132*gamma221 + 
      h21L*gamma122*gamma231 + h22L*gamma222*gamma231 + 
      h11L*gamma122*gamma232 + h22L*gamma221*gamma232 + 
      h21L*gamma222*gamma232 + h32L*gamma131*gamma321 + 
      2*h31L*gamma132*gamma321 + h32L*gamma232*gamma321 + 
      h32L*gamma231*gamma322 + h31L*gamma232*gamma322 + 
      h21L*gamma132*gamma331 + h22L*gamma232*gamma331 + 
      h11L*gamma132*gamma332 + h32L*gamma221*gamma332 + 
      h21L*gamma232*gamma332 + h33L*gamma321*gamma332 + 
      h32L*gamma331*gamma332 + gamma121*(h21L*gamma131 + 2*h11L*gamma132 + 
      h21L*gamma232 + h31L*gamma332) + h31L*pow(gamma332,2))));
    
    CCTK_REAL R31 CCTK_ATTRIBUTE_UNUSED = 0.5*(Gam1L*(h31L*gamma111 + 
      h11L*gamma131 + h32L*gamma211 + h21L*gamma231 + h33L*gamma311 + 
      h31L*gamma331) + Gam2L*(h31L*gamma121 + h11L*gamma132 + h32L*gamma221 + 
      h21L*gamma232 + h33L*gamma321 + h31L*gamma332) + Gam3L*(h31L*gamma131 + 
      h11L*gamma133 + h32L*gamma231 + h21L*gamma233 + h33L*gamma331 + 
      h31L*gamma333) - hInv11*PDstandardNth11h31 - 
      2*hInv12*PDstandardNth12h31 - 2*hInv13*PDstandardNth13h31 + 
      h31L*PDstandardNth1Gam1 + h32L*PDstandardNth1Gam2 + 
      h33L*PDstandardNth1Gam3 - hInv22*PDstandardNth22h31 - 
      2*hInv23*PDstandardNth23h31 - hInv33*PDstandardNth33h31 + 
      h11L*PDstandardNth3Gam1 + h21L*PDstandardNth3Gam2 + 
      h31L*PDstandardNth3Gam3 + 2*((h31L*gamma122*gamma211 + 
      h32L*gamma111*gamma221 + 2*h21L*gamma131*gamma221 + 
      h32L*gamma211*gamma222 + h11L*gamma122*gamma231 + 
      h22L*gamma221*gamma231 + h21L*gamma222*gamma231 + 
      h31L*gamma132*gamma311 + h32L*gamma232*gamma311 + 
      h33L*gamma111*gamma321 + 2*h31L*gamma131*gamma321 + 
      h32L*gamma231*gamma321 + h33L*gamma211*gamma322 + 
      h31L*gamma231*gamma322 + h11L*gamma132*gamma331 + 
      h32L*gamma221*gamma331 + h21L*gamma232*gamma331 + 
      h33L*gamma321*gamma331 + gamma121*(2*h11L*gamma131 + h21L*gamma231 + 
      h31L*(gamma111 + gamma331)) + h33L*gamma311*gamma332 + 
      h31L*gamma331*gamma332)*hInv12 + (h32L*gamma131*gamma221 + 
      2*h21L*gamma133*gamma221 + h31L*gamma122*gamma231 + 
      h32L*gamma222*gamma231 + h11L*gamma122*gamma233 + 
      h22L*gamma221*gamma233 + h21L*gamma222*gamma233 + 
      h33L*gamma131*gamma321 + 2*h31L*gamma133*gamma321 + 
      h32L*gamma233*gamma321 + h33L*gamma231*gamma322 + 
      h31L*gamma233*gamma322 + h31L*gamma132*gamma331 + 
      h32L*gamma232*gamma331 + h33L*gamma331*gamma332 + 
      h11L*gamma132*gamma333 + h32L*gamma221*gamma333 + 
      h21L*gamma232*gamma333 + h33L*gamma321*gamma333 + 
      h31L*gamma332*gamma333 + gamma121*(2*h11L*gamma133 + h21L*gamma233 + 
      h31L*(gamma131 + gamma333)))*hInv23 + hInv12*(2*h21L*gamma132*gamma211 
      + h22L*gamma211*gamma232 + h21L*gamma221*gamma232 + 
      2*h31L*gamma132*gamma311 + h32L*gamma232*gamma311 + 
      gamma121*(h32L*gamma211 + h31L*gamma221 + h11L*gamma232 + 
      h33L*gamma311) + h31L*gamma131*gamma321 + h33L*gamma221*gamma321 + 
      h32L*gamma231*gamma321 + h31L*gamma232*gamma321 + 
      h33L*gamma321*gamma331 + h11L*gamma131*gamma332 + 
      h32L*gamma211*gamma332 + h21L*gamma231*gamma332 + 
      h33L*gamma311*gamma332 + h31L*gamma331*gamma332 + 
      gamma111*(2*h11L*gamma132 + h21L*gamma232 + h31L*(gamma121 + gamma332)) 
      + h32L*pow(gamma221,2)) + hInv23*(h32L*gamma121*gamma231 + 
      h32L*gamma221*gamma232 + h22L*gamma231*gamma232 + 
      h31L*gamma133*gamma321 + h32L*gamma233*gamma321 + 
      h33L*gamma121*gamma331 + h32L*gamma232*gamma331 + 
      gamma132*(2*h21L*gamma231 + h11L*gamma232 + h31L*(gamma221 + 
      2*gamma331)) + h11L*gamma133*gamma332 + h33L*gamma221*gamma332 + 
      h32L*gamma231*gamma332 + h31L*gamma232*gamma332 + 
      h21L*gamma233*gamma332 + h33L*gamma331*gamma332 + 
      gamma131*(2*h11L*gamma132 + h21L*gamma232 + h31L*(gamma121 + gamma332)) 
      + h33L*gamma321*gamma333 + h31L*gamma332*gamma333 + 
      h21L*pow(gamma232,2)) + hInv11*(h31L*gamma121*gamma211 + 
      2*h21L*gamma131*gamma211 + h32L*gamma211*gamma221 + 
      h11L*gamma121*gamma231 + h22L*gamma211*gamma231 + 
      h21L*gamma221*gamma231 + 3*h31L*gamma131*gamma311 + 
      2*h32L*gamma231*gamma311 + h33L*gamma211*gamma321 + 
      h31L*gamma231*gamma321 + h11L*gamma131*gamma331 + 
      h32L*gamma211*gamma331 + h21L*gamma231*gamma331 + 
      2*h33L*gamma311*gamma331 + gamma111*(2*h11L*gamma131 + h32L*gamma211 + 
      h21L*gamma231 + h33L*gamma311 + h31L*gamma331) + h31L*pow(gamma111,2) + 
      h31L*pow(gamma331,2)) + hInv13*(h32L*gamma131*gamma211 + 
      2*h21L*gamma133*gamma211 + h31L*gamma121*gamma231 + 
      h32L*gamma221*gamma231 + h11L*gamma121*gamma233 + 
      h22L*gamma211*gamma233 + h21L*gamma221*gamma233 + 
      h33L*gamma131*gamma311 + 2*h31L*gamma133*gamma311 + 
      h32L*gamma233*gamma311 + h33L*gamma231*gamma321 + 
      h31L*gamma233*gamma321 + h31L*gamma131*gamma331 + 
      h32L*gamma231*gamma331 + h11L*gamma131*gamma333 + 
      h32L*gamma211*gamma333 + h21L*gamma231*gamma333 + 
      h33L*gamma311*gamma333 + h31L*gamma331*gamma333 + 
      gamma111*(2*h11L*gamma133 + h21L*gamma233 + h31L*(gamma131 + gamma333)) 
      + h33L*pow(gamma331,2)) + hInv13*(h31L*gamma132*gamma211 + 
      h32L*gamma111*gamma231 + h11L*gamma132*gamma231 + 
      h32L*gamma211*gamma232 + h21L*gamma231*gamma232 + 
      h31L*gamma133*gamma311 + h32L*gamma233*gamma311 + 
      h33L*gamma111*gamma331 + h11L*gamma133*gamma331 + 
      2*h32L*gamma231*gamma331 + h21L*gamma233*gamma331 + 
      gamma131*(3*h21L*gamma231 + h31L*(gamma111 + 3*gamma331)) + 
      h33L*gamma211*gamma332 + h31L*gamma231*gamma332 + 
      h33L*gamma311*gamma333 + h31L*gamma331*gamma333 + 
      2*h11L*pow(gamma131,2) + h22L*pow(gamma231,2) + h33L*pow(gamma331,2)) + 
      hInv22*(h31L*gamma122*gamma221 + 2*h21L*gamma132*gamma221 + 
      h32L*gamma221*gamma222 + h11L*gamma122*gamma232 + 
      h22L*gamma221*gamma232 + h21L*gamma222*gamma232 + 
      3*h31L*gamma132*gamma321 + 2*h32L*gamma232*gamma321 + 
      h33L*gamma221*gamma322 + h31L*gamma232*gamma322 + 
      h11L*gamma132*gamma332 + h32L*gamma221*gamma332 + 
      h21L*gamma232*gamma332 + 2*h33L*gamma321*gamma332 + 
      gamma121*(2*h11L*gamma132 + h32L*gamma221 + h21L*gamma232 + 
      h33L*gamma321 + h31L*gamma332) + h31L*pow(gamma121,2) + 
      h31L*pow(gamma332,2)) + hInv33*(h31L*gamma132*gamma231 + 
      2*h21L*gamma133*gamma231 + h32L*gamma231*gamma232 + 
      h11L*gamma132*gamma233 + h22L*gamma231*gamma233 + 
      h21L*gamma232*gamma233 + 3*h31L*gamma133*gamma331 + 
      2*h32L*gamma233*gamma331 + h33L*gamma231*gamma332 + 
      h31L*gamma233*gamma332 + h11L*gamma133*gamma333 + 
      h32L*gamma231*gamma333 + h21L*gamma233*gamma333 + 
      2*h33L*gamma331*gamma333 + gamma131*(2*h11L*gamma133 + h32L*gamma231 + 
      h21L*gamma233 + h33L*gamma331 + h31L*gamma333) + h31L*pow(gamma131,2) + 
      h31L*pow(gamma333,2))));
    
    CCTK_REAL R22 CCTK_ATTRIBUTE_UNUSED = gamma222*(4*h21L*gamma121*hInv12 
      + 6*h22L*gamma221*hInv12 + 4*h21L*gamma122*hInv22 + 
      4*h21L*gamma132*hInv23) + gamma221*(2*h22L*gamma121*hInv12 + 
      4*h21L*gamma122*hInv12 + 4*h21L*gamma132*hInv13 + 
      2*h22L*gamma122*hInv22 + 2*h22L*gamma132*hInv23) + 
      gamma322*(4*h32L*gamma221*hInv12 + 2*h22L*gamma231*hInv12 + 
      2*h31L*gamma122*hInv22 + 2*h31L*gamma132*hInv23) + 
      gamma321*(2*h21L*gamma131*hInv11 + 4*h32L*gamma221*hInv11 + 
      2*h32L*gamma121*hInv12 + 4*h32L*gamma222*hInv12 + 
      2*h22L*gamma232*hInv12 + 2*h32L*gamma332*hInv12 + 
      2*h31L*gamma132*hInv13 + 2*h21L*gamma133*hInv13 + 
      4*h32L*gamma232*hInv13 + 2*h32L*gamma333*hInv13 + 
      2*h32L*gamma132*hInv23) + gamma232*(4*h21L*gamma121*hInv13 + 
      6*h22L*gamma221*hInv13 + 2*h22L*gamma322*hInv22 + 
      4*h21L*gamma122*hInv23 + 4*h21L*gamma132*hInv33) + 
      gamma332*(4*h32L*gamma221*hInv13 + 2*h22L*gamma231*hInv13 + 
      2*h32L*gamma322*hInv22 + 2*h31L*gamma122*hInv23 + 
      4*h32L*gamma222*hInv23 + 2*h22L*gamma232*hInv23 + 
      2*h31L*gamma132*hInv33) - hInv12*PDstandardNth12h22 - 
      hInv13*PDstandardNth13h22 + hInv23*(2*h11L*gamma122*gamma132 + 
      6*h22L*gamma222*gamma232 + 2*h21L*gamma133*gamma322 + 
      2*h22L*gamma233*gamma322 + 2*h32L*gamma122*gamma331 + 
      2*h21L*gamma132*gamma332 + 2*h32L*gamma322*gamma333 - 
      PDstandardNth23h22) + h32L*(Gam1L*gamma321 + Gam2L*gamma322 + 
      Gam3L*gamma332 + 2*gamma121*gamma311*hInv11 + 
      2*gamma122*gamma311*hInv12 + 2*gamma132*gamma311*hInv13 + 
      2*gamma122*gamma321*hInv22 + 4*gamma232*gamma322*hInv23 + 
      2*gamma132*gamma331*hInv33 + PDstandardNth2Gam3) + h21L*(Gam1L*gamma121 
      + Gam2L*gamma122 + Gam3L*gamma132 + 2*gamma111*gamma121*hInv11 + 
      2*gamma121*gamma131*hInv13 + 2*gamma121*gamma122*hInv22 + 
      2*gamma122*gamma131*hInv23 + 2*gamma131*gamma132*hInv33 + 
      PDstandardNth2Gam1 + 2*hInv12*pow(gamma121,2)) + h22L*(Gam1L*gamma221 + 
      Gam2L*gamma222 + Gam3L*gamma232 + 2*gamma121*gamma211*hInv11 + 
      2*gamma122*gamma211*hInv12 + 2*gamma132*gamma211*hInv13 + 
      2*gamma122*gamma231*hInv23 + 2*gamma132*gamma231*hInv33 + 
      PDstandardNth2Gam2 + 3*hInv22*pow(gamma222,2)) + 
      hInv11*(4*h21L*gamma121*gamma221 + 2*h31L*gamma121*gamma321 + 
      2*h22L*gamma231*gamma321 + 2*h32L*gamma321*gamma331 - 
      0.5*PDstandardNth11h22 + h11L*pow(gamma121,2) + 3*h22L*pow(gamma221,2) 
      + h33L*pow(gamma321,2)) + hInv22*(4*h32L*gamma222*gamma322 - 
      0.5*PDstandardNth22h22 + h11L*pow(gamma122,2) + h33L*pow(gamma322,2)) + 
      hInv33*(4*h32L*gamma232*gamma332 + 2*h22L*gamma233*gamma332 + 
      2*h32L*gamma332*gamma333 - 0.5*PDstandardNth33h22 + 
      h11L*pow(gamma132,2) + 3*h22L*pow(gamma232,2) + h33L*pow(gamma332,2)) + 
      2*(h21L*gamma111*gamma122*hInv12 + h11L*gamma121*gamma122*hInv12 + 
      h31L*gamma122*gamma321*hInv12 + h21L*gamma132*gamma321*hInv12 + 
      h31L*gamma121*gamma322*hInv12 + h21L*gamma131*gamma322*hInv12 + 
      h33L*gamma321*gamma322*hInv12 + h32L*gamma322*gamma331*hInv12 + 
      h21L*gamma111*gamma132*hInv13 + h11L*gamma121*gamma132*hInv13 + 
      h22L*gamma121*gamma231*hInv13 + h22L*gamma233*gamma321*hInv13 + 
      h32L*gamma121*gamma331*hInv13 + h31L*gamma121*gamma332*hInv13 + 
      h21L*gamma131*gamma332*hInv13 + h33L*gamma321*gamma332*hInv13 + 
      h32L*gamma331*gamma332*hInv13 + h21L*gamma132*gamma322*hInv22 + 
      h21L*gamma121*gamma132*hInv23 + h33L*gamma322*gamma332*hInv23 + 
      h21L*gamma133*gamma332*hInv33 + h32L*hInv23*pow(gamma332,2));
    
    CCTK_REAL R32 CCTK_ATTRIBUTE_UNUSED = 0.5*(Gam1L*(h31L*gamma121 + 
      h21L*gamma131 + h32L*gamma221 + h22L*gamma231 + h33L*gamma321 + 
      h32L*gamma331) + Gam2L*(h31L*gamma122 + h21L*gamma132 + h32L*gamma222 + 
      h22L*gamma232 + h33L*gamma322 + h32L*gamma332) + Gam3L*(h31L*gamma132 + 
      h21L*gamma133 + h32L*gamma232 + h22L*gamma233 + h33L*gamma332 + 
      h32L*gamma333) - hInv11*PDstandardNth11h32 - 
      2*hInv12*PDstandardNth12h32 - 2*hInv13*PDstandardNth13h32 - 
      hInv22*PDstandardNth22h32 - 2*hInv23*PDstandardNth23h32 + 
      h31L*PDstandardNth2Gam1 + h32L*PDstandardNth2Gam2 + 
      h33L*PDstandardNth2Gam3 - hInv33*PDstandardNth33h32 + 
      h21L*PDstandardNth3Gam1 + h22L*PDstandardNth3Gam2 + 
      h32L*PDstandardNth3Gam3 + 2*((h31L*gamma111*gamma122 + 
      h21L*gamma111*gamma132 + h32L*gamma122*gamma211 + 
      h22L*gamma132*gamma211 + h21L*gamma132*gamma221 + 
      h32L*gamma221*gamma222 + 2*h22L*gamma221*gamma232 + 
      h33L*gamma122*gamma311 + h32L*gamma132*gamma311 + 
      h31L*gamma132*gamma321 + h33L*gamma222*gamma321 + 
      2*h32L*gamma232*gamma321 + h31L*gamma131*gamma322 + 
      h32L*gamma231*gamma322 + h33L*gamma322*gamma331 + 
      h21L*gamma131*gamma332 + h32L*gamma221*gamma332 + 
      h22L*gamma231*gamma332 + h33L*gamma321*gamma332 + 
      h32L*gamma331*gamma332 + gamma121*(h11L*gamma132 + 2*h21L*gamma232 + 
      h31L*(gamma222 + gamma332)))*hInv12 + (h31L*gamma111*gamma132 + 
      h21L*gamma111*gamma133 + h32L*gamma132*gamma211 + 
      h22L*gamma133*gamma211 + h21L*gamma133*gamma221 + 
      h32L*gamma221*gamma232 + 2*h22L*gamma221*gamma233 + 
      h33L*gamma132*gamma311 + h32L*gamma133*gamma311 + 
      h31L*gamma133*gamma321 + h33L*gamma232*gamma321 + 
      2*h32L*gamma233*gamma321 + h31L*gamma131*gamma332 + 
      h32L*gamma231*gamma332 + h33L*gamma331*gamma332 + 
      h21L*gamma131*gamma333 + h32L*gamma221*gamma333 + 
      h22L*gamma231*gamma333 + h33L*gamma321*gamma333 + 
      h32L*gamma331*gamma333 + gamma121*(h11L*gamma133 + 2*h21L*gamma233 + 
      h31L*(gamma232 + gamma333)))*hInv13 + hInv12*(h22L*gamma131*gamma221 + 
      h21L*gamma131*gamma222 + h32L*gamma221*gamma222 + 
      2*h22L*gamma222*gamma231 + h32L*gamma131*gamma321 + 
      h31L*gamma132*gamma321 + h32L*gamma232*gamma321 + 
      gamma121*(h21L*gamma131 + h32L*gamma221 + h33L*gamma321) + 
      h31L*gamma131*gamma322 + h33L*gamma221*gamma322 + 
      2*h32L*gamma231*gamma322 + h21L*gamma132*gamma331 + 
      h32L*gamma222*gamma331 + h22L*gamma232*gamma331 + 
      h33L*gamma322*gamma331 + gamma122*(h11L*gamma131 + 2*h21L*gamma231 + 
      h31L*(gamma221 + gamma331)) + h33L*gamma321*gamma332 + 
      h32L*gamma331*gamma332 + h31L*pow(gamma121,2)) + 
      hInv13*(h32L*gamma121*gamma231 + h32L*gamma221*gamma232 + 
      2*h22L*gamma231*gamma232 + h31L*gamma133*gamma321 + 
      h32L*gamma233*gamma321 + h33L*gamma121*gamma331 + 
      h21L*gamma133*gamma331 + h32L*gamma232*gamma331 + 
      h22L*gamma233*gamma331 + gamma132*(2*h21L*gamma231 + h31L*(gamma221 + 
      gamma331)) + h33L*gamma221*gamma332 + 2*h32L*gamma231*gamma332 + 
      h33L*gamma331*gamma332 + gamma131*(h31L*gamma121 + h11L*gamma132 + 
      h22L*gamma231 + h21L*gamma232 + h32L*gamma331 + h31L*gamma332) + 
      h33L*gamma321*gamma333 + h32L*gamma331*gamma333 + h21L*pow(gamma131,2)) 
      + hInv11*(h21L*gamma111*gamma131 + h22L*gamma131*gamma211 + 
      h21L*gamma131*gamma221 + 2*h22L*gamma221*gamma231 + 
      h32L*gamma131*gamma311 + 2*h31L*gamma131*gamma321 + 
      h33L*gamma221*gamma321 + 3*h32L*gamma231*gamma321 + 
      h21L*gamma131*gamma331 + h32L*gamma221*gamma331 + 
      h22L*gamma231*gamma331 + 2*h33L*gamma321*gamma331 + 
      gamma121*(h31L*gamma111 + h11L*gamma131 + h32L*gamma211 + h31L*gamma221 
      + 2*h21L*gamma231 + h33L*gamma311 + h31L*gamma331) + 
      h32L*pow(gamma221,2) + h32L*pow(gamma331,2)) + 
      hInv22*(h21L*gamma121*gamma132 + h22L*gamma132*gamma221 + 
      h21L*gamma132*gamma222 + 2*h22L*gamma222*gamma232 + 
      h32L*gamma132*gamma321 + 2*h31L*gamma132*gamma322 + 
      h33L*gamma222*gamma322 + 3*h32L*gamma232*gamma322 + 
      h21L*gamma132*gamma332 + h32L*gamma222*gamma332 + 
      h22L*gamma232*gamma332 + 2*h33L*gamma322*gamma332 + 
      gamma122*(h31L*gamma121 + h11L*gamma132 + h32L*gamma221 + h31L*gamma222 
      + 2*h21L*gamma232 + h33L*gamma321 + h31L*gamma332) + 
      h32L*pow(gamma222,2) + h32L*pow(gamma332,2)) + 
      hInv23*(h31L*gamma121*gamma132 + h21L*gamma121*gamma133 + 
      h32L*gamma132*gamma221 + h22L*gamma133*gamma221 + 
      h21L*gamma133*gamma222 + h32L*gamma222*gamma232 + 
      2*h22L*gamma222*gamma233 + h33L*gamma132*gamma321 + 
      h32L*gamma133*gamma321 + h31L*gamma133*gamma322 + 
      h33L*gamma232*gamma322 + 2*h32L*gamma233*gamma322 + 
      h31L*gamma132*gamma332 + h32L*gamma232*gamma332 + 
      h21L*gamma132*gamma333 + h32L*gamma222*gamma333 + 
      h22L*gamma232*gamma333 + h33L*gamma322*gamma333 + 
      h32L*gamma332*gamma333 + gamma122*(h11L*gamma133 + 2*h21L*gamma233 + 
      h31L*(gamma232 + gamma333)) + h33L*pow(gamma332,2)) + 
      hInv23*(h31L*gamma122*gamma131 + h32L*gamma122*gamma231 + 
      h32L*gamma222*gamma232 + h31L*gamma133*gamma322 + 
      h32L*gamma233*gamma322 + h33L*gamma122*gamma331 + 
      h21L*gamma133*gamma332 + h33L*gamma222*gamma332 + 
      3*h32L*gamma232*gamma332 + h22L*gamma233*gamma332 + 
      gamma132*(h21L*gamma131 + h31L*gamma222 + h22L*gamma231 + 
      3*h21L*gamma232 + h32L*gamma331 + 2*h31L*gamma332) + 
      h33L*gamma322*gamma333 + h32L*gamma332*gamma333 + h11L*pow(gamma132,2) 
      + 2*h22L*pow(gamma232,2) + h33L*pow(gamma332,2)) + 
      hInv33*(h21L*gamma131*gamma133 + h22L*gamma133*gamma231 + 
      h21L*gamma133*gamma232 + 2*h22L*gamma232*gamma233 + 
      h32L*gamma133*gamma331 + 2*h31L*gamma133*gamma332 + 
      h33L*gamma232*gamma332 + 3*h32L*gamma233*gamma332 + 
      h21L*gamma133*gamma333 + h32L*gamma232*gamma333 + 
      h22L*gamma233*gamma333 + 2*h33L*gamma332*gamma333 + 
      gamma132*(h31L*gamma131 + h11L*gamma133 + h32L*gamma231 + h31L*gamma232 
      + 2*h21L*gamma233 + h33L*gamma331 + h31L*gamma333) + 
      h32L*pow(gamma232,2) + h32L*pow(gamma333,2))));
    
    CCTK_REAL R33 CCTK_ATTRIBUTE_UNUSED = gamma333*(4*h32L*gamma231*hInv13 
      + 6*h33L*gamma331*hInv13 + 4*h31L*gamma132*hInv23 + 
      4*h32L*gamma232*hInv23 + 6*h33L*gamma332*hInv23 + 
      4*h31L*gamma133*hInv33) + gamma331*(2*h33L*gamma131*hInv13 + 
      4*h32L*gamma233*hInv13 + 2*h33L*gamma132*hInv23 + 
      2*h33L*gamma133*hInv33) + gamma332*(2*h33L*gamma232*hInv23 + 
      4*h32L*gamma233*hInv23 + 2*h33L*gamma233*hInv33) + 
      4*(h31L*gamma132*gamma331*hInv12 + h32L*gamma232*gamma331*hInv12 + 
      h31L*gamma131*gamma332*hInv12 + h32L*gamma231*gamma332*hInv12 + 
      h31L*gamma133*gamma331*hInv13 + h31L*gamma132*gamma332*hInv22 + 
      h31L*gamma133*gamma332*hInv23 + h32L*gamma233*gamma333*hInv33) - 
      hInv12*PDstandardNth12h33 + hInv13*(2*h33L*gamma231*gamma332 + 
      4*h31L*gamma131*gamma333 - PDstandardNth13h33) - 
      hInv23*PDstandardNth23h33 + h32L*(Gam1L*gamma231 + Gam2L*gamma232 + 
      Gam3L*gamma233 + 2*gamma131*gamma211*hInv11 + 
      2*gamma132*gamma211*hInv12 + 2*gamma133*gamma211*hInv13 + 
      2*gamma132*gamma221*hInv22 + 2*gamma133*gamma221*hInv23 + 
      2*gamma133*gamma231*hInv33 + PDstandardNth3Gam2) + h31L*(Gam1L*gamma131 
      + Gam2L*gamma132 + Gam3L*gamma133 + 2*gamma111*gamma131*hInv11 + 
      2*gamma121*gamma131*hInv12 + 2*gamma121*gamma132*hInv22 + 
      2*gamma131*gamma132*hInv23 + 2*gamma131*gamma133*hInv33 + 
      PDstandardNth3Gam1 + 2*hInv13*pow(gamma131,2)) + 
      2*(h21L*gamma131*gamma231*hInv11 + h32L*gamma221*gamma231*hInv11 + 
      h31L*gamma111*gamma132*hInv12 + h11L*gamma131*gamma132*hInv12 + 
      h32L*gamma131*gamma221*hInv12 + h31L*gamma122*gamma231*hInv12 + 
      h21L*gamma132*gamma231*hInv12 + h32L*gamma222*gamma231*hInv12 + 
      h31L*gamma121*gamma232*hInv12 + h21L*gamma131*gamma232*hInv12 + 
      h32L*gamma221*gamma232*hInv12 + h22L*gamma231*gamma232*hInv12 + 
      h33L*gamma131*gamma321*hInv12 + h33L*gamma232*gamma321*hInv12 + 
      h33L*gamma231*gamma322*hInv12 + h31L*gamma111*gamma133*hInv13 + 
      h11L*gamma131*gamma133*hInv13 + h32L*gamma131*gamma231*hInv13 + 
      h31L*gamma132*gamma231*hInv13 + h21L*gamma133*gamma231*hInv13 + 
      h32L*gamma231*gamma232*hInv13 + h31L*gamma121*gamma233*hInv13 + 
      h21L*gamma131*gamma233*hInv13 + h32L*gamma221*gamma233*hInv13 + 
      h22L*gamma231*gamma233*hInv13 + h33L*gamma233*gamma321*hInv13 + 
      h31L*gamma122*gamma232*hInv22 + h31L*gamma121*gamma133*hInv23 + 
      h11L*gamma132*gamma133*hInv23 + h32L*gamma132*gamma231*hInv23 + 
      h31L*gamma132*gamma232*hInv23 + h21L*gamma133*gamma232*hInv23 + 
      h31L*gamma122*gamma233*hInv23 + h21L*gamma132*gamma233*hInv23 + 
      h32L*gamma222*gamma233*hInv23 + h22L*gamma232*gamma233*hInv23 + 
      h33L*gamma233*gamma322*hInv23 + h31L*gamma132*gamma233*hInv33 + 
      h32L*hInv23*pow(gamma232,2)) + hInv33*(2*h21L*gamma133*gamma233 + 
      2*h32L*gamma232*gamma233 - 0.5*PDstandardNth33h33 + 
      h11L*pow(gamma133,2) + h22L*pow(gamma233,2)) + 
      hInv11*(2*h31L*gamma121*gamma231 + 2*h33L*gamma231*gamma321 + 
      4*h31L*gamma131*gamma331 + 4*h32L*gamma231*gamma331 - 
      0.5*PDstandardNth11h33 + h11L*pow(gamma131,2) + h22L*pow(gamma231,2) + 
      3*h33L*pow(gamma331,2)) + hInv22*(2*h21L*gamma132*gamma232 + 
      2*h32L*gamma222*gamma232 + 2*h33L*gamma232*gamma322 + 
      4*h32L*gamma232*gamma332 - 0.5*PDstandardNth22h33 + 
      h11L*pow(gamma132,2) + h22L*pow(gamma232,2) + 3*h33L*pow(gamma332,2)) + 
      h33L*(Gam1L*gamma331 + Gam2L*gamma332 + Gam3L*gamma333 + 
      2*gamma131*gamma311*hInv11 + 2*gamma132*gamma311*hInv12 + 
      6*gamma331*gamma332*hInv12 + 2*gamma133*gamma311*hInv13 + 
      2*gamma132*gamma321*hInv22 + 2*gamma133*gamma321*hInv23 + 
      PDstandardNth3Gam3 + 3*hInv33*pow(gamma333,2));
    
    CCTK_REAL Rphi11 CCTK_ATTRIBUTE_UNUSED = -2*(DDphi11 + 
      h11L*(DDphi11*hInv11 + (DDphi12 + DDphi21)*hInv12 + (DDphi13 + 
      DDphi31)*hInv13 + DDphi22*hInv22 + (DDphi23 + DDphi32)*hInv23 + 
      DDphi33*hInv33)) + 4*pow(Dphi1,2) - 4*h11L*(2*(Dphi1*(Dphi2*hInv12 + 
      Dphi3*hInv13) + Dphi2*Dphi3*hInv23) + hInv11*pow(Dphi1,2) + 
      hInv22*pow(Dphi2,2) + hInv33*pow(Dphi3,2));
    
    CCTK_REAL Rphi21 CCTK_ATTRIBUTE_UNUSED = 4*Dphi1*Dphi2 - 2*(DDphi21 + 
      h21L*(DDphi11*hInv11 + (DDphi12 + DDphi21)*hInv12 + (DDphi13 + 
      DDphi31)*hInv13 + DDphi22*hInv22 + (DDphi23 + DDphi32)*hInv23 + 
      DDphi33*hInv33)) - 4*h21L*(2*(Dphi1*(Dphi2*hInv12 + Dphi3*hInv13) + 
      Dphi2*Dphi3*hInv23) + hInv11*pow(Dphi1,2) + hInv22*pow(Dphi2,2) + 
      hInv33*pow(Dphi3,2));
    
    CCTK_REAL Rphi31 CCTK_ATTRIBUTE_UNUSED = 4*Dphi1*Dphi3 - 2*(DDphi31 + 
      h31L*(DDphi11*hInv11 + (DDphi12 + DDphi21)*hInv12 + (DDphi13 + 
      DDphi31)*hInv13 + DDphi22*hInv22 + (DDphi23 + DDphi32)*hInv23 + 
      DDphi33*hInv33)) - 4*h31L*(2*(Dphi1*(Dphi2*hInv12 + Dphi3*hInv13) + 
      Dphi2*Dphi3*hInv23) + hInv11*pow(Dphi1,2) + hInv22*pow(Dphi2,2) + 
      hInv33*pow(Dphi3,2));
    
    CCTK_REAL Rphi22 CCTK_ATTRIBUTE_UNUSED = -2*(DDphi22 + 
      h22L*(DDphi11*hInv11 + (DDphi12 + DDphi21)*hInv12 + (DDphi13 + 
      DDphi31)*hInv13 + DDphi22*hInv22 + (DDphi23 + DDphi32)*hInv23 + 
      DDphi33*hInv33)) + 4*pow(Dphi2,2) - 4*h22L*(2*(Dphi1*(Dphi2*hInv12 + 
      Dphi3*hInv13) + Dphi2*Dphi3*hInv23) + hInv11*pow(Dphi1,2) + 
      hInv22*pow(Dphi2,2) + hInv33*pow(Dphi3,2));
    
    CCTK_REAL Rphi32 CCTK_ATTRIBUTE_UNUSED = 4*Dphi2*Dphi3 - 2*(DDphi32 + 
      h32L*(DDphi11*hInv11 + (DDphi12 + DDphi21)*hInv12 + (DDphi13 + 
      DDphi31)*hInv13 + DDphi22*hInv22 + (DDphi23 + DDphi32)*hInv23 + 
      DDphi33*hInv33)) - 4*h32L*(2*(Dphi1*(Dphi2*hInv12 + Dphi3*hInv13) + 
      Dphi2*Dphi3*hInv23) + hInv11*pow(Dphi1,2) + hInv22*pow(Dphi2,2) + 
      hInv33*pow(Dphi3,2));
    
    CCTK_REAL Rphi33 CCTK_ATTRIBUTE_UNUSED = -2*(DDphi33 + 
      h33L*(DDphi11*hInv11 + (DDphi12 + DDphi21)*hInv12 + (DDphi13 + 
      DDphi31)*hInv13 + DDphi22*hInv22 + (DDphi23 + DDphi32)*hInv23 + 
      DDphi33*hInv33)) + 4*pow(Dphi3,2) - 4*h33L*(2*(Dphi1*(Dphi2*hInv12 + 
      Dphi3*hInv13) + Dphi2*Dphi3*hInv23) + hInv11*pow(Dphi1,2) + 
      hInv22*pow(Dphi2,2) + hInv33*pow(Dphi3,2));
    
    bssnhamL = -(A11L*AInv11) - 2*A21L*AInv12 - 2*A31L*AInv13 - 
      A22L*AInv22 - 2*A32L*AInv23 - A33L*AInv33 + em4phi*(hInv11*(R11 + 
      Rphi11) + 2*hInv12*(R21 + Rphi21) + hInv22*(R22 + Rphi22) + 
      2*hInv13*(R31 + Rphi31) + 2*hInv23*(R32 + Rphi32) + hInv33*(R33 + 
      Rphi33)) + 0.666666666666666666666666666667*pow(KL,2);
    
    bssnmom1L = (-(A22L*gamma221) + A21L*(-gamma121 - gamma222) - 
      A32L*gamma321)*hInv22 + (6*A31L*Dphi2 - A32L*gamma221 - 
      A33L*gamma321)*hInv23 + (6*A21L*Dphi3 - A22L*gamma231 - 
      A32L*gamma331)*hInv23 + A21L*((6*Dphi1 - gamma111)*hInv12 - 
      3*gamma231*hInv13 + 6*Dphi2*hInv22 - gamma131*hInv23) - 
      2*(A21L*gamma211*hInv11 + A31L*gamma311*hInv11 + A11L*gamma132*hInv23) 
      + (6*A31L*Dphi3 - A32L*gamma231 - A33L*gamma331)*hInv33 + 
      A31L*(-2*gamma332*hInv23 - gamma131*hInv33) + A11L*(6*Dphi1*hInv11 - 
      2*gamma111*hInv11 + 6*Dphi2*hInv12 + 6*Dphi3*hInv13 - gamma122*hInv22 - 
      gamma133*hInv33) + A21L*(-2*gamma232*hInv23 - gamma233*hInv33) + 
      A31L*((6*Dphi1 - gamma111)*hInv13 - gamma322*hInv22 - gamma121*hInv23 - 
      gamma333*hInv33) + hInv11*PDstandardNth1A11 + hInv12*(-3*A11L*gamma121 
      - A22L*gamma211 - 3*A21L*gamma221 - A32L*gamma311 - 3*A31L*gamma321 + 
      PDstandardNth1A21) + hInv13*(-3*A11L*gamma131 - A32L*gamma211 - 
      A33L*gamma311 - 3*A31L*gamma331 + PDstandardNth1A31) - 
      0.666666666666666666666666666667*PDstandardNth1K + 
      hInv12*PDstandardNth2A11 + hInv22*PDstandardNth2A21 + 
      hInv23*PDstandardNth2A31 + hInv13*PDstandardNth3A11 + 
      hInv23*PDstandardNth3A21 + hInv33*PDstandardNth3A31;
    
    bssnmom2L = (-(A11L*gamma121) - A22L*gamma211 - A21L*gamma221 - 
      A32L*gamma311 - A31L*gamma321)*hInv11 + (6*A22L*Dphi1 - 
      A11L*gamma122)*hInv12 + (6*A21L*Dphi2 - A31L*gamma322)*hInv12 + 
      (6*A32L*Dphi1 - A11L*gamma132)*hInv13 + (6*A21L*Dphi3 - 
      A31L*gamma332)*hInv13 + gamma121*(-3*A21L*hInv12 - A31L*hInv13) + 
      gamma221*(-3*A22L*hInv12 - A32L*hInv13) + gamma321*(-3*A32L*hInv12 - 
      A33L*hInv13) + A21L*((6*Dphi1 - gamma111)*hInv11 - gamma222*hInv12 - 
      gamma232*hInv13) + 6*A22L*Dphi2*hInv22 - 2*A22L*gamma222*hInv22 + 
      6*A22L*Dphi3*hInv23 - 3*A22L*gamma232*hInv23 + gamma122*(-2*A21L*hInv22 
      - A31L*hInv23) + gamma322*(-2*A32L*hInv22 - A33L*hInv23) + 
      6*A32L*Dphi3*hInv33 + gamma132*(-3*A21L*hInv23 - A31L*hInv33) + 
      gamma332*(-3*A32L*hInv23 - A33L*hInv33) + A21L*(-2*gamma131*hInv13 - 
      gamma133*hInv33) + A32L*(6*Dphi2*hInv23 - gamma232*hInv33) + 
      A22L*(-2*gamma231*hInv13 - gamma233*hInv33) + A32L*(-2*gamma331*hInv13 
      - gamma222*hInv23 - gamma333*hInv33) + hInv11*PDstandardNth1A21 + 
      hInv12*PDstandardNth1A22 + hInv13*PDstandardNth1A32 + 
      hInv12*PDstandardNth2A21 + hInv22*PDstandardNth2A22 + 
      hInv23*PDstandardNth2A32 - 
      0.666666666666666666666666666667*PDstandardNth2K + 
      hInv13*PDstandardNth3A21 + hInv23*PDstandardNth3A22 + 
      hInv33*PDstandardNth3A32;
    
    bssnmom3L = (-(A11L*gamma131) - A32L*gamma211 - A21L*gamma231 - 
      A33L*gamma311 - A31L*gamma331)*hInv11 + (-2*A31L*gamma121 - 
      A22L*gamma231)*hInv12 + (6*A32L*Dphi1 + A21L*(-gamma131 - 
      gamma232))*hInv12 + (6*A31L*Dphi2 - A11L*gamma132 - 
      A32L*gamma331)*hInv12 + (6*A33L*Dphi1 - A11L*gamma133)*hInv13 + 
      (6*A31L*Dphi3 - A21L*gamma233)*hInv13 - 3*A33L*gamma331*hInv13 + 
      A31L*((6*Dphi1 - gamma111)*hInv11 - gamma332*hInv12 - gamma333*hInv13) 
      + (6*A32L*Dphi2 - A21L*gamma132 - A22L*gamma232)*hInv22 + 
      A31L*(-3*gamma131*hInv13 - gamma122*hInv22) + A33L*(-2*gamma321*hInv12 
      - gamma322*hInv22) + A32L*(-3*gamma231*hInv13 - gamma332*hInv22) - 
      3*A31L*gamma132*hInv23 + (6*A33L*Dphi2 - A21L*gamma133)*hInv23 - 
      3*A32L*gamma232*hInv23 + (6*A32L*Dphi3 - A22L*gamma233)*hInv23 - 
      3*A33L*gamma332*hInv23 + A32L*(-2*gamma221*hInv12 - gamma222*hInv22 - 
      gamma333*hInv23) + 6*A33L*Dphi3*hInv33 - 2*A31L*gamma133*hInv33 - 
      2*A32L*gamma233*hInv33 - 2*A33L*gamma333*hInv33 + 
      hInv11*PDstandardNth1A31 + hInv12*PDstandardNth1A32 + 
      hInv13*PDstandardNth1A33 + hInv12*PDstandardNth2A31 + 
      hInv22*PDstandardNth2A32 + hInv23*PDstandardNth2A33 + 
      hInv13*PDstandardNth3A31 + hInv23*PDstandardNth3A32 + 
      hInv33*PDstandardNth3A33 - 
      0.666666666666666666666666666667*PDstandardNth3K;
    
    CCTK_REAL bssnCcons1L CCTK_ATTRIBUTE_UNUSED = Gam1L - gamma111*hInv11 
      - 2*gamma121*hInv12 - 2*gamma131*hInv13 - gamma122*hInv22 - 
      2*gamma132*hInv23 - gamma133*hInv33;
    
    CCTK_REAL bssnCcons2L CCTK_ATTRIBUTE_UNUSED = Gam2L - gamma211*hInv11 
      - 2*gamma221*hInv12 - 2*gamma231*hInv13 - gamma222*hInv22 - 
      2*gamma232*hInv23 - gamma233*hInv33;
    
    CCTK_REAL bssnCcons3L CCTK_ATTRIBUTE_UNUSED = Gam3L - gamma311*hInv11 
      - 2*gamma321*hInv12 - 2*gamma331*hInv13 - gamma322*hInv22 - 
      2*gamma332*hInv23 - gamma333*hInv33;
    
    CCTK_REAL bssnDL CCTK_ATTRIBUTE_UNUSED = log(deth);
    
    CCTK_REAL bssnTL CCTK_ATTRIBUTE_UNUSED = A11L*hInv11 + A22L*hInv22 + 
      2*(A21L*hInv12 + A31L*hInv13 + A32L*hInv23) + A33L*hInv33;
    
    bssnhamL = bssnhamL - 16*pi*rho;
    
    bssnmom1L = bssnmom1L - 8*pi*Si1;
    
    bssnmom2L = bssnmom2L - 8*pi*Si2;
    
    bssnmom3L = bssnmom3L - 8*pi*Si3;
    /* Copy local copies back to grid functions */
    bssnCcons1[index] = bssnCcons1L;
    bssnCcons2[index] = bssnCcons2L;
    bssnCcons3[index] = bssnCcons3L;
    bssnD[index] = bssnDL;
    bssnham[index] = bssnhamL;
    bssnmom1[index] = bssnmom1L;
    bssnmom2[index] = bssnmom2L;
    bssnmom3[index] = bssnmom3L;
    bssnT[index] = bssnTL;
  }
  CCTK_ENDLOOP3(bssnchimatter_calc_constraints_matter);
}
extern "C" void bssnchimatter_calc_constraints_matter(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering bssnchimatter_calc_constraints_matter_Body");
  }
  if (cctk_iteration % bssnchimatter_calc_constraints_matter_calc_every != bssnchimatter_calc_constraints_matter_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "Kranc2BSSNChiMatter::A_group",
    "Kranc2BSSNChiMatter::alpha_group",
    "Kranc2BSSNChiMatter::beta_group",
    "Kranc2BSSNChiMatter::bssnCcons_group",
    "Kranc2BSSNChiMatter::bssnmom_group",
    "Kranc2BSSNChiMatter::chi_group",
    "Kranc2BSSNChiMatter::Gam_group",
    "Kranc2BSSNChiMatter::h_group",
    "Kranc2BSSNChiMatter::K_group",
    "Kranc2BSSNChiMatter::scalarconstraints"};
  AssertGroupStorage(cctkGH, "bssnchimatter_calc_constraints_matter", 10, groups);
  
  switch (constraintsfdorder)
  {
    case 2:
    {
      EnsureStencilFits(cctkGH, "bssnchimatter_calc_constraints_matter", 1, 1, 1);
      break;
    }
    
    case 4:
    {
      EnsureStencilFits(cctkGH, "bssnchimatter_calc_constraints_matter", 2, 2, 2);
      break;
    }
    
    case 6:
    {
      EnsureStencilFits(cctkGH, "bssnchimatter_calc_constraints_matter", 3, 3, 3);
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  LoopOverInterior(cctkGH, bssnchimatter_calc_constraints_matter_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving bssnchimatter_calc_constraints_matter_Body");
  }
}

} // namespace Kranc2BSSNChiMatter
