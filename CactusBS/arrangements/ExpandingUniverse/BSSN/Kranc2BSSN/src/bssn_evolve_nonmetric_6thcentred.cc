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

extern "C" void bssn_evolve_nonmetric_6thcentred_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % bssn_evolve_nonmetric_6thcentred_calc_every != bssn_evolve_nonmetric_6thcentred_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "Kranc2BSSN::A_grouprhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for Kranc2BSSN::A_grouprhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "Kranc2BSSN::Gam_grouprhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for Kranc2BSSN::Gam_grouprhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "Kranc2BSSN::K_grouprhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for Kranc2BSSN::K_grouprhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "Kranc2BSSN::scalarconstraints","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for Kranc2BSSN::scalarconstraints.");
  return;
}

static void bssn_evolve_nonmetric_6thcentred_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  CCTK_LOOP3(bssn_evolve_nonmetric_6thcentred,
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
    CCTK_REAL phiL CCTK_ATTRIBUTE_UNUSED = phi[index];
    
    /* Include user supplied include files */
    /* Precompute derivatives */
    CCTK_REAL PDstandard6th1A11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th2A11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th3A11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th1A21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th2A21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th3A21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th1A22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th2A22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th3A22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th1A31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th2A31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th3A31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th1A32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th2A32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th3A32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th1A33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th2A33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th3A33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th1alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th2alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th3alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th11alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th22alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th33alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th12alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th13alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th23alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th1beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th2beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th3beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th11beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th22beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th33beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th12beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th13beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th23beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th1beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th2beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th3beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th11beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th22beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th33beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th12beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th13beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th23beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th1beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th2beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th3beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th11beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th22beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th33beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th12beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th13beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th23beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th1Gam1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th2Gam1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th3Gam1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th1Gam2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th2Gam2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th3Gam2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th1Gam3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th2Gam3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th3Gam3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th1h11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th2h11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th3h11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th11h11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th22h11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th33h11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th12h11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th13h11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th23h11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th1h21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th2h21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th3h21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th11h21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th22h21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th33h21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th12h21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th13h21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th23h21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th1h22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th2h22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th3h22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th11h22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th22h22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th33h22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th12h22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th13h22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th23h22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th1h31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th2h31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th3h31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th11h31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th22h31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th33h31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th12h31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th13h31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th23h31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th1h32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th2h32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th3h32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th11h32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th22h32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th33h32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th12h32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th13h32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th23h32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th1h33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th2h33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th3h33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th11h33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th22h33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th33h33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th12h33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th13h33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th23h33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th1K CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th2K CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th3K CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th1phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th2phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th3phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th11phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th22phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th33phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th12phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th13phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard6th23phi CCTK_ATTRIBUTE_UNUSED;
    
    switch (constraintsfdorder)
    {
      case 2:
      {
        PDstandard6th1A11 = PDstandard6th1(&A11[index]);
        PDstandard6th2A11 = PDstandard6th2(&A11[index]);
        PDstandard6th3A11 = PDstandard6th3(&A11[index]);
        PDstandard6th1A21 = PDstandard6th1(&A21[index]);
        PDstandard6th2A21 = PDstandard6th2(&A21[index]);
        PDstandard6th3A21 = PDstandard6th3(&A21[index]);
        PDstandard6th1A22 = PDstandard6th1(&A22[index]);
        PDstandard6th2A22 = PDstandard6th2(&A22[index]);
        PDstandard6th3A22 = PDstandard6th3(&A22[index]);
        PDstandard6th1A31 = PDstandard6th1(&A31[index]);
        PDstandard6th2A31 = PDstandard6th2(&A31[index]);
        PDstandard6th3A31 = PDstandard6th3(&A31[index]);
        PDstandard6th1A32 = PDstandard6th1(&A32[index]);
        PDstandard6th2A32 = PDstandard6th2(&A32[index]);
        PDstandard6th3A32 = PDstandard6th3(&A32[index]);
        PDstandard6th1A33 = PDstandard6th1(&A33[index]);
        PDstandard6th2A33 = PDstandard6th2(&A33[index]);
        PDstandard6th3A33 = PDstandard6th3(&A33[index]);
        PDstandard6th1alpha = PDstandard6th1(&alpha[index]);
        PDstandard6th2alpha = PDstandard6th2(&alpha[index]);
        PDstandard6th3alpha = PDstandard6th3(&alpha[index]);
        PDstandard6th11alpha = PDstandard6th11(&alpha[index]);
        PDstandard6th22alpha = PDstandard6th22(&alpha[index]);
        PDstandard6th33alpha = PDstandard6th33(&alpha[index]);
        PDstandard6th12alpha = PDstandard6th12(&alpha[index]);
        PDstandard6th13alpha = PDstandard6th13(&alpha[index]);
        PDstandard6th23alpha = PDstandard6th23(&alpha[index]);
        PDstandard6th1beta1 = PDstandard6th1(&beta1[index]);
        PDstandard6th2beta1 = PDstandard6th2(&beta1[index]);
        PDstandard6th3beta1 = PDstandard6th3(&beta1[index]);
        PDstandard6th11beta1 = PDstandard6th11(&beta1[index]);
        PDstandard6th22beta1 = PDstandard6th22(&beta1[index]);
        PDstandard6th33beta1 = PDstandard6th33(&beta1[index]);
        PDstandard6th12beta1 = PDstandard6th12(&beta1[index]);
        PDstandard6th13beta1 = PDstandard6th13(&beta1[index]);
        PDstandard6th23beta1 = PDstandard6th23(&beta1[index]);
        PDstandard6th1beta2 = PDstandard6th1(&beta2[index]);
        PDstandard6th2beta2 = PDstandard6th2(&beta2[index]);
        PDstandard6th3beta2 = PDstandard6th3(&beta2[index]);
        PDstandard6th11beta2 = PDstandard6th11(&beta2[index]);
        PDstandard6th22beta2 = PDstandard6th22(&beta2[index]);
        PDstandard6th33beta2 = PDstandard6th33(&beta2[index]);
        PDstandard6th12beta2 = PDstandard6th12(&beta2[index]);
        PDstandard6th13beta2 = PDstandard6th13(&beta2[index]);
        PDstandard6th23beta2 = PDstandard6th23(&beta2[index]);
        PDstandard6th1beta3 = PDstandard6th1(&beta3[index]);
        PDstandard6th2beta3 = PDstandard6th2(&beta3[index]);
        PDstandard6th3beta3 = PDstandard6th3(&beta3[index]);
        PDstandard6th11beta3 = PDstandard6th11(&beta3[index]);
        PDstandard6th22beta3 = PDstandard6th22(&beta3[index]);
        PDstandard6th33beta3 = PDstandard6th33(&beta3[index]);
        PDstandard6th12beta3 = PDstandard6th12(&beta3[index]);
        PDstandard6th13beta3 = PDstandard6th13(&beta3[index]);
        PDstandard6th23beta3 = PDstandard6th23(&beta3[index]);
        PDstandard6th1Gam1 = PDstandard6th1(&Gam1[index]);
        PDstandard6th2Gam1 = PDstandard6th2(&Gam1[index]);
        PDstandard6th3Gam1 = PDstandard6th3(&Gam1[index]);
        PDstandard6th1Gam2 = PDstandard6th1(&Gam2[index]);
        PDstandard6th2Gam2 = PDstandard6th2(&Gam2[index]);
        PDstandard6th3Gam2 = PDstandard6th3(&Gam2[index]);
        PDstandard6th1Gam3 = PDstandard6th1(&Gam3[index]);
        PDstandard6th2Gam3 = PDstandard6th2(&Gam3[index]);
        PDstandard6th3Gam3 = PDstandard6th3(&Gam3[index]);
        PDstandard6th1h11 = PDstandard6th1(&h11[index]);
        PDstandard6th2h11 = PDstandard6th2(&h11[index]);
        PDstandard6th3h11 = PDstandard6th3(&h11[index]);
        PDstandard6th11h11 = PDstandard6th11(&h11[index]);
        PDstandard6th22h11 = PDstandard6th22(&h11[index]);
        PDstandard6th33h11 = PDstandard6th33(&h11[index]);
        PDstandard6th12h11 = PDstandard6th12(&h11[index]);
        PDstandard6th13h11 = PDstandard6th13(&h11[index]);
        PDstandard6th23h11 = PDstandard6th23(&h11[index]);
        PDstandard6th1h21 = PDstandard6th1(&h21[index]);
        PDstandard6th2h21 = PDstandard6th2(&h21[index]);
        PDstandard6th3h21 = PDstandard6th3(&h21[index]);
        PDstandard6th11h21 = PDstandard6th11(&h21[index]);
        PDstandard6th22h21 = PDstandard6th22(&h21[index]);
        PDstandard6th33h21 = PDstandard6th33(&h21[index]);
        PDstandard6th12h21 = PDstandard6th12(&h21[index]);
        PDstandard6th13h21 = PDstandard6th13(&h21[index]);
        PDstandard6th23h21 = PDstandard6th23(&h21[index]);
        PDstandard6th1h22 = PDstandard6th1(&h22[index]);
        PDstandard6th2h22 = PDstandard6th2(&h22[index]);
        PDstandard6th3h22 = PDstandard6th3(&h22[index]);
        PDstandard6th11h22 = PDstandard6th11(&h22[index]);
        PDstandard6th22h22 = PDstandard6th22(&h22[index]);
        PDstandard6th33h22 = PDstandard6th33(&h22[index]);
        PDstandard6th12h22 = PDstandard6th12(&h22[index]);
        PDstandard6th13h22 = PDstandard6th13(&h22[index]);
        PDstandard6th23h22 = PDstandard6th23(&h22[index]);
        PDstandard6th1h31 = PDstandard6th1(&h31[index]);
        PDstandard6th2h31 = PDstandard6th2(&h31[index]);
        PDstandard6th3h31 = PDstandard6th3(&h31[index]);
        PDstandard6th11h31 = PDstandard6th11(&h31[index]);
        PDstandard6th22h31 = PDstandard6th22(&h31[index]);
        PDstandard6th33h31 = PDstandard6th33(&h31[index]);
        PDstandard6th12h31 = PDstandard6th12(&h31[index]);
        PDstandard6th13h31 = PDstandard6th13(&h31[index]);
        PDstandard6th23h31 = PDstandard6th23(&h31[index]);
        PDstandard6th1h32 = PDstandard6th1(&h32[index]);
        PDstandard6th2h32 = PDstandard6th2(&h32[index]);
        PDstandard6th3h32 = PDstandard6th3(&h32[index]);
        PDstandard6th11h32 = PDstandard6th11(&h32[index]);
        PDstandard6th22h32 = PDstandard6th22(&h32[index]);
        PDstandard6th33h32 = PDstandard6th33(&h32[index]);
        PDstandard6th12h32 = PDstandard6th12(&h32[index]);
        PDstandard6th13h32 = PDstandard6th13(&h32[index]);
        PDstandard6th23h32 = PDstandard6th23(&h32[index]);
        PDstandard6th1h33 = PDstandard6th1(&h33[index]);
        PDstandard6th2h33 = PDstandard6th2(&h33[index]);
        PDstandard6th3h33 = PDstandard6th3(&h33[index]);
        PDstandard6th11h33 = PDstandard6th11(&h33[index]);
        PDstandard6th22h33 = PDstandard6th22(&h33[index]);
        PDstandard6th33h33 = PDstandard6th33(&h33[index]);
        PDstandard6th12h33 = PDstandard6th12(&h33[index]);
        PDstandard6th13h33 = PDstandard6th13(&h33[index]);
        PDstandard6th23h33 = PDstandard6th23(&h33[index]);
        PDstandard6th1K = PDstandard6th1(&K[index]);
        PDstandard6th2K = PDstandard6th2(&K[index]);
        PDstandard6th3K = PDstandard6th3(&K[index]);
        PDstandard6th1phi = PDstandard6th1(&phi[index]);
        PDstandard6th2phi = PDstandard6th2(&phi[index]);
        PDstandard6th3phi = PDstandard6th3(&phi[index]);
        PDstandard6th11phi = PDstandard6th11(&phi[index]);
        PDstandard6th22phi = PDstandard6th22(&phi[index]);
        PDstandard6th33phi = PDstandard6th33(&phi[index]);
        PDstandard6th12phi = PDstandard6th12(&phi[index]);
        PDstandard6th13phi = PDstandard6th13(&phi[index]);
        PDstandard6th23phi = PDstandard6th23(&phi[index]);
        break;
      }
      
      case 4:
      {
        PDstandard6th1A11 = PDstandard6th1(&A11[index]);
        PDstandard6th2A11 = PDstandard6th2(&A11[index]);
        PDstandard6th3A11 = PDstandard6th3(&A11[index]);
        PDstandard6th1A21 = PDstandard6th1(&A21[index]);
        PDstandard6th2A21 = PDstandard6th2(&A21[index]);
        PDstandard6th3A21 = PDstandard6th3(&A21[index]);
        PDstandard6th1A22 = PDstandard6th1(&A22[index]);
        PDstandard6th2A22 = PDstandard6th2(&A22[index]);
        PDstandard6th3A22 = PDstandard6th3(&A22[index]);
        PDstandard6th1A31 = PDstandard6th1(&A31[index]);
        PDstandard6th2A31 = PDstandard6th2(&A31[index]);
        PDstandard6th3A31 = PDstandard6th3(&A31[index]);
        PDstandard6th1A32 = PDstandard6th1(&A32[index]);
        PDstandard6th2A32 = PDstandard6th2(&A32[index]);
        PDstandard6th3A32 = PDstandard6th3(&A32[index]);
        PDstandard6th1A33 = PDstandard6th1(&A33[index]);
        PDstandard6th2A33 = PDstandard6th2(&A33[index]);
        PDstandard6th3A33 = PDstandard6th3(&A33[index]);
        PDstandard6th1alpha = PDstandard6th1(&alpha[index]);
        PDstandard6th2alpha = PDstandard6th2(&alpha[index]);
        PDstandard6th3alpha = PDstandard6th3(&alpha[index]);
        PDstandard6th11alpha = PDstandard6th11(&alpha[index]);
        PDstandard6th22alpha = PDstandard6th22(&alpha[index]);
        PDstandard6th33alpha = PDstandard6th33(&alpha[index]);
        PDstandard6th12alpha = PDstandard6th12(&alpha[index]);
        PDstandard6th13alpha = PDstandard6th13(&alpha[index]);
        PDstandard6th23alpha = PDstandard6th23(&alpha[index]);
        PDstandard6th1beta1 = PDstandard6th1(&beta1[index]);
        PDstandard6th2beta1 = PDstandard6th2(&beta1[index]);
        PDstandard6th3beta1 = PDstandard6th3(&beta1[index]);
        PDstandard6th11beta1 = PDstandard6th11(&beta1[index]);
        PDstandard6th22beta1 = PDstandard6th22(&beta1[index]);
        PDstandard6th33beta1 = PDstandard6th33(&beta1[index]);
        PDstandard6th12beta1 = PDstandard6th12(&beta1[index]);
        PDstandard6th13beta1 = PDstandard6th13(&beta1[index]);
        PDstandard6th23beta1 = PDstandard6th23(&beta1[index]);
        PDstandard6th1beta2 = PDstandard6th1(&beta2[index]);
        PDstandard6th2beta2 = PDstandard6th2(&beta2[index]);
        PDstandard6th3beta2 = PDstandard6th3(&beta2[index]);
        PDstandard6th11beta2 = PDstandard6th11(&beta2[index]);
        PDstandard6th22beta2 = PDstandard6th22(&beta2[index]);
        PDstandard6th33beta2 = PDstandard6th33(&beta2[index]);
        PDstandard6th12beta2 = PDstandard6th12(&beta2[index]);
        PDstandard6th13beta2 = PDstandard6th13(&beta2[index]);
        PDstandard6th23beta2 = PDstandard6th23(&beta2[index]);
        PDstandard6th1beta3 = PDstandard6th1(&beta3[index]);
        PDstandard6th2beta3 = PDstandard6th2(&beta3[index]);
        PDstandard6th3beta3 = PDstandard6th3(&beta3[index]);
        PDstandard6th11beta3 = PDstandard6th11(&beta3[index]);
        PDstandard6th22beta3 = PDstandard6th22(&beta3[index]);
        PDstandard6th33beta3 = PDstandard6th33(&beta3[index]);
        PDstandard6th12beta3 = PDstandard6th12(&beta3[index]);
        PDstandard6th13beta3 = PDstandard6th13(&beta3[index]);
        PDstandard6th23beta3 = PDstandard6th23(&beta3[index]);
        PDstandard6th1Gam1 = PDstandard6th1(&Gam1[index]);
        PDstandard6th2Gam1 = PDstandard6th2(&Gam1[index]);
        PDstandard6th3Gam1 = PDstandard6th3(&Gam1[index]);
        PDstandard6th1Gam2 = PDstandard6th1(&Gam2[index]);
        PDstandard6th2Gam2 = PDstandard6th2(&Gam2[index]);
        PDstandard6th3Gam2 = PDstandard6th3(&Gam2[index]);
        PDstandard6th1Gam3 = PDstandard6th1(&Gam3[index]);
        PDstandard6th2Gam3 = PDstandard6th2(&Gam3[index]);
        PDstandard6th3Gam3 = PDstandard6th3(&Gam3[index]);
        PDstandard6th1h11 = PDstandard6th1(&h11[index]);
        PDstandard6th2h11 = PDstandard6th2(&h11[index]);
        PDstandard6th3h11 = PDstandard6th3(&h11[index]);
        PDstandard6th11h11 = PDstandard6th11(&h11[index]);
        PDstandard6th22h11 = PDstandard6th22(&h11[index]);
        PDstandard6th33h11 = PDstandard6th33(&h11[index]);
        PDstandard6th12h11 = PDstandard6th12(&h11[index]);
        PDstandard6th13h11 = PDstandard6th13(&h11[index]);
        PDstandard6th23h11 = PDstandard6th23(&h11[index]);
        PDstandard6th1h21 = PDstandard6th1(&h21[index]);
        PDstandard6th2h21 = PDstandard6th2(&h21[index]);
        PDstandard6th3h21 = PDstandard6th3(&h21[index]);
        PDstandard6th11h21 = PDstandard6th11(&h21[index]);
        PDstandard6th22h21 = PDstandard6th22(&h21[index]);
        PDstandard6th33h21 = PDstandard6th33(&h21[index]);
        PDstandard6th12h21 = PDstandard6th12(&h21[index]);
        PDstandard6th13h21 = PDstandard6th13(&h21[index]);
        PDstandard6th23h21 = PDstandard6th23(&h21[index]);
        PDstandard6th1h22 = PDstandard6th1(&h22[index]);
        PDstandard6th2h22 = PDstandard6th2(&h22[index]);
        PDstandard6th3h22 = PDstandard6th3(&h22[index]);
        PDstandard6th11h22 = PDstandard6th11(&h22[index]);
        PDstandard6th22h22 = PDstandard6th22(&h22[index]);
        PDstandard6th33h22 = PDstandard6th33(&h22[index]);
        PDstandard6th12h22 = PDstandard6th12(&h22[index]);
        PDstandard6th13h22 = PDstandard6th13(&h22[index]);
        PDstandard6th23h22 = PDstandard6th23(&h22[index]);
        PDstandard6th1h31 = PDstandard6th1(&h31[index]);
        PDstandard6th2h31 = PDstandard6th2(&h31[index]);
        PDstandard6th3h31 = PDstandard6th3(&h31[index]);
        PDstandard6th11h31 = PDstandard6th11(&h31[index]);
        PDstandard6th22h31 = PDstandard6th22(&h31[index]);
        PDstandard6th33h31 = PDstandard6th33(&h31[index]);
        PDstandard6th12h31 = PDstandard6th12(&h31[index]);
        PDstandard6th13h31 = PDstandard6th13(&h31[index]);
        PDstandard6th23h31 = PDstandard6th23(&h31[index]);
        PDstandard6th1h32 = PDstandard6th1(&h32[index]);
        PDstandard6th2h32 = PDstandard6th2(&h32[index]);
        PDstandard6th3h32 = PDstandard6th3(&h32[index]);
        PDstandard6th11h32 = PDstandard6th11(&h32[index]);
        PDstandard6th22h32 = PDstandard6th22(&h32[index]);
        PDstandard6th33h32 = PDstandard6th33(&h32[index]);
        PDstandard6th12h32 = PDstandard6th12(&h32[index]);
        PDstandard6th13h32 = PDstandard6th13(&h32[index]);
        PDstandard6th23h32 = PDstandard6th23(&h32[index]);
        PDstandard6th1h33 = PDstandard6th1(&h33[index]);
        PDstandard6th2h33 = PDstandard6th2(&h33[index]);
        PDstandard6th3h33 = PDstandard6th3(&h33[index]);
        PDstandard6th11h33 = PDstandard6th11(&h33[index]);
        PDstandard6th22h33 = PDstandard6th22(&h33[index]);
        PDstandard6th33h33 = PDstandard6th33(&h33[index]);
        PDstandard6th12h33 = PDstandard6th12(&h33[index]);
        PDstandard6th13h33 = PDstandard6th13(&h33[index]);
        PDstandard6th23h33 = PDstandard6th23(&h33[index]);
        PDstandard6th1K = PDstandard6th1(&K[index]);
        PDstandard6th2K = PDstandard6th2(&K[index]);
        PDstandard6th3K = PDstandard6th3(&K[index]);
        PDstandard6th1phi = PDstandard6th1(&phi[index]);
        PDstandard6th2phi = PDstandard6th2(&phi[index]);
        PDstandard6th3phi = PDstandard6th3(&phi[index]);
        PDstandard6th11phi = PDstandard6th11(&phi[index]);
        PDstandard6th22phi = PDstandard6th22(&phi[index]);
        PDstandard6th33phi = PDstandard6th33(&phi[index]);
        PDstandard6th12phi = PDstandard6th12(&phi[index]);
        PDstandard6th13phi = PDstandard6th13(&phi[index]);
        PDstandard6th23phi = PDstandard6th23(&phi[index]);
        break;
      }
      
      case 6:
      {
        PDstandard6th1A11 = PDstandard6th1(&A11[index]);
        PDstandard6th2A11 = PDstandard6th2(&A11[index]);
        PDstandard6th3A11 = PDstandard6th3(&A11[index]);
        PDstandard6th1A21 = PDstandard6th1(&A21[index]);
        PDstandard6th2A21 = PDstandard6th2(&A21[index]);
        PDstandard6th3A21 = PDstandard6th3(&A21[index]);
        PDstandard6th1A22 = PDstandard6th1(&A22[index]);
        PDstandard6th2A22 = PDstandard6th2(&A22[index]);
        PDstandard6th3A22 = PDstandard6th3(&A22[index]);
        PDstandard6th1A31 = PDstandard6th1(&A31[index]);
        PDstandard6th2A31 = PDstandard6th2(&A31[index]);
        PDstandard6th3A31 = PDstandard6th3(&A31[index]);
        PDstandard6th1A32 = PDstandard6th1(&A32[index]);
        PDstandard6th2A32 = PDstandard6th2(&A32[index]);
        PDstandard6th3A32 = PDstandard6th3(&A32[index]);
        PDstandard6th1A33 = PDstandard6th1(&A33[index]);
        PDstandard6th2A33 = PDstandard6th2(&A33[index]);
        PDstandard6th3A33 = PDstandard6th3(&A33[index]);
        PDstandard6th1alpha = PDstandard6th1(&alpha[index]);
        PDstandard6th2alpha = PDstandard6th2(&alpha[index]);
        PDstandard6th3alpha = PDstandard6th3(&alpha[index]);
        PDstandard6th11alpha = PDstandard6th11(&alpha[index]);
        PDstandard6th22alpha = PDstandard6th22(&alpha[index]);
        PDstandard6th33alpha = PDstandard6th33(&alpha[index]);
        PDstandard6th12alpha = PDstandard6th12(&alpha[index]);
        PDstandard6th13alpha = PDstandard6th13(&alpha[index]);
        PDstandard6th23alpha = PDstandard6th23(&alpha[index]);
        PDstandard6th1beta1 = PDstandard6th1(&beta1[index]);
        PDstandard6th2beta1 = PDstandard6th2(&beta1[index]);
        PDstandard6th3beta1 = PDstandard6th3(&beta1[index]);
        PDstandard6th11beta1 = PDstandard6th11(&beta1[index]);
        PDstandard6th22beta1 = PDstandard6th22(&beta1[index]);
        PDstandard6th33beta1 = PDstandard6th33(&beta1[index]);
        PDstandard6th12beta1 = PDstandard6th12(&beta1[index]);
        PDstandard6th13beta1 = PDstandard6th13(&beta1[index]);
        PDstandard6th23beta1 = PDstandard6th23(&beta1[index]);
        PDstandard6th1beta2 = PDstandard6th1(&beta2[index]);
        PDstandard6th2beta2 = PDstandard6th2(&beta2[index]);
        PDstandard6th3beta2 = PDstandard6th3(&beta2[index]);
        PDstandard6th11beta2 = PDstandard6th11(&beta2[index]);
        PDstandard6th22beta2 = PDstandard6th22(&beta2[index]);
        PDstandard6th33beta2 = PDstandard6th33(&beta2[index]);
        PDstandard6th12beta2 = PDstandard6th12(&beta2[index]);
        PDstandard6th13beta2 = PDstandard6th13(&beta2[index]);
        PDstandard6th23beta2 = PDstandard6th23(&beta2[index]);
        PDstandard6th1beta3 = PDstandard6th1(&beta3[index]);
        PDstandard6th2beta3 = PDstandard6th2(&beta3[index]);
        PDstandard6th3beta3 = PDstandard6th3(&beta3[index]);
        PDstandard6th11beta3 = PDstandard6th11(&beta3[index]);
        PDstandard6th22beta3 = PDstandard6th22(&beta3[index]);
        PDstandard6th33beta3 = PDstandard6th33(&beta3[index]);
        PDstandard6th12beta3 = PDstandard6th12(&beta3[index]);
        PDstandard6th13beta3 = PDstandard6th13(&beta3[index]);
        PDstandard6th23beta3 = PDstandard6th23(&beta3[index]);
        PDstandard6th1Gam1 = PDstandard6th1(&Gam1[index]);
        PDstandard6th2Gam1 = PDstandard6th2(&Gam1[index]);
        PDstandard6th3Gam1 = PDstandard6th3(&Gam1[index]);
        PDstandard6th1Gam2 = PDstandard6th1(&Gam2[index]);
        PDstandard6th2Gam2 = PDstandard6th2(&Gam2[index]);
        PDstandard6th3Gam2 = PDstandard6th3(&Gam2[index]);
        PDstandard6th1Gam3 = PDstandard6th1(&Gam3[index]);
        PDstandard6th2Gam3 = PDstandard6th2(&Gam3[index]);
        PDstandard6th3Gam3 = PDstandard6th3(&Gam3[index]);
        PDstandard6th1h11 = PDstandard6th1(&h11[index]);
        PDstandard6th2h11 = PDstandard6th2(&h11[index]);
        PDstandard6th3h11 = PDstandard6th3(&h11[index]);
        PDstandard6th11h11 = PDstandard6th11(&h11[index]);
        PDstandard6th22h11 = PDstandard6th22(&h11[index]);
        PDstandard6th33h11 = PDstandard6th33(&h11[index]);
        PDstandard6th12h11 = PDstandard6th12(&h11[index]);
        PDstandard6th13h11 = PDstandard6th13(&h11[index]);
        PDstandard6th23h11 = PDstandard6th23(&h11[index]);
        PDstandard6th1h21 = PDstandard6th1(&h21[index]);
        PDstandard6th2h21 = PDstandard6th2(&h21[index]);
        PDstandard6th3h21 = PDstandard6th3(&h21[index]);
        PDstandard6th11h21 = PDstandard6th11(&h21[index]);
        PDstandard6th22h21 = PDstandard6th22(&h21[index]);
        PDstandard6th33h21 = PDstandard6th33(&h21[index]);
        PDstandard6th12h21 = PDstandard6th12(&h21[index]);
        PDstandard6th13h21 = PDstandard6th13(&h21[index]);
        PDstandard6th23h21 = PDstandard6th23(&h21[index]);
        PDstandard6th1h22 = PDstandard6th1(&h22[index]);
        PDstandard6th2h22 = PDstandard6th2(&h22[index]);
        PDstandard6th3h22 = PDstandard6th3(&h22[index]);
        PDstandard6th11h22 = PDstandard6th11(&h22[index]);
        PDstandard6th22h22 = PDstandard6th22(&h22[index]);
        PDstandard6th33h22 = PDstandard6th33(&h22[index]);
        PDstandard6th12h22 = PDstandard6th12(&h22[index]);
        PDstandard6th13h22 = PDstandard6th13(&h22[index]);
        PDstandard6th23h22 = PDstandard6th23(&h22[index]);
        PDstandard6th1h31 = PDstandard6th1(&h31[index]);
        PDstandard6th2h31 = PDstandard6th2(&h31[index]);
        PDstandard6th3h31 = PDstandard6th3(&h31[index]);
        PDstandard6th11h31 = PDstandard6th11(&h31[index]);
        PDstandard6th22h31 = PDstandard6th22(&h31[index]);
        PDstandard6th33h31 = PDstandard6th33(&h31[index]);
        PDstandard6th12h31 = PDstandard6th12(&h31[index]);
        PDstandard6th13h31 = PDstandard6th13(&h31[index]);
        PDstandard6th23h31 = PDstandard6th23(&h31[index]);
        PDstandard6th1h32 = PDstandard6th1(&h32[index]);
        PDstandard6th2h32 = PDstandard6th2(&h32[index]);
        PDstandard6th3h32 = PDstandard6th3(&h32[index]);
        PDstandard6th11h32 = PDstandard6th11(&h32[index]);
        PDstandard6th22h32 = PDstandard6th22(&h32[index]);
        PDstandard6th33h32 = PDstandard6th33(&h32[index]);
        PDstandard6th12h32 = PDstandard6th12(&h32[index]);
        PDstandard6th13h32 = PDstandard6th13(&h32[index]);
        PDstandard6th23h32 = PDstandard6th23(&h32[index]);
        PDstandard6th1h33 = PDstandard6th1(&h33[index]);
        PDstandard6th2h33 = PDstandard6th2(&h33[index]);
        PDstandard6th3h33 = PDstandard6th3(&h33[index]);
        PDstandard6th11h33 = PDstandard6th11(&h33[index]);
        PDstandard6th22h33 = PDstandard6th22(&h33[index]);
        PDstandard6th33h33 = PDstandard6th33(&h33[index]);
        PDstandard6th12h33 = PDstandard6th12(&h33[index]);
        PDstandard6th13h33 = PDstandard6th13(&h33[index]);
        PDstandard6th23h33 = PDstandard6th23(&h33[index]);
        PDstandard6th1K = PDstandard6th1(&K[index]);
        PDstandard6th2K = PDstandard6th2(&K[index]);
        PDstandard6th3K = PDstandard6th3(&K[index]);
        PDstandard6th1phi = PDstandard6th1(&phi[index]);
        PDstandard6th2phi = PDstandard6th2(&phi[index]);
        PDstandard6th3phi = PDstandard6th3(&phi[index]);
        PDstandard6th11phi = PDstandard6th11(&phi[index]);
        PDstandard6th22phi = PDstandard6th22(&phi[index]);
        PDstandard6th33phi = PDstandard6th33(&phi[index]);
        PDstandard6th12phi = PDstandard6th12(&phi[index]);
        PDstandard6th13phi = PDstandard6th13(&phi[index]);
        PDstandard6th23phi = PDstandard6th23(&phi[index]);
        break;
      }
      default:
        CCTK_BUILTIN_UNREACHABLE();
    }
    /* Calculate temporaries and grid functions */
    CCTK_REAL deth CCTK_ATTRIBUTE_UNUSED = 2*h21L*h31L*h32L - 
      h33L*pow(h21L,2) + h22L*(h11L*h33L - pow(h31L,2)) - h11L*pow(h32L,2);
    
    CCTK_REAL invdeth CCTK_ATTRIBUTE_UNUSED = pow(deth,-1);
    
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
      0.5*hInv11*PDstandard6th1h11 + hInv12*(PDstandard6th1h21 - 
      0.5*PDstandard6th2h11) + hInv13*(PDstandard6th1h31 - 
      0.5*PDstandard6th3h11);
    
    CCTK_REAL gamma211 CCTK_ATTRIBUTE_UNUSED = 
      0.5*hInv12*PDstandard6th1h11 + hInv22*(PDstandard6th1h21 - 
      0.5*PDstandard6th2h11) + hInv23*(PDstandard6th1h31 - 
      0.5*PDstandard6th3h11);
    
    CCTK_REAL gamma311 CCTK_ATTRIBUTE_UNUSED = 
      0.5*hInv13*PDstandard6th1h11 + hInv23*(PDstandard6th1h21 - 
      0.5*PDstandard6th2h11) + hInv33*(PDstandard6th1h31 - 
      0.5*PDstandard6th3h11);
    
    CCTK_REAL gamma121 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv12*PDstandard6th1h22 + hInv11*PDstandard6th2h11 + 
      hInv13*(PDstandard6th1h32 + PDstandard6th2h31 - PDstandard6th3h21));
    
    CCTK_REAL gamma221 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv22*PDstandard6th1h22 + hInv12*PDstandard6th2h11 + 
      hInv23*(PDstandard6th1h32 + PDstandard6th2h31 - PDstandard6th3h21));
    
    CCTK_REAL gamma321 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv23*PDstandard6th1h22 + hInv13*PDstandard6th2h11 + 
      hInv33*(PDstandard6th1h32 + PDstandard6th2h31 - PDstandard6th3h21));
    
    CCTK_REAL gamma131 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv13*PDstandard6th1h33 + hInv11*PDstandard6th3h11 + 
      hInv12*(PDstandard6th1h32 - PDstandard6th2h31 + PDstandard6th3h21));
    
    CCTK_REAL gamma231 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv23*PDstandard6th1h33 + hInv12*PDstandard6th3h11 + 
      hInv22*(PDstandard6th1h32 - PDstandard6th2h31 + PDstandard6th3h21));
    
    CCTK_REAL gamma331 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv33*PDstandard6th1h33 + hInv13*PDstandard6th3h11 + 
      hInv23*(PDstandard6th1h32 - PDstandard6th2h31 + PDstandard6th3h21));
    
    CCTK_REAL gamma122 CCTK_ATTRIBUTE_UNUSED = 
      hInv11*(-0.5*PDstandard6th1h22 + PDstandard6th2h21) + 
      0.5*hInv12*PDstandard6th2h22 + hInv13*(PDstandard6th2h32 - 
      0.5*PDstandard6th3h22);
    
    CCTK_REAL gamma222 CCTK_ATTRIBUTE_UNUSED = 
      hInv12*(-0.5*PDstandard6th1h22 + PDstandard6th2h21) + 
      0.5*hInv22*PDstandard6th2h22 + hInv23*(PDstandard6th2h32 - 
      0.5*PDstandard6th3h22);
    
    CCTK_REAL gamma322 CCTK_ATTRIBUTE_UNUSED = 
      hInv13*(-0.5*PDstandard6th1h22 + PDstandard6th2h21) + 
      0.5*hInv23*PDstandard6th2h22 + hInv33*(PDstandard6th2h32 - 
      0.5*PDstandard6th3h22);
    
    CCTK_REAL gamma132 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv13*PDstandard6th2h33 + hInv11*(-PDstandard6th1h32 + 
      PDstandard6th2h31 + PDstandard6th3h21) + hInv12*PDstandard6th3h22);
    
    CCTK_REAL gamma232 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv23*PDstandard6th2h33 + hInv12*(-PDstandard6th1h32 + 
      PDstandard6th2h31 + PDstandard6th3h21) + hInv22*PDstandard6th3h22);
    
    CCTK_REAL gamma332 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv33*PDstandard6th2h33 + hInv13*(-PDstandard6th1h32 + 
      PDstandard6th2h31 + PDstandard6th3h21) + hInv23*PDstandard6th3h22);
    
    CCTK_REAL gamma133 CCTK_ATTRIBUTE_UNUSED = 
      hInv11*(-0.5*PDstandard6th1h33 + PDstandard6th3h31) + 
      hInv12*(-0.5*PDstandard6th2h33 + PDstandard6th3h32) + 
      0.5*hInv13*PDstandard6th3h33;
    
    CCTK_REAL gamma233 CCTK_ATTRIBUTE_UNUSED = 
      hInv12*(-0.5*PDstandard6th1h33 + PDstandard6th3h31) + 
      hInv22*(-0.5*PDstandard6th2h33 + PDstandard6th3h32) + 
      0.5*hInv23*PDstandard6th3h33;
    
    CCTK_REAL gamma333 CCTK_ATTRIBUTE_UNUSED = 
      hInv13*(-0.5*PDstandard6th1h33 + PDstandard6th3h31) + 
      hInv23*(-0.5*PDstandard6th2h33 + PDstandard6th3h32) + 
      0.5*hInv33*PDstandard6th3h33;
    
    CCTK_REAL AInv11 CCTK_ATTRIBUTE_UNUSED = 2*(A32L*hInv12*hInv13 + 
      hInv11*(A21L*hInv12 + A31L*hInv13)) + A11L*pow(hInv11,2) + 
      A22L*pow(hInv12,2) + A33L*pow(hInv13,2);
    
    CCTK_REAL AInv12 CCTK_ATTRIBUTE_UNUSED = hInv11*(A11L*hInv12 + 
      A21L*hInv22 + A31L*hInv23) + hInv12*(A31L*hInv13 + A22L*hInv22 + 
      A32L*hInv23) + hInv13*(A32L*hInv22 + A33L*hInv23) + A21L*pow(hInv12,2);
    
    CCTK_REAL AInv13 CCTK_ATTRIBUTE_UNUSED = hInv11*(A11L*hInv13 + 
      A21L*hInv23 + A31L*hInv33) + hInv12*(A21L*hInv13 + A22L*hInv23 + 
      A32L*hInv33) + hInv13*(A32L*hInv23 + A33L*hInv33) + A31L*pow(hInv13,2);
    
    CCTK_REAL AInv22 CCTK_ATTRIBUTE_UNUSED = 2*(A32L*hInv22*hInv23 + 
      hInv12*(A21L*hInv22 + A31L*hInv23)) + A11L*pow(hInv12,2) + 
      A22L*pow(hInv22,2) + A33L*pow(hInv23,2);
    
    CCTK_REAL AInv23 CCTK_ATTRIBUTE_UNUSED = hInv13*(A21L*hInv22 + 
      A31L*hInv23) + A33L*hInv23*hInv33 + hInv12*(A11L*hInv13 + A21L*hInv23 + 
      A31L*hInv33) + hInv22*(A22L*hInv23 + A32L*hInv33) + A32L*pow(hInv23,2);
    
    CCTK_REAL AInv33 CCTK_ATTRIBUTE_UNUSED = 2*(A32L*hInv23*hInv33 + 
      hInv13*(A21L*hInv23 + A31L*hInv33)) + A11L*pow(hInv13,2) + 
      A22L*pow(hInv23,2) + A33L*pow(hInv33,2);
    
    CCTK_REAL em4phi CCTK_ATTRIBUTE_UNUSED = exp(-4*phiL);
    
    CCTK_REAL R11 CCTK_ATTRIBUTE_UNUSED = h11L*(Gam1L*gamma111 + 
      Gam2L*gamma121 + Gam3L*gamma131 + PDstandard6th1Gam1) + 
      h21L*(Gam1L*gamma211 + Gam2L*gamma221 + Gam3L*gamma231 + 
      PDstandard6th1Gam2) + h31L*(Gam1L*gamma311 + Gam2L*gamma321 + 
      Gam3L*gamma331 + PDstandard6th1Gam3) + 
      hInv23*(2*(gamma231*(h11L*gamma122 + h21L*gamma222 + h32L*gamma321 + 
      h31L*gamma322) + h11L*gamma132*(gamma221 + gamma331) + 
      gamma121*(3*h11L*gamma131 + 2*(h21L*gamma231 + h31L*gamma331)) + 
      gamma331*(h21L*gamma232 + h31L*gamma332) + gamma221*(h22L*gamma231 + 
      h21L*(2*gamma131 + gamma232) + h32L*gamma331 + h31L*gamma332) + 
      gamma321*(h11L*gamma133 + h21L*gamma233 + h33L*gamma331 + 
      h31L*(2*gamma131 + gamma333))) - PDstandard6th23h11) + 
      hInv12*(-PDstandard6th12h11 + 2*(gamma121*(3*h11L*gamma111 + 
      2*(h21L*gamma211 + h31L*gamma311)) + gamma221*(2*h21L*gamma111 + 
      h22L*gamma211 + h32L*gamma311 + h31L*gamma321) + 
      h11L*(gamma122*gamma211 + gamma121*gamma221 + gamma132*gamma311 + 
      gamma131*gamma321) + gamma321*(h32L*gamma211 + h33L*gamma311 + 
      h31L*(2*gamma111 + gamma331)) + h31L*(gamma211*gamma322 + 
      gamma311*gamma332) + h21L*(gamma211*gamma222 + gamma232*gamma311 + 
      gamma231*gamma321 + pow(gamma221,2)))) + 
      hInv11*(4*gamma111*(h21L*gamma211 + h31L*gamma311) + 
      2*(gamma211*(h21L*gamma221 + h32L*gamma311 + h31L*gamma321) + 
      gamma311*(h21L*gamma231 + h31L*gamma331)) - 0.5*PDstandard6th11h11 + 
      h11L*(2*(gamma121*gamma211 + gamma131*gamma311) + 3*pow(gamma111,2)) + 
      h22L*pow(gamma211,2) + h33L*pow(gamma311,2)) + 
      hInv22*(4*gamma121*(h21L*gamma221 + h31L*gamma321) + 
      2*(gamma221*(h21L*gamma222 + h32L*gamma321 + h31L*gamma322) + 
      gamma321*(h21L*gamma232 + h31L*gamma332)) - 0.5*PDstandard6th22h11 + 
      h11L*(2*(gamma122*gamma221 + gamma132*gamma321) + 3*pow(gamma121,2)) + 
      h22L*pow(gamma221,2) + h33L*pow(gamma321,2)) + 
      hInv33*(4*gamma131*(h21L*gamma231 + h31L*gamma331) + 
      2*(gamma231*(h21L*gamma232 + h32L*gamma331 + h31L*gamma332) + 
      gamma331*(h21L*gamma233 + h31L*gamma333)) - 0.5*PDstandard6th33h11 + 
      h11L*(2*(gamma132*gamma231 + gamma133*gamma331) + 3*pow(gamma131,2)) + 
      h22L*pow(gamma231,2) + h33L*pow(gamma331,2)) + 
      hInv13*(-PDstandard6th13h11 + 2*(h11L*(gamma132*gamma211 + 
      gamma121*gamma231 + gamma133*gamma311) + gamma231*(h22L*gamma211 + 
      h21L*(2*gamma111 + gamma221) + h32L*gamma311 + h31L*gamma321) + 
      (2*h31L*gamma111 + h32L*gamma211 + h33L*gamma311)*gamma331 + 
      h21L*(gamma211*gamma232 + gamma233*gamma311 + gamma231*gamma331) + 
      gamma131*(2*(h21L*gamma211 + h31L*gamma311) + h11L*(3*gamma111 + 
      gamma331)) + h31L*(gamma211*gamma332 + gamma311*gamma333 + 
      pow(gamma331,2))));
    
    CCTK_REAL R21 CCTK_ATTRIBUTE_UNUSED = 0.5*(Gam1L*(h11L*gamma121 + 
      h22L*gamma211 + h21L*(gamma111 + gamma221) + h32L*gamma311 + 
      h31L*gamma321) + Gam2L*(h11L*gamma122 + h22L*gamma221 + h21L*(gamma121 
      + gamma222) + h32L*gamma321 + h31L*gamma322) + Gam3L*(h11L*gamma132 + 
      h22L*gamma231 + h21L*(gamma131 + gamma232) + h32L*gamma331 + 
      h31L*gamma332) + h22L*PDstandard6th1Gam2 + h32L*PDstandard6th1Gam3 + 
      h11L*PDstandard6th2Gam1 + h21L*(PDstandard6th1Gam1 + 
      PDstandard6th2Gam2) + h31L*PDstandard6th2Gam3) + 
      hInv11*(gamma221*(h11L*gamma121 + 2*h22L*gamma211 + h32L*gamma311 + 
      h31L*gamma321) + gamma111*(2*h11L*gamma121 + h22L*gamma211 + 
      h21L*gamma221 + h32L*gamma311 + h31L*gamma321) + 
      gamma321*(h11L*gamma131 + 2*h32L*gamma211 + h31L*gamma331) + 
      gamma311*(2*h31L*gamma121 + h22L*gamma231 + h33L*gamma321 + 
      h32L*gamma331) - 0.5*PDstandard6th11h21 + h21L*(3*gamma121*gamma211 + 
      gamma131*gamma311 + gamma231*gamma321 + pow(gamma111,2) + 
      pow(gamma221,2))) + hInv22*(gamma222*(h11L*gamma122 + 2*h22L*gamma221 + 
      h32L*gamma321 + h31L*gamma322) + gamma121*(2*h11L*gamma122 + 
      h22L*gamma221 + h21L*gamma222 + h32L*gamma321 + h31L*gamma322) + 
      gamma322*(h11L*gamma132 + 2*h32L*gamma221 + h31L*gamma332) + 
      gamma321*(2*h31L*gamma122 + h22L*gamma232 + h33L*gamma322 + 
      h32L*gamma332) - 0.5*PDstandard6th22h21 + h21L*(3*gamma122*gamma221 + 
      gamma132*gamma321 + gamma232*gamma322 + pow(gamma121,2) + 
      pow(gamma222,2))) + hInv33*(gamma232*(h11L*gamma132 + 2*h22L*gamma231 + 
      h32L*gamma331 + h31L*gamma332) + gamma131*(2*h11L*gamma132 + 
      h22L*gamma231 + h21L*gamma232 + h32L*gamma331 + h31L*gamma332) + 
      gamma332*(h11L*gamma133 + 2*h32L*gamma231 + h31L*gamma333) + 
      gamma331*(2*h31L*gamma132 + h22L*gamma233 + h33L*gamma332 + 
      h32L*gamma333) - 0.5*PDstandard6th33h21 + h21L*(3*gamma132*gamma231 + 
      gamma133*gamma331 + gamma233*gamma332 + pow(gamma131,2) + 
      pow(gamma232,2))) + hInv12*(gamma222*(h21L*gamma111 + h11L*gamma121 + 
      h32L*gamma311 + h31L*gamma321) + gamma221*(h22L*gamma111 + 
      h11L*gamma122 + h21L*(4*gamma121 + 2*gamma222) + 3*h32L*gamma321 + 
      h31L*gamma322) + gamma211*(3*h21L*gamma122 + h22L*(gamma121 + 
      2*gamma222) + 2*h32L*gamma322) + gamma322*(h11L*gamma131 + 
      h21L*gamma231 + h31L*(gamma111 + gamma331)) + gamma311*(h21L*gamma132 + 
      h22L*gamma232 + h33L*gamma322 + h32L*(gamma121 + gamma332)) + 
      gamma321*(h11L*gamma132 + h22L*gamma231 + h21L*(gamma131 + gamma232) + 
      h32L*(gamma111 + gamma331) + h31L*(3*gamma121 + gamma332)) - 
      PDstandard6th12h21 + 2*(gamma111*(h21L*gamma121 + h11L*gamma122) + 
      h31L*gamma122*gamma311 + h11L*pow(gamma121,2) + h22L*pow(gamma221,2)) + 
      h33L*pow(gamma321,2)) + hInv13*(gamma232*(h21L*gamma111 + h11L*gamma121 
      + h32L*gamma311 + h31L*gamma321) + 2*((h21L*gamma111 + 
      h11L*gamma121)*gamma131 + h22L*gamma221*gamma231 + 
      gamma132*(h11L*gamma111 + h31L*gamma311) + h31L*gamma121*gamma331) + 
      (h31L*gamma111 + h11L*gamma131)*gamma332 + gamma331*(h32L*gamma111 + 
      h21L*gamma131 + h31L*gamma332) + gamma221*(h11L*gamma132 + 
      h21L*(gamma131 + 2*gamma232) + h32L*gamma331 + h31L*gamma332) + 
      gamma211*(3*h21L*gamma132 + h22L*(gamma131 + 2*gamma232) + 
      2*h32L*gamma332) + gamma311*(h32L*gamma131 + h21L*gamma133 + 
      h22L*gamma233 + h33L*gamma332) + gamma231*(2*h32L*gamma321 + 
      h22L*(gamma111 + gamma331) + h21L*(3*gamma121 + gamma332)) + 
      gamma321*(h11L*gamma133 + h21L*gamma233 + h33L*gamma331 + 
      h31L*(gamma131 + gamma333)) - PDstandard6th13h21 + 
      h32L*(gamma311*gamma333 + pow(gamma331,2))) + 
      hInv23*(gamma231*(h22L*gamma121 + 3*h21L*gamma122 + 2*h32L*gamma322) + 
      2*((h21L*gamma121 + h11L*gamma122)*gamma131 + h22L*gamma222*gamma231 + 
      gamma132*(h11L*gamma121 + h31L*gamma321) + h31L*gamma122*gamma331) + 
      (h31L*gamma121 + h11L*gamma132)*gamma332 + gamma222*(h11L*gamma132 + 
      h21L*(gamma131 + 2*gamma232) + h32L*gamma331 + h31L*gamma332) + 
      gamma221*(3*h21L*gamma132 + h22L*(gamma131 + 2*gamma232) + 
      2*h32L*gamma332) + gamma232*(h11L*gamma122 + h32L*gamma321 + 
      h31L*gamma322 + h22L*gamma331 + h21L*(gamma121 + gamma332)) + 
      gamma331*(h21L*gamma132 + h32L*(gamma121 + gamma332)) + 
      gamma322*(h11L*gamma133 + h21L*gamma233 + h33L*gamma331 + 
      h31L*(gamma131 + gamma333)) + gamma321*(h21L*gamma133 + h22L*gamma233 + 
      h33L*gamma332 + h32L*(gamma131 + gamma333)) - PDstandard6th23h21 + 
      h31L*pow(gamma332,2));
    
    CCTK_REAL R31 CCTK_ATTRIBUTE_UNUSED = 0.5*(h32L*(Gam1L*gamma211 + 
      Gam2L*gamma221 + Gam3L*gamma231 + PDstandard6th1Gam2) + 
      h33L*(Gam1L*gamma311 + Gam2L*gamma321 + Gam3L*gamma331 + 
      PDstandard6th1Gam3) + h11L*(Gam1L*gamma131 + Gam2L*gamma132 + 
      Gam3L*gamma133 + PDstandard6th3Gam1) + h21L*(Gam1L*gamma231 + 
      Gam2L*gamma232 + Gam3L*gamma233 + PDstandard6th3Gam2) + 
      h31L*(Gam1L*(gamma111 + gamma331) + Gam2L*(gamma121 + gamma332) + 
      Gam3L*(gamma131 + gamma333) + PDstandard6th1Gam1 + PDstandard6th3Gam3)) 
      + hInv12*(gamma221*(h32L*gamma111 + h31L*gamma121 + h22L*gamma231 + 
      h21L*gamma232 + h33L*gamma321) + (h31L*gamma121 + 
      h11L*gamma132)*gamma331 + gamma232*(h11L*gamma121 + h31L*gamma321 + 
      h21L*(gamma111 + gamma331)) + gamma321*(3*h31L*gamma131 + 
      h33L*(gamma111 + 2*gamma331)) + (h31L*gamma111 + 
      h11L*gamma131)*gamma332 + 2*(gamma121*(h31L*gamma111 + h11L*gamma131) + 
      gamma132*(h11L*gamma111 + h21L*gamma211) + h21L*gamma131*gamma221 + 
      h32L*gamma231*gamma321 + h31L*gamma331*gamma332) + 
      gamma231*(h11L*gamma122 + h31L*gamma322 + h21L*(gamma121 + gamma222 + 
      gamma332)) + gamma211*(h31L*gamma122 + h22L*gamma232 + h33L*gamma322 + 
      h32L*(gamma121 + gamma222 + gamma332)) + gamma311*(3*h31L*gamma132 + 
      h33L*(gamma121 + 2*gamma332)) - PDstandard6th12h31 + 
      h32L*(2*gamma232*gamma311 + gamma221*gamma331 + pow(gamma221,2))) + 
      hInv23*(gamma233*(h11L*gamma122 + h21L*gamma222 + h31L*gamma322) + 
      (h33L*gamma121 + 3*h31L*gamma132)*gamma331 + 2*(gamma131*(h31L*gamma121 
      + h11L*gamma132) + gamma133*(h11L*gamma121 + h21L*gamma221) + 
      h21L*gamma132*gamma231 + h32L*gamma233*gamma321 + 
      h33L*gamma331*gamma332) + gamma231*(h31L*gamma122 + h22L*gamma232 + 
      h33L*gamma322 + h32L*(gamma121 + gamma222 + gamma332)) + (h31L*gamma121 
      + h11L*gamma132)*gamma333 + gamma232*(h11L*gamma132 + 2*h32L*gamma331 + 
      h31L*gamma332 + h21L*(gamma131 + gamma333)) + gamma221*(h31L*gamma132 + 
      h22L*gamma233 + h33L*gamma332 + h32L*(gamma131 + gamma232 + gamma333)) 
      + gamma332*(h11L*gamma133 + h31L*(gamma131 + 2*gamma333)) + 
      gamma321*(3*h31L*gamma133 + h33L*(gamma131 + 2*gamma333)) - 
      PDstandard6th23h31 + h21L*(gamma233*(gamma121 + gamma332) + 
      pow(gamma232,2))) + hInv11*((h11L*gamma131 + 2*h33L*gamma311)*gamma331 
      + gamma111*(2*h11L*gamma131 + h32L*gamma211 + h21L*gamma231 + 
      h33L*gamma311 + h31L*gamma331) + gamma231*(h11L*gamma121 + 
      2*h32L*gamma311 + h31L*gamma321 + h21L*(gamma221 + gamma331)) + 
      gamma211*(2*h21L*gamma131 + h22L*gamma231 + h33L*gamma321 + 
      h32L*(gamma221 + gamma331)) - 0.5*PDstandard6th11h31 + 
      h31L*(gamma121*gamma211 + 3*gamma131*gamma311 + pow(gamma111,2) + 
      pow(gamma331,2))) + hInv13*(gamma233*(h11L*gamma121 + h31L*gamma321 + 
      h21L*(gamma111 + gamma221 + gamma331)) + (h31L*gamma111 + 
      h11L*gamma131)*gamma333 + gamma211*(h31L*gamma132 + h33L*gamma332 + 
      h32L*(gamma131 + gamma232 + gamma333)) + gamma231*(h11L*gamma132 + 
      h33L*gamma321 + h32L*(gamma111 + gamma221 + 3*gamma331) + 
      h31L*(gamma121 + gamma332) + h21L*(3*gamma131 + gamma232 + gamma333)) + 
      gamma311*(3*h31L*gamma133 + h33L*(gamma131 + 2*gamma333)) + 
      gamma331*(h33L*gamma111 + h11L*gamma133 + h31L*(4*gamma131 + 
      2*gamma333)) - PDstandard6th13h31 + h22L*(gamma211*gamma233 + 
      pow(gamma231,2)) + 2*(gamma111*(h31L*gamma131 + h11L*gamma133) + 
      h21L*gamma133*gamma211 + h32L*gamma233*gamma311 + h11L*pow(gamma131,2) 
      + h33L*pow(gamma331,2))) + hInv22*((h11L*gamma132 + 
      2*h33L*gamma321)*gamma332 + gamma121*(2*h11L*gamma132 + h32L*gamma221 + 
      h21L*gamma232 + h33L*gamma321 + h31L*gamma332) + 
      gamma232*(h11L*gamma122 + 2*h32L*gamma321 + h31L*gamma322 + 
      h21L*(gamma222 + gamma332)) + gamma221*(2*h21L*gamma132 + h22L*gamma232 
      + h33L*gamma322 + h32L*(gamma222 + gamma332)) - 0.5*PDstandard6th22h31 
      + h31L*(gamma122*gamma221 + 3*gamma132*gamma321 + pow(gamma121,2) + 
      pow(gamma332,2))) + hInv33*((h11L*gamma133 + 2*h33L*gamma331)*gamma333 
      + gamma131*(2*h11L*gamma133 + h32L*gamma231 + h21L*gamma233 + 
      h33L*gamma331 + h31L*gamma333) + gamma233*(h11L*gamma132 + 
      2*h32L*gamma331 + h31L*gamma332 + h21L*(gamma232 + gamma333)) + 
      gamma231*(2*h21L*gamma133 + h22L*gamma233 + h33L*gamma332 + 
      h32L*(gamma232 + gamma333)) - 0.5*PDstandard6th33h31 + 
      h31L*(gamma132*gamma231 + 3*gamma133*gamma331 + pow(gamma131,2) + 
      pow(gamma333,2)));
    
    CCTK_REAL R22 CCTK_ATTRIBUTE_UNUSED = 
      hInv13*(2*(gamma132*(h22L*gamma211 + h21L*(gamma111 + 2*gamma221) + 
      h32L*gamma311 + h31L*gamma321) + gamma232*(3*h22L*gamma221 + 
      2*h32L*gamma321) + (h21L*gamma131 + h22L*gamma231 + h32L*(2*gamma221 + 
      gamma331))*gamma332 + gamma121*(h11L*gamma132 + h22L*gamma231 + 
      h21L*(gamma131 + 2*gamma232) + h32L*gamma331 + h31L*gamma332) + 
      gamma321*(h21L*gamma133 + h22L*gamma233 + h33L*gamma332 + 
      h32L*gamma333)) - PDstandard6th13h22) + h21L*(Gam1L*gamma121 + 
      Gam2L*gamma122 + Gam3L*gamma132 + PDstandard6th2Gam1) + 
      h22L*(Gam1L*gamma221 + Gam2L*gamma222 + Gam3L*gamma232 + 
      PDstandard6th2Gam2) + h32L*(Gam1L*gamma321 + Gam2L*gamma322 + 
      Gam3L*gamma332 + PDstandard6th2Gam3) + hInv12*(-PDstandard6th12h22 + 
      2*(gamma122*(h32L*gamma311 + h31L*gamma321) + gamma121*(h11L*gamma122 + 
      h22L*gamma221 + 2*h21L*gamma222 + h32L*gamma321 + h31L*gamma322) + 
      h22L*(gamma122*gamma211 + 3*gamma221*gamma222 + gamma232*gamma321 + 
      gamma231*gamma322) + gamma322*(h33L*gamma321 + h32L*gamma331) + 
      h32L*(2*(gamma222*gamma321 + gamma221*gamma322) + gamma321*gamma332) + 
      h21L*(gamma122*(gamma111 + 2*gamma221) + gamma132*gamma321 + 
      gamma131*gamma322 + pow(gamma121,2)))) + 
      hInv11*(2*gamma121*(h22L*gamma211 + h21L*(gamma111 + 2*gamma221) + 
      h32L*gamma311 + h31L*gamma321) + gamma321*(4*h32L*gamma221 + 
      2*(h21L*gamma131 + h32L*gamma331)) - 0.5*PDstandard6th11h22 + 
      h11L*pow(gamma121,2) + h22L*(2*gamma231*gamma321 + 3*pow(gamma221,2)) + 
      h33L*pow(gamma321,2)) + hInv22*(2*gamma122*(h22L*gamma221 + 
      h21L*(gamma121 + 2*gamma222) + h32L*gamma321 + h31L*gamma322) + 
      gamma322*(4*h32L*gamma222 + 2*(h21L*gamma132 + h32L*gamma332)) - 
      0.5*PDstandard6th22h22 + h11L*pow(gamma122,2) + 
      h22L*(2*gamma232*gamma322 + 3*pow(gamma222,2)) + h33L*pow(gamma322,2)) 
      + hInv33*(2*gamma132*(h22L*gamma231 + h21L*(gamma131 + 2*gamma232) + 
      h32L*gamma331 + h31L*gamma332) + gamma332*(4*h32L*gamma232 + 
      2*(h21L*gamma133 + h32L*gamma333)) - 0.5*PDstandard6th33h22 + 
      h11L*pow(gamma132,2) + h22L*(2*gamma233*gamma332 + 3*pow(gamma232,2)) + 
      h33L*pow(gamma332,2)) + hInv23*(-PDstandard6th23h22 + 
      2*(gamma122*(h11L*gamma132 + h22L*gamma231 + h21L*(gamma131 + 
      2*gamma232) + h32L*gamma331 + h31L*gamma332) + gamma132*(h22L*gamma221 
      + h32L*gamma321 + h31L*gamma322 + h21L*(gamma121 + 2*gamma222 + 
      gamma332)) + gamma232*(2*h32L*gamma322 + h22L*(3*gamma222 + gamma332)) 
      + gamma322*(h21L*gamma133 + h22L*gamma233 + h33L*gamma332 + 
      h32L*gamma333) + h32L*(2*gamma222*gamma332 + pow(gamma332,2))));
    
    CCTK_REAL R32 CCTK_ATTRIBUTE_UNUSED = 0.5*(h31L*(Gam1L*gamma121 + 
      Gam2L*gamma122 + Gam3L*gamma132 + PDstandard6th2Gam1) + 
      h33L*(Gam1L*gamma321 + Gam2L*gamma322 + Gam3L*gamma332 + 
      PDstandard6th2Gam3) + h21L*(Gam1L*gamma131 + Gam2L*gamma132 + 
      Gam3L*gamma133 + PDstandard6th3Gam1) + h22L*(Gam1L*gamma231 + 
      Gam2L*gamma232 + Gam3L*gamma233 + PDstandard6th3Gam2) + 
      h32L*(Gam1L*(gamma221 + gamma331) + Gam2L*(gamma222 + gamma332) + 
      Gam3L*(gamma232 + gamma333) + PDstandard6th2Gam2 + PDstandard6th3Gam3)) 
      + hInv12*(gamma132*(h11L*gamma121 + h22L*gamma211 + h32L*gamma311 + 
      h21L*(gamma111 + gamma221 + gamma331)) + h22L*(gamma232*gamma331 + 
      gamma231*gamma332) + 2*((h21L*gamma122 + h22L*gamma222)*gamma231 + 
      (h21L*gamma121 + h22L*gamma221)*gamma232 + gamma322*(h31L*gamma131 + 
      h33L*gamma331) + h32L*gamma331*gamma332) + gamma131*(h11L*gamma122 + 
      h22L*gamma221 + h32L*gamma321 + h21L*(gamma121 + gamma222 + gamma332)) 
      + h32L*(gamma122*gamma211 + 3*(gamma232*gamma321 + gamma231*gamma322) + 
      gamma222*gamma331 + gamma221*(gamma121 + 2*gamma222 + gamma332)) + 
      h33L*(gamma122*gamma311 + gamma221*gamma322 + gamma321*(gamma121 + 
      gamma222 + 2*gamma332)) - PDstandard6th12h32 + 
      h31L*(2*gamma132*gamma321 + gamma122*(gamma111 + gamma221 + gamma331) + 
      gamma121*(gamma222 + gamma332) + pow(gamma121,2))) + 
      hInv13*(h11L*(gamma131*gamma132 + gamma121*gamma133) + (h32L*gamma132 + 
      h22L*gamma133)*gamma211 + (h33L*gamma132 + h32L*gamma133)*gamma311 + 
      2*(h21L*gamma121*gamma233 + gamma221*(h32L*gamma232 + h22L*gamma233) + 
      h31L*gamma131*gamma332) + gamma221*(h33L*gamma332 + h32L*gamma333) + 
      h31L*(2*gamma133*gamma321 + gamma132*(gamma111 + gamma221 + gamma331) + 
      gamma121*(gamma131 + gamma232 + gamma333)) + gamma231*(h32L*(gamma121 + 
      3*gamma332) + h22L*(gamma131 + 2*gamma232 + gamma333)) + 
      gamma321*(3*h32L*gamma233 + h33L*(gamma232 + 2*gamma333)) + 
      gamma331*(h22L*gamma233 + h33L*(gamma121 + 2*gamma332) + h32L*(gamma131 
      + gamma232 + 2*gamma333)) - PDstandard6th13h32 + 
      h21L*(2*gamma132*gamma231 + gamma133*(gamma111 + gamma221 + gamma331) + 
      gamma131*(gamma232 + gamma333) + pow(gamma131,2))) + 
      hInv11*((h22L*gamma231 + 2*h33L*gamma321)*gamma331 + 
      gamma221*(2*h22L*gamma231 + h33L*gamma321 + h32L*gamma331) + 
      gamma131*(h22L*gamma211 + 2*h31L*gamma321 + h21L*(gamma111 + gamma221 + 
      gamma331)) + gamma121*(h11L*gamma131 + h32L*gamma211 + 2*h21L*gamma231 
      + h33L*gamma311 + h31L*(gamma111 + gamma221 + gamma331)) - 
      0.5*PDstandard6th11h32 + h32L*(gamma131*gamma311 + 3*gamma231*gamma321 
      + pow(gamma221,2) + pow(gamma331,2))) + hInv22*((h22L*gamma232 + 
      2*h33L*gamma322)*gamma332 + gamma222*(2*h22L*gamma232 + h33L*gamma322 + 
      h32L*gamma332) + gamma132*(h22L*gamma221 + 2*h31L*gamma322 + 
      h21L*(gamma121 + gamma222 + gamma332)) + gamma122*(h11L*gamma132 + 
      h32L*gamma221 + 2*h21L*gamma232 + h33L*gamma321 + h31L*(gamma121 + 
      gamma222 + gamma332)) - 0.5*PDstandard6th22h32 + 
      h32L*(gamma132*gamma321 + 3*gamma232*gamma322 + pow(gamma222,2) + 
      pow(gamma332,2))) + hInv23*((h32L*gamma132 + h22L*gamma133)*gamma221 + 
      (h32L*gamma122 + h22L*gamma132)*gamma231 + (h33L*gamma132 + 
      h32L*gamma133)*gamma321 + (h33L*gamma122 + h32L*gamma132)*gamma331 + 
      (h32L*gamma222 + h22L*gamma232)*gamma333 + h31L*(2*gamma133*gamma322 + 
      gamma132*(gamma121 + gamma222 + 3*gamma332) + gamma122*(gamma131 + 
      gamma232 + gamma333)) + h21L*(2*gamma122*gamma233 + gamma133*(gamma121 
      + gamma222 + gamma332) + gamma132*(gamma131 + 3*gamma232 + gamma333)) + 
      gamma322*(3*h32L*gamma233 + h33L*(gamma232 + 2*gamma333)) + 
      gamma332*(h33L*gamma222 + h22L*gamma233 + h32L*(4*gamma232 + 
      2*gamma333)) - PDstandard6th23h32 + h11L*(gamma122*gamma133 + 
      pow(gamma132,2)) + 2*(gamma222*(h32L*gamma232 + h22L*gamma233) + 
      h22L*pow(gamma232,2) + h33L*pow(gamma332,2))) + hInv33*((h22L*gamma233 
      + 2*h33L*gamma332)*gamma333 + gamma232*(2*h22L*gamma233 + h33L*gamma332 
      + h32L*gamma333) + gamma133*(h22L*gamma231 + 2*h31L*gamma332 + 
      h21L*(gamma131 + gamma232 + gamma333)) + gamma132*(h11L*gamma133 + 
      h32L*gamma231 + 2*h21L*gamma233 + h33L*gamma331 + h31L*(gamma131 + 
      gamma232 + gamma333)) - 0.5*PDstandard6th33h32 + 
      h32L*(gamma133*gamma331 + 3*gamma233*gamma332 + pow(gamma232,2) + 
      pow(gamma333,2)));
    
    CCTK_REAL R33 CCTK_ATTRIBUTE_UNUSED = 
      hInv12*(2*(gamma132*(h32L*gamma211 + h21L*gamma231 + h33L*gamma311 + 
      h31L*(gamma111 + 2*gamma331)) + gamma232*(h31L*gamma121 + h33L*gamma321 
      + h32L*(gamma221 + 2*gamma331)) + 3*h33L*gamma331*gamma332 + 
      gamma131*(h11L*gamma132 + h32L*gamma221 + h21L*gamma232 + h33L*gamma321 
      + h31L*(gamma121 + 2*gamma332)) + gamma231*(h31L*gamma122 + 
      h22L*gamma232 + h33L*gamma322 + h32L*(gamma222 + 2*gamma332))) - 
      PDstandard6th12h33) + h31L*(Gam1L*gamma131 + Gam2L*gamma132 + 
      Gam3L*gamma133 + PDstandard6th3Gam1) + h32L*(Gam1L*gamma231 + 
      Gam2L*gamma232 + Gam3L*gamma233 + PDstandard6th3Gam2) + 
      h33L*(Gam1L*gamma331 + Gam2L*gamma332 + Gam3L*gamma333 + 
      PDstandard6th3Gam3) + hInv13*(-PDstandard6th13h33 + 
      2*(gamma233*(h22L*gamma231 + h33L*gamma321) + gamma133*(h32L*gamma211 + 
      h21L*gamma231 + h33L*gamma311 + h31L*(gamma111 + 2*gamma331)) + 
      gamma131*(h11L*gamma133 + h32L*gamma231 + h21L*gamma233 + h33L*gamma331 
      + 2*h31L*gamma333) + h33L*(gamma231*gamma332 + 3*gamma331*gamma333) + 
      h32L*(gamma233*(gamma221 + 2*gamma331) + gamma231*(gamma232 + 
      2*gamma333)) + h31L*(gamma132*gamma231 + gamma121*gamma233 + 
      pow(gamma131,2)))) + hInv23*(-PDstandard6th23h33 + 
      2*(gamma233*(h31L*gamma122 + h22L*gamma232 + h33L*gamma322) + 
      gamma133*(h32L*gamma221 + h21L*gamma232 + h33L*gamma321 + 
      h31L*(gamma121 + 2*gamma332)) + 3*h33L*gamma332*gamma333 + 
      gamma232*(h33L*gamma332 + 2*h32L*gamma333) + gamma132*(h11L*gamma133 + 
      h32L*gamma231 + h21L*gamma233 + h33L*gamma331 + h31L*(gamma131 + 
      gamma232 + 2*gamma333)) + h32L*(gamma233*(gamma222 + 2*gamma332) + 
      pow(gamma232,2)))) + hInv11*(4*h32L*gamma231*gamma331 + 
      2*(gamma231*(h31L*gamma121 + h32L*gamma221 + h33L*gamma321) + 
      gamma131*(h32L*gamma211 + h21L*gamma231 + h33L*gamma311 + 
      h31L*(gamma111 + 2*gamma331))) - 0.5*PDstandard6th11h33 + 
      h11L*pow(gamma131,2) + h22L*pow(gamma231,2) + 3*h33L*pow(gamma331,2)) + 
      hInv22*(4*h32L*gamma232*gamma332 + 2*(gamma232*(h31L*gamma122 + 
      h32L*gamma222 + h33L*gamma322) + gamma132*(h32L*gamma221 + 
      h21L*gamma232 + h33L*gamma321 + h31L*(gamma121 + 2*gamma332))) - 
      0.5*PDstandard6th22h33 + h11L*pow(gamma132,2) + h22L*pow(gamma232,2) + 
      3*h33L*pow(gamma332,2)) + hInv33*(4*h32L*gamma233*gamma333 + 
      2*(gamma233*(h31L*gamma132 + h32L*gamma232 + h33L*gamma332) + 
      gamma133*(h32L*gamma231 + h21L*gamma233 + h33L*gamma331 + 
      h31L*(gamma131 + 2*gamma333))) - 0.5*PDstandard6th33h33 + 
      h11L*pow(gamma133,2) + h22L*pow(gamma233,2) + 3*h33L*pow(gamma333,2));
    
    CCTK_REAL DDphi11 CCTK_ATTRIBUTE_UNUSED = PDstandard6th11phi - 
      gamma111*PDstandard6th1phi - gamma211*PDstandard6th2phi - 
      gamma311*PDstandard6th3phi;
    
    CCTK_REAL DDphi21 CCTK_ATTRIBUTE_UNUSED = PDstandard6th12phi - 
      gamma121*PDstandard6th1phi - gamma221*PDstandard6th2phi - 
      gamma321*PDstandard6th3phi;
    
    CCTK_REAL DDphi31 CCTK_ATTRIBUTE_UNUSED = PDstandard6th13phi - 
      gamma131*PDstandard6th1phi - gamma231*PDstandard6th2phi - 
      gamma331*PDstandard6th3phi;
    
    CCTK_REAL DDphi12 CCTK_ATTRIBUTE_UNUSED = PDstandard6th12phi - 
      gamma121*PDstandard6th1phi - gamma221*PDstandard6th2phi - 
      gamma321*PDstandard6th3phi;
    
    CCTK_REAL DDphi22 CCTK_ATTRIBUTE_UNUSED = 
      -(gamma122*PDstandard6th1phi) + PDstandard6th22phi - 
      gamma222*PDstandard6th2phi - gamma322*PDstandard6th3phi;
    
    CCTK_REAL DDphi32 CCTK_ATTRIBUTE_UNUSED = 
      -(gamma132*PDstandard6th1phi) + PDstandard6th23phi - 
      gamma232*PDstandard6th2phi - gamma332*PDstandard6th3phi;
    
    CCTK_REAL DDphi13 CCTK_ATTRIBUTE_UNUSED = PDstandard6th13phi - 
      gamma131*PDstandard6th1phi - gamma231*PDstandard6th2phi - 
      gamma331*PDstandard6th3phi;
    
    CCTK_REAL DDphi23 CCTK_ATTRIBUTE_UNUSED = 
      -(gamma132*PDstandard6th1phi) + PDstandard6th23phi - 
      gamma232*PDstandard6th2phi - gamma332*PDstandard6th3phi;
    
    CCTK_REAL DDphi33 CCTK_ATTRIBUTE_UNUSED = 
      -(gamma133*PDstandard6th1phi) - gamma233*PDstandard6th2phi + 
      PDstandard6th33phi - gamma333*PDstandard6th3phi;
    
    CCTK_REAL Dphi1 CCTK_ATTRIBUTE_UNUSED = PDstandard6th1phi;
    
    CCTK_REAL Dphi2 CCTK_ATTRIBUTE_UNUSED = PDstandard6th2phi;
    
    CCTK_REAL Dphi3 CCTK_ATTRIBUTE_UNUSED = PDstandard6th3phi;
    
    CCTK_REAL Rphi11 CCTK_ATTRIBUTE_UNUSED = 4*pow(Dphi1,2) - 2*(DDphi11 + 
      h11L*((DDphi12 + DDphi21 + 4*Dphi1*Dphi2)*hInv12 + (DDphi13 + DDphi31 + 
      4*Dphi1*Dphi3)*hInv13 + (DDphi23 + DDphi32 + 4*Dphi2*Dphi3)*hInv23 + 
      hInv11*(DDphi11 + 2*pow(Dphi1,2)) + hInv22*(DDphi22 + 2*pow(Dphi2,2)) + 
      hInv33*(DDphi33 + 2*pow(Dphi3,2))));
    
    CCTK_REAL Rphi21 CCTK_ATTRIBUTE_UNUSED = 4*Dphi1*Dphi2 - 2*(DDphi21 + 
      h21L*((DDphi12 + DDphi21 + 4*Dphi1*Dphi2)*hInv12 + (DDphi13 + DDphi31 + 
      4*Dphi1*Dphi3)*hInv13 + (DDphi23 + DDphi32 + 4*Dphi2*Dphi3)*hInv23 + 
      hInv11*(DDphi11 + 2*pow(Dphi1,2)) + hInv22*(DDphi22 + 2*pow(Dphi2,2)) + 
      hInv33*(DDphi33 + 2*pow(Dphi3,2))));
    
    CCTK_REAL Rphi31 CCTK_ATTRIBUTE_UNUSED = 4*Dphi1*Dphi3 - 2*(DDphi31 + 
      h31L*((DDphi12 + DDphi21 + 4*Dphi1*Dphi2)*hInv12 + (DDphi13 + DDphi31 + 
      4*Dphi1*Dphi3)*hInv13 + (DDphi23 + DDphi32 + 4*Dphi2*Dphi3)*hInv23 + 
      hInv11*(DDphi11 + 2*pow(Dphi1,2)) + hInv22*(DDphi22 + 2*pow(Dphi2,2)) + 
      hInv33*(DDphi33 + 2*pow(Dphi3,2))));
    
    CCTK_REAL Rphi22 CCTK_ATTRIBUTE_UNUSED = 4*pow(Dphi2,2) - 2*(DDphi22 + 
      h22L*((DDphi12 + DDphi21 + 4*Dphi1*Dphi2)*hInv12 + (DDphi13 + DDphi31 + 
      4*Dphi1*Dphi3)*hInv13 + (DDphi23 + DDphi32 + 4*Dphi2*Dphi3)*hInv23 + 
      hInv11*(DDphi11 + 2*pow(Dphi1,2)) + hInv22*(DDphi22 + 2*pow(Dphi2,2)) + 
      hInv33*(DDphi33 + 2*pow(Dphi3,2))));
    
    CCTK_REAL Rphi32 CCTK_ATTRIBUTE_UNUSED = 4*Dphi2*Dphi3 - 2*(DDphi32 + 
      h32L*((DDphi12 + DDphi21 + 4*Dphi1*Dphi2)*hInv12 + (DDphi13 + DDphi31 + 
      4*Dphi1*Dphi3)*hInv13 + (DDphi23 + DDphi32 + 4*Dphi2*Dphi3)*hInv23 + 
      hInv11*(DDphi11 + 2*pow(Dphi1,2)) + hInv22*(DDphi22 + 2*pow(Dphi2,2)) + 
      hInv33*(DDphi33 + 2*pow(Dphi3,2))));
    
    CCTK_REAL Rphi33 CCTK_ATTRIBUTE_UNUSED = 4*pow(Dphi3,2) - 2*(DDphi33 + 
      h33L*((DDphi12 + DDphi21 + 4*Dphi1*Dphi2)*hInv12 + (DDphi13 + DDphi31 + 
      4*Dphi1*Dphi3)*hInv13 + (DDphi23 + DDphi32 + 4*Dphi2*Dphi3)*hInv23 + 
      hInv11*(DDphi11 + 2*pow(Dphi1,2)) + hInv22*(DDphi22 + 2*pow(Dphi2,2)) + 
      hInv33*(DDphi33 + 2*pow(Dphi3,2))));
    
    CCTK_REAL DDalpha11 CCTK_ATTRIBUTE_UNUSED = PDstandard6th11alpha - 
      gamma111*PDstandard6th1alpha - gamma211*PDstandard6th2alpha - 
      gamma311*PDstandard6th3alpha;
    
    CCTK_REAL DDalpha21 CCTK_ATTRIBUTE_UNUSED = PDstandard6th12alpha - 
      gamma121*PDstandard6th1alpha - gamma221*PDstandard6th2alpha - 
      gamma321*PDstandard6th3alpha;
    
    CCTK_REAL DDalpha31 CCTK_ATTRIBUTE_UNUSED = PDstandard6th13alpha - 
      gamma131*PDstandard6th1alpha - gamma231*PDstandard6th2alpha - 
      gamma331*PDstandard6th3alpha;
    
    CCTK_REAL DDalpha12 CCTK_ATTRIBUTE_UNUSED = PDstandard6th12alpha - 
      gamma121*PDstandard6th1alpha - gamma221*PDstandard6th2alpha - 
      gamma321*PDstandard6th3alpha;
    
    CCTK_REAL DDalpha22 CCTK_ATTRIBUTE_UNUSED = 
      -(gamma122*PDstandard6th1alpha) + PDstandard6th22alpha - 
      gamma222*PDstandard6th2alpha - gamma322*PDstandard6th3alpha;
    
    CCTK_REAL DDalpha32 CCTK_ATTRIBUTE_UNUSED = 
      -(gamma132*PDstandard6th1alpha) + PDstandard6th23alpha - 
      gamma232*PDstandard6th2alpha - gamma332*PDstandard6th3alpha;
    
    CCTK_REAL DDalpha13 CCTK_ATTRIBUTE_UNUSED = PDstandard6th13alpha - 
      gamma131*PDstandard6th1alpha - gamma231*PDstandard6th2alpha - 
      gamma331*PDstandard6th3alpha;
    
    CCTK_REAL DDalpha23 CCTK_ATTRIBUTE_UNUSED = 
      -(gamma132*PDstandard6th1alpha) + PDstandard6th23alpha - 
      gamma232*PDstandard6th2alpha - gamma332*PDstandard6th3alpha;
    
    CCTK_REAL DDalpha33 CCTK_ATTRIBUTE_UNUSED = 
      -(gamma133*PDstandard6th1alpha) - gamma233*PDstandard6th2alpha + 
      PDstandard6th33alpha - gamma333*PDstandard6th3alpha;
    
    CCTK_REAL B11 CCTK_ATTRIBUTE_UNUSED = -DDalpha11 + 
      4*Dphi1*PDstandard6th1alpha + alphaL*(R11 + Rphi11);
    
    CCTK_REAL B12 CCTK_ATTRIBUTE_UNUSED = -DDalpha21 + 
      2*(Dphi2*PDstandard6th1alpha + Dphi1*PDstandard6th2alpha) + alphaL*(R21 
      + Rphi21);
    
    CCTK_REAL B13 CCTK_ATTRIBUTE_UNUSED = -DDalpha31 + 
      2*(Dphi3*PDstandard6th1alpha + Dphi1*PDstandard6th3alpha) + alphaL*(R31 
      + Rphi31);
    
    CCTK_REAL B21 CCTK_ATTRIBUTE_UNUSED = -DDalpha12 + 
      2*(Dphi2*PDstandard6th1alpha + Dphi1*PDstandard6th2alpha) + alphaL*(R21 
      + Rphi21);
    
    CCTK_REAL B22 CCTK_ATTRIBUTE_UNUSED = -DDalpha22 + 
      4*Dphi2*PDstandard6th2alpha + alphaL*(R22 + Rphi22);
    
    CCTK_REAL B23 CCTK_ATTRIBUTE_UNUSED = -DDalpha32 + 
      2*(Dphi3*PDstandard6th2alpha + Dphi2*PDstandard6th3alpha) + alphaL*(R32 
      + Rphi32);
    
    CCTK_REAL B31 CCTK_ATTRIBUTE_UNUSED = -DDalpha13 + 
      2*(Dphi3*PDstandard6th1alpha + Dphi1*PDstandard6th3alpha) + alphaL*(R31 
      + Rphi31);
    
    CCTK_REAL B32 CCTK_ATTRIBUTE_UNUSED = -DDalpha23 + 
      2*(Dphi3*PDstandard6th2alpha + Dphi2*PDstandard6th3alpha) + alphaL*(R32 
      + Rphi32);
    
    CCTK_REAL B33 CCTK_ATTRIBUTE_UNUSED = -DDalpha33 + 
      4*Dphi3*PDstandard6th3alpha + alphaL*(R33 + Rphi33);
    
    CCTK_REAL trB CCTK_ATTRIBUTE_UNUSED = B11*hInv11 + (B12 + B21)*hInv12 
      + (B13 + B31)*hInv13 + B22*hInv22 + (B23 + B32)*hInv23 + B33*hInv33;
    
    CCTK_REAL bssnhamL CCTK_ATTRIBUTE_UNUSED = -(A11L*AInv11) - 
      2*A21L*AInv12 - 2*A31L*AInv13 - A22L*AInv22 - 2*A32L*AInv23 - 
      A33L*AInv33 + em4phi*hInv11*(R11 + Rphi11) + 2*em4phi*hInv12*(R21 + 
      Rphi21) + em4phi*hInv22*(R22 + Rphi22) + 2*em4phi*hInv13*(R31 + Rphi31) 
      + 2*em4phi*hInv23*(R32 + Rphi32) + em4phi*hInv33*(R33 + Rphi33) + 
      0.666666666666666666666666666667*pow(KL,2);
    
    ptrdiff_t dir1 CCTK_ATTRIBUTE_UNUSED = isgn(beta1L);
    
    ptrdiff_t dir2 CCTK_ATTRIBUTE_UNUSED = isgn(beta2L);
    
    ptrdiff_t dir3 CCTK_ATTRIBUTE_UNUSED = isgn(beta3L);
    
    CCTK_REAL KrhsL CCTK_ATTRIBUTE_UNUSED = beta1L*PDstandard6th1K + 
      beta2L*PDstandard6th2K + em4phi*(-(hInv11*(DDalpha11 + 
      2*Dphi1*PDstandard6th1alpha)) - hInv12*(DDalpha12 + DDalpha21 + 
      2*Dphi2*PDstandard6th1alpha + 2*Dphi1*PDstandard6th2alpha) - 
      hInv22*(DDalpha22 + 2*Dphi2*PDstandard6th2alpha) - hInv13*(DDalpha13 + 
      DDalpha31 + 2*Dphi3*PDstandard6th1alpha + 2*Dphi1*PDstandard6th3alpha)) 
      + em4phi*(-(hInv23*(DDalpha23 + DDalpha32 + 2*Dphi3*PDstandard6th2alpha 
      + 2*Dphi2*PDstandard6th3alpha)) - hInv33*(DDalpha33 + 
      2*Dphi3*PDstandard6th3alpha)) + beta3L*PDstandard6th3K + 
      alphaL*(A11L*AInv11 + 2*A21L*AInv12 + 2*A31L*AInv13 + A22L*AInv22 + 
      2*A32L*AInv23 + A33L*AInv33 + 
      0.333333333333333333333333333333*pow(KL,2));
    
    CCTK_REAL A11rhsL CCTK_ATTRIBUTE_UNUSED = beta1L*PDstandard6th1A11 + 
      2*(A11L*PDstandard6th1beta1 + A21L*PDstandard6th1beta2 + 
      A31L*PDstandard6th1beta3) + beta2L*PDstandard6th2A11 + 
      beta3L*PDstandard6th3A11 - 
      0.666666666666666666666666666667*A11L*(PDstandard6th1beta1 + 
      PDstandard6th2beta2 + PDstandard6th3beta3) + em4phi*(B11 - 
      0.333333333333333333333333333333*h11L*trB) + alphaL*(A11L*(KL - 
      4*(A21L*hInv12 + A31L*hInv13)) - 4*A21L*A31L*hInv23 - 
      2*(hInv11*pow(A11L,2) + hInv22*pow(A21L,2) + hInv33*pow(A31L,2)));
    
    CCTK_REAL A21rhsL CCTK_ATTRIBUTE_UNUSED = beta1L*PDstandard6th1A21 + 
      A22L*PDstandard6th1beta2 + A32L*PDstandard6th1beta3 + 
      beta2L*PDstandard6th2A21 + A11L*PDstandard6th2beta1 + 
      A31L*PDstandard6th2beta3 + beta3L*PDstandard6th3A21 + A21L*(alphaL*(KL 
      - 2*(A11L*hInv11 + A22L*hInv22)) + 
      0.333333333333333333333333333333*(PDstandard6th1beta1 + 
      PDstandard6th2beta2) - 
      0.666666666666666666666666666667*PDstandard6th3beta3) + em4phi*(B12 - 
      0.333333333333333333333333333333*h21L*trB) - 2*alphaL*((A21L*A31L + 
      A11L*A32L)*hInv13 + (A22L*A31L + A21L*A32L)*hInv23 + A31L*A32L*hInv33 + 
      hInv12*(A11L*A22L + pow(A21L,2)));
    
    CCTK_REAL A31rhsL CCTK_ATTRIBUTE_UNUSED = beta1L*PDstandard6th1A31 + 
      A32L*PDstandard6th1beta2 + A33L*PDstandard6th1beta3 + 
      beta2L*PDstandard6th2A31 + beta3L*PDstandard6th3A31 + 
      A11L*PDstandard6th3beta1 + A21L*PDstandard6th3beta2 + A31L*(alphaL*(KL 
      - 2*(A11L*hInv11 + A33L*hInv33)) - 
      0.666666666666666666666666666667*PDstandard6th2beta2 + 
      0.333333333333333333333333333333*(PDstandard6th1beta1 + 
      PDstandard6th3beta3)) + em4phi*(B13 - 
      0.333333333333333333333333333333*h31L*trB) - 2*alphaL*((A21L*A31L + 
      A11L*A32L)*hInv12 + A21L*A32L*hInv22 + (A31L*A32L + A21L*A33L)*hInv23 + 
      hInv13*(A11L*A33L + pow(A31L,2)));
    
    CCTK_REAL A22rhsL CCTK_ATTRIBUTE_UNUSED = beta1L*PDstandard6th1A22 + 
      beta2L*PDstandard6th2A22 + 2*(A21L*PDstandard6th2beta1 + 
      A22L*PDstandard6th2beta2 + A32L*PDstandard6th2beta3) + 
      beta3L*PDstandard6th3A22 - 
      0.666666666666666666666666666667*A22L*(PDstandard6th1beta1 + 
      PDstandard6th2beta2 + PDstandard6th3beta3) + em4phi*(B22 - 
      0.333333333333333333333333333333*h22L*trB) + 
      alphaL*(-4*A21L*A32L*hInv13 + A22L*(KL - 4*(A21L*hInv12 + A32L*hInv23)) 
      - 2*(hInv11*pow(A21L,2) + hInv22*pow(A22L,2) + hInv33*pow(A32L,2)));
    
    CCTK_REAL A32rhsL CCTK_ATTRIBUTE_UNUSED = beta1L*PDstandard6th1A32 + 
      beta2L*PDstandard6th2A32 + A31L*PDstandard6th2beta1 + 
      A33L*PDstandard6th2beta3 + beta3L*PDstandard6th3A32 + 
      A21L*PDstandard6th3beta1 + A22L*PDstandard6th3beta2 + 
      A32L*(-0.666666666666666666666666666667*PDstandard6th1beta1 + 
      0.333333333333333333333333333333*(PDstandard6th2beta2 + 
      PDstandard6th3beta3)) + em4phi*(B23 - 
      0.333333333333333333333333333333*h32L*trB) + alphaL*(A32L*(KL - 
      2*(A22L*hInv22 + A33L*hInv33)) - 2*(A21L*A31L*hInv11 + (A22L*A31L + 
      A21L*A32L)*hInv12 + (A31L*A32L + A21L*A33L)*hInv13 + hInv23*(A22L*A33L 
      + pow(A32L,2))));
    
    CCTK_REAL A33rhsL CCTK_ATTRIBUTE_UNUSED = beta1L*PDstandard6th1A33 + 
      beta2L*PDstandard6th2A33 + beta3L*PDstandard6th3A33 - 
      0.666666666666666666666666666667*A33L*(PDstandard6th1beta1 + 
      PDstandard6th2beta2 + PDstandard6th3beta3) + 
      2*(A31L*PDstandard6th3beta1 + A32L*PDstandard6th3beta2 + 
      A33L*PDstandard6th3beta3) + em4phi*(B33 - 
      0.333333333333333333333333333333*h33L*trB) + alphaL*(A33L*(KL - 
      4*A31L*hInv13) - 4*A32L*(A31L*hInv12 + A33L*hInv23) - 
      2*(hInv11*pow(A31L,2) + hInv22*pow(A32L,2) + hInv33*pow(A33L,2)));
    
    CCTK_REAL Gam1rhsL CCTK_ATTRIBUTE_UNUSED = alphaL*(12*(AInv11*Dphi1 + 
      AInv12*Dphi2 + AInv13*Dphi3) + 2*AInv11*gamma111 + 4*AInv12*gamma121 + 
      2*AInv22*gamma122 + 4*AInv13*gamma131 + 4*AInv23*gamma132 + 
      2*AInv33*gamma133) + beta1L*PDstandard6th1Gam1 - 
      Gam2L*PDstandard6th2beta1 + beta2L*PDstandard6th2Gam1 - 
      2*(AInv11*PDstandard6th1alpha + AInv12*PDstandard6th2alpha + 
      AInv13*PDstandard6th3alpha) - Gam3L*PDstandard6th3beta1 - 
      1.33333333333333333333333333333*Gam1L*addYoTerm*(PDstandard6th1beta1 + 
      PDstandard6th2beta2 + PDstandard6th3beta3) + 
      Gam1L*(-PDstandard6th1beta1 + 
      0.666666666666666666666666666667*(PDstandard6th1beta1 + 
      PDstandard6th2beta2 + PDstandard6th3beta3)) + 
      0.333333333333333333333333333333*hInv11*(4*PDstandard6th11beta1 + 
      PDstandard6th12beta2 + PDstandard6th13beta3 - 4*alphaL*PDstandard6th1K 
      + 4*addYoTerm*gamma111*(PDstandard6th1beta1 + PDstandard6th2beta2 + 
      PDstandard6th3beta3)) + 
      0.333333333333333333333333333333*hInv12*(7*PDstandard6th12beta1 + 
      PDstandard6th22beta2 + PDstandard6th23beta3 - 4*alphaL*PDstandard6th2K 
      + 8*addYoTerm*gamma121*(PDstandard6th1beta1 + PDstandard6th2beta2 + 
      PDstandard6th3beta3)) + hInv22*(PDstandard6th22beta1 + 
      1.33333333333333333333333333333*addYoTerm*gamma122*(PDstandard6th1beta1 
      + PDstandard6th2beta2 + PDstandard6th3beta3)) + 
      hInv23*(2*PDstandard6th23beta1 + 
      2.66666666666666666666666666667*addYoTerm*gamma132*(PDstandard6th1beta1 
      + PDstandard6th2beta2 + PDstandard6th3beta3)) + 
      hInv33*(PDstandard6th33beta1 + 
      1.33333333333333333333333333333*addYoTerm*gamma133*(PDstandard6th1beta1 
      + PDstandard6th2beta2 + PDstandard6th3beta3)) + 
      beta3L*PDstandard6th3Gam1 + 
      0.333333333333333333333333333333*hInv13*(7*PDstandard6th13beta1 + 
      PDstandard6th23beta2 + PDstandard6th33beta3 + 
      8*addYoTerm*gamma131*(PDstandard6th1beta1 + PDstandard6th2beta2 + 
      PDstandard6th3beta3) - 4*alphaL*PDstandard6th3K);
    
    CCTK_REAL Gam2rhsL CCTK_ATTRIBUTE_UNUSED = alphaL*(12*(AInv12*Dphi1 + 
      AInv22*Dphi2 + AInv23*Dphi3) + 2*AInv11*gamma211 + 4*AInv12*gamma221 + 
      2*AInv22*gamma222 + 4*AInv13*gamma231 + 4*AInv23*gamma232 + 
      2*AInv33*gamma233) - Gam1L*PDstandard6th1beta2 + 
      beta1L*PDstandard6th1Gam2 + beta2L*PDstandard6th2Gam2 - 
      2*(AInv12*PDstandard6th1alpha + AInv22*PDstandard6th2alpha + 
      AInv23*PDstandard6th3alpha) - Gam3L*PDstandard6th3beta2 - 
      1.33333333333333333333333333333*Gam2L*addYoTerm*(PDstandard6th1beta1 + 
      PDstandard6th2beta2 + PDstandard6th3beta3) + 
      Gam2L*(-PDstandard6th2beta2 + 
      0.666666666666666666666666666667*(PDstandard6th1beta1 + 
      PDstandard6th2beta2 + PDstandard6th3beta3)) + 
      hInv11*(PDstandard6th11beta2 + 
      1.33333333333333333333333333333*addYoTerm*gamma211*(PDstandard6th1beta1 
      + PDstandard6th2beta2 + PDstandard6th3beta3)) + 
      0.333333333333333333333333333333*hInv12*(PDstandard6th11beta1 + 
      7*PDstandard6th12beta2 + PDstandard6th13beta3 - 
      4*alphaL*PDstandard6th1K + 8*addYoTerm*gamma221*(PDstandard6th1beta1 + 
      PDstandard6th2beta2 + PDstandard6th3beta3)) + 
      0.333333333333333333333333333333*hInv22*(PDstandard6th12beta1 + 
      4*PDstandard6th22beta2 + PDstandard6th23beta3 - 
      4*alphaL*PDstandard6th2K + 4*addYoTerm*gamma222*(PDstandard6th1beta1 + 
      PDstandard6th2beta2 + PDstandard6th3beta3)) + 
      hInv13*(2*PDstandard6th13beta2 + 
      2.66666666666666666666666666667*addYoTerm*gamma231*(PDstandard6th1beta1 
      + PDstandard6th2beta2 + PDstandard6th3beta3)) + 
      hInv33*(PDstandard6th33beta2 + 
      1.33333333333333333333333333333*addYoTerm*gamma233*(PDstandard6th1beta1 
      + PDstandard6th2beta2 + PDstandard6th3beta3)) + 
      beta3L*PDstandard6th3Gam2 + 
      0.333333333333333333333333333333*hInv23*(PDstandard6th13beta1 + 
      7*PDstandard6th23beta2 + PDstandard6th33beta3 + 
      8*addYoTerm*gamma232*(PDstandard6th1beta1 + PDstandard6th2beta2 + 
      PDstandard6th3beta3) - 4*alphaL*PDstandard6th3K);
    
    CCTK_REAL Gam3rhsL CCTK_ATTRIBUTE_UNUSED = alphaL*(12*(AInv13*Dphi1 + 
      AInv23*Dphi2 + AInv33*Dphi3) + 2*AInv11*gamma311 + 4*AInv12*gamma321 + 
      2*AInv22*gamma322 + 4*AInv13*gamma331 + 4*AInv23*gamma332 + 
      2*AInv33*gamma333) - Gam1L*PDstandard6th1beta3 + 
      beta1L*PDstandard6th1Gam3 - Gam2L*PDstandard6th2beta3 + 
      beta2L*PDstandard6th2Gam3 - 2*(AInv13*PDstandard6th1alpha + 
      AInv23*PDstandard6th2alpha + AInv33*PDstandard6th3alpha) - 
      1.33333333333333333333333333333*Gam3L*addYoTerm*(PDstandard6th1beta1 + 
      PDstandard6th2beta2 + PDstandard6th3beta3) + 
      Gam3L*(-PDstandard6th3beta3 + 
      0.666666666666666666666666666667*(PDstandard6th1beta1 + 
      PDstandard6th2beta2 + PDstandard6th3beta3)) + 
      hInv11*(PDstandard6th11beta3 + 
      1.33333333333333333333333333333*addYoTerm*gamma311*(PDstandard6th1beta1 
      + PDstandard6th2beta2 + PDstandard6th3beta3)) + 
      hInv12*(2*PDstandard6th12beta3 + 
      2.66666666666666666666666666667*addYoTerm*gamma321*(PDstandard6th1beta1 
      + PDstandard6th2beta2 + PDstandard6th3beta3)) + 
      hInv22*(PDstandard6th22beta3 + 
      1.33333333333333333333333333333*addYoTerm*gamma322*(PDstandard6th1beta1 
      + PDstandard6th2beta2 + PDstandard6th3beta3)) + 
      0.333333333333333333333333333333*hInv13*(PDstandard6th11beta1 + 
      PDstandard6th12beta2 + 7*PDstandard6th13beta3 - 
      4*alphaL*PDstandard6th1K + 8*addYoTerm*gamma331*(PDstandard6th1beta1 + 
      PDstandard6th2beta2 + PDstandard6th3beta3)) + 
      0.333333333333333333333333333333*hInv23*(PDstandard6th12beta1 + 
      PDstandard6th22beta2 + 7*PDstandard6th23beta3 - 
      4*alphaL*PDstandard6th2K + 8*addYoTerm*gamma332*(PDstandard6th1beta1 + 
      PDstandard6th2beta2 + PDstandard6th3beta3)) + beta3L*PDstandard6th3Gam3 
      + 0.333333333333333333333333333333*hInv33*(PDstandard6th13beta1 + 
      PDstandard6th23beta2 + 4*PDstandard6th33beta3 + 
      4*addYoTerm*gamma333*(PDstandard6th1beta1 + PDstandard6th2beta2 + 
      PDstandard6th3beta3) - 4*alphaL*PDstandard6th3K);
    /* Copy local copies back to grid functions */
    A11rhs[index] = A11rhsL;
    A21rhs[index] = A21rhsL;
    A22rhs[index] = A22rhsL;
    A31rhs[index] = A31rhsL;
    A32rhs[index] = A32rhsL;
    A33rhs[index] = A33rhsL;
    bssnham[index] = bssnhamL;
    Gam1rhs[index] = Gam1rhsL;
    Gam2rhs[index] = Gam2rhsL;
    Gam3rhs[index] = Gam3rhsL;
    Krhs[index] = KrhsL;
  }
  CCTK_ENDLOOP3(bssn_evolve_nonmetric_6thcentred);
}
extern "C" void bssn_evolve_nonmetric_6thcentred(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering bssn_evolve_nonmetric_6thcentred_Body");
  }
  if (cctk_iteration % bssn_evolve_nonmetric_6thcentred_calc_every != bssn_evolve_nonmetric_6thcentred_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "Kranc2BSSN::A_group",
    "Kranc2BSSN::A_grouprhs",
    "Kranc2BSSN::alpha_group",
    "Kranc2BSSN::beta_group",
    "Kranc2BSSN::Gam_group",
    "Kranc2BSSN::Gam_grouprhs",
    "Kranc2BSSN::h_group",
    "Kranc2BSSN::K_group",
    "Kranc2BSSN::K_grouprhs",
    "Kranc2BSSN::phi_group",
    "Kranc2BSSN::scalarconstraints"};
  AssertGroupStorage(cctkGH, "bssn_evolve_nonmetric_6thcentred", 11, groups);
  
  switch (constraintsfdorder)
  {
    case 2:
    {
      EnsureStencilFits(cctkGH, "bssn_evolve_nonmetric_6thcentred", 3, 3, 3);
      break;
    }
    
    case 4:
    {
      EnsureStencilFits(cctkGH, "bssn_evolve_nonmetric_6thcentred", 3, 3, 3);
      break;
    }
    
    case 6:
    {
      EnsureStencilFits(cctkGH, "bssn_evolve_nonmetric_6thcentred", 3, 3, 3);
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  LoopOverInterior(cctkGH, bssn_evolve_nonmetric_6thcentred_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving bssn_evolve_nonmetric_6thcentred_Body");
  }
}

} // namespace Kranc2BSSN
