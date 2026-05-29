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

namespace Kranc2BSSNMatter {

extern "C" void bssnmatter_mueller_shift_4thcentred_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % bssnmatter_mueller_shift_4thcentred_calc_every != bssnmatter_mueller_shift_4thcentred_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "Kranc2BSSNMatter::beta_grouprhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for Kranc2BSSNMatter::beta_grouprhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "Kranc2BSSNMatter::betat_grouprhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for Kranc2BSSNMatter::betat_grouprhs.");
  return;
}

static void bssnmatter_mueller_shift_4thcentred_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  CCTK_LOOP3(bssnmatter_mueller_shift_4thcentred,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL alphaL CCTK_ATTRIBUTE_UNUSED = alpha[index];
    CCTK_REAL beta1L CCTK_ATTRIBUTE_UNUSED = beta1[index];
    CCTK_REAL beta2L CCTK_ATTRIBUTE_UNUSED = beta2[index];
    CCTK_REAL beta3L CCTK_ATTRIBUTE_UNUSED = beta3[index];
    CCTK_REAL betat1L CCTK_ATTRIBUTE_UNUSED = betat1[index];
    CCTK_REAL betat2L CCTK_ATTRIBUTE_UNUSED = betat2[index];
    CCTK_REAL betat3L CCTK_ATTRIBUTE_UNUSED = betat3[index];
    CCTK_REAL Gam1L CCTK_ATTRIBUTE_UNUSED = Gam1[index];
    CCTK_REAL Gam1rhsL CCTK_ATTRIBUTE_UNUSED = Gam1rhs[index];
    CCTK_REAL Gam2L CCTK_ATTRIBUTE_UNUSED = Gam2[index];
    CCTK_REAL Gam2rhsL CCTK_ATTRIBUTE_UNUSED = Gam2rhs[index];
    CCTK_REAL Gam3L CCTK_ATTRIBUTE_UNUSED = Gam3[index];
    CCTK_REAL Gam3rhsL CCTK_ATTRIBUTE_UNUSED = Gam3rhs[index];
    CCTK_REAL h11L CCTK_ATTRIBUTE_UNUSED = h11[index];
    CCTK_REAL h21L CCTK_ATTRIBUTE_UNUSED = h21[index];
    CCTK_REAL h22L CCTK_ATTRIBUTE_UNUSED = h22[index];
    CCTK_REAL h31L CCTK_ATTRIBUTE_UNUSED = h31[index];
    CCTK_REAL h32L CCTK_ATTRIBUTE_UNUSED = h32[index];
    CCTK_REAL h33L CCTK_ATTRIBUTE_UNUSED = h33[index];
    CCTK_REAL phiL CCTK_ATTRIBUTE_UNUSED = phi[index];
    
    /* Include user supplied include files */
    /* Precompute derivatives */
    CCTK_REAL PDstandard4th1beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th2beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th3beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th1beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th2beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th3beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th1beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th2beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th3beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th1betat1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th2betat1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th3betat1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th1betat2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th2betat2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th3betat2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th1betat3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th2betat3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th3betat3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th1Gam1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th2Gam1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th3Gam1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th1Gam2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th2Gam2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th3Gam2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th1Gam3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th2Gam3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th3Gam3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th1phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th2phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th3phi CCTK_ATTRIBUTE_UNUSED;
    
    switch (constraintsfdorder)
    {
      case 2:
      {
        PDstandard4th1beta1 = PDstandard4th1(&beta1[index]);
        PDstandard4th2beta1 = PDstandard4th2(&beta1[index]);
        PDstandard4th3beta1 = PDstandard4th3(&beta1[index]);
        PDstandard4th1beta2 = PDstandard4th1(&beta2[index]);
        PDstandard4th2beta2 = PDstandard4th2(&beta2[index]);
        PDstandard4th3beta2 = PDstandard4th3(&beta2[index]);
        PDstandard4th1beta3 = PDstandard4th1(&beta3[index]);
        PDstandard4th2beta3 = PDstandard4th2(&beta3[index]);
        PDstandard4th3beta3 = PDstandard4th3(&beta3[index]);
        PDstandard4th1betat1 = PDstandard4th1(&betat1[index]);
        PDstandard4th2betat1 = PDstandard4th2(&betat1[index]);
        PDstandard4th3betat1 = PDstandard4th3(&betat1[index]);
        PDstandard4th1betat2 = PDstandard4th1(&betat2[index]);
        PDstandard4th2betat2 = PDstandard4th2(&betat2[index]);
        PDstandard4th3betat2 = PDstandard4th3(&betat2[index]);
        PDstandard4th1betat3 = PDstandard4th1(&betat3[index]);
        PDstandard4th2betat3 = PDstandard4th2(&betat3[index]);
        PDstandard4th3betat3 = PDstandard4th3(&betat3[index]);
        PDstandard4th1Gam1 = PDstandard4th1(&Gam1[index]);
        PDstandard4th2Gam1 = PDstandard4th2(&Gam1[index]);
        PDstandard4th3Gam1 = PDstandard4th3(&Gam1[index]);
        PDstandard4th1Gam2 = PDstandard4th1(&Gam2[index]);
        PDstandard4th2Gam2 = PDstandard4th2(&Gam2[index]);
        PDstandard4th3Gam2 = PDstandard4th3(&Gam2[index]);
        PDstandard4th1Gam3 = PDstandard4th1(&Gam3[index]);
        PDstandard4th2Gam3 = PDstandard4th2(&Gam3[index]);
        PDstandard4th3Gam3 = PDstandard4th3(&Gam3[index]);
        PDstandard4th1phi = PDstandard4th1(&phi[index]);
        PDstandard4th2phi = PDstandard4th2(&phi[index]);
        PDstandard4th3phi = PDstandard4th3(&phi[index]);
        break;
      }
      
      case 4:
      {
        PDstandard4th1beta1 = PDstandard4th1(&beta1[index]);
        PDstandard4th2beta1 = PDstandard4th2(&beta1[index]);
        PDstandard4th3beta1 = PDstandard4th3(&beta1[index]);
        PDstandard4th1beta2 = PDstandard4th1(&beta2[index]);
        PDstandard4th2beta2 = PDstandard4th2(&beta2[index]);
        PDstandard4th3beta2 = PDstandard4th3(&beta2[index]);
        PDstandard4th1beta3 = PDstandard4th1(&beta3[index]);
        PDstandard4th2beta3 = PDstandard4th2(&beta3[index]);
        PDstandard4th3beta3 = PDstandard4th3(&beta3[index]);
        PDstandard4th1betat1 = PDstandard4th1(&betat1[index]);
        PDstandard4th2betat1 = PDstandard4th2(&betat1[index]);
        PDstandard4th3betat1 = PDstandard4th3(&betat1[index]);
        PDstandard4th1betat2 = PDstandard4th1(&betat2[index]);
        PDstandard4th2betat2 = PDstandard4th2(&betat2[index]);
        PDstandard4th3betat2 = PDstandard4th3(&betat2[index]);
        PDstandard4th1betat3 = PDstandard4th1(&betat3[index]);
        PDstandard4th2betat3 = PDstandard4th2(&betat3[index]);
        PDstandard4th3betat3 = PDstandard4th3(&betat3[index]);
        PDstandard4th1Gam1 = PDstandard4th1(&Gam1[index]);
        PDstandard4th2Gam1 = PDstandard4th2(&Gam1[index]);
        PDstandard4th3Gam1 = PDstandard4th3(&Gam1[index]);
        PDstandard4th1Gam2 = PDstandard4th1(&Gam2[index]);
        PDstandard4th2Gam2 = PDstandard4th2(&Gam2[index]);
        PDstandard4th3Gam2 = PDstandard4th3(&Gam2[index]);
        PDstandard4th1Gam3 = PDstandard4th1(&Gam3[index]);
        PDstandard4th2Gam3 = PDstandard4th2(&Gam3[index]);
        PDstandard4th3Gam3 = PDstandard4th3(&Gam3[index]);
        PDstandard4th1phi = PDstandard4th1(&phi[index]);
        PDstandard4th2phi = PDstandard4th2(&phi[index]);
        PDstandard4th3phi = PDstandard4th3(&phi[index]);
        break;
      }
      
      case 6:
      {
        PDstandard4th1beta1 = PDstandard4th1(&beta1[index]);
        PDstandard4th2beta1 = PDstandard4th2(&beta1[index]);
        PDstandard4th3beta1 = PDstandard4th3(&beta1[index]);
        PDstandard4th1beta2 = PDstandard4th1(&beta2[index]);
        PDstandard4th2beta2 = PDstandard4th2(&beta2[index]);
        PDstandard4th3beta2 = PDstandard4th3(&beta2[index]);
        PDstandard4th1beta3 = PDstandard4th1(&beta3[index]);
        PDstandard4th2beta3 = PDstandard4th2(&beta3[index]);
        PDstandard4th3beta3 = PDstandard4th3(&beta3[index]);
        PDstandard4th1betat1 = PDstandard4th1(&betat1[index]);
        PDstandard4th2betat1 = PDstandard4th2(&betat1[index]);
        PDstandard4th3betat1 = PDstandard4th3(&betat1[index]);
        PDstandard4th1betat2 = PDstandard4th1(&betat2[index]);
        PDstandard4th2betat2 = PDstandard4th2(&betat2[index]);
        PDstandard4th3betat2 = PDstandard4th3(&betat2[index]);
        PDstandard4th1betat3 = PDstandard4th1(&betat3[index]);
        PDstandard4th2betat3 = PDstandard4th2(&betat3[index]);
        PDstandard4th3betat3 = PDstandard4th3(&betat3[index]);
        PDstandard4th1Gam1 = PDstandard4th1(&Gam1[index]);
        PDstandard4th2Gam1 = PDstandard4th2(&Gam1[index]);
        PDstandard4th3Gam1 = PDstandard4th3(&Gam1[index]);
        PDstandard4th1Gam2 = PDstandard4th1(&Gam2[index]);
        PDstandard4th2Gam2 = PDstandard4th2(&Gam2[index]);
        PDstandard4th3Gam2 = PDstandard4th3(&Gam2[index]);
        PDstandard4th1Gam3 = PDstandard4th1(&Gam3[index]);
        PDstandard4th2Gam3 = PDstandard4th2(&Gam3[index]);
        PDstandard4th3Gam3 = PDstandard4th3(&Gam3[index]);
        PDstandard4th1phi = PDstandard4th1(&phi[index]);
        PDstandard4th2phi = PDstandard4th2(&phi[index]);
        PDstandard4th3phi = PDstandard4th3(&phi[index]);
        break;
      }
      default:
        CCTK_BUILTIN_UNREACHABLE();
    }
    /* Calculate temporaries and grid functions */
    ptrdiff_t dir1 CCTK_ATTRIBUTE_UNUSED = isgn(beta1L);
    
    ptrdiff_t dir2 CCTK_ATTRIBUTE_UNUSED = isgn(beta2L);
    
    ptrdiff_t dir3 CCTK_ATTRIBUTE_UNUSED = isgn(beta3L);
    
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
    
    CCTK_REAL e2phi CCTK_ATTRIBUTE_UNUSED = exp(2*phiL);
    
    CCTK_REAL muellerTop CCTK_ATTRIBUTE_UNUSED = 
      2*(hInv23*PDstandard4th2phi*PDstandard4th3phi + 
      PDstandard4th1phi*(hInv12*PDstandard4th2phi + 
      hInv13*PDstandard4th3phi)) + hInv11*pow(PDstandard4th1phi,2) + 
      hInv22*pow(PDstandard4th2phi,2) + hInv33*pow(PDstandard4th3phi,2);
    
    CCTK_REAL muellerBottom CCTK_ATTRIBUTE_UNUSED = 0.5*e2phi*pow(1 - 
      pow(e2phi,-ShiftA),ShiftB);
    
    CCTK_REAL mueller CCTK_ATTRIBUTE_UNUSED = 
      etaBeta*pow(muellerBottom,-1)*pow(fabs(muellerTop),0.5);
    
    CCTK_REAL beta1rhsL CCTK_ATTRIBUTE_UNUSED = 
      newNASAAdvection*(beta1L*PDstandard4th1beta1 + 
      beta2L*PDstandard4th2beta1 + beta3L*PDstandard4th3beta1) + 
      betat1L*betaDotAlphaFactor*pow(alphaL,gammaDriverLambda);
    
    CCTK_REAL beta2rhsL CCTK_ATTRIBUTE_UNUSED = 
      newNASAAdvection*(beta1L*PDstandard4th1beta2 + 
      beta2L*PDstandard4th2beta2 + beta3L*PDstandard4th3beta2) + 
      betat2L*betaDotAlphaFactor*pow(alphaL,gammaDriverLambda);
    
    CCTK_REAL beta3rhsL CCTK_ATTRIBUTE_UNUSED = 
      newNASAAdvection*(beta1L*PDstandard4th1beta3 + 
      beta2L*PDstandard4th2beta3 + beta3L*PDstandard4th3beta3) + 
      betat3L*betaDotAlphaFactor*pow(alphaL,gammaDriverLambda);
    
    CCTK_REAL betat1rhsL CCTK_ATTRIBUTE_UNUSED = -(betat1L*mueller) + 
      betatAdvection*(beta1L*PDstandard4th1betat1 + 
      beta2L*PDstandard4th2betat1 + beta3L*PDstandard4th3betat1) - 
      NASAAdvection*(beta1L*PDstandard4th1Gam1 + beta2L*PDstandard4th2Gam1 + 
      beta3L*PDstandard4th3Gam1) + 
      Gam1rhsL*chiBeta*pow(alphaL,gammaDriverLapsePower);
    
    CCTK_REAL betat2rhsL CCTK_ATTRIBUTE_UNUSED = -(betat2L*mueller) + 
      betatAdvection*(beta1L*PDstandard4th1betat2 + 
      beta2L*PDstandard4th2betat2 + beta3L*PDstandard4th3betat2) - 
      NASAAdvection*(beta1L*PDstandard4th1Gam2 + beta2L*PDstandard4th2Gam2 + 
      beta3L*PDstandard4th3Gam2) + 
      Gam2rhsL*chiBeta*pow(alphaL,gammaDriverLapsePower);
    
    CCTK_REAL betat3rhsL CCTK_ATTRIBUTE_UNUSED = -(betat3L*mueller) + 
      betatAdvection*(beta1L*PDstandard4th1betat3 + 
      beta2L*PDstandard4th2betat3 + beta3L*PDstandard4th3betat3) - 
      NASAAdvection*(beta1L*PDstandard4th1Gam3 + beta2L*PDstandard4th2Gam3 + 
      beta3L*PDstandard4th3Gam3) + 
      Gam3rhsL*chiBeta*pow(alphaL,gammaDriverLapsePower);
    /* Copy local copies back to grid functions */
    beta1rhs[index] = beta1rhsL;
    beta2rhs[index] = beta2rhsL;
    beta3rhs[index] = beta3rhsL;
    betat1rhs[index] = betat1rhsL;
    betat2rhs[index] = betat2rhsL;
    betat3rhs[index] = betat3rhsL;
  }
  CCTK_ENDLOOP3(bssnmatter_mueller_shift_4thcentred);
}
extern "C" void bssnmatter_mueller_shift_4thcentred(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering bssnmatter_mueller_shift_4thcentred_Body");
  }
  if (cctk_iteration % bssnmatter_mueller_shift_4thcentred_calc_every != bssnmatter_mueller_shift_4thcentred_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "Kranc2BSSNMatter::alpha_group",
    "Kranc2BSSNMatter::beta_group",
    "Kranc2BSSNMatter::beta_grouprhs",
    "Kranc2BSSNMatter::betat_group",
    "Kranc2BSSNMatter::betat_grouprhs",
    "Kranc2BSSNMatter::Gam_group",
    "Kranc2BSSNMatter::Gam_grouprhs",
    "Kranc2BSSNMatter::h_group",
    "Kranc2BSSNMatter::phi_group"};
  AssertGroupStorage(cctkGH, "bssnmatter_mueller_shift_4thcentred", 9, groups);
  
  switch (constraintsfdorder)
  {
    case 2:
    {
      EnsureStencilFits(cctkGH, "bssnmatter_mueller_shift_4thcentred", 2, 2, 2);
      break;
    }
    
    case 4:
    {
      EnsureStencilFits(cctkGH, "bssnmatter_mueller_shift_4thcentred", 2, 2, 2);
      break;
    }
    
    case 6:
    {
      EnsureStencilFits(cctkGH, "bssnmatter_mueller_shift_4thcentred", 2, 2, 2);
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  LoopOverInterior(cctkGH, bssnmatter_mueller_shift_4thcentred_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving bssnmatter_mueller_shift_4thcentred_Body");
  }
}

} // namespace Kranc2BSSNMatter
