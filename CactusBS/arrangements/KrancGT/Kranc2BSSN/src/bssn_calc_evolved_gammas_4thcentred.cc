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

extern "C" void bssn_calc_evolved_gammas_4thcentred_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % bssn_calc_evolved_gammas_4thcentred_calc_every != bssn_calc_evolved_gammas_4thcentred_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "Kranc2BSSN::Gam_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for Kranc2BSSN::Gam_group.");
  return;
}

static void bssn_calc_evolved_gammas_4thcentred_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  CCTK_LOOP3(bssn_calc_evolved_gammas_4thcentred,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL beta1L CCTK_ATTRIBUTE_UNUSED = beta1[index];
    CCTK_REAL beta2L CCTK_ATTRIBUTE_UNUSED = beta2[index];
    CCTK_REAL beta3L CCTK_ATTRIBUTE_UNUSED = beta3[index];
    CCTK_REAL h11L CCTK_ATTRIBUTE_UNUSED = h11[index];
    CCTK_REAL h21L CCTK_ATTRIBUTE_UNUSED = h21[index];
    CCTK_REAL h22L CCTK_ATTRIBUTE_UNUSED = h22[index];
    CCTK_REAL h31L CCTK_ATTRIBUTE_UNUSED = h31[index];
    CCTK_REAL h32L CCTK_ATTRIBUTE_UNUSED = h32[index];
    CCTK_REAL h33L CCTK_ATTRIBUTE_UNUSED = h33[index];
    
    /* Include user supplied include files */
    /* Precompute derivatives */
    CCTK_REAL PDstandard4th1h11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th2h11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th3h11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th1h21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th2h21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th3h21 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th1h22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th2h22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th3h22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th1h31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th2h31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th3h31 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th1h32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th2h32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th3h32 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th1h33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th2h33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard4th3h33 CCTK_ATTRIBUTE_UNUSED;
    
    switch (constraintsfdorder)
    {
      case 2:
      {
        PDstandard4th1h11 = PDstandard4th1(&h11[index]);
        PDstandard4th2h11 = PDstandard4th2(&h11[index]);
        PDstandard4th3h11 = PDstandard4th3(&h11[index]);
        PDstandard4th1h21 = PDstandard4th1(&h21[index]);
        PDstandard4th2h21 = PDstandard4th2(&h21[index]);
        PDstandard4th3h21 = PDstandard4th3(&h21[index]);
        PDstandard4th1h22 = PDstandard4th1(&h22[index]);
        PDstandard4th2h22 = PDstandard4th2(&h22[index]);
        PDstandard4th3h22 = PDstandard4th3(&h22[index]);
        PDstandard4th1h31 = PDstandard4th1(&h31[index]);
        PDstandard4th2h31 = PDstandard4th2(&h31[index]);
        PDstandard4th3h31 = PDstandard4th3(&h31[index]);
        PDstandard4th1h32 = PDstandard4th1(&h32[index]);
        PDstandard4th2h32 = PDstandard4th2(&h32[index]);
        PDstandard4th3h32 = PDstandard4th3(&h32[index]);
        PDstandard4th1h33 = PDstandard4th1(&h33[index]);
        PDstandard4th2h33 = PDstandard4th2(&h33[index]);
        PDstandard4th3h33 = PDstandard4th3(&h33[index]);
        break;
      }
      
      case 4:
      {
        PDstandard4th1h11 = PDstandard4th1(&h11[index]);
        PDstandard4th2h11 = PDstandard4th2(&h11[index]);
        PDstandard4th3h11 = PDstandard4th3(&h11[index]);
        PDstandard4th1h21 = PDstandard4th1(&h21[index]);
        PDstandard4th2h21 = PDstandard4th2(&h21[index]);
        PDstandard4th3h21 = PDstandard4th3(&h21[index]);
        PDstandard4th1h22 = PDstandard4th1(&h22[index]);
        PDstandard4th2h22 = PDstandard4th2(&h22[index]);
        PDstandard4th3h22 = PDstandard4th3(&h22[index]);
        PDstandard4th1h31 = PDstandard4th1(&h31[index]);
        PDstandard4th2h31 = PDstandard4th2(&h31[index]);
        PDstandard4th3h31 = PDstandard4th3(&h31[index]);
        PDstandard4th1h32 = PDstandard4th1(&h32[index]);
        PDstandard4th2h32 = PDstandard4th2(&h32[index]);
        PDstandard4th3h32 = PDstandard4th3(&h32[index]);
        PDstandard4th1h33 = PDstandard4th1(&h33[index]);
        PDstandard4th2h33 = PDstandard4th2(&h33[index]);
        PDstandard4th3h33 = PDstandard4th3(&h33[index]);
        break;
      }
      
      case 6:
      {
        PDstandard4th1h11 = PDstandard4th1(&h11[index]);
        PDstandard4th2h11 = PDstandard4th2(&h11[index]);
        PDstandard4th3h11 = PDstandard4th3(&h11[index]);
        PDstandard4th1h21 = PDstandard4th1(&h21[index]);
        PDstandard4th2h21 = PDstandard4th2(&h21[index]);
        PDstandard4th3h21 = PDstandard4th3(&h21[index]);
        PDstandard4th1h22 = PDstandard4th1(&h22[index]);
        PDstandard4th2h22 = PDstandard4th2(&h22[index]);
        PDstandard4th3h22 = PDstandard4th3(&h22[index]);
        PDstandard4th1h31 = PDstandard4th1(&h31[index]);
        PDstandard4th2h31 = PDstandard4th2(&h31[index]);
        PDstandard4th3h31 = PDstandard4th3(&h31[index]);
        PDstandard4th1h32 = PDstandard4th1(&h32[index]);
        PDstandard4th2h32 = PDstandard4th2(&h32[index]);
        PDstandard4th3h32 = PDstandard4th3(&h32[index]);
        PDstandard4th1h33 = PDstandard4th1(&h33[index]);
        PDstandard4th2h33 = PDstandard4th2(&h33[index]);
        PDstandard4th3h33 = PDstandard4th3(&h33[index]);
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
    
    CCTK_REAL gamma111 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv11*PDstandard4th1h11 + 2*(hInv12*PDstandard4th1h21 + 
      hInv13*PDstandard4th1h31) - hInv12*PDstandard4th2h11 - 
      hInv13*PDstandard4th3h11);
    
    CCTK_REAL gamma211 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv12*PDstandard4th1h11 + 2*(hInv22*PDstandard4th1h21 + 
      hInv23*PDstandard4th1h31) - hInv22*PDstandard4th2h11 - 
      hInv23*PDstandard4th3h11);
    
    CCTK_REAL gamma311 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv13*PDstandard4th1h11 + 2*(hInv23*PDstandard4th1h21 + 
      hInv33*PDstandard4th1h31) - hInv23*PDstandard4th2h11 - 
      hInv33*PDstandard4th3h11);
    
    CCTK_REAL gamma121 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv12*PDstandard4th1h22 + hInv11*PDstandard4th2h11 + 
      hInv13*(PDstandard4th1h32 + PDstandard4th2h31 - PDstandard4th3h21));
    
    CCTK_REAL gamma221 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv22*PDstandard4th1h22 + hInv12*PDstandard4th2h11 + 
      hInv23*(PDstandard4th1h32 + PDstandard4th2h31 - PDstandard4th3h21));
    
    CCTK_REAL gamma321 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv23*PDstandard4th1h22 + hInv13*PDstandard4th2h11 + 
      hInv33*(PDstandard4th1h32 + PDstandard4th2h31 - PDstandard4th3h21));
    
    CCTK_REAL gamma131 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv13*PDstandard4th1h33 + hInv11*PDstandard4th3h11 + 
      hInv12*(PDstandard4th1h32 - PDstandard4th2h31 + PDstandard4th3h21));
    
    CCTK_REAL gamma231 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv23*PDstandard4th1h33 + hInv12*PDstandard4th3h11 + 
      hInv22*(PDstandard4th1h32 - PDstandard4th2h31 + PDstandard4th3h21));
    
    CCTK_REAL gamma331 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv33*PDstandard4th1h33 + hInv13*PDstandard4th3h11 + 
      hInv23*(PDstandard4th1h32 - PDstandard4th2h31 + PDstandard4th3h21));
    
    CCTK_REAL gamma122 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv11*(-PDstandard4th1h22 + 2*PDstandard4th2h21) + 
      hInv12*PDstandard4th2h22 + hInv13*(2*PDstandard4th2h32 - 
      PDstandard4th3h22));
    
    CCTK_REAL gamma222 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv12*(-PDstandard4th1h22 + 2*PDstandard4th2h21) + 
      hInv22*PDstandard4th2h22 + hInv23*(2*PDstandard4th2h32 - 
      PDstandard4th3h22));
    
    CCTK_REAL gamma322 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv13*(-PDstandard4th1h22 + 2*PDstandard4th2h21) + 
      hInv23*PDstandard4th2h22 + hInv33*(2*PDstandard4th2h32 - 
      PDstandard4th3h22));
    
    CCTK_REAL gamma132 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv13*PDstandard4th2h33 + hInv11*(-PDstandard4th1h32 + 
      PDstandard4th2h31 + PDstandard4th3h21) + hInv12*PDstandard4th3h22);
    
    CCTK_REAL gamma232 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv23*PDstandard4th2h33 + hInv12*(-PDstandard4th1h32 + 
      PDstandard4th2h31 + PDstandard4th3h21) + hInv22*PDstandard4th3h22);
    
    CCTK_REAL gamma332 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv33*PDstandard4th2h33 + hInv13*(-PDstandard4th1h32 + 
      PDstandard4th2h31 + PDstandard4th3h21) + hInv23*PDstandard4th3h22);
    
    CCTK_REAL gamma133 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv11*(-PDstandard4th1h33 + 2*PDstandard4th3h31) + 
      hInv12*(-PDstandard4th2h33 + 2*PDstandard4th3h32) + 
      hInv13*PDstandard4th3h33);
    
    CCTK_REAL gamma233 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv12*(-PDstandard4th1h33 + 2*PDstandard4th3h31) + 
      hInv22*(-PDstandard4th2h33 + 2*PDstandard4th3h32) + 
      hInv23*PDstandard4th3h33);
    
    CCTK_REAL gamma333 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(hInv13*(-PDstandard4th1h33 + 2*PDstandard4th3h31) + 
      hInv23*(-PDstandard4th2h33 + 2*PDstandard4th3h32) + 
      hInv33*PDstandard4th3h33);
    
    CCTK_REAL Gam1L CCTK_ATTRIBUTE_UNUSED = gamma111*hInv11 + 
      gamma122*hInv22 + 2*(gamma121*hInv12 + gamma131*hInv13 + 
      gamma132*hInv23) + gamma133*hInv33;
    
    CCTK_REAL Gam2L CCTK_ATTRIBUTE_UNUSED = gamma211*hInv11 + 
      gamma222*hInv22 + 2*(gamma221*hInv12 + gamma231*hInv13 + 
      gamma232*hInv23) + gamma233*hInv33;
    
    CCTK_REAL Gam3L CCTK_ATTRIBUTE_UNUSED = gamma311*hInv11 + 
      gamma322*hInv22 + 2*(gamma321*hInv12 + gamma331*hInv13 + 
      gamma332*hInv23) + gamma333*hInv33;
    /* Copy local copies back to grid functions */
    Gam1[index] = Gam1L;
    Gam2[index] = Gam2L;
    Gam3[index] = Gam3L;
  }
  CCTK_ENDLOOP3(bssn_calc_evolved_gammas_4thcentred);
}
extern "C" void bssn_calc_evolved_gammas_4thcentred(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering bssn_calc_evolved_gammas_4thcentred_Body");
  }
  if (cctk_iteration % bssn_calc_evolved_gammas_4thcentred_calc_every != bssn_calc_evolved_gammas_4thcentred_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "Kranc2BSSN::beta_group",
    "Kranc2BSSN::Gam_group",
    "Kranc2BSSN::h_group"};
  AssertGroupStorage(cctkGH, "bssn_calc_evolved_gammas_4thcentred", 3, groups);
  
  switch (constraintsfdorder)
  {
    case 2:
    {
      EnsureStencilFits(cctkGH, "bssn_calc_evolved_gammas_4thcentred", 2, 2, 2);
      break;
    }
    
    case 4:
    {
      EnsureStencilFits(cctkGH, "bssn_calc_evolved_gammas_4thcentred", 2, 2, 2);
      break;
    }
    
    case 6:
    {
      EnsureStencilFits(cctkGH, "bssn_calc_evolved_gammas_4thcentred", 2, 2, 2);
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  LoopOverInterior(cctkGH, bssn_calc_evolved_gammas_4thcentred_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving bssn_calc_evolved_gammas_4thcentred_Body");
  }
}

} // namespace Kranc2BSSN
