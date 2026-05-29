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

extern "C" void bssnmatter_evolve_nonmetric_matter_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % bssnmatter_evolve_nonmetric_matter_calc_every != bssnmatter_evolve_nonmetric_matter_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "Kranc2BSSNMatter::A_grouprhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for Kranc2BSSNMatter::A_grouprhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "Kranc2BSSNMatter::Gam_grouprhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for Kranc2BSSNMatter::Gam_grouprhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "Kranc2BSSNMatter::K_grouprhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for Kranc2BSSNMatter::K_grouprhs.");
  return;
}

static void bssnmatter_evolve_nonmetric_matter_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  CCTK_LOOP3(bssnmatter_evolve_nonmetric_matter,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL A11rhsL CCTK_ATTRIBUTE_UNUSED = A11rhs[index];
    CCTK_REAL A21rhsL CCTK_ATTRIBUTE_UNUSED = A21rhs[index];
    CCTK_REAL A22rhsL CCTK_ATTRIBUTE_UNUSED = A22rhs[index];
    CCTK_REAL A31rhsL CCTK_ATTRIBUTE_UNUSED = A31rhs[index];
    CCTK_REAL A32rhsL CCTK_ATTRIBUTE_UNUSED = A32rhs[index];
    CCTK_REAL A33rhsL CCTK_ATTRIBUTE_UNUSED = A33rhs[index];
    CCTK_REAL alphaL CCTK_ATTRIBUTE_UNUSED = alpha[index];
    CCTK_REAL beta1L CCTK_ATTRIBUTE_UNUSED = beta1[index];
    CCTK_REAL beta2L CCTK_ATTRIBUTE_UNUSED = beta2[index];
    CCTK_REAL beta3L CCTK_ATTRIBUTE_UNUSED = beta3[index];
    CCTK_REAL Gam1rhsL CCTK_ATTRIBUTE_UNUSED = Gam1rhs[index];
    CCTK_REAL Gam2rhsL CCTK_ATTRIBUTE_UNUSED = Gam2rhs[index];
    CCTK_REAL Gam3rhsL CCTK_ATTRIBUTE_UNUSED = Gam3rhs[index];
    CCTK_REAL h11L CCTK_ATTRIBUTE_UNUSED = h11[index];
    CCTK_REAL h21L CCTK_ATTRIBUTE_UNUSED = h21[index];
    CCTK_REAL h22L CCTK_ATTRIBUTE_UNUSED = h22[index];
    CCTK_REAL h31L CCTK_ATTRIBUTE_UNUSED = h31[index];
    CCTK_REAL h32L CCTK_ATTRIBUTE_UNUSED = h32[index];
    CCTK_REAL h33L CCTK_ATTRIBUTE_UNUSED = h33[index];
    CCTK_REAL KrhsL CCTK_ATTRIBUTE_UNUSED = Krhs[index];
    CCTK_REAL phiL CCTK_ATTRIBUTE_UNUSED = phi[index];
    
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
    
    CCTK_REAL em4phi CCTK_ATTRIBUTE_UNUSED = exp(-4*phiL);
    
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
    
    CCTK_REAL S CCTK_ATTRIBUTE_UNUSED = em4phi*(hInv11*Sij11 + 
      hInv22*Sij22 + 2*(hInv12*Sij21 + hInv13*Sij31 + hInv23*Sij32) + 
      hInv33*Sij33);
    
    CCTK_REAL rho CCTK_ATTRIBUTE_UNUSED = pow(alphaL,-2)*(eTttL + 
      2*(beta2L*beta3L*Sij32 + beta1L*(beta2L*Sij21 + beta3L*Sij31 - Tj1)) - 
      2*(beta2L*Tj2 + beta3L*Tj3) + Sij11*pow(beta1L,2) + Sij22*pow(beta2L,2) 
      + Sij33*pow(beta3L,2));
    
    CCTK_REAL pi CCTK_ATTRIBUTE_UNUSED = 3.1415926535897932385;
    
    KrhsL = KrhsL + 4*alphaL*pi*(rho + S);
    
    A11rhsL = A11rhsL + 2.66666666666666666666666666667*alphaL*pi*(h11L*S 
      - 3*em4phi*Sij11);
    
    A21rhsL = A21rhsL + 2.66666666666666666666666666667*alphaL*pi*(h21L*S 
      - 3*em4phi*Sij21);
    
    A31rhsL = A31rhsL + 2.66666666666666666666666666667*alphaL*pi*(h31L*S 
      - 3*em4phi*Sij31);
    
    A22rhsL = A22rhsL + 2.66666666666666666666666666667*alphaL*pi*(h22L*S 
      - 3*em4phi*Sij22);
    
    A32rhsL = A32rhsL + 2.66666666666666666666666666667*alphaL*pi*(h32L*S 
      - 3*em4phi*Sij32);
    
    A33rhsL = A33rhsL + 2.66666666666666666666666666667*alphaL*pi*(h33L*S 
      - 3*em4phi*Sij33);
    
    Gam1rhsL = Gam1rhsL - 16*alphaL*pi*(hInv11*Si1 + hInv12*Si2 + 
      hInv13*Si3);
    
    Gam2rhsL = Gam2rhsL - 16*alphaL*pi*(hInv12*Si1 + hInv22*Si2 + 
      hInv23*Si3);
    
    Gam3rhsL = Gam3rhsL - 16*alphaL*pi*(hInv13*Si1 + hInv23*Si2 + 
      hInv33*Si3);
    /* Copy local copies back to grid functions */
    A11rhs[index] = A11rhsL;
    A21rhs[index] = A21rhsL;
    A22rhs[index] = A22rhsL;
    A31rhs[index] = A31rhsL;
    A32rhs[index] = A32rhsL;
    A33rhs[index] = A33rhsL;
    Gam1rhs[index] = Gam1rhsL;
    Gam2rhs[index] = Gam2rhsL;
    Gam3rhs[index] = Gam3rhsL;
    Krhs[index] = KrhsL;
  }
  CCTK_ENDLOOP3(bssnmatter_evolve_nonmetric_matter);
}
extern "C" void bssnmatter_evolve_nonmetric_matter(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering bssnmatter_evolve_nonmetric_matter_Body");
  }
  if (cctk_iteration % bssnmatter_evolve_nonmetric_matter_calc_every != bssnmatter_evolve_nonmetric_matter_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "Kranc2BSSNMatter::A_grouprhs",
    "Kranc2BSSNMatter::alpha_group",
    "Kranc2BSSNMatter::beta_group",
    "Kranc2BSSNMatter::Gam_grouprhs",
    "Kranc2BSSNMatter::h_group",
    "Kranc2BSSNMatter::K_grouprhs",
    "Kranc2BSSNMatter::phi_group"};
  AssertGroupStorage(cctkGH, "bssnmatter_evolve_nonmetric_matter", 7, groups);
  
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
  
  LoopOverInterior(cctkGH, bssnmatter_evolve_nonmetric_matter_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving bssnmatter_evolve_nonmetric_matter_Body");
  }
}

} // namespace Kranc2BSSNMatter
