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

extern "C" void bssnchimatter_boundary_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % bssnchimatter_boundary_calc_every != bssnchimatter_boundary_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "Kranc2BSSNChiMatter::A_grouprhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for Kranc2BSSNChiMatter::A_grouprhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "Kranc2BSSNChiMatter::alpha_grouprhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for Kranc2BSSNChiMatter::alpha_grouprhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "Kranc2BSSNChiMatter::beta_grouprhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for Kranc2BSSNChiMatter::beta_grouprhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "Kranc2BSSNChiMatter::betat_grouprhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for Kranc2BSSNChiMatter::betat_grouprhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "Kranc2BSSNChiMatter::bssnmom_group","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for Kranc2BSSNChiMatter::bssnmom_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "Kranc2BSSNChiMatter::chi_grouprhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for Kranc2BSSNChiMatter::chi_grouprhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "Kranc2BSSNChiMatter::Gam_grouprhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for Kranc2BSSNChiMatter::Gam_grouprhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "Kranc2BSSNChiMatter::h_grouprhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for Kranc2BSSNChiMatter::h_grouprhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "Kranc2BSSNChiMatter::K_grouprhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for Kranc2BSSNChiMatter::K_grouprhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "Kranc2BSSNChiMatter::scalarconstraints","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for Kranc2BSSNChiMatter::scalarconstraints.");
  return;
}

static void bssnchimatter_boundary_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  CCTK_LOOP3(bssnchimatter_boundary,
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
    CCTK_REAL betat1L CCTK_ATTRIBUTE_UNUSED = betat1[index];
    CCTK_REAL betat2L CCTK_ATTRIBUTE_UNUSED = betat2[index];
    CCTK_REAL betat3L CCTK_ATTRIBUTE_UNUSED = betat3[index];
    CCTK_REAL chiL CCTK_ATTRIBUTE_UNUSED = chi[index];
    CCTK_REAL h11L CCTK_ATTRIBUTE_UNUSED = h11[index];
    CCTK_REAL h21L CCTK_ATTRIBUTE_UNUSED = h21[index];
    CCTK_REAL h22L CCTK_ATTRIBUTE_UNUSED = h22[index];
    CCTK_REAL h31L CCTK_ATTRIBUTE_UNUSED = h31[index];
    CCTK_REAL h32L CCTK_ATTRIBUTE_UNUSED = h32[index];
    CCTK_REAL h33L CCTK_ATTRIBUTE_UNUSED = h33[index];
    CCTK_REAL KL CCTK_ATTRIBUTE_UNUSED = K[index];
    CCTK_REAL rL CCTK_ATTRIBUTE_UNUSED = r[index];
    CCTK_REAL xL CCTK_ATTRIBUTE_UNUSED = x[index];
    CCTK_REAL yL CCTK_ATTRIBUTE_UNUSED = y[index];
    CCTK_REAL zL CCTK_ATTRIBUTE_UNUSED = z[index];
    
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
    CCTK_REAL n1 CCTK_ATTRIBUTE_UNUSED = -(xL*pow(rL,-1));
    
    CCTK_REAL n2 CCTK_ATTRIBUTE_UNUSED = -(yL*pow(rL,-1));
    
    CCTK_REAL n3 CCTK_ATTRIBUTE_UNUSED = -(zL*pow(rL,-1));
    
    ptrdiff_t dir1 CCTK_ATTRIBUTE_UNUSED = isgn(n1);
    
    ptrdiff_t dir2 CCTK_ATTRIBUTE_UNUSED = isgn(n2);
    
    ptrdiff_t dir3 CCTK_ATTRIBUTE_UNUSED = isgn(n3);
    
    CCTK_REAL h11rhsL CCTK_ATTRIBUTE_UNUSED = 
      PDonesided2nd1(&h11[index])*n1 + PDonesided2nd2(&h11[index])*n2 + 
      PDonesided2nd3(&h11[index])*n3 + (1 - h11L)*pow(rL,-1);
    
    CCTK_REAL h21rhsL CCTK_ATTRIBUTE_UNUSED = 
      PDonesided2nd1(&h21[index])*n1 + PDonesided2nd2(&h21[index])*n2 + 
      PDonesided2nd3(&h21[index])*n3 - h21L*pow(rL,-1);
    
    CCTK_REAL h31rhsL CCTK_ATTRIBUTE_UNUSED = 
      PDonesided2nd1(&h31[index])*n1 + PDonesided2nd2(&h31[index])*n2 + 
      PDonesided2nd3(&h31[index])*n3 - h31L*pow(rL,-1);
    
    CCTK_REAL h22rhsL CCTK_ATTRIBUTE_UNUSED = 
      PDonesided2nd1(&h22[index])*n1 + PDonesided2nd2(&h22[index])*n2 + 
      PDonesided2nd3(&h22[index])*n3 + (1 - h22L)*pow(rL,-1);
    
    CCTK_REAL h32rhsL CCTK_ATTRIBUTE_UNUSED = 
      PDonesided2nd1(&h32[index])*n1 + PDonesided2nd2(&h32[index])*n2 + 
      PDonesided2nd3(&h32[index])*n3 - h32L*pow(rL,-1);
    
    CCTK_REAL h33rhsL CCTK_ATTRIBUTE_UNUSED = 
      PDonesided2nd1(&h33[index])*n1 + PDonesided2nd2(&h33[index])*n2 + 
      PDonesided2nd3(&h33[index])*n3 + (1 - h33L)*pow(rL,-1);
    
    CCTK_REAL A11rhsL CCTK_ATTRIBUTE_UNUSED = 
      PDonesided2nd1(&A11[index])*n1 + PDonesided2nd2(&A11[index])*n2 + 
      PDonesided2nd3(&A11[index])*n3 - A11L*pow(rL,-1);
    
    CCTK_REAL A21rhsL CCTK_ATTRIBUTE_UNUSED = 
      PDonesided2nd1(&A21[index])*n1 + PDonesided2nd2(&A21[index])*n2 + 
      PDonesided2nd3(&A21[index])*n3 - A21L*pow(rL,-1);
    
    CCTK_REAL A31rhsL CCTK_ATTRIBUTE_UNUSED = 
      PDonesided2nd1(&A31[index])*n1 + PDonesided2nd2(&A31[index])*n2 + 
      PDonesided2nd3(&A31[index])*n3 - A31L*pow(rL,-1);
    
    CCTK_REAL A22rhsL CCTK_ATTRIBUTE_UNUSED = 
      PDonesided2nd1(&A22[index])*n1 + PDonesided2nd2(&A22[index])*n2 + 
      PDonesided2nd3(&A22[index])*n3 - A22L*pow(rL,-1);
    
    CCTK_REAL A32rhsL CCTK_ATTRIBUTE_UNUSED = 
      PDonesided2nd1(&A32[index])*n1 + PDonesided2nd2(&A32[index])*n2 + 
      PDonesided2nd3(&A32[index])*n3 - A32L*pow(rL,-1);
    
    CCTK_REAL A33rhsL CCTK_ATTRIBUTE_UNUSED = 
      PDonesided2nd1(&A33[index])*n1 + PDonesided2nd2(&A33[index])*n2 + 
      PDonesided2nd3(&A33[index])*n3 - A33L*pow(rL,-1);
    
    CCTK_REAL KrhsL CCTK_ATTRIBUTE_UNUSED = PDonesided2nd1(&K[index])*n1 + 
      PDonesided2nd2(&K[index])*n2 + PDonesided2nd3(&K[index])*n3 - 
      KL*pow(rL,-1);
    
    CCTK_REAL Gam1rhsL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL Gam2rhsL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL Gam3rhsL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL chirhsL CCTK_ATTRIBUTE_UNUSED = 
      PDonesided2nd1(&chi[index])*n1 + PDonesided2nd2(&chi[index])*n2 + 
      PDonesided2nd3(&chi[index])*n3 + (1 - chiL)*pow(rL,-1);
    
    CCTK_REAL alpharhsL CCTK_ATTRIBUTE_UNUSED = 
      PDonesided2nd1(&alpha[index])*n1 + PDonesided2nd2(&alpha[index])*n2 + 
      PDonesided2nd3(&alpha[index])*n3 + (1 - alphaL)*pow(rL,-1);
    
    CCTK_REAL beta1rhsL CCTK_ATTRIBUTE_UNUSED = 
      PDonesided2nd1(&beta1[index])*n1 + PDonesided2nd2(&beta1[index])*n2 + 
      PDonesided2nd3(&beta1[index])*n3 - beta1L*pow(rL,-1);
    
    CCTK_REAL beta2rhsL CCTK_ATTRIBUTE_UNUSED = 
      PDonesided2nd1(&beta2[index])*n1 + PDonesided2nd2(&beta2[index])*n2 + 
      PDonesided2nd3(&beta2[index])*n3 - beta2L*pow(rL,-1);
    
    CCTK_REAL beta3rhsL CCTK_ATTRIBUTE_UNUSED = 
      PDonesided2nd1(&beta3[index])*n1 + PDonesided2nd2(&beta3[index])*n2 + 
      PDonesided2nd3(&beta3[index])*n3 - beta3L*pow(rL,-1);
    
    CCTK_REAL betat1rhsL CCTK_ATTRIBUTE_UNUSED = 
      PDonesided2nd1(&betat1[index])*n1 + PDonesided2nd2(&betat1[index])*n2 + 
      PDonesided2nd3(&betat1[index])*n3 - betat1L*pow(rL,-1);
    
    CCTK_REAL betat2rhsL CCTK_ATTRIBUTE_UNUSED = 
      PDonesided2nd1(&betat2[index])*n1 + PDonesided2nd2(&betat2[index])*n2 + 
      PDonesided2nd3(&betat2[index])*n3 - betat2L*pow(rL,-1);
    
    CCTK_REAL betat3rhsL CCTK_ATTRIBUTE_UNUSED = 
      PDonesided2nd1(&betat3[index])*n1 + PDonesided2nd2(&betat3[index])*n2 + 
      PDonesided2nd3(&betat3[index])*n3 - betat3L*pow(rL,-1);
    
    CCTK_REAL bssnhamL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL bssnDL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL bssnTL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL bssnmom1L CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL bssnmom2L CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL bssnmom3L CCTK_ATTRIBUTE_UNUSED = 0;
    /* Copy local copies back to grid functions */
    A11rhs[index] = A11rhsL;
    A21rhs[index] = A21rhsL;
    A22rhs[index] = A22rhsL;
    A31rhs[index] = A31rhsL;
    A32rhs[index] = A32rhsL;
    A33rhs[index] = A33rhsL;
    alpharhs[index] = alpharhsL;
    beta1rhs[index] = beta1rhsL;
    beta2rhs[index] = beta2rhsL;
    beta3rhs[index] = beta3rhsL;
    betat1rhs[index] = betat1rhsL;
    betat2rhs[index] = betat2rhsL;
    betat3rhs[index] = betat3rhsL;
    bssnD[index] = bssnDL;
    bssnham[index] = bssnhamL;
    bssnmom1[index] = bssnmom1L;
    bssnmom2[index] = bssnmom2L;
    bssnmom3[index] = bssnmom3L;
    bssnT[index] = bssnTL;
    chirhs[index] = chirhsL;
    Gam1rhs[index] = Gam1rhsL;
    Gam2rhs[index] = Gam2rhsL;
    Gam3rhs[index] = Gam3rhsL;
    h11rhs[index] = h11rhsL;
    h21rhs[index] = h21rhsL;
    h22rhs[index] = h22rhsL;
    h31rhs[index] = h31rhsL;
    h32rhs[index] = h32rhsL;
    h33rhs[index] = h33rhsL;
    Krhs[index] = KrhsL;
  }
  CCTK_ENDLOOP3(bssnchimatter_boundary);
}
extern "C" void bssnchimatter_boundary(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering bssnchimatter_boundary_Body");
  }
  if (cctk_iteration % bssnchimatter_boundary_calc_every != bssnchimatter_boundary_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "Kranc2BSSNChiMatter::A_group",
    "Kranc2BSSNChiMatter::A_grouprhs",
    "Kranc2BSSNChiMatter::alpha_group",
    "Kranc2BSSNChiMatter::alpha_grouprhs",
    "Kranc2BSSNChiMatter::beta_group",
    "Kranc2BSSNChiMatter::beta_grouprhs",
    "Kranc2BSSNChiMatter::betat_group",
    "Kranc2BSSNChiMatter::betat_grouprhs",
    "Kranc2BSSNChiMatter::bssnmom_group",
    "Kranc2BSSNChiMatter::chi_group",
    "Kranc2BSSNChiMatter::chi_grouprhs",
    "Kranc2BSSNChiMatter::Gam_grouprhs",
    "grid::coordinates",
    "Kranc2BSSNChiMatter::h_group",
    "Kranc2BSSNChiMatter::h_grouprhs",
    "Kranc2BSSNChiMatter::K_group",
    "Kranc2BSSNChiMatter::K_grouprhs",
    "Kranc2BSSNChiMatter::scalarconstraints"};
  AssertGroupStorage(cctkGH, "bssnchimatter_boundary", 18, groups);
  
  switch (constraintsfdorder)
  {
    case 2:
    {
      EnsureStencilFits(cctkGH, "bssnchimatter_boundary", 2, 2, 2);
      break;
    }
    
    case 4:
    {
      EnsureStencilFits(cctkGH, "bssnchimatter_boundary", 2, 2, 2);
      break;
    }
    
    case 6:
    {
      EnsureStencilFits(cctkGH, "bssnchimatter_boundary", 2, 2, 2);
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  LoopOverBoundary(cctkGH, bssnchimatter_boundary_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving bssnchimatter_boundary_Body");
  }
}

} // namespace Kranc2BSSNChiMatter
