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

extern "C" void EvolveScalarFields_boundaryPhiOfT_SelectBCs(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_EvolveScalarFields_boundaryPhiOfT_SelectBCs
  DECLARE_CCTK_ARGUMENTS_CHECKED(EvolveScalarFields_boundaryPhiOfT_SelectBCs);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % EvolveScalarFields_boundaryPhiOfT_calc_every != EvolveScalarFields_boundaryPhiOfT_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "EvolveScalarFields::phia_grouprhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for EvolveScalarFields::phia_grouprhs.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "EvolveScalarFields::phib_grouprhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for EvolveScalarFields::phib_grouprhs.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "EvolveScalarFields::pia_grouprhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for EvolveScalarFields::pia_grouprhs.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "EvolveScalarFields::pib_grouprhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for EvolveScalarFields::pib_grouprhs.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "EvolveScalarFields::rhoE_group","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for EvolveScalarFields::rhoE_group.");
  return;
}

static void EvolveScalarFields_boundaryPhiOfT_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  CCTK_LOOP3(EvolveScalarFields_boundaryPhiOfT,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL phiaL CCTK_ATTRIBUTE_UNUSED = phia[index];
    CCTK_REAL phibL CCTK_ATTRIBUTE_UNUSED = phib[index];
    CCTK_REAL piaL CCTK_ATTRIBUTE_UNUSED = pia[index];
    CCTK_REAL pibL CCTK_ATTRIBUTE_UNUSED = pib[index];
    CCTK_REAL rL CCTK_ATTRIBUTE_UNUSED = r[index];
    CCTK_REAL xL CCTK_ATTRIBUTE_UNUSED = x[index];
    CCTK_REAL yL CCTK_ATTRIBUTE_UNUSED = y[index];
    CCTK_REAL zL CCTK_ATTRIBUTE_UNUSED = z[index];
    
    /* Include user supplied include files */
    /* Precompute derivatives */
    /* Calculate temporaries and grid functions */
    CCTK_REAL n1 CCTK_ATTRIBUTE_UNUSED = -(xL*pow(rL,-1));
    
    CCTK_REAL n2 CCTK_ATTRIBUTE_UNUSED = -(yL*pow(rL,-1));
    
    CCTK_REAL n3 CCTK_ATTRIBUTE_UNUSED = -(zL*pow(rL,-1));
    
    ptrdiff_t dir1 CCTK_ATTRIBUTE_UNUSED = isgn(n1);
    
    ptrdiff_t dir2 CCTK_ATTRIBUTE_UNUSED = isgn(n2);
    
    ptrdiff_t dir3 CCTK_ATTRIBUTE_UNUSED = isgn(n3);
    
    CCTK_REAL piarhsL CCTK_ATTRIBUTE_UNUSED = 
      PDonesided2nd1(&pia[index])*n1 + PDonesided2nd2(&pia[index])*n2 + 
      PDonesided2nd3(&pia[index])*n3 + (-piaL + piaBG)*pow(rL,-1);
    
    CCTK_REAL phiarhsL CCTK_ATTRIBUTE_UNUSED = 
      PDonesided2nd1(&phia[index])*n1 + PDonesided2nd2(&phia[index])*n2 + 
      PDonesided2nd3(&phia[index])*n3 + (-phiaL + piaBG*t)*pow(rL,-1);
    
    CCTK_REAL pibrhsL CCTK_ATTRIBUTE_UNUSED = 
      PDonesided2nd1(&pib[index])*n1 + PDonesided2nd2(&pib[index])*n2 + 
      PDonesided2nd3(&pib[index])*n3 + (-pibL + pibBG)*pow(rL,-1);
    
    CCTK_REAL phibrhsL CCTK_ATTRIBUTE_UNUSED = 
      PDonesided2nd1(&phib[index])*n1 + PDonesided2nd2(&phib[index])*n2 + 
      PDonesided2nd3(&phib[index])*n3 + (-phibL + pibBG*t)*pow(rL,-1);
    
    CCTK_REAL rhoEL CCTK_ATTRIBUTE_UNUSED = 0;
    /* Copy local copies back to grid functions */
    phiarhs[index] = phiarhsL;
    phibrhs[index] = phibrhsL;
    piarhs[index] = piarhsL;
    pibrhs[index] = pibrhsL;
    rhoE[index] = rhoEL;
  }
  CCTK_ENDLOOP3(EvolveScalarFields_boundaryPhiOfT);
}
extern "C" void EvolveScalarFields_boundaryPhiOfT(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_EvolveScalarFields_boundaryPhiOfT
  DECLARE_CCTK_ARGUMENTS_CHECKED(EvolveScalarFields_boundaryPhiOfT);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering EvolveScalarFields_boundaryPhiOfT_Body");
  }
  if (cctk_iteration % EvolveScalarFields_boundaryPhiOfT_calc_every != EvolveScalarFields_boundaryPhiOfT_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "grid::coordinates",
    "EvolveScalarFields::phia_group",
    "EvolveScalarFields::phia_grouprhs",
    "EvolveScalarFields::phib_group",
    "EvolveScalarFields::phib_grouprhs",
    "EvolveScalarFields::pia_group",
    "EvolveScalarFields::pia_grouprhs",
    "EvolveScalarFields::pib_group",
    "EvolveScalarFields::pib_grouprhs",
    "EvolveScalarFields::rhoE_group"};
  AssertGroupStorage(cctkGH, "EvolveScalarFields_boundaryPhiOfT", 10, groups);
  
  EnsureStencilFits(cctkGH, "EvolveScalarFields_boundaryPhiOfT", 2, 2, 2);
  
  LoopOverBoundary(cctkGH, EvolveScalarFields_boundaryPhiOfT_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving EvolveScalarFields_boundaryPhiOfT_Body");
  }
}

} // namespace EvolveScalarFields
