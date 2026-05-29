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

extern "C" void EvolveScalarFields_EvolveScalarField_6th_potential_onlymass_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % EvolveScalarFields_EvolveScalarField_6th_potential_onlymass_calc_every != EvolveScalarFields_EvolveScalarField_6th_potential_onlymass_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "EvolveScalarFields::phia_grouprhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for EvolveScalarFields::phia_grouprhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "EvolveScalarFields::phib_grouprhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for EvolveScalarFields::phib_grouprhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "EvolveScalarFields::pia_grouprhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for EvolveScalarFields::pia_grouprhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "EvolveScalarFields::pib_grouprhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for EvolveScalarFields::pib_grouprhs.");
  return;
}

static void EvolveScalarFields_EvolveScalarField_6th_potential_onlymass_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  CCTK_LOOP3(EvolveScalarFields_EvolveScalarField_6th_potential_onlymass,
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
    CCTK_REAL phiaL CCTK_ATTRIBUTE_UNUSED = phia[index];
    CCTK_REAL phibL CCTK_ATTRIBUTE_UNUSED = phib[index];
    CCTK_REAL piaL CCTK_ATTRIBUTE_UNUSED = pia[index];
    CCTK_REAL pibL CCTK_ATTRIBUTE_UNUSED = pib[index];
    
    /* Include user supplied include files */
    /* Precompute derivatives */
    const CCTK_REAL PDstandard6th1alp CCTK_ATTRIBUTE_UNUSED = PDstandard6th1(&alp[index]);
    const CCTK_REAL PDstandard6th2alp CCTK_ATTRIBUTE_UNUSED = PDstandard6th2(&alp[index]);
    const CCTK_REAL PDstandard6th3alp CCTK_ATTRIBUTE_UNUSED = PDstandard6th3(&alp[index]);
    const CCTK_REAL PDstandard6th1gxx CCTK_ATTRIBUTE_UNUSED = PDstandard6th1(&gxx[index]);
    const CCTK_REAL PDstandard6th2gxx CCTK_ATTRIBUTE_UNUSED = PDstandard6th2(&gxx[index]);
    const CCTK_REAL PDstandard6th3gxx CCTK_ATTRIBUTE_UNUSED = PDstandard6th3(&gxx[index]);
    const CCTK_REAL PDstandard6th1gxy CCTK_ATTRIBUTE_UNUSED = PDstandard6th1(&gxy[index]);
    const CCTK_REAL PDstandard6th2gxy CCTK_ATTRIBUTE_UNUSED = PDstandard6th2(&gxy[index]);
    const CCTK_REAL PDstandard6th3gxy CCTK_ATTRIBUTE_UNUSED = PDstandard6th3(&gxy[index]);
    const CCTK_REAL PDstandard6th1gxz CCTK_ATTRIBUTE_UNUSED = PDstandard6th1(&gxz[index]);
    const CCTK_REAL PDstandard6th2gxz CCTK_ATTRIBUTE_UNUSED = PDstandard6th2(&gxz[index]);
    const CCTK_REAL PDstandard6th3gxz CCTK_ATTRIBUTE_UNUSED = PDstandard6th3(&gxz[index]);
    const CCTK_REAL PDstandard6th1gyy CCTK_ATTRIBUTE_UNUSED = PDstandard6th1(&gyy[index]);
    const CCTK_REAL PDstandard6th2gyy CCTK_ATTRIBUTE_UNUSED = PDstandard6th2(&gyy[index]);
    const CCTK_REAL PDstandard6th3gyy CCTK_ATTRIBUTE_UNUSED = PDstandard6th3(&gyy[index]);
    const CCTK_REAL PDstandard6th1gyz CCTK_ATTRIBUTE_UNUSED = PDstandard6th1(&gyz[index]);
    const CCTK_REAL PDstandard6th2gyz CCTK_ATTRIBUTE_UNUSED = PDstandard6th2(&gyz[index]);
    const CCTK_REAL PDstandard6th3gyz CCTK_ATTRIBUTE_UNUSED = PDstandard6th3(&gyz[index]);
    const CCTK_REAL PDstandard6th1gzz CCTK_ATTRIBUTE_UNUSED = PDstandard6th1(&gzz[index]);
    const CCTK_REAL PDstandard6th2gzz CCTK_ATTRIBUTE_UNUSED = PDstandard6th2(&gzz[index]);
    const CCTK_REAL PDstandard6th3gzz CCTK_ATTRIBUTE_UNUSED = PDstandard6th3(&gzz[index]);
    const CCTK_REAL PDstandard6th1phia CCTK_ATTRIBUTE_UNUSED = PDstandard6th1(&phia[index]);
    const CCTK_REAL PDstandard6th2phia CCTK_ATTRIBUTE_UNUSED = PDstandard6th2(&phia[index]);
    const CCTK_REAL PDstandard6th3phia CCTK_ATTRIBUTE_UNUSED = PDstandard6th3(&phia[index]);
    const CCTK_REAL PDstandard6th11phia CCTK_ATTRIBUTE_UNUSED = PDstandard6th11(&phia[index]);
    const CCTK_REAL PDstandard6th22phia CCTK_ATTRIBUTE_UNUSED = PDstandard6th22(&phia[index]);
    const CCTK_REAL PDstandard6th33phia CCTK_ATTRIBUTE_UNUSED = PDstandard6th33(&phia[index]);
    const CCTK_REAL PDstandard6th12phia CCTK_ATTRIBUTE_UNUSED = PDstandard6th12(&phia[index]);
    const CCTK_REAL PDstandard6th13phia CCTK_ATTRIBUTE_UNUSED = PDstandard6th13(&phia[index]);
    const CCTK_REAL PDstandard6th23phia CCTK_ATTRIBUTE_UNUSED = PDstandard6th23(&phia[index]);
    const CCTK_REAL PDstandard6th1phib CCTK_ATTRIBUTE_UNUSED = PDstandard6th1(&phib[index]);
    const CCTK_REAL PDstandard6th2phib CCTK_ATTRIBUTE_UNUSED = PDstandard6th2(&phib[index]);
    const CCTK_REAL PDstandard6th3phib CCTK_ATTRIBUTE_UNUSED = PDstandard6th3(&phib[index]);
    const CCTK_REAL PDstandard6th11phib CCTK_ATTRIBUTE_UNUSED = PDstandard6th11(&phib[index]);
    const CCTK_REAL PDstandard6th22phib CCTK_ATTRIBUTE_UNUSED = PDstandard6th22(&phib[index]);
    const CCTK_REAL PDstandard6th33phib CCTK_ATTRIBUTE_UNUSED = PDstandard6th33(&phib[index]);
    const CCTK_REAL PDstandard6th12phib CCTK_ATTRIBUTE_UNUSED = PDstandard6th12(&phib[index]);
    const CCTK_REAL PDstandard6th13phib CCTK_ATTRIBUTE_UNUSED = PDstandard6th13(&phib[index]);
    const CCTK_REAL PDstandard6th23phib CCTK_ATTRIBUTE_UNUSED = PDstandard6th23(&phib[index]);
    /* Calculate temporaries and grid functions */
    ptrdiff_t dir1 CCTK_ATTRIBUTE_UNUSED = isgn(betaxL);
    
    ptrdiff_t dir2 CCTK_ATTRIBUTE_UNUSED = isgn(betayL);
    
    ptrdiff_t dir3 CCTK_ATTRIBUTE_UNUSED = isgn(betazL);
    
    CCTK_REAL gInv11 CCTK_ATTRIBUTE_UNUSED = (-(gyyL*gzzL) + 
      pow(gyzL,2))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*(-(gxxL*gzzL) + pow(gxzL,2)) + gxxL*pow(gyzL,2),-1);
    
    CCTK_REAL gInv12 CCTK_ATTRIBUTE_UNUSED = (-(gxzL*gyzL) + 
      gxyL*gzzL)*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + gyyL*pow(gxzL,2) 
      + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-1);
    
    CCTK_REAL gInv13 CCTK_ATTRIBUTE_UNUSED = (gxzL*gyyL - 
      gxyL*gyzL)*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*(-(gxxL*gzzL) + pow(gxzL,2)) + gxxL*pow(gyzL,2),-1);
    
    CCTK_REAL gInv22 CCTK_ATTRIBUTE_UNUSED = (-(gxxL*gzzL) + 
      pow(gxzL,2))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*(-(gxxL*gzzL) + pow(gxzL,2)) + gxxL*pow(gyzL,2),-1);
    
    CCTK_REAL gInv23 CCTK_ATTRIBUTE_UNUSED = (-(gxyL*gxzL) + 
      gxxL*gyzL)*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + gyyL*pow(gxzL,2) 
      + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-1);
    
    CCTK_REAL gInv33 CCTK_ATTRIBUTE_UNUSED = (-(gxxL*gyyL) + 
      pow(gxyL,2))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*(-(gxxL*gzzL) + pow(gxzL,2)) + gxxL*pow(gyzL,2),-1);
    
    CCTK_REAL ktrace CCTK_ATTRIBUTE_UNUSED = kxxL*gInv11 + kyyL*gInv22 + 
      2*(kxyL*gInv12 + kxzL*gInv13 + kyzL*gInv23) + kzzL*gInv33;
    
    CCTK_REAL Gam3111 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv11*PDstandard6th1gxx + 2*(gInv12*PDstandard6th1gxy + 
      gInv13*PDstandard6th1gxz) - gInv12*PDstandard6th2gxx - 
      gInv13*PDstandard6th3gxx);
    
    CCTK_REAL Gam3211 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv12*PDstandard6th1gxx + 2*(gInv22*PDstandard6th1gxy + 
      gInv23*PDstandard6th1gxz) - gInv22*PDstandard6th2gxx - 
      gInv23*PDstandard6th3gxx);
    
    CCTK_REAL Gam3311 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv13*PDstandard6th1gxx + 2*(gInv23*PDstandard6th1gxy + 
      gInv33*PDstandard6th1gxz) - gInv23*PDstandard6th2gxx - 
      gInv33*PDstandard6th3gxx);
    
    CCTK_REAL Gam3121 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv12*PDstandard6th1gyy + gInv11*PDstandard6th2gxx + 
      gInv13*(PDstandard6th1gyz + PDstandard6th2gxz - PDstandard6th3gxy));
    
    CCTK_REAL Gam3221 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv22*PDstandard6th1gyy + gInv12*PDstandard6th2gxx + 
      gInv23*(PDstandard6th1gyz + PDstandard6th2gxz - PDstandard6th3gxy));
    
    CCTK_REAL Gam3321 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv23*PDstandard6th1gyy + gInv13*PDstandard6th2gxx + 
      gInv33*(PDstandard6th1gyz + PDstandard6th2gxz - PDstandard6th3gxy));
    
    CCTK_REAL Gam3131 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv13*PDstandard6th1gzz + gInv11*PDstandard6th3gxx + 
      gInv12*(PDstandard6th1gyz - PDstandard6th2gxz + PDstandard6th3gxy));
    
    CCTK_REAL Gam3231 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv23*PDstandard6th1gzz + gInv12*PDstandard6th3gxx + 
      gInv22*(PDstandard6th1gyz - PDstandard6th2gxz + PDstandard6th3gxy));
    
    CCTK_REAL Gam3331 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv33*PDstandard6th1gzz + gInv13*PDstandard6th3gxx + 
      gInv23*(PDstandard6th1gyz - PDstandard6th2gxz + PDstandard6th3gxy));
    
    CCTK_REAL Gam3122 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv11*(-PDstandard6th1gyy + 2*PDstandard6th2gxy) + 
      gInv12*PDstandard6th2gyy + gInv13*(2*PDstandard6th2gyz - 
      PDstandard6th3gyy));
    
    CCTK_REAL Gam3222 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv12*(-PDstandard6th1gyy + 2*PDstandard6th2gxy) + 
      gInv22*PDstandard6th2gyy + gInv23*(2*PDstandard6th2gyz - 
      PDstandard6th3gyy));
    
    CCTK_REAL Gam3322 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv13*(-PDstandard6th1gyy + 2*PDstandard6th2gxy) + 
      gInv23*PDstandard6th2gyy + gInv33*(2*PDstandard6th2gyz - 
      PDstandard6th3gyy));
    
    CCTK_REAL Gam3132 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv13*PDstandard6th2gzz + gInv11*(-PDstandard6th1gyz + 
      PDstandard6th2gxz + PDstandard6th3gxy) + gInv12*PDstandard6th3gyy);
    
    CCTK_REAL Gam3232 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv23*PDstandard6th2gzz + gInv12*(-PDstandard6th1gyz + 
      PDstandard6th2gxz + PDstandard6th3gxy) + gInv22*PDstandard6th3gyy);
    
    CCTK_REAL Gam3332 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv33*PDstandard6th2gzz + gInv13*(-PDstandard6th1gyz + 
      PDstandard6th2gxz + PDstandard6th3gxy) + gInv23*PDstandard6th3gyy);
    
    CCTK_REAL Gam3133 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv11*(-PDstandard6th1gzz + 2*PDstandard6th3gxz) + 
      gInv12*(-PDstandard6th2gzz + 2*PDstandard6th3gyz) + 
      gInv13*PDstandard6th3gzz);
    
    CCTK_REAL Gam3233 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv12*(-PDstandard6th1gzz + 2*PDstandard6th3gxz) + 
      gInv22*(-PDstandard6th2gzz + 2*PDstandard6th3gyz) + 
      gInv23*PDstandard6th3gzz);
    
    CCTK_REAL Gam3333 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv13*(-PDstandard6th1gzz + 2*PDstandard6th3gxz) + 
      gInv23*(-PDstandard6th2gzz + 2*PDstandard6th3gyz) + 
      gInv33*PDstandard6th3gzz);
    
    CCTK_REAL derivVa CCTK_ATTRIBUTE_UNUSED = 
      phiaL*controlmasspositiveornegative*pow(phiamass,2);
    
    CCTK_REAL derivVb CCTK_ATTRIBUTE_UNUSED = 
      phibL*controlmasspositiveornegative*pow(phibmass,2);
    
    CCTK_REAL phiarhsL CCTK_ATTRIBUTE_UNUSED = 
      betaxL*PDlopsided6th1(&phia[index]) + 
      betayL*PDlopsided6th2(&phia[index]) + 
      betazL*PDlopsided6th3(&phia[index]) - alpL*piaL;
    
    CCTK_REAL piarhsL CCTK_ATTRIBUTE_UNUSED = 
      betaxL*PDlopsided6th1(&pia[index]) + betayL*PDlopsided6th2(&pia[index]) 
      + betazL*PDlopsided6th3(&pia[index]) - 
      PDstandard6th1alp*(gInv11*PDstandard6th1phia + 
      gInv12*PDstandard6th2phia + gInv13*PDstandard6th3phia) - 
      PDstandard6th2alp*(gInv12*PDstandard6th1phia + 
      gInv22*PDstandard6th2phia + gInv23*PDstandard6th3phia) - 
      PDstandard6th3alp*(gInv13*PDstandard6th1phia + 
      gInv23*PDstandard6th2phia + gInv33*PDstandard6th3phia) + alpL*(derivVa 
      + piaL*ktrace - gInv11*PDstandard6th11phia - 
      2*gInv12*PDstandard6th12phia - 2*gInv13*PDstandard6th13phia + 
      (Gam3111*gInv11 + 2*Gam3121*gInv12 + 2*Gam3131*gInv13 + Gam3122*gInv22 
      + 2*Gam3132*gInv23 + Gam3133*gInv33)*PDstandard6th1phia - 
      gInv22*PDstandard6th22phia - 2*gInv23*PDstandard6th23phia + 
      (Gam3211*gInv11 + 2*Gam3221*gInv12 + 2*Gam3231*gInv13 + Gam3222*gInv22 
      + 2*Gam3232*gInv23 + Gam3233*gInv33)*PDstandard6th2phia - 
      gInv33*PDstandard6th33phia + (Gam3311*gInv11 + 2*Gam3321*gInv12 + 
      2*Gam3331*gInv13 + Gam3322*gInv22 + 2*Gam3332*gInv23 + 
      Gam3333*gInv33)*PDstandard6th3phia);
    
    CCTK_REAL phibrhsL CCTK_ATTRIBUTE_UNUSED = 
      betaxL*PDlopsided6th1(&phib[index]) + 
      betayL*PDlopsided6th2(&phib[index]) + 
      betazL*PDlopsided6th3(&phib[index]) - alpL*pibL;
    
    CCTK_REAL pibrhsL CCTK_ATTRIBUTE_UNUSED = 
      betaxL*PDlopsided6th1(&pib[index]) + betayL*PDlopsided6th2(&pib[index]) 
      + betazL*PDlopsided6th3(&pib[index]) - 
      PDstandard6th1alp*(gInv11*PDstandard6th1phib + 
      gInv12*PDstandard6th2phib + gInv13*PDstandard6th3phib) - 
      PDstandard6th2alp*(gInv12*PDstandard6th1phib + 
      gInv22*PDstandard6th2phib + gInv23*PDstandard6th3phib) - 
      PDstandard6th3alp*(gInv13*PDstandard6th1phib + 
      gInv23*PDstandard6th2phib + gInv33*PDstandard6th3phib) + alpL*(derivVb 
      + pibL*ktrace - gInv11*PDstandard6th11phib - 
      2*gInv12*PDstandard6th12phib - 2*gInv13*PDstandard6th13phib + 
      (Gam3111*gInv11 + 2*Gam3121*gInv12 + 2*Gam3131*gInv13 + Gam3122*gInv22 
      + 2*Gam3132*gInv23 + Gam3133*gInv33)*PDstandard6th1phib - 
      gInv22*PDstandard6th22phib - 2*gInv23*PDstandard6th23phib + 
      (Gam3211*gInv11 + 2*Gam3221*gInv12 + 2*Gam3231*gInv13 + Gam3222*gInv22 
      + 2*Gam3232*gInv23 + Gam3233*gInv33)*PDstandard6th2phib - 
      gInv33*PDstandard6th33phib + (Gam3311*gInv11 + 2*Gam3321*gInv12 + 
      2*Gam3331*gInv13 + Gam3322*gInv22 + 2*Gam3332*gInv23 + 
      Gam3333*gInv33)*PDstandard6th3phib);
    /* Copy local copies back to grid functions */
    phiarhs[index] = phiarhsL;
    phibrhs[index] = phibrhsL;
    piarhs[index] = piarhsL;
    pibrhs[index] = pibrhsL;
  }
  CCTK_ENDLOOP3(EvolveScalarFields_EvolveScalarField_6th_potential_onlymass);
}
extern "C" void EvolveScalarFields_EvolveScalarField_6th_potential_onlymass(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering EvolveScalarFields_EvolveScalarField_6th_potential_onlymass_Body");
  }
  if (cctk_iteration % EvolveScalarFields_EvolveScalarField_6th_potential_onlymass_calc_every != EvolveScalarFields_EvolveScalarField_6th_potential_onlymass_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "admbase::curv",
    "admbase::lapse",
    "admbase::metric",
    "admbase::shift",
    "EvolveScalarFields::phia_group",
    "EvolveScalarFields::phia_grouprhs",
    "EvolveScalarFields::phib_group",
    "EvolveScalarFields::phib_grouprhs",
    "EvolveScalarFields::pia_group",
    "EvolveScalarFields::pia_grouprhs",
    "EvolveScalarFields::pib_group",
    "EvolveScalarFields::pib_grouprhs"};
  AssertGroupStorage(cctkGH, "EvolveScalarFields_EvolveScalarField_6th_potential_onlymass", 12, groups);
  
  EnsureStencilFits(cctkGH, "EvolveScalarFields_EvolveScalarField_6th_potential_onlymass", 4, 4, 4);
  
  LoopOverInterior(cctkGH, EvolveScalarFields_EvolveScalarField_6th_potential_onlymass_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving EvolveScalarFields_EvolveScalarField_6th_potential_onlymass_Body");
  }
}

} // namespace EvolveScalarFields
