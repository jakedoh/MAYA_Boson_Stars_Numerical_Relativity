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

namespace ScalarTensor {

extern "C" void ScalarTensor_EvolveScalarField_6th_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % ScalarTensor_EvolveScalarField_6th_calc_every != ScalarTensor_EvolveScalarField_6th_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ScalarTensor::STphi_grouprhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ScalarTensor::STphi_grouprhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ScalarTensor::STpi_grouprhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ScalarTensor::STpi_grouprhs.");
  return;
}

static void ScalarTensor_EvolveScalarField_6th_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  CCTK_LOOP3(ScalarTensor_EvolveScalarField_6th,
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
    CCTK_REAL STphiL CCTK_ATTRIBUTE_UNUSED = STphi[index];
    CCTK_REAL STpiL CCTK_ATTRIBUTE_UNUSED = STpi[index];
    
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
    const CCTK_REAL PDstandard6th1STphi CCTK_ATTRIBUTE_UNUSED = PDstandard6th1(&STphi[index]);
    const CCTK_REAL PDstandard6th2STphi CCTK_ATTRIBUTE_UNUSED = PDstandard6th2(&STphi[index]);
    const CCTK_REAL PDstandard6th3STphi CCTK_ATTRIBUTE_UNUSED = PDstandard6th3(&STphi[index]);
    const CCTK_REAL PDstandard6th11STphi CCTK_ATTRIBUTE_UNUSED = PDstandard6th11(&STphi[index]);
    const CCTK_REAL PDstandard6th22STphi CCTK_ATTRIBUTE_UNUSED = PDstandard6th22(&STphi[index]);
    const CCTK_REAL PDstandard6th33STphi CCTK_ATTRIBUTE_UNUSED = PDstandard6th33(&STphi[index]);
    const CCTK_REAL PDstandard6th12STphi CCTK_ATTRIBUTE_UNUSED = PDstandard6th12(&STphi[index]);
    const CCTK_REAL PDstandard6th13STphi CCTK_ATTRIBUTE_UNUSED = PDstandard6th13(&STphi[index]);
    const CCTK_REAL PDstandard6th23STphi CCTK_ATTRIBUTE_UNUSED = PDstandard6th23(&STphi[index]);
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
    
    CCTK_REAL ConfFactorA2 CCTK_ATTRIBUTE_UNUSED = exp((2*STphiL - 
      2*phiBG)*STalpha0 + STBeta*(1.*pow(STphiL,2) - 1.*pow(phiBG,2)));
    
    CCTK_REAL TttLocal CCTK_ATTRIBUTE_UNUSED = eTttL*ConfFactorA2;
    
    CCTK_REAL TtxLocal CCTK_ATTRIBUTE_UNUSED = eTtxL*ConfFactorA2;
    
    CCTK_REAL TtyLocal CCTK_ATTRIBUTE_UNUSED = eTtyL*ConfFactorA2;
    
    CCTK_REAL TtzLocal CCTK_ATTRIBUTE_UNUSED = eTtzL*ConfFactorA2;
    
    CCTK_REAL TmnLocal11 CCTK_ATTRIBUTE_UNUSED = eTxxL*ConfFactorA2;
    
    CCTK_REAL TmnLocal21 CCTK_ATTRIBUTE_UNUSED = eTxyL*ConfFactorA2;
    
    CCTK_REAL TmnLocal31 CCTK_ATTRIBUTE_UNUSED = eTxzL*ConfFactorA2;
    
    CCTK_REAL TmnLocal22 CCTK_ATTRIBUTE_UNUSED = eTyyL*ConfFactorA2;
    
    CCTK_REAL TmnLocal32 CCTK_ATTRIBUTE_UNUSED = eTyzL*ConfFactorA2;
    
    CCTK_REAL TmnLocal33 CCTK_ATTRIBUTE_UNUSED = eTzzL*ConfFactorA2;
    
    CCTK_REAL ttrace CCTK_ATTRIBUTE_UNUSED = -(pow(alpL,-2)*(TttLocal + 
      2*(betayL*betazL*TmnLocal32 + betaxL*(betayL*TmnLocal21 + 
      betazL*TmnLocal31 - TtxLocal)) - 2*(betayL*TtyLocal + betazL*TtzLocal) 
      - (gInv11*TmnLocal11 + gInv22*TmnLocal22 + 2*(gInv12*TmnLocal21 + 
      gInv13*TmnLocal31 + gInv23*TmnLocal32) + gInv33*TmnLocal33)*pow(alpL,2) 
      + TmnLocal11*pow(betaxL,2) + TmnLocal22*pow(betayL,2) + 
      TmnLocal33*pow(betazL,2)));
    
    CCTK_REAL STphirhsL CCTK_ATTRIBUTE_UNUSED = 
      betaxL*PDlopsided6th1(&STphi[index]) + 
      betayL*PDlopsided6th2(&STphi[index]) + 
      betazL*PDlopsided6th3(&STphi[index]) + alpL*STpiL;
    
    CCTK_REAL STpirhsL CCTK_ATTRIBUTE_UNUSED = 
      betaxL*PDlopsided6th1(&STpi[index]) + 
      betayL*PDlopsided6th2(&STpi[index]) + 
      betazL*PDlopsided6th3(&STpi[index]) + 
      PDstandard6th1alp*(gInv11*PDstandard6th1STphi + 
      gInv12*PDstandard6th2STphi + gInv13*PDstandard6th3STphi) + 
      PDstandard6th2alp*(gInv12*PDstandard6th1STphi + 
      gInv22*PDstandard6th2STphi + gInv23*PDstandard6th3STphi) + 
      PDstandard6th3alp*(gInv13*PDstandard6th1STphi + 
      gInv23*PDstandard6th2STphi + gInv33*PDstandard6th3STphi) + 
      alpL*(STpiL*ktrace + gInv11*PDstandard6th11STphi - (Gam3111*gInv11 + 
      2*Gam3121*gInv12 + 2*Gam3131*gInv13 + Gam3122*gInv22 + 2*Gam3132*gInv23 
      + Gam3133*gInv33)*PDstandard6th1STphi + gInv22*PDstandard6th22STphi + 
      2*(gInv12*PDstandard6th12STphi + gInv13*PDstandard6th13STphi + 
      gInv23*PDstandard6th23STphi) - (Gam3211*gInv11 + 2*Gam3221*gInv12 + 
      2*Gam3231*gInv13 + Gam3222*gInv22 + 2*Gam3232*gInv23 + 
      Gam3233*gInv33)*PDstandard6th2STphi + gInv33*PDstandard6th33STphi - 
      (Gam3311*gInv11 + 2*Gam3321*gInv12 + 2*Gam3331*gInv13 + Gam3322*gInv22 
      + 2*Gam3332*gInv23 + Gam3333*gInv33)*PDstandard6th3STphi + (STalpha0 + 
      1.*STphiL*STBeta)*ttrace);
    /* Copy local copies back to grid functions */
    STphirhs[index] = STphirhsL;
    STpirhs[index] = STpirhsL;
  }
  CCTK_ENDLOOP3(ScalarTensor_EvolveScalarField_6th);
}
extern "C" void ScalarTensor_EvolveScalarField_6th(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ScalarTensor_EvolveScalarField_6th_Body");
  }
  if (cctk_iteration % ScalarTensor_EvolveScalarField_6th_calc_every != ScalarTensor_EvolveScalarField_6th_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "admbase::curv",
    "admbase::lapse",
    "admbase::metric",
    "admbase::shift",
    "ScalarTensor::STphi_group",
    "ScalarTensor::STphi_grouprhs",
    "ScalarTensor::STpi_group",
    "ScalarTensor::STpi_grouprhs"};
  AssertGroupStorage(cctkGH, "ScalarTensor_EvolveScalarField_6th", 8, groups);
  
  EnsureStencilFits(cctkGH, "ScalarTensor_EvolveScalarField_6th", 4, 4, 4);
  
  LoopOverInterior(cctkGH, ScalarTensor_EvolveScalarField_6th_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ScalarTensor_EvolveScalarField_6th_Body");
  }
}

} // namespace ScalarTensor
