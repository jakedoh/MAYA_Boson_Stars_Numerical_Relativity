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

namespace CT_Analytic {


static void CT_Analytic_Exact_Calc_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  CCTK_LOOP3(CT_Analytic_Exact_Calc,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL testc1L CCTK_ATTRIBUTE_UNUSED = testc1[index];
    CCTK_REAL testc2L CCTK_ATTRIBUTE_UNUSED = testc2[index];
    CCTK_REAL xL CCTK_ATTRIBUTE_UNUSED = x[index];
    CCTK_REAL yL CCTK_ATTRIBUTE_UNUSED = y[index];
    CCTK_REAL zL CCTK_ATTRIBUTE_UNUSED = z[index];
    
    /* Include user supplied include files */
    /* Precompute derivatives */
    /* Calculate temporaries and grid functions */
    CCTK_REAL d2soln CCTK_ATTRIBUTE_UNUSED = 
      ampG*((amp[0])*(-(exp(-0.5*pow((sigmax[0]),-2)*pow(-(x0[0]) + xL,2) - 
      0.5*pow((sigmay[0]),-2)*pow(-(y0[0]) + yL,2) - 
      0.5*pow((sigmaz[0]),-2)*pow(-(z0[0]) + zL,2))*pow((sigmax[0]),-2)) + 
      exp(-0.5*pow((sigmax[0]),-2)*pow(-(x0[0]) + xL,2) - 
      0.5*pow((sigmay[0]),-2)*pow(-(y0[0]) + yL,2) - 
      0.5*pow((sigmaz[0]),-2)*pow(-(z0[0]) + 
      zL,2))*pow((sigmax[0]),-4)*pow(-(x0[0]) + xL,2)) + 
      (amp[1])*(-(exp(-0.5*pow((sigmax[1]),-2)*pow(-(x0[1]) + xL,2) - 
      0.5*pow((sigmay[1]),-2)*pow(-(y0[1]) + yL,2) - 
      0.5*pow((sigmaz[1]),-2)*pow(-(z0[1]) + zL,2))*pow((sigmax[1]),-2)) + 
      exp(-0.5*pow((sigmax[1]),-2)*pow(-(x0[1]) + xL,2) - 
      0.5*pow((sigmay[1]),-2)*pow(-(y0[1]) + yL,2) - 
      0.5*pow((sigmaz[1]),-2)*pow(-(z0[1]) + 
      zL,2))*pow((sigmax[1]),-4)*pow(-(x0[1]) + xL,2)) + 
      (amp[2])*(-(exp(-0.5*pow((sigmax[2]),-2)*pow(-(x0[2]) + xL,2) - 
      0.5*pow((sigmay[2]),-2)*pow(-(y0[2]) + yL,2) - 
      0.5*pow((sigmaz[2]),-2)*pow(-(z0[2]) + zL,2))*pow((sigmax[2]),-2)) + 
      exp(-0.5*pow((sigmax[2]),-2)*pow(-(x0[2]) + xL,2) - 
      0.5*pow((sigmay[2]),-2)*pow(-(y0[2]) + yL,2) - 
      0.5*pow((sigmaz[2]),-2)*pow(-(z0[2]) + 
      zL,2))*pow((sigmax[2]),-4)*pow(-(x0[2]) + xL,2)) + 
      (amp[3])*(-(exp(-0.5*pow((sigmax[3]),-2)*pow(-(x0[3]) + xL,2) - 
      0.5*pow((sigmay[3]),-2)*pow(-(y0[3]) + yL,2) - 
      0.5*pow((sigmaz[3]),-2)*pow(-(z0[3]) + zL,2))*pow((sigmax[3]),-2)) + 
      exp(-0.5*pow((sigmax[3]),-2)*pow(-(x0[3]) + xL,2) - 
      0.5*pow((sigmay[3]),-2)*pow(-(y0[3]) + yL,2) - 
      0.5*pow((sigmaz[3]),-2)*pow(-(z0[3]) + 
      zL,2))*pow((sigmax[3]),-4)*pow(-(x0[3]) + xL,2)) + 
      (amp[4])*(-(exp(-0.5*pow((sigmax[4]),-2)*pow(-(x0[4]) + xL,2) - 
      0.5*pow((sigmay[4]),-2)*pow(-(y0[4]) + yL,2) - 
      0.5*pow((sigmaz[4]),-2)*pow(-(z0[4]) + zL,2))*pow((sigmax[4]),-2)) + 
      exp(-0.5*pow((sigmax[4]),-2)*pow(-(x0[4]) + xL,2) - 
      0.5*pow((sigmay[4]),-2)*pow(-(y0[4]) + yL,2) - 
      0.5*pow((sigmaz[4]),-2)*pow(-(z0[4]) + 
      zL,2))*pow((sigmax[4]),-4)*pow(-(x0[4]) + xL,2))) + 
      ampG*((amp[0])*(-(exp(-0.5*pow((sigmax[0]),-2)*pow(-(x0[0]) + xL,2) - 
      0.5*pow((sigmay[0]),-2)*pow(-(y0[0]) + yL,2) - 
      0.5*pow((sigmaz[0]),-2)*pow(-(z0[0]) + zL,2))*pow((sigmay[0]),-2)) + 
      exp(-0.5*pow((sigmax[0]),-2)*pow(-(x0[0]) + xL,2) - 
      0.5*pow((sigmay[0]),-2)*pow(-(y0[0]) + yL,2) - 
      0.5*pow((sigmaz[0]),-2)*pow(-(z0[0]) + 
      zL,2))*pow((sigmay[0]),-4)*pow(-(y0[0]) + yL,2)) + 
      (amp[1])*(-(exp(-0.5*pow((sigmax[1]),-2)*pow(-(x0[1]) + xL,2) - 
      0.5*pow((sigmay[1]),-2)*pow(-(y0[1]) + yL,2) - 
      0.5*pow((sigmaz[1]),-2)*pow(-(z0[1]) + zL,2))*pow((sigmay[1]),-2)) + 
      exp(-0.5*pow((sigmax[1]),-2)*pow(-(x0[1]) + xL,2) - 
      0.5*pow((sigmay[1]),-2)*pow(-(y0[1]) + yL,2) - 
      0.5*pow((sigmaz[1]),-2)*pow(-(z0[1]) + 
      zL,2))*pow((sigmay[1]),-4)*pow(-(y0[1]) + yL,2)) + 
      (amp[2])*(-(exp(-0.5*pow((sigmax[2]),-2)*pow(-(x0[2]) + xL,2) - 
      0.5*pow((sigmay[2]),-2)*pow(-(y0[2]) + yL,2) - 
      0.5*pow((sigmaz[2]),-2)*pow(-(z0[2]) + zL,2))*pow((sigmay[2]),-2)) + 
      exp(-0.5*pow((sigmax[2]),-2)*pow(-(x0[2]) + xL,2) - 
      0.5*pow((sigmay[2]),-2)*pow(-(y0[2]) + yL,2) - 
      0.5*pow((sigmaz[2]),-2)*pow(-(z0[2]) + 
      zL,2))*pow((sigmay[2]),-4)*pow(-(y0[2]) + yL,2)) + 
      (amp[3])*(-(exp(-0.5*pow((sigmax[3]),-2)*pow(-(x0[3]) + xL,2) - 
      0.5*pow((sigmay[3]),-2)*pow(-(y0[3]) + yL,2) - 
      0.5*pow((sigmaz[3]),-2)*pow(-(z0[3]) + zL,2))*pow((sigmay[3]),-2)) + 
      exp(-0.5*pow((sigmax[3]),-2)*pow(-(x0[3]) + xL,2) - 
      0.5*pow((sigmay[3]),-2)*pow(-(y0[3]) + yL,2) - 
      0.5*pow((sigmaz[3]),-2)*pow(-(z0[3]) + 
      zL,2))*pow((sigmay[3]),-4)*pow(-(y0[3]) + yL,2)) + 
      (amp[4])*(-(exp(-0.5*pow((sigmax[4]),-2)*pow(-(x0[4]) + xL,2) - 
      0.5*pow((sigmay[4]),-2)*pow(-(y0[4]) + yL,2) - 
      0.5*pow((sigmaz[4]),-2)*pow(-(z0[4]) + zL,2))*pow((sigmay[4]),-2)) + 
      exp(-0.5*pow((sigmax[4]),-2)*pow(-(x0[4]) + xL,2) - 
      0.5*pow((sigmay[4]),-2)*pow(-(y0[4]) + yL,2) - 
      0.5*pow((sigmaz[4]),-2)*pow(-(z0[4]) + 
      zL,2))*pow((sigmay[4]),-4)*pow(-(y0[4]) + yL,2))) + 
      ampG*((amp[0])*(-(exp(-0.5*pow((sigmax[0]),-2)*pow(-(x0[0]) + xL,2) - 
      0.5*pow((sigmay[0]),-2)*pow(-(y0[0]) + yL,2) - 
      0.5*pow((sigmaz[0]),-2)*pow(-(z0[0]) + zL,2))*pow((sigmaz[0]),-2)) + 
      exp(-0.5*pow((sigmax[0]),-2)*pow(-(x0[0]) + xL,2) - 
      0.5*pow((sigmay[0]),-2)*pow(-(y0[0]) + yL,2) - 
      0.5*pow((sigmaz[0]),-2)*pow(-(z0[0]) + 
      zL,2))*pow((sigmaz[0]),-4)*pow(-(z0[0]) + zL,2)) + 
      (amp[1])*(-(exp(-0.5*pow((sigmax[1]),-2)*pow(-(x0[1]) + xL,2) - 
      0.5*pow((sigmay[1]),-2)*pow(-(y0[1]) + yL,2) - 
      0.5*pow((sigmaz[1]),-2)*pow(-(z0[1]) + zL,2))*pow((sigmaz[1]),-2)) + 
      exp(-0.5*pow((sigmax[1]),-2)*pow(-(x0[1]) + xL,2) - 
      0.5*pow((sigmay[1]),-2)*pow(-(y0[1]) + yL,2) - 
      0.5*pow((sigmaz[1]),-2)*pow(-(z0[1]) + 
      zL,2))*pow((sigmaz[1]),-4)*pow(-(z0[1]) + zL,2)) + 
      (amp[2])*(-(exp(-0.5*pow((sigmax[2]),-2)*pow(-(x0[2]) + xL,2) - 
      0.5*pow((sigmay[2]),-2)*pow(-(y0[2]) + yL,2) - 
      0.5*pow((sigmaz[2]),-2)*pow(-(z0[2]) + zL,2))*pow((sigmaz[2]),-2)) + 
      exp(-0.5*pow((sigmax[2]),-2)*pow(-(x0[2]) + xL,2) - 
      0.5*pow((sigmay[2]),-2)*pow(-(y0[2]) + yL,2) - 
      0.5*pow((sigmaz[2]),-2)*pow(-(z0[2]) + 
      zL,2))*pow((sigmaz[2]),-4)*pow(-(z0[2]) + zL,2)) + 
      (amp[3])*(-(exp(-0.5*pow((sigmax[3]),-2)*pow(-(x0[3]) + xL,2) - 
      0.5*pow((sigmay[3]),-2)*pow(-(y0[3]) + yL,2) - 
      0.5*pow((sigmaz[3]),-2)*pow(-(z0[3]) + zL,2))*pow((sigmaz[3]),-2)) + 
      exp(-0.5*pow((sigmax[3]),-2)*pow(-(x0[3]) + xL,2) - 
      0.5*pow((sigmay[3]),-2)*pow(-(y0[3]) + yL,2) - 
      0.5*pow((sigmaz[3]),-2)*pow(-(z0[3]) + 
      zL,2))*pow((sigmaz[3]),-4)*pow(-(z0[3]) + zL,2)) + 
      (amp[4])*(-(exp(-0.5*pow((sigmax[4]),-2)*pow(-(x0[4]) + xL,2) - 
      0.5*pow((sigmay[4]),-2)*pow(-(y0[4]) + yL,2) - 
      0.5*pow((sigmaz[4]),-2)*pow(-(z0[4]) + zL,2))*pow((sigmaz[4]),-2)) + 
      exp(-0.5*pow((sigmax[4]),-2)*pow(-(x0[4]) + xL,2) - 
      0.5*pow((sigmay[4]),-2)*pow(-(y0[4]) + yL,2) - 
      0.5*pow((sigmaz[4]),-2)*pow(-(z0[4]) + 
      zL,2))*pow((sigmaz[4]),-4)*pow(-(z0[4]) + zL,2))) + 
      0.5*ampSg*(massa*(3*pow(xL - xa,2)*pow(pow(eps,2) + pow(xL - xa,2) + 
      pow(yL - ya,2) + pow(zL - za,2),-2.5) - pow(pow(eps,2) + pow(xL - xa,2) 
      + pow(yL - ya,2) + pow(zL - za,2),-1.5)) + massb*(3*pow(xL - 
      xb,2)*pow(pow(eps,2) + pow(xL - xb,2) + pow(yL - yb,2) + pow(zL - 
      zb,2),-2.5) - pow(pow(eps,2) + pow(xL - xb,2) + pow(yL - yb,2) + pow(zL 
      - zb,2),-1.5))) + 0.5*ampSg*(massa*(3*pow(yL - ya,2)*pow(pow(eps,2) + 
      pow(xL - xa,2) + pow(yL - ya,2) + pow(zL - za,2),-2.5) - pow(pow(eps,2) 
      + pow(xL - xa,2) + pow(yL - ya,2) + pow(zL - za,2),-1.5)) + 
      massb*(3*pow(yL - yb,2)*pow(pow(eps,2) + pow(xL - xb,2) + pow(yL - 
      yb,2) + pow(zL - zb,2),-2.5) - pow(pow(eps,2) + pow(xL - xb,2) + pow(yL 
      - yb,2) + pow(zL - zb,2),-1.5))) + 0.5*ampSg*(massa*(3*pow(zL - 
      za,2)*pow(pow(eps,2) + pow(xL - xa,2) + pow(yL - ya,2) + pow(zL - 
      za,2),-2.5) - pow(pow(eps,2) + pow(xL - xa,2) + pow(yL - ya,2) + pow(zL 
      - za,2),-1.5)) + massb*(3*pow(zL - zb,2)*pow(pow(eps,2) + pow(xL - 
      xb,2) + pow(yL - yb,2) + pow(zL - zb,2),-2.5) - pow(pow(eps,2) + pow(xL 
      - xb,2) + pow(yL - yb,2) + pow(zL - zb,2),-1.5))) - 
      39.47841760435743447533796399950460454125*ampS*pow(kx,2)*sin(6.283185307179586476925286766559005768394*(xL*kx 
      + phasex))*sin(6.283185307179586476925286766559005768394*(yL*ky + 
      phasey))*sin(6.283185307179586476925286766559005768394*(zL*kz + 
      phasez)) - 
      39.47841760435743447533796399950460454125*ampS*pow(ky,2)*sin(6.283185307179586476925286766559005768394*(xL*kx 
      + phasex))*sin(6.283185307179586476925286766559005768394*(yL*ky + 
      phasey))*sin(6.283185307179586476925286766559005768394*(zL*kz + 
      phasez)) - 
      39.47841760435743447533796399950460454125*ampS*pow(kz,2)*sin(6.283185307179586476925286766559005768394*(xL*kx 
      + phasex))*sin(6.283185307179586476925286766559005768394*(yL*ky + 
      phasey))*sin(6.283185307179586476925286766559005768394*(zL*kz + 
      phasez));
    
    CCTK_REAL epsiL CCTK_ATTRIBUTE_UNUSED = ampC + 
      ampG*((amp[0])*exp(-0.5*pow((sigmax[0]),-2)*pow(-(x0[0]) + xL,2) - 
      0.5*pow((sigmay[0]),-2)*pow(-(y0[0]) + yL,2) - 
      0.5*pow((sigmaz[0]),-2)*pow(-(z0[0]) + zL,2)) + 
      (amp[1])*exp(-0.5*pow((sigmax[1]),-2)*pow(-(x0[1]) + xL,2) - 
      0.5*pow((sigmay[1]),-2)*pow(-(y0[1]) + yL,2) - 
      0.5*pow((sigmaz[1]),-2)*pow(-(z0[1]) + zL,2)) + 
      (amp[2])*exp(-0.5*pow((sigmax[2]),-2)*pow(-(x0[2]) + xL,2) - 
      0.5*pow((sigmay[2]),-2)*pow(-(y0[2]) + yL,2) - 
      0.5*pow((sigmaz[2]),-2)*pow(-(z0[2]) + zL,2)) + 
      (amp[3])*exp(-0.5*pow((sigmax[3]),-2)*pow(-(x0[3]) + xL,2) - 
      0.5*pow((sigmay[3]),-2)*pow(-(y0[3]) + yL,2) - 
      0.5*pow((sigmaz[3]),-2)*pow(-(z0[3]) + zL,2)) + 
      (amp[4])*exp(-0.5*pow((sigmax[4]),-2)*pow(-(x0[4]) + xL,2) - 
      0.5*pow((sigmay[4]),-2)*pow(-(y0[4]) + yL,2) - 
      0.5*pow((sigmaz[4]),-2)*pow(-(z0[4]) + zL,2))) + 
      0.5*ampSg*(massa*pow(pow(eps,2) + pow(xL - xa,2) + pow(yL - ya,2) + 
      pow(zL - za,2),-0.5) + massb*pow(pow(eps,2) + pow(xL - xb,2) + pow(yL - 
      yb,2) + pow(zL - zb,2),-0.5)) + 
      ampS*sin(6.283185307179586476925286766559005768394*(xL*kx + 
      phasex))*sin(6.283185307179586476925286766559005768394*(yL*ky + 
      phasey))*sin(6.283185307179586476925286766559005768394*(zL*kz + 
      phasez));
    
    CCTK_REAL elaplacianL CCTK_ATTRIBUTE_UNUSED = d2soln;
    
    CCTK_REAL testinipsiL CCTK_ATTRIBUTE_UNUSED = ampI*(ampC + 
      ampG*((amp[0])*exp(-0.5*pow((sigmax[0]),-2)*pow(-(x0[0]) + xL,2) - 
      0.5*pow((sigmay[0]),-2)*pow(-(y0[0]) + yL,2) - 
      0.5*pow((sigmaz[0]),-2)*pow(-(z0[0]) + zL,2)) + 
      (amp[1])*exp(-0.5*pow((sigmax[1]),-2)*pow(-(x0[1]) + xL,2) - 
      0.5*pow((sigmay[1]),-2)*pow(-(y0[1]) + yL,2) - 
      0.5*pow((sigmaz[1]),-2)*pow(-(z0[1]) + zL,2)) + 
      (amp[2])*exp(-0.5*pow((sigmax[2]),-2)*pow(-(x0[2]) + xL,2) - 
      0.5*pow((sigmay[2]),-2)*pow(-(y0[2]) + yL,2) - 
      0.5*pow((sigmaz[2]),-2)*pow(-(z0[2]) + zL,2)) + 
      (amp[3])*exp(-0.5*pow((sigmax[3]),-2)*pow(-(x0[3]) + xL,2) - 
      0.5*pow((sigmay[3]),-2)*pow(-(y0[3]) + yL,2) - 
      0.5*pow((sigmaz[3]),-2)*pow(-(z0[3]) + zL,2)) + 
      (amp[4])*exp(-0.5*pow((sigmax[4]),-2)*pow(-(x0[4]) + xL,2) - 
      0.5*pow((sigmay[4]),-2)*pow(-(y0[4]) + yL,2) - 
      0.5*pow((sigmaz[4]),-2)*pow(-(z0[4]) + zL,2))) + 
      0.5*ampSg*(massa*pow(pow(eps,2) + pow(xL - xa,2) + pow(yL - ya,2) + 
      pow(zL - za,2),-0.5) + massb*pow(pow(eps,2) + pow(xL - xb,2) + pow(yL - 
      yb,2) + pow(zL - zb,2),-0.5)) + 
      ampS*sin(6.283185307179586476925286766559005768394*(xL*kx + 
      phasex))*sin(6.283185307179586476925286766559005768394*(yL*ky + 
      phasey))*sin(6.283185307179586476925286766559005768394*(zL*kz + 
      phasez)));
    
    CCTK_REAL testcxxL CCTK_ATTRIBUTE_UNUSED = 1;
    
    CCTK_REAL testcxyL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL testcxzL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL testcyyL CCTK_ATTRIBUTE_UNUSED = 1;
    
    CCTK_REAL testcyzL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL testczzL CCTK_ATTRIBUTE_UNUSED = 1;
    
    CCTK_REAL testcxL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL testcyL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL testczL CCTK_ATTRIBUTE_UNUSED = 0;
    
    testc1L = ampC1;
    
    testc2L = 0;
    
    CCTK_REAL testc3L CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL testc4L CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL testc0L CCTK_ATTRIBUTE_UNUSED = -d2soln - testc2L*pow(ampC + 
      ampG*((amp[0])*exp(-0.5*pow((sigmax[0]),-2)*pow(-(x0[0]) + xL,2) - 
      0.5*pow((sigmay[0]),-2)*pow(-(y0[0]) + yL,2) - 
      0.5*pow((sigmaz[0]),-2)*pow(-(z0[0]) + zL,2)) + 
      (amp[1])*exp(-0.5*pow((sigmax[1]),-2)*pow(-(x0[1]) + xL,2) - 
      0.5*pow((sigmay[1]),-2)*pow(-(y0[1]) + yL,2) - 
      0.5*pow((sigmaz[1]),-2)*pow(-(z0[1]) + zL,2)) + 
      (amp[2])*exp(-0.5*pow((sigmax[2]),-2)*pow(-(x0[2]) + xL,2) - 
      0.5*pow((sigmay[2]),-2)*pow(-(y0[2]) + yL,2) - 
      0.5*pow((sigmaz[2]),-2)*pow(-(z0[2]) + zL,2)) + 
      (amp[3])*exp(-0.5*pow((sigmax[3]),-2)*pow(-(x0[3]) + xL,2) - 
      0.5*pow((sigmay[3]),-2)*pow(-(y0[3]) + yL,2) - 
      0.5*pow((sigmaz[3]),-2)*pow(-(z0[3]) + zL,2)) + 
      (amp[4])*exp(-0.5*pow((sigmax[4]),-2)*pow(-(x0[4]) + xL,2) - 
      0.5*pow((sigmay[4]),-2)*pow(-(y0[4]) + yL,2) - 
      0.5*pow((sigmaz[4]),-2)*pow(-(z0[4]) + zL,2))) + 
      0.5*ampSg*(massa*pow(pow(eps,2) + pow(xL - xa,2) + pow(yL - ya,2) + 
      pow(zL - za,2),-0.5) + massb*pow(pow(eps,2) + pow(xL - xb,2) + pow(yL - 
      yb,2) + pow(zL - zb,2),-0.5)) + 
      ampS*sin(6.283185307179586476925286766559005768394*(xL*kx + 
      phasex))*sin(6.283185307179586476925286766559005768394*(yL*ky + 
      phasey))*sin(6.283185307179586476925286766559005768394*(zL*kz + 
      phasez)),2) - testc1L*(ampC + 
      ampG*((amp[0])*exp(-0.5*pow((sigmax[0]),-2)*pow(-(x0[0]) + xL,2) - 
      0.5*pow((sigmay[0]),-2)*pow(-(y0[0]) + yL,2) - 
      0.5*pow((sigmaz[0]),-2)*pow(-(z0[0]) + zL,2)) + 
      (amp[1])*exp(-0.5*pow((sigmax[1]),-2)*pow(-(x0[1]) + xL,2) - 
      0.5*pow((sigmay[1]),-2)*pow(-(y0[1]) + yL,2) - 
      0.5*pow((sigmaz[1]),-2)*pow(-(z0[1]) + zL,2)) + 
      (amp[2])*exp(-0.5*pow((sigmax[2]),-2)*pow(-(x0[2]) + xL,2) - 
      0.5*pow((sigmay[2]),-2)*pow(-(y0[2]) + yL,2) - 
      0.5*pow((sigmaz[2]),-2)*pow(-(z0[2]) + zL,2)) + 
      (amp[3])*exp(-0.5*pow((sigmax[3]),-2)*pow(-(x0[3]) + xL,2) - 
      0.5*pow((sigmay[3]),-2)*pow(-(y0[3]) + yL,2) - 
      0.5*pow((sigmaz[3]),-2)*pow(-(z0[3]) + zL,2)) + 
      (amp[4])*exp(-0.5*pow((sigmax[4]),-2)*pow(-(x0[4]) + xL,2) - 
      0.5*pow((sigmay[4]),-2)*pow(-(y0[4]) + yL,2) - 
      0.5*pow((sigmaz[4]),-2)*pow(-(z0[4]) + zL,2))) + 
      0.5*ampSg*(massa*pow(pow(eps,2) + pow(xL - xa,2) + pow(yL - ya,2) + 
      pow(zL - za,2),-0.5) + massb*pow(pow(eps,2) + pow(xL - xb,2) + pow(yL - 
      yb,2) + pow(zL - zb,2),-0.5)) + 
      ampS*sin(6.283185307179586476925286766559005768394*(xL*kx + 
      phasex))*sin(6.283185307179586476925286766559005768394*(yL*ky + 
      phasey))*sin(6.283185307179586476925286766559005768394*(zL*kz + 
      phasez)));
    
    testc1L = ampC1;
    /* Copy local copies back to grid functions */
    elaplacian[index] = elaplacianL;
    epsi[index] = epsiL;
    testc0[index] = testc0L;
    testc1[index] = testc1L;
    testc2[index] = testc2L;
    testc3[index] = testc3L;
    testc4[index] = testc4L;
    testcx[index] = testcxL;
    testcxx[index] = testcxxL;
    testcxy[index] = testcxyL;
    testcxz[index] = testcxzL;
    testcy[index] = testcyL;
    testcyy[index] = testcyyL;
    testcyz[index] = testcyzL;
    testcz[index] = testczL;
    testczz[index] = testczzL;
    testinipsi[index] = testinipsiL;
  }
  CCTK_ENDLOOP3(CT_Analytic_Exact_Calc);
}
extern "C" void CT_Analytic_Exact_Calc(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_CT_Analytic_Exact_Calc
  DECLARE_CCTK_ARGUMENTS_CHECKED(CT_Analytic_Exact_Calc);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering CT_Analytic_Exact_Calc_Body");
  }
  if (cctk_iteration % CT_Analytic_Exact_Calc_calc_every != CT_Analytic_Exact_Calc_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "CT_Analytic::CT_elaplacian",
    "CT_Analytic::CT_epsi",
    "CT_Analytic::CT_testc0",
    "CT_Analytic::CT_testc1",
    "CT_Analytic::CT_testc2",
    "CT_Analytic::CT_testc3",
    "CT_Analytic::CT_testc4",
    "CT_Analytic::CT_testcx",
    "CT_Analytic::CT_testcxx",
    "CT_Analytic::CT_testcxy",
    "CT_Analytic::CT_testcxz",
    "CT_Analytic::CT_testcy",
    "CT_Analytic::CT_testcyy",
    "CT_Analytic::CT_testcyz",
    "CT_Analytic::CT_testcz",
    "CT_Analytic::CT_testczz",
    "CT_Analytic::CT_testinipsi",
    "grid::coordinates"};
  AssertGroupStorage(cctkGH, "CT_Analytic_Exact_Calc", 18, groups);
  
  
  LoopOverEverything(cctkGH, CT_Analytic_Exact_Calc_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving CT_Analytic_Exact_Calc_Body");
  }
}

} // namespace CT_Analytic
