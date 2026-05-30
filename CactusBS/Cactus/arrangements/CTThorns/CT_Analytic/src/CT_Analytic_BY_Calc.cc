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


static void CT_Analytic_BY_Calc_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  CCTK_LOOP3(CT_Analytic_BY_Calc,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL xL CCTK_ATTRIBUTE_UNUSED = x[index];
    CCTK_REAL yL CCTK_ATTRIBUTE_UNUSED = y[index];
    CCTK_REAL zL CCTK_ATTRIBUTE_UNUSED = z[index];
    
    /* Include user supplied include files */
    /* Precompute derivatives */
    /* Calculate temporaries and grid functions */
    CCTK_REAL ra CCTK_ATTRIBUTE_UNUSED = pow(pow(eps,2) + pow(xL - xa,2) + 
      pow(yL - ya,2) + pow(zL - za,2),0.5);
    
    CCTK_REAL nax CCTK_ATTRIBUTE_UNUSED = (xL - xa)*pow(ra,-1);
    
    CCTK_REAL nay CCTK_ATTRIBUTE_UNUSED = (yL - ya)*pow(ra,-1);
    
    CCTK_REAL naz CCTK_ATTRIBUTE_UNUSED = (zL - za)*pow(ra,-1);
    
    CCTK_REAL rb CCTK_ATTRIBUTE_UNUSED = pow(pow(eps,2) + pow(xL - xb,2) + 
      pow(yL - yb,2) + pow(zL - zb,2),0.5);
    
    CCTK_REAL nbx CCTK_ATTRIBUTE_UNUSED = (xL - xb)*pow(rb,-1);
    
    CCTK_REAL nby CCTK_ATTRIBUTE_UNUSED = (yL - yb)*pow(rb,-1);
    
    CCTK_REAL nbz CCTK_ATTRIBUTE_UNUSED = (zL - zb)*pow(rb,-1);
    
    CCTK_REAL nPa CCTK_ATTRIBUTE_UNUSED = nax*Pax + nay*Pay + naz*Paz;
    
    CCTK_REAL nPb CCTK_ATTRIBUTE_UNUSED = nbx*Pbx + nby*Pby + nbz*Pbz;
    
    CCTK_REAL epsiL CCTK_ATTRIBUTE_UNUSED = 1;
    
    CCTK_REAL testinipsiL CCTK_ATTRIBUTE_UNUSED = ampI;
    
    CCTK_REAL testcxxL CCTK_ATTRIBUTE_UNUSED = 1;
    
    CCTK_REAL testcyyL CCTK_ATTRIBUTE_UNUSED = 1;
    
    CCTK_REAL testczzL CCTK_ATTRIBUTE_UNUSED = 1;
    
    CCTK_REAL testa0L CCTK_ATTRIBUTE_UNUSED = 0.5*(massa*pow(pow(eps,2) + 
      pow(xL - xa,2) + pow(yL - ya,2) + pow(zL - za,2),-0.5) + 
      massb*pow(pow(eps,2) + pow(xL - xb,2) + pow(yL - yb,2) + pow(zL - 
      zb,2),-0.5));
    
    CCTK_REAL testa1L CCTK_ATTRIBUTE_UNUSED = 
      (amp[0])*exp(0.5*(-(pow((sigmax[0]),-2)*pow((x0[0]) - xL,2)) - 
      pow((sigmay[0]),-2)*pow((y0[0]) - yL,2) - 
      pow((sigmaz[0]),-2)*pow((z0[0]) - zL,2))) + 
      (amp[1])*exp(0.5*(-(pow((sigmax[1]),-2)*pow((x0[1]) - xL,2)) - 
      pow((sigmay[1]),-2)*pow((y0[1]) - yL,2) - 
      pow((sigmaz[1]),-2)*pow((z0[1]) - zL,2))) + 
      (amp[2])*exp(0.5*(-(pow((sigmax[2]),-2)*pow((x0[2]) - xL,2)) - 
      pow((sigmay[2]),-2)*pow((y0[2]) - yL,2) - 
      pow((sigmaz[2]),-2)*pow((z0[2]) - zL,2))) + 
      (amp[3])*exp(0.5*(-(pow((sigmax[3]),-2)*pow((x0[3]) - xL,2)) - 
      pow((sigmay[3]),-2)*pow((y0[3]) - yL,2) - 
      pow((sigmaz[3]),-2)*pow((z0[3]) - zL,2))) + 
      (amp[4])*exp(0.5*(-(pow((sigmax[4]),-2)*pow((x0[4]) - xL,2)) - 
      pow((sigmay[4]),-2)*pow((y0[4]) - yL,2) - 
      pow((sigmaz[4]),-2)*pow((z0[4]) - zL,2)));
    
    CCTK_REAL testXxL CCTK_ATTRIBUTE_UNUSED = -0.25*((nax*nPa + 
      7*Pax)*pow(ra,-1) + (nbx*nPb + 7*Pbx)*pow(rb,-1));
    
    CCTK_REAL testXyL CCTK_ATTRIBUTE_UNUSED = -0.25*((nay*nPa + 
      7*Pay)*pow(ra,-1) + (nby*nPb + 7*Pby)*pow(rb,-1));
    
    CCTK_REAL testXzL CCTK_ATTRIBUTE_UNUSED = -0.25*((naz*nPa + 
      7*Paz)*pow(ra,-1) + (nbz*nPb + 7*Pbz)*pow(rb,-1));
    
    CCTK_REAL testAxxL CCTK_ATTRIBUTE_UNUSED = (3.*nax*Pax + nPa*(-1.5 + 
      1.5*pow(nax,2)))*pow(ra,-2) + (3.*nbx*Pbx + nPb*(-1.5 + 
      1.5*pow(nbx,2)))*pow(rb,-2);
    
    CCTK_REAL testAxyL CCTK_ATTRIBUTE_UNUSED = 1.5*((nay*Pax + 
      nax*(nay*nPa + Pay))*pow(ra,-2) + (nby*Pbx + nbx*(nby*nPb + 
      Pby))*pow(rb,-2));
    
    CCTK_REAL testAxzL CCTK_ATTRIBUTE_UNUSED = 1.5*((naz*Pax + 
      nax*(naz*nPa + Paz))*pow(ra,-2) + (nbz*Pbx + nbx*(nbz*nPb + 
      Pbz))*pow(rb,-2));
    
    CCTK_REAL testAyyL CCTK_ATTRIBUTE_UNUSED = (3.*nay*Pay + nPa*(-1.5 + 
      1.5*pow(nay,2)))*pow(ra,-2) + (3.*nby*Pby + nPb*(-1.5 + 
      1.5*pow(nby,2)))*pow(rb,-2);
    
    CCTK_REAL testAyzL CCTK_ATTRIBUTE_UNUSED = 1.5*((naz*Pay + 
      nay*(naz*nPa + Paz))*pow(ra,-2) + (nbz*Pby + nby*(nbz*nPb + 
      Pbz))*pow(rb,-2));
    
    CCTK_REAL testAzzL CCTK_ATTRIBUTE_UNUSED = (3.*naz*Paz + nPa*(-1.5 + 
      1.5*pow(naz,2)))*pow(ra,-2) + (3.*nbz*Pbz + nPb*(-1.5 + 
      1.5*pow(nbz,2)))*pow(rb,-2);
    
    CCTK_REAL testZL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL testKL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL testdxKL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL testdyKL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL testdzKL CCTK_ATTRIBUTE_UNUSED = 0;
    /* Copy local copies back to grid functions */
    epsi[index] = epsiL;
    testa0[index] = testa0L;
    testa1[index] = testa1L;
    testAxx[index] = testAxxL;
    testAxy[index] = testAxyL;
    testAxz[index] = testAxzL;
    testAyy[index] = testAyyL;
    testAyz[index] = testAyzL;
    testAzz[index] = testAzzL;
    testcxx[index] = testcxxL;
    testcyy[index] = testcyyL;
    testczz[index] = testczzL;
    testdxK[index] = testdxKL;
    testdyK[index] = testdyKL;
    testdzK[index] = testdzKL;
    testinipsi[index] = testinipsiL;
    testK[index] = testKL;
    testXx[index] = testXxL;
    testXy[index] = testXyL;
    testXz[index] = testXzL;
    testZ[index] = testZL;
  }
  CCTK_ENDLOOP3(CT_Analytic_BY_Calc);
}
extern "C" void CT_Analytic_BY_Calc(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_CT_Analytic_BY_Calc
  DECLARE_CCTK_ARGUMENTS_CHECKED(CT_Analytic_BY_Calc);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering CT_Analytic_BY_Calc_Body");
  }
  if (cctk_iteration % CT_Analytic_BY_Calc_calc_every != CT_Analytic_BY_Calc_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "CT_Analytic::CT_epsi",
    "CT_Analytic::CT_testa0",
    "CT_Analytic::CT_testa1",
    "CT_Analytic::CT_testAxx",
    "CT_Analytic::CT_testAxy",
    "CT_Analytic::CT_testAxz",
    "CT_Analytic::CT_testAyy",
    "CT_Analytic::CT_testAyz",
    "CT_Analytic::CT_testAzz",
    "CT_Analytic::CT_testcxx",
    "CT_Analytic::CT_testcyy",
    "CT_Analytic::CT_testczz",
    "CT_Analytic::CT_testdxK",
    "CT_Analytic::CT_testdyK",
    "CT_Analytic::CT_testdzK",
    "CT_Analytic::CT_testinipsi",
    "CT_Analytic::CT_testK",
    "CT_Analytic::CT_testXx",
    "CT_Analytic::CT_testXy",
    "CT_Analytic::CT_testXz",
    "CT_Analytic::CT_testZ",
    "grid::coordinates"};
  AssertGroupStorage(cctkGH, "CT_Analytic_BY_Calc", 22, groups);
  
  
  LoopOverEverything(cctkGH, CT_Analytic_BY_Calc_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving CT_Analytic_BY_Calc_Body");
  }
}

} // namespace CT_Analytic
