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
#include "loopcontrol.h"
#include "vectors.h"
#include "FLRWBackground.h"

#define KRANC_C

namespace FLRWBackground {

static void FLRWBackground_MetricBoundaryConformal_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /* Loop over the grid points */
  const int imin0=imin[0];
  const int imin1=imin[1];
  const int imin2=imin[2];
  const int imax0=imax[0];
  const int imax1=imax[1];
  const int imax2=imax[2];

  const CCTK_REAL t = cctk_time;
  double expansion = a0*t/t0;
  double Hubble = H0*pow(expansion , -2.0);

  double g0 = pow(expansion , 2);
  double k0 = -Hubble;

  #pragma omp parallel
  CCTK_LOOP3STR(FLRWBackground_MetricBoundary,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2],
    vecimin,vecimax, CCTK_REAL_VEC_SIZE)
  {
    CCTK_INT ind = CCTK_GFINDEX3D(cctkGH, i , j , k);

    gxx[ind] = g0;
    gxy[ind] = 0.;
    gxz[ind] = 0.;
    gyy[ind] = g0;
    gyz[ind] = 0.;
    gzz[ind] = g0;

    kxx[ind] = k0;
    kxy[ind] = 0.;
    kxz[ind] = 0.;
    kyy[ind] = k0;
    kyz[ind] = 0.;
    kzz[ind] = k0;

    alp[ind] = expansion;
    betax[ind] = 0.;
    betay[ind] = 0.;
    betaz[ind] = 0.;
  }
  CCTK_ENDLOOP3STR(FLRWBackground_MetricBoundary);
}



extern "C" void FLRWBackground_MetricBoundaryConformal(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering FLRWBackground_MetricBoundaryConformal_Body");
  }

  LoopOverBoundary(cctkGH, FLRWBackground_MetricBoundaryConformal_Body);

  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving FLRWBackground_MetricBoundaryConformal_Body");
  }
}

} // namespace FLRWBackground
