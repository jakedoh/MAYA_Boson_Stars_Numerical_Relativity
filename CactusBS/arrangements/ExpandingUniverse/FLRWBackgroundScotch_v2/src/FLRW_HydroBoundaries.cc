
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


extern "C" void FLRWBackground_SelectHydroBoundConds(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;

  /* BCs for HydroBase::rho */
  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                   "HydroBase::rho", "flat");
  if (ierr < 0)
    CCTK_WARN(0, "Failed to register Flat BC for HydroBase::rho!");

  /* BCs for HydroBase::press */
  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                   "HydroBase::press", "flat");
  if (ierr < 0)
    CCTK_WARN(0, "Failed to register Flat BC for HydroBase::press!");

  /* BCs for HydroBase::eps */

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                   "HydroBase::eps", "flat");
  if (ierr < 0)
    CCTK_WARN(0, "Failed to register Flat BC for HydroBase::eps!");


  /* BCs for HydroBase::vel */
  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                   "HydroBase::velx", "flat");
  ierr = ierr + Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                   "HydroBase::vely", "flat");
  ierr = ierr + Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                   "HydroBase::velz", "flat");
  if (ierr < 0)
    CCTK_WARN(0, "Failed to register Flat BC for HydroBase::vel!");

  return;
}



static void FLRWBackground_HydroBoundary_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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

  double expansion = a0 * pow(t/t0 , 1.0/2.0);
  double Hubble = H0*t0/t;

  double rhocrit = 3.0*H0*H0/(8.0*PI);

  double epst = getEps( rhocrit , K_intrinsic , expansion , tol);

  double rhot = rhocrit*pow(expansion , -4)/(1+epst);
  double presst = rhot*epst/3.0;
  double velt = 0.;


  #pragma omp parallel
  CCTK_LOOP3STR(FLRWBackground_HydroBoundary,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2],
    vecimin,vecimax, CCTK_REAL_VEC_SIZE)
  {
    CCTK_INT ind = CCTK_GFINDEX3D(cctkGH, i , j , k);

    CCTK_INT indx = CCTK_GFINDEX3D(cctkGH, i , j , k );
    CCTK_INT indy = CCTK_GFINDEX3D(cctkGH, i , j , k );
    CCTK_INT indz = CCTK_GFINDEX3D(cctkGH, i , j , k );

    /* Copy local copies back to grid functions */
    rho[ind] = rhot;
    press[ind] = presst;
    eps[ind] = epst;
    velx[indx] = velt;
    vely[indy] = velt;
    velz[indz] = velt;
    rho[ind] = 1e-10;

// I commented because it is defined in the Old HydroBase interface used by Scotch    w_lorentz[ind] = 1.0;

  }
  CCTK_ENDLOOP3STR(FLRWBackground_HydroBoundary);
}



extern "C" void FLRWBackground_HydroBoundary(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering FLRWBackground_HydroBoundary_Body");
  }

  const char* const groups[] = {
    "HydroBase::rho",
    "HydroBase::press",
    "HydroBase::eps",
    "HydroBase::vel"};

  AssertGroupStorage(cctkGH, "FLRWBackground_HydroBoundary", 4, groups);

  LoopOverBoundaryWithGhosts(cctkGH, FLRWBackground_HydroBoundary_Body);

  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving FLRWBackground_HydroBoundary_Body");
  }
}

} // namespace FLRWBackground
