// Please consult GiRaFFEfood/README for the list of authors and current maintainers //
/*
 * This program implements the initial configuration called FFE Breakdown, a 1D solution
 * that initially obeys the FFE conditions, but evolves to a state where it does not.
 * The vector potential is calculated with formula (121), arXiv:1310.3274v2.
 * We calculate the magnetic and electric fields with eqs. (120), arXiv:1310.3274v2.
 * Then we use eq. (99) from the same paper to obtain the velocity.
 * */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstdlib>

#ifndef M_PI
#define M_PI 3.141592653589793238463
#endif

using namespace std;

extern "C" void GiRaFFEfood_FFEBreakdown(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT imin = 0, jmin = 0, kmin = 0;
  CCTK_INT imax = cctk_lsh[0], jmax = cctk_lsh[1], kmax = cctk_lsh[2];

  CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  CCTK_REAL dy = CCTK_DELTA_SPACE(1);

#pragma omp parallel for
  for(int k=kmin;k<kmax;k++)
    for(int j=jmin;j<jmax;j++)
      for(int i=imin;i<imax;i++) {
        int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        //The vector potential is staggered. Getting the staggered coordinates.
        CCTK_REAL xL_plus_half = x[index] + 0.5*dx;
        CCTK_REAL yL_plus_half = y[index] + 0.5*dy;

        //Equation (121), arXiv:1310.3274v2. Calculate the vector potential.
        //Calculate Ay(i+0.5,j,k+0.5).
        CCTK_REAL AyL;
        //$A_y=x-0.2$ if $x<0$
        if (xL_plus_half <= 0.0) {
          AyL = xL_plus_half - 0.2;
        //$A_y=-5.0x^2+x-0.2$ if $0<x<0.2$
        } else if ((xL_plus_half > 0.0) and (xL_plus_half < 0.2)) {
          AyL = -5.0*xL_plus_half*xL_plus_half + xL_plus_half - 0.2;
        //$A_y=-x$ if $x>0.2$
        } else { // if (xL_plus_half >= 0.2) {
          AyL = -xL_plus_half;
        }

        //Calculate Az(i+0.5,j+0.5,k)
        //$A_z=y-A_y$
        CCTK_REAL AzL = yL_plus_half - AyL;

        //Fill in the interface variables:
        Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = 0.0;
        Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = AyL;
        Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] = AzL;

        // Initialize the velocity:
        // We want unstaggered v, so we want unstaggered B and E.
        CCTK_REAL xL = x[index];
        CCTK_REAL zx = -10.0*xL + 1.0;
        CCTK_REAL BxL, ByL, BzL;
        //Equation (120), arXiv:1310.3274v2. Calculate E and B.
        if (xL <= 0.0) {
          BxL = 1.0, ByL = 1.0, BzL = 1.0;
        } else if ((xL > 0.0) and (xL < 0.2)) {
          BxL = 1.0, ByL = zx, BzL = zx;
        } else { // if (xL >=  0.2) {
          BxL = 1.0, ByL = -1.0, BzL = -1.0;
        }

        CCTK_REAL ExL = 0.0;
        CCTK_REAL EyL = 0.5;
        CCTK_REAL EzL = -0.5;

// Equation (99), arXiv:1310.3274v2. Calculate the velocity.
        CCTK_REAL B2L = BxL*BxL + ByL*ByL + BzL*BzL;

        vel[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = (EyL*BzL - EzL*ByL)/B2L;
        vel[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = (EzL*BxL - ExL*BzL)/B2L;
        vel[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] = (ExL*ByL - EyL*BxL)/B2L;

      }
}


// Diagnostic for outputting the difference B^2 - E^2
extern "C" void GiRaFFEfood_OutputB2mE2(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int imin=0,jmin=0,kmin=0;
  int imax=cctk_lsh[0],jmax=cctk_lsh[1],kmax=cctk_lsh[2];

  const CCTK_REAL fourPI = 4.0*M_PI;

#pragma omp parallel for
  for(int k=kmin;k<kmax;k++)
    for(int j=jmin;j<jmax;j++)
      for(int i=imin;i<imax;i++) {
        int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        CCTK_REAL BxL = Bx[index];
        CCTK_REAL ByL = By[index];
        CCTK_REAL BzL = Bz[index];

        CCTK_REAL B2 = BxL*BxL + ByL*ByL + BzL*BzL;

        CCTK_REAL SxL = mhd_st_x[index];
        CCTK_REAL SyL = mhd_st_y[index];
        CCTK_REAL SzL = mhd_st_z[index];

        CCTK_REAL ExL = fourPI*(ByL*SzL-BzL*SyL)/B2;
        CCTK_REAL EyL = fourPI*(BzL*SxL-BxL*SzL)/B2;
        CCTK_REAL EzL = fourPI*(BxL*SyL-ByL*SxL)/B2;

        CCTK_REAL E2 = ExL*ExL + EyL*EyL + EzL*EzL;

        B2mE2[index] = B2 - E2;

        Ex[index] = ExL;
        Ey[index] = EyL;
        Ez[index] = EzL;

      }
}

