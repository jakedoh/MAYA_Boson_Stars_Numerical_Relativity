// Please consult GiRaFFEfood/README for the list of authors and current maintainers //
/*
 * This program implements the initial configuration called Fast Wave, which is a
 * 1D characteristic wave of force-free electrodynamics, propagating at the speed of light.
 * The vector potential is calculated with formula (100), arXiv:1310.3274v2.
 * We calculate the magnetic and electric fields with eqs. (97, 98), arXiv:1310.3274v2.
 * Then we use eq. (99) from the same paper to obtain the velocity.
 **/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstdlib>

using namespace std;

extern "C" void GiRaFFEfood_FastWave(CCTK_ARGUMENTS) {

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

        //Only Az is set in this test!
        //Calculate Az(i+0.5,j+0.5,k) - the vector potential is staggered
        CCTK_REAL xL_plus_half = x[index] + 0.5*dx;
        CCTK_REAL yL_plus_half = y[index] + 0.5*dy;

        //Equation (100), arXiv:1310.3274v2
        CCTK_REAL AzL = yL_plus_half;
        if (xL_plus_half <= -0.1) {
          AzL += (-xL_plus_half - 0.0075);
        } else if ((xL_plus_half > -0.1) and (xL_plus_half < 0.1)) {
          AzL += (0.75*xL_plus_half*xL_plus_half - 0.85*xL_plus_half);
        } else if (xL_plus_half >= 0.1) {
          AzL += (-0.7*xL_plus_half - 0.0075);
        }

        //Fill in the interface variables:
        Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = 0.0;
        Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = 0.0;
        Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] = AzL;

        // Initialize the velocity:
        //Equations (97), arXiv:1310.3274v2. Calculate B.
        CCTK_REAL BxL = 1.0;

        CCTK_REAL ByL;
        CCTK_REAL xL = x[index];
        if (xL <= -0.1) {
          ByL = 1.0;
        } else if ((xL > -0.1) and (xL < 0.1)) {
          ByL = 1.0 - 1.5*(xL + 0.1);
        } else { // if (xL>=0.1) {
          ByL = 0.7;
        }

        CCTK_REAL BzL = 0.0;

        //Equations (98), arXiv:1310.3274v2. Calculate E.
        CCTK_REAL ExL = 0.0;
        CCTK_REAL EyL = 0.0;
        CCTK_REAL EzL = - ByL;

        // Equation (99), arXiv:1310.3274v2.
        CCTK_REAL B2L = BxL*BxL + ByL*ByL + BzL*BzL;

        vel[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = (EyL*BzL - EzL*ByL)/B2L;
        vel[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = (EzL*BxL - ExL*BzL)/B2L;
        vel[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] = (ExL*ByL - EyL*BxL)/B2L;

      }
}

