// Please consult GiRaFFEfood/README for the list of authors and current maintainers //
/*
 * This program implements the initial configuration called Alfven Wave,
 * which is a characteristic wave of an FFE system, and represents
 * an exact 1D test, first proposed by S.S. Komissarov in arXiv:0202447.
 * The vector potential is calculated with formula (106), arXiv:1310.3274v2.
 * We calculate the magnetic and electric fields with eqs. (102, 103, 104, 105),
 * arXiv:1310.3274v2. Then we use eq. (99), arxiv:1310.3274v2 to obtain the velocity.
 */

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

extern "C" void GiRaFFEfood_AlfvenWave(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  using namespace std;

  CCTK_INT imin = 0, jmin = 0, kmin = 0;
  CCTK_INT imax = cctk_lsh[0], jmax = cctk_lsh[1], kmax = cctk_lsh[2];

  CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  CCTK_REAL dy = CCTK_DELTA_SPACE(1);

  //Lorentz factor,\gamma_\mu=1/\sqrt{1-\mu^2}
  CCTK_REAL gmu = 1.0/sqrt(1 - wave_speed*wave_speed);

#pragma omp parallel for
  for(int k=kmin;k<kmax;k++)
    for(int j=jmin;j<jmax;j++)
      for(int i=imin;i<imax;i++) {
        int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        //Getting the staggered coordinates. The vector potential is staggered.
        CCTK_REAL xL_plus_half = x[index] + 0.5*dx;
        CCTK_REAL yL_plus_half = y[index] + 0.5*dy;

        //Equation (106), arXiv:1310.3274v2. Only Ay and Az are nonzero here.
        //Calculating Ay:
        CCTK_REAL AyL;
        //$g(x)=\cos(5\pi\gamma_{\mu}x)/\pi$
        CCTK_REAL gx = cos(5.0*M_PI*gmu*xL_plus_half)/M_PI;

        //$A_y=\gamma_\mu x-0.015$ if $x \leq -0.1/\gamma_\mu$
        if ( xL_plus_half <= -0.1/gmu) {
          AyL = gmu*xL_plus_half - 0.015;
          //$A_y=1.15\gamma_\mu x-0.03 g(x)$ if $-0.1/\gamma_\mu \leq x \leq 0.1/\gamma_\mu$
        } else if ((xL_plus_half > -0.1/gmu) and (xL_plus_half < 0.1/gmu)) {
          AyL = 1.15*gmu*xL_plus_half - 0.03*gx;
          //$A_y=1.3\gamma_\mu x-0.015$ if $x \geq 0.1/\gamma_\mu$
        } else { // if ( xL_plus_half >=  0.1/gmu) {
          AyL = 1.3*gmu*xL_plus_half - 0.015;
        }

        //Calculating Az:
        CCTK_REAL AzL;
        //$A_z=y-\gamma_\mu(1-\mu)x$
        AzL = yL_plus_half - gmu*(1 -wave_speed)*xL_plus_half;

        //Fill in the interface variables:
        Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = 0.0;
        Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = AyL;
        Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] = AzL;

        //Initialize the velocity
        //Equation (102) arXiv:1310.3274v2. First calculate B'
        CCTK_REAL Bxp = 1.0;
        CCTK_REAL Byp = 1.0;
        CCTK_REAL Bzp;

        // We want unstaggered v, so we want unstaggered B and E.
        CCTK_REAL xp = x[index]*gmu;
        CCTK_REAL f_xp = 1.0 + sin(5.0*M_PI*xp);

        if (xp <= -0.1) {
          Bzp = 1.0;
        } else if ((xp > -0.1) and (xp < 0.1)) {
          Bzp = 1.0 + 0.15*f_xp;
        } else { // if (xp >=  0.1) {
          Bzp = 1.3;
        }

        //Equation (103) arXiv:1310.3274v2. Next calculate E'
        CCTK_REAL Exp = -Bzp;
        CCTK_REAL Eyp = 0.0;
        CCTK_REAL Ezp = 1.0;

        //Equation (104) arXiv:1310.3274v2. Calculate B
        CCTK_REAL BxL = Bxp;
        CCTK_REAL ByL = gmu*(Byp - wave_speed*Ezp);
        CCTK_REAL BzL = gmu*(Bzp + wave_speed*Eyp);

        //Equation (105) arXiv:1310.3274v2. Calculate E
        CCTK_REAL ExL = Exp;
        CCTK_REAL EyL = gmu*(Eyp + wave_speed*Bzp);
        CCTK_REAL EzL = gmu*(Ezp - wave_speed*Byp);

        //Equation (99) arXiv:1310.3274v2. Now calculate the velocity for flat space.
        CCTK_REAL B2L = BxL*BxL + ByL*ByL + BzL*BzL;

        vel[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = (EyL*BzL - EzL*ByL)/B2L;
        vel[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = (EzL*BxL - ExL*BzL)/B2L;
        vel[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] = (ExL*ByL - EyL*BxL)/B2L;

      }
}

