// Please consult GiRaFFEfood/README for the list of authors and current maintainers //
/*
 * This program implements the initial configuration called Degenerate Alfven Wave,
 * a 1D flat space test consisting of the superposition of right and left-going
 * Alfven waves, possessing the same wave speed.
 * The vector potential is calculated with formulas (109) and (110), arXiv:1310.3274v2.
 * We calculate the magnetic and electric fields with eqs. (108, 107, 103, 104, 105),
 * arXiv:1310.3274v2. Then we use eq. (99) from the same paper, to obtain the velocity.
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

using namespace std;

extern "C" void GiRaFFEfood_DegenAlfvenWave(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT imin = 0, jmin = 0, kmin = 0;
  CCTK_INT imax = cctk_lsh[0], jmax = cctk_lsh[1], kmax = cctk_lsh[2];

  CCTK_REAL dx = CCTK_DELTA_SPACE(0);

  //Lorentz factor, $\gamma_\mu=1/\sqrt{1-\mu^2}$
  CCTK_REAL gmu = 1.0/sqrt(1.0 - wave_speed*wave_speed);

#pragma omp parallel for
  for(int k=kmin;k<kmax;k++)
    for(int j=jmin;j<jmax;j++)
      for(int i=imin;i<imax;i++) {
        int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        //The vector potential is staggered. Getting the staggered x-coordinate.
        CCTK_REAL xL_plus_half  = x[index] + 0.5*dx;

        //Equation (109), arXiv:1310.3274v2. Calculate Ay(i+0.5,j,k+0.5).
        CCTK_REAL AyL;
        //$h_1(x)=\cos[2.5\pi(\gamma_\mu x+0.1)]$
        CCTK_REAL h1x = cos(2.5*M_PI*(gmu*xL_plus_half+0.1));

        //$A_y=-0.8/\pi$ if $x \leq -0.1/\gamma_\mu$
        if (xL_plus_half <= -0.1/gmu) {
          AyL = - 0.8/M_PI;
          //$A_y=-(0.8/\pi)h_1(x)$ if $-0.1/\gamma_\mu \leq x \leq 0.1/\gamma_\mu$
          //} else if ((xL_plus_half > -0.1/gmu) and (xL_plus_half = 0.1/gmu)) {
        } else if ((xL_plus_half > -0.1/gmu) and (xL_plus_half <= 0.1/gmu)) {
          AyL = - (0.8/M_PI)*h1x;
          //$A_y=2(\gamma_\mu x-0.1)$ if $x \geq 0.1/\gamma_\mu$
        } else { // if (xL_plus_half > 0.1/gmu) {
          AyL = 2.0*(gmu*xL_plus_half - 0.1);
        }

        //Equation (110), arXiv:1310.3274v2. Calculate Az(i+0.5,j+0.5,k).
        CCTK_REAL AzL;
        //$h_2(x)=\sin[2.5\pi(\gamma_\mu x+0.1)]$
        CCTK_REAL h2x = sin(2.5*M_PI*(gmu*xL_plus_half+0.1));

        //$A_z=-2(\gamma_\mu x+0.1)$ if $x \leq -0.1/\gamma_\mu$
        if (xL_plus_half <= -0.1/gmu) {
          AzL = - 2.0*(gmu*xL_plus_half + 0.1);
          //$A_z=-(0.8/\pi)h_2(x)$ if $-0.1/\gamma_\mu \leq x \leq 0.1/\gamma_\mu$
          //} else if ((xL_plus_half > -0.1/gmu) and (xL_plus_half < 0.1/gmu)) {
        } else if ((xL_plus_half > -0.1/gmu) and (xL_plus_half <= 0.1/gmu)) {
          AzL = - (0.8/M_PI)*h2x;
          //$A_z=-0.8/\pi$ if $x \geq 0.1/\gamma_\mu$
        } else { // if ( xL_plus_half >  0.1/gmu) {
          AzL = - 0.8/M_PI;
        }

        //Fill in the interface variables:
        Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = 0.0;
        Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = AyL;
        Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] = AzL;

        // Initialize the velocity:
        // We want unstaggered v, so we want unstaggered B and E.
        CCTK_REAL xp = x[index]*gmu;
        CCTK_REAL phip;
        //Equation (108), arXiv:1310.3274v2. Calculate \phi(x')
        if (xp <= -0.1) {
          phip = 0.0;
        } else if ((xp > -0.1) and (xp <= 0.1)) {
          phip = 2.5*M_PI*(xp + 0.1);
        } else { // if (xp >  0.1) {
          phip = 0.5*M_PI;
        }

        //Equation (107), arXiv:1310.3274v2. Calculate E' and B'
        CCTK_REAL Exp = 0.0;
        CCTK_REAL Bxp = 0.0;
        CCTK_REAL Byp = 2.0*cos(phip);
        CCTK_REAL Bzp = 2.0*sin(phip);

        //Equation (103), arXiv:1310.3274v2. <---should we follow Equation (107) here?
        CCTK_REAL Eyp = 0.0;
        CCTK_REAL Ezp = 0.0;

        //Equation (104), arXiv:1310.3274v2. Calculate B
        CCTK_REAL BxL = Bxp;
        CCTK_REAL ByL = gmu*(Byp - wave_speed*Ezp);
        CCTK_REAL BzL = gmu*(Bzp + wave_speed*Eyp);

        //Equation (105), arXiv:1310.3274v2. Calculate E
        CCTK_REAL ExL = Exp;
        CCTK_REAL EyL = gmu*(Eyp + wave_speed*Bzp);
        CCTK_REAL EzL = gmu*(Ezp - wave_speed*Byp);

        // Equation (99), arXiv:1310.3274v2. Calculate the velocity v
        CCTK_REAL B2L = BxL*BxL + ByL*ByL + BzL*BzL;

        vel[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = (EyL*BzL - EzL*ByL)/B2L;
        vel[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = (EzL*BxL - ExL*BzL)/B2L;
        vel[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] = (ExL*ByL - EyL*BxL)/B2L;

      }
}

// Diagnostic for outputting flat-space electric field

extern "C" void GiRaFFEfood_OutputFlatSpaceEfield(CCTK_ARGUMENTS) {

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

        Ex[index] = fourPI*(ByL*SzL-BzL*SyL)/B2;
        Ey[index] = fourPI*(BzL*SxL-BxL*SzL)/B2;
        Ez[index] = fourPI*(BxL*SyL-ByL*SxL)/B2;

      }
}


