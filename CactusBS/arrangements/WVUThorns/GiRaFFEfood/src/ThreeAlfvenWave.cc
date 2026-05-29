// Please consult GiRaFFEfood/README for the list of authors and current maintainers //
/*
 * This program implements the initial configuration called Three Alfven Wave,
 * a 1D test which consists of a superposition of a stationary Alfven wave and
 * two fast discontinuities: a right and left-going fast wave.
 * The vector potential is calculated with formula (114), arXiv:1310.3274v2
 * We calculate the magnetic and electric fields with eqs. (116, 117, 118),
 * arXiv:1310.3274v2. Then we use eq. (99) from the same paper, to obtain the velocity.
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

using namespace std;

void GiRaFFEfood_ThreeAlfvenWave_ExactSolution_StationaryWave_B_and_E(CCTK_REAL x, CCTK_REAL &Bx,CCTK_REAL &By,CCTK_REAL &Bz,CCTK_REAL &Ex,CCTK_REAL &Ey,CCTK_REAL &Ez) {
  if(x<=0.0) {
    //Equation (116) stationary Alfven wave
    Bx = 1.0, By = 1.0, Bz = 2.0;
    Ex = -1.0; Ey = 1.0; Ez = 0.0;
  } else {
    //Equation (116) stationary Alfven wave
    Bx = 1.0, By = 1.5, Bz = 2.0;
    Ex = -1.5; Ey = 1.0; Ez = 0.0;
  }
}

void GiRaFFEfood_ThreeAlfvenWave_ExactSolution_RightMovingWave_B_and_E(CCTK_REAL x, CCTK_REAL &Bx,CCTK_REAL &By,CCTK_REAL &Bz,CCTK_REAL &Ex,CCTK_REAL &Ey,CCTK_REAL &Ez) {
  if(x<=0.0) {
    //Equation (117) right-going Alfven wave
    Bx = 0.0; By = 0.0; Bz = 0.0;
    Ex = 0.0; Ey = 0.0; Ez = 0.0;
  } else {
    //Equation (117) right-going Alfven wave
    Bx = 0.0; By = 1.5; Bz = 1.0;
    Ex = 0.0; Ey = 1.0; Ez = -1.5;
  }
}

void GiRaFFEfood_ThreeAlfvenWave_ExactSolution_LeftMovingWave_B_and_E(CCTK_REAL x, CCTK_REAL &Bx,CCTK_REAL &By,CCTK_REAL &Bz,CCTK_REAL &Ex,CCTK_REAL &Ey,CCTK_REAL &Ez) {
  if(x<=0.0) {
    //Equation (118) left-going Alfven wave
    Bx = 0.0; By = 0.5; Bz = 1.5;
    Ex = 0.0; Ey = -1.5; Ez = 0.5;
  } else {
    //Equation (118) left-going Alfven wave
    Bx = 0.0; By = 0.0; Bz = 0.0;
    Ex = 0.0; Ey = 0.0; Ez = 0.0;
  }
}

extern "C" void GiRaFFEfood_ThreeAlfvenWave(CCTK_ARGUMENTS) {

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

        //Heaviside Step Function
        CCTK_REAL Hmx,Hpx;

        //Getting the staggered x-coordinate.
        CCTK_REAL xL_plus_half = x[index] + 0.5*dx;
        if ( xL_plus_half < 0.0) {
          Hpx = 0.0, Hmx = 1.0;
        } else if ( xL_plus_half == 0.0) {
          Hpx = 0.5, Hmx = 0.5;
        } else { // if ( xL_plus_half > 0.0) {
          Hpx = 1.0, Hmx = 0.0;
        }

        //Calculate Ay(i+0.5,j,k+0.5)
        //$A_y=3.5x H(-x)+3.0x H(x)$
        CCTK_REAL AyL = 3.5*xL_plus_half*Hmx + 3.0*xL_plus_half*Hpx;

        //Calculate Az(i+0.5,j+0.5,k)
        //Getting the staggered y-coordinate.
        CCTK_REAL yL_plus_half = y[index] + 0.5*dy;
        //$A_z=y-1.5x H(-x)-3.0x H(x)$
        CCTK_REAL AzL = yL_plus_half - 1.5*xL_plus_half*Hmx - 3.0*xL_plus_half*Hpx;

        //The vector potential is staggered.
        //Only Ay and Az are nonzero here.
        Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = 0.0;
        Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = AyL;
        Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] = AzL;

        CCTK_REAL xL = x[index];

        CCTK_REAL Bxa,Bya,Bza,Exa,Eya,Eza;
        GiRaFFEfood_ThreeAlfvenWave_ExactSolution_StationaryWave_B_and_E (xL, Bxa,Bya,Bza,Exa,Eya,Eza);
        CCTK_REAL Bxp,Byp,Bzp,Exp,Eyp,Ezp;
        GiRaFFEfood_ThreeAlfvenWave_ExactSolution_RightMovingWave_B_and_E(xL, Bxp,Byp,Bzp,Exp,Eyp,Ezp);
        CCTK_REAL Bxm,Bym,Bzm,Exm,Eym,Ezm;
        GiRaFFEfood_ThreeAlfvenWave_ExactSolution_LeftMovingWave_B_and_E (xL, Bxm,Bym,Bzm,Exm,Eym,Ezm);

        // Equation(115)
        CCTK_REAL BxL = Bxa + Bxp + Bxm;
        CCTK_REAL ExL = Exa + Exp + Exm;
        CCTK_REAL ByL = Bya + Byp + Bym;
        CCTK_REAL EyL = Eya + Eyp + Eym;
        CCTK_REAL BzL = Bza + Bzp + Bzm;
        CCTK_REAL EzL = Eza + Ezp + Ezm;

        // Flat-space B^2
        CCTK_REAL B2 = BxL*BxL + ByL*ByL + BzL*BzL;

        // Equation (99)
        CCTK_REAL vxL = (EyL*BzL - EzL*ByL)/B2;
        CCTK_REAL vyL = (EzL*BxL - ExL*BzL)/B2;
        CCTK_REAL vzL = (ExL*ByL - EyL*BxL)/B2;
        vel[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = vxL;
        vel[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = vyL;
        vel[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] = vzL;
      }
}

/*  Following the SplitMonopole routine to calculate the error between the exact Bi and the primitive Bi. */
extern "C" void GiRaFFEfood_ErrorBThreeWaves(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT imin = 0, jmin = 0, kmin = 0;
  CCTK_INT imax = cctk_lsh[0]-0, jmax = cctk_lsh[1]-0, kmax = cctk_lsh[2]-0;

  CCTK_REAL errorsumx = 0.0;
  CCTK_REAL errorsumy = 0.0;
  CCTK_REAL errorsumz = 0.0;

#pragma omp parallel for reduction(+:errorsumx,errorsumy,errorsumz) schedule(static)
  for(int k=kmin;k<kmax;k++)
    for(int j=jmin;j<jmax;j++)
      for(int i=imin;i<imax;i++) {
        int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        CCTK_REAL xL = x[index];


        // Notice the time shifts added!
        CCTK_REAL Bxa,Bya,Bza,Exa,Eya,Eza;
        GiRaFFEfood_ThreeAlfvenWave_ExactSolution_StationaryWave_B_and_E (xL,           Bxa,Bya,Bza,Exa,Eya,Eza);
        CCTK_REAL Bxp,Byp,Bzp,Exp,Eyp,Ezp;
        GiRaFFEfood_ThreeAlfvenWave_ExactSolution_RightMovingWave_B_and_E(xL-cctk_time, Bxp,Byp,Bzp,Exp,Eyp,Ezp);
        CCTK_REAL Bxm,Bym,Bzm,Exm,Eym,Ezm;
        GiRaFFEfood_ThreeAlfvenWave_ExactSolution_LeftMovingWave_B_and_E (xL+cctk_time, Bxm,Bym,Bzm,Exm,Eym,Ezm);

        // Equation(115)
        CCTK_REAL BxL = Bxa + Bxp + Bxm;
        CCTK_REAL ByL = Bya + Byp + Bym;
        CCTK_REAL BzL = Bza + Bzp + Bzm;

        exactBx_ThreeWaves[index] = BxL;
        exactBy_ThreeWaves[index] = ByL;
        exactBz_ThreeWaves[index] = BzL;

        //Fill in the interface variables:
        delBx_ThreeWaves[index] = Bx[index] - BxL;
        delBy_ThreeWaves[index] = By[index] - ByL;
        delBz_ThreeWaves[index] = Bz[index] - BzL;

        errorsumx += (delBx_ThreeWaves[index])*(delBy_ThreeWaves[index]);
        errorsumy += (delBy_ThreeWaves[index])*(delBy_ThreeWaves[index]);
        errorsumz += (delBz_ThreeWaves[index])*(delBz_ThreeWaves[index]);
      }

  // Useful diagnostic:
  /*
    CCTK_REAL dx = CCTK_DELTA_SPACE(0);
    CCTK_REAL dy = CCTK_DELTA_SPACE(1);
    CCTK_REAL dz = CCTK_DELTA_SPACE(2);
    CCTK_REAL norm2Bx = sqrt(errorsumx*dx*dy*dz);
    CCTK_REAL norm2By = sqrt(errorsumy*dx*dy*dz);
    CCTK_REAL norm2Bz = sqrt(errorsumz*dx*dy*dz);
  */
  CCTK_VInfo(CCTK_THORNSTRING,
             "norm2BxByBz: %e %e %e | %e\n",errorsumx,errorsumy,errorsumz,delBx_ThreeWaves[CCTK_GFINDEX3D(cctkGH,3,3,3)]);
}
