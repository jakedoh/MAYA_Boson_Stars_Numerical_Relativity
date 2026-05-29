// Please consult GiRaFFEfood/README for the list of authors and current maintainers //
/*
 * This program implements the initial configuration called the Aligned Rotator,
 * a time-dependent toy model of a pulsar magnetosphere, on flat space background.
 * Both the star and the magnetic field rotate with constant angular velocity, and
 * the dipolar magnetic field threads the surface of the star, extending to infinity.
 * The light cylinder is at 5 radii away from the star's center, as in arXiv:astro-ph/0603147.
 * The vector potential is calculated with equation (158) arXiv:1310.3274v2.
 * The velocity within the light cylinder is given by equation (157), arXiv:1310.3274v2.
 * Outside the light cylinder the velocity is set to zero.
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

void GiRaFFEfood_alignA_i(CCTK_REAL &AxL, CCTK_REAL &AyL, CCTK_REAL &AzL,
                          CCTK_REAL &xL, CCTK_REAL &yL, CCTK_REAL &zL);

extern "C" void GiRaFFEfood_AlignedRotator(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  using namespace std;

  int imin = 0;
  int jmin = 0;
  int kmin = 0;

  int imax = cctk_lsh[0];
  int jmax = cctk_lsh[1];
  int kmax = cctk_lsh[2];

  CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  CCTK_REAL dz = CCTK_DELTA_SPACE(2);

#pragma omp parallel for
  for(int k=kmin;k<kmax;k++)
    for(int j=jmin;j<jmax;j++)
      for(int i=imin;i<imax;i++) {
        int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        //Initialize the vector potential
        CCTK_REAL AxL = 0.0;
        CCTK_REAL AyL = 0.0;
        CCTK_REAL AzL = 0.0;

        //Initialize the velocity:
        CCTK_REAL VxL = 0.0;
        CCTK_REAL VyL = 0.0;
        CCTK_REAL VzL = 0.0;

        //Cartesian coordinates
        CCTK_REAL xL = x[index];
        CCTK_REAL yL = y[index];
        CCTK_REAL zL = z[index];
        //Staggered coordinates
        CCTK_REAL xh = xL + 0.5*dx;
        CCTK_REAL yh = yL + 0.5*dy;
        CCTK_REAL zh = zL + 0.5*dz;

        //Avoid NaNs
        CCTK_REAL fudge_factor = 1e-5;
        if(xL == 0) {xL = fudge_factor*dx;}
        if(yL == 0) {yL = fudge_factor*dy;}
        if(zL == 0) {zL = fudge_factor*dz;}

        //Call the function for vector potential. Attention: the vector potential is staggered.
        CCTK_REAL dummy0,dummy1,dummy2;
        //Ax(i,j+0.5,k+0.5)
        GiRaFFEfood_alignA_i(AxL,dummy1,dummy2, xL,yh,zh);
        //Ay(i+0.5,j,k+0.5)
        GiRaFFEfood_alignA_i(dummy0,AyL,dummy2, xh,yL,zh);
        //Az(i+0.5,j+0.5,k)
        GiRaFFEfood_alignA_i(dummy0,dummy1,AzL, xh,yh,zL);

        //Set the vector potential:
        Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = AxL;
        Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = AyL;
        Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] = AzL;

        //Equation (157) arXiv:1310.3274v2. Set the velocity field:
        //$\textbf{v} = \Omega\textbf{e_z}\times\textbf{r}$
        VxL = -Omega_aligned_rotator*yL;
        VyL =  Omega_aligned_rotator*xL;
        VzL =  0.0;
        // "In the exterior, the 3-velocity is set to 0." <- Statement just below Eq. 157.
        if(sqrt(xL*xL + yL*yL + zL*zL) > R_NS_aligned_rotator) {
          VxL = 0.0;
          VyL = 0.0;
          VzL = 0.0;
        }

        vel[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = VxL;
        vel[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = VyL;
        vel[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] = VzL;

      }
}

//Equation (158) arXiv:1310.3274v2. Calculate the vector potential
void GiRaFFEfood_alignA_i(CCTK_REAL &AxL, CCTK_REAL &AyL, CCTK_REAL &AzL,
                          CCTK_REAL &xL, CCTK_REAL &yL, CCTK_REAL &zL) {

  DECLARE_CCTK_PARAMETERS;
  //$A_phi=\mu*\omega^2/r^3$ where $\mu=B_p*R_{NS}^3/2$ (Eq.158)
  // Just paste the output from Mathematica, to minimize chance of human error.
  /*
    omegaL = Sqrt[xL*xL + yL*yL];
    rL = Sqrt[xL*xL + yL*yL + zL*zL];

    mu = Bpalignedrotator*RNSalignedrotator^3/2.0;
    Aph = mu*omegaL^2/rL^3;

    dphdx = -yL/(xL*xL + yL*yL);
    dphdy = xL/(xL*xL + yL*yL);
    dphdz = 0.0;

    FullSimplify[Aph]

    //Transform into Cartesian coordinates:
    AxL = Aph*dphdx;
    AyL = Aph*dphdy;
    AzL = Aph*dphdz;
    FullSimplify[AxL]
    FullSimplify[AyL]
    FullSimplify[AzL]

    // Do these one at a time to get lines of code.
    CForm[FullSimplify[AxL]]

  */

  CCTK_REAL Bpalignedrotator = B_p_aligned_rotator;
  CCTK_REAL RNSalignedrotator3 = R_NS_aligned_rotator*R_NS_aligned_rotator*R_NS_aligned_rotator;

  CCTK_REAL rL = sqrt(xL*xL + yL*yL + zL*zL);
  CCTK_REAL rL3 = rL*rL*rL;

  AxL = (-0.5*Bpalignedrotator*RNSalignedrotator3*yL)/rL3;
  AyL = (0.5*Bpalignedrotator*RNSalignedrotator3*xL)/rL3;
  AzL = 0;

}
