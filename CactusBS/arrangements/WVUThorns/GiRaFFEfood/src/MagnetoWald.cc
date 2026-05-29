// Please consult GiRaFFEfood/README for the list of authors and current maintainers //
/*
 * This program implements the initial configuration called Magnetic Wald,
 * describing a spinning black hole immersed in uniform magnetic field.
 * It does not have an analytic solution, but a steady state is reached,
 * with a current sheet forming in the equatorial plane.
 * The vector potential is calculated with formulas (156), arXiv:1310.3274v2
 * The velocity is simply given by v^i = -\beta^i, then transformed to the
 * Valencia formalism, arXiv:1501.07276v2, eq.(11)
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

void GiRaFFEfood_MagnW_V_i(CCTK_REAL &VxL, CCTK_REAL &VyL, CCTK_REAL &VzL,
                           CCTK_REAL &betarL, CCTK_REAL &betath, CCTK_REAL &betaph,
                           CCTK_REAL xL, CCTK_REAL yL, CCTK_REAL zL);

void GiRaFFEfood_MagnW_A_i(CCTK_REAL &AxL, CCTK_REAL &AyL, CCTK_REAL &AzL,
                           CCTK_REAL xL, CCTK_REAL yL, CCTK_REAL zL);

void GiRaFFEfood_MagnW_Sph_to_Cart_XformVector(CCTK_REAL &VxL, CCTK_REAL &VyL, CCTK_REAL &VzL,
                                               CCTK_REAL &VrL, CCTK_REAL &Vth, CCTK_REAL &Vph,
                                               CCTK_REAL xL, CCTK_REAL yL, CCTK_REAL zL);

extern "C" void GiRaFFEfood_MagnetoWald(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT imin = 0, jmin = 0, kmin = 0;
  CCTK_INT imax = cctk_lsh[0], jmax = cctk_lsh[1], kmax = cctk_lsh[2];

  CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  CCTK_REAL dz = CCTK_DELTA_SPACE(2);

  //#pragma omp parallel for
  for(int k=kmin;k<kmax;k++)
    for(int j=jmin;j<jmax;j++)
      for(int i=imin;i<imax;i++) {
        int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        //Initialize the vector potential
        CCTK_REAL AxL = 0.0;
        CCTK_REAL AyL = 0.0;
        CCTK_REAL AzL = 0.0;

        //Initialize the velocity
        CCTK_REAL VxL = 0.0;
        CCTK_REAL VyL = 0.0;
        CCTK_REAL VzL = 0.0;

        //Cartesian coordinates
        CCTK_REAL xL = x[index];
        CCTK_REAL yL = y[index];
        CCTK_REAL zL = z[index];

        //Avoid NaNs
        CCTK_REAL fudge_factor = 1e-5;
        if(fabs(xL) < fudge_factor*dx) {xL = fudge_factor*dx;}
        if(fabs(yL) < fudge_factor*dy) {yL = fudge_factor*dy;}
        if(fabs(zL) < fudge_factor*dz) {zL = fudge_factor*dz;}

        //Staggered coordinates
        CCTK_REAL xh = xL + 0.5*dx;
        CCTK_REAL yh = yL + 0.5*dy;
        CCTK_REAL zh = zL + 0.5*dz;

        //I need the metric at staggered coordinates!
        //Calculate the potentials
        CCTK_REAL dummy0,dummy1,dummy2;
        //The vector potential is staggered. We can't use the outside metric
        // grhrhL,grhthL,grhphL,gththL,gthphL,gphphL,alphaL,betarL,
        //Ax(i,j+0.5,k+0.5)
        GiRaFFEfood_MagnW_A_i(AxL, dummy1, dummy2, xL, yh, zh);
        //Ay(i+0.5,j,k+0.5)
        GiRaFFEfood_MagnW_A_i(dummy0, AyL, dummy2, xh, yL, zh);
        //Az(i+0.5,j+0.5,k)
        GiRaFFEfood_MagnW_A_i(dummy0, dummy1, AzL, xh, yh, zL);

        //Fill in the interface variables:
        Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = AxL;
        Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = AyL;
        Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] = AzL;

        // Transform from spherical to cartesian, the velocity is not staggered
        CCTK_REAL betarL = SKSbetar[index];
        CCTK_REAL betath = SKSbetath[index];
        CCTK_REAL betaph = SKSbetaph[index];
        GiRaFFEfood_MagnW_V_i(VxL,VyL,VzL,
                              betarL,betath,betaph,
                              xL,yL,zL);

        CCTK_REAL alpL = alp[index];
        CCTK_REAL betaxL = betax[index];
        CCTK_REAL betayL = betay[index];
        CCTK_REAL betazL = betaz[index];

        //Bypass transformation:
        //VxL =-betaxL;
        //VyL =-betayL;
        //VzL =-betazL;

        CCTK_REAL ETVx = (VxL + betaxL)/alpL;
        CCTK_REAL ETVy = (VyL + betayL)/alpL;
        CCTK_REAL ETVz = (VzL + betazL)/alpL;

        //Fill in the interface variables:
        vel[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = ETVx; //0.0;
        vel[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = ETVy; //0.0;
        vel[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] = ETVz; //0.0;
        //Should the velocity be zero: v^i = (-beta^i + beta^i)/alpha???

      }
}

void GiRaFFEfood_MagnW_A_i(CCTK_REAL &AxL, CCTK_REAL &AyL, CCTK_REAL &AzL,
                           CCTK_REAL xL, CCTK_REAL yL, CCTK_REAL zL) {
  DECLARE_CCTK_PARAMETERS;

  // Input parameters
  CCTK_REAL aL = BH_spin;
  CCTK_REAL B0 = Wald_B0;
  CCTK_REAL mL = BH_mass;

  // Spherical coordinates
  CCTK_REAL rL = sqrt(xL*xL + yL*yL + zL*zL);
  CCTK_REAL costh = zL/rL;
  CCTK_REAL costh2 = costh*costh;
  CCTK_REAL sinth2 = 1.0 - costh2;

  // Kerr-Schild spherical coordinates
  CCTK_REAL r_KS = rL + KerrSchild_radial_shift;
  CCTK_REAL rho2 = r_KS*r_KS + aL*aL*costh2;

  //I. Regular spherical Kerr-Schild Metric. Should really read this from ShiftedKerrSchild thorn!

  //http://relativity.livingreviews.org/open?pubNo=lrr-2000-5&page=articlesu8.html
  CCTK_REAL grhrhL = 1.0 + 2.0*mL*r_KS/rho2;
  // unused; kept for completeness: CCTK_REAL grhthL = 0.0;
  CCTK_REAL grhphL =-grhrhL*aL*sinth2;
  // unused; kept for completeness: CCTK_REAL gththL = rho2;
  CCTK_REAL gthphL = 0.0;
  CCTK_REAL gphphL = (r_KS*r_KS + aL*aL + 2.0*mL*r_KS/rho2*aL*aL*sinth2)*sinth2;

  //https://arxiv.org/pdf/1310.3575.pdf, eq. (2a), (2c)
  CCTK_REAL gtrhL = 2.0*mL*r_KS/rho2;
  CCTK_REAL gtthL = 0.0;
  CCTK_REAL gtphL =-gtrhL*aL*sinth2;

  //II. Regular spherical Kerr-Schild Metric, in Boyer-Lindquist coordinates

  //The 4-vector potential in spherical coordinates: eq.(156), 1310.3274.pdf
  CCTK_REAL Arh = 0.5*B0*(grhphL + 2.0*aL*gtrhL);
  CCTK_REAL Ath = 0.5*B0*(gthphL + 2.0*aL*gtthL);
  CCTK_REAL Aph = 0.5*B0*(gphphL + 2.0*aL*gtphL);

  //The Jacobian from the Spherical to Cartesian:
  CCTK_REAL drL_dx = xL/rL;
  CCTK_REAL drL_dy = yL/rL;
  CCTK_REAL drL_dz = zL/rL;

  CCTK_REAL dth_dx = xL*zL/(rL*rL*sqrt(xL*xL + yL*yL));
  CCTK_REAL dth_dy = yL*zL/(rL*rL*sqrt(xL*xL + yL*yL));
  CCTK_REAL dth_dz =-sqrt(xL*xL + yL*yL)/(rL*rL);

  CCTK_REAL dph_dx =-yL/(xL*xL + yL*yL);
  CCTK_REAL dph_dy = xL/(xL*xL + yL*yL);
  CCTK_REAL dph_dz = 0.0;

  //Transform to Cartesian coordinates:
  AxL = Arh*drL_dx + Ath*dth_dx + Aph*dph_dx;
  AyL = Arh*drL_dy + Ath*dth_dy + Aph*dph_dy;
  AzL = Arh*drL_dz + Ath*dth_dz + Aph*dph_dz;

}

void GiRaFFEfood_MagnW_Sph_to_Cart_XformVector(CCTK_REAL &VxL, CCTK_REAL &VyL, CCTK_REAL &VzL,
                                               CCTK_REAL &VrL, CCTK_REAL &Vth, CCTK_REAL &Vph,
                                               CCTK_REAL xL, CCTK_REAL yL, CCTK_REAL zL) {

  //Coordinate transformation from Cartesian to spherical:
  CCTK_REAL rL = sqrt(xL*xL + yL*yL + zL*zL);
  CCTK_REAL th = atan2(sqrt(xL*xL + yL*yL),zL);
  CCTK_REAL ph = atan2(yL,xL);

  CCTK_REAL sinth = sin(th);
  CCTK_REAL costh = cos(th);
  CCTK_REAL sinph = sin(ph);
  CCTK_REAL cosph = cos(ph);

  // E.g., can find transformation matrix here: https://en.wikipedia.org/wiki/List_of_common_coordinate_transformations#From_spherical_coordinates

  CCTK_REAL dx_drL = sinth*cosph;
  CCTK_REAL dy_drL = sinth*sinph;
  CCTK_REAL dz_drL = costh;

  CCTK_REAL dx_dth = rL*costh*cosph;
  CCTK_REAL dy_dth = rL*costh*sinph;
  CCTK_REAL dz_dth =-rL*sinth;

  CCTK_REAL dx_dph =-rL*sinth*sinph;
  CCTK_REAL dy_dph = rL*sinth*cosph;
  CCTK_REAL dz_dph = 0.0;

  //Transform. Note that we must be careful since we define v^i (the vector):
  VxL = VrL*dx_drL + Vth*dx_dth + Vph*dx_dph;
  VyL = VrL*dy_drL + Vth*dy_dth + Vph*dy_dph;
  VzL = VrL*dz_drL + Vth*dz_dth + Vph*dz_dph;

}

void GiRaFFEfood_MagnW_V_i(CCTK_REAL &VxL, CCTK_REAL &VyL, CCTK_REAL &VzL,
                           CCTK_REAL &betarL, CCTK_REAL &betathL, CCTK_REAL &betaphL,
                           CCTK_REAL xL, CCTK_REAL yL, CCTK_REAL zL) {
  DECLARE_CCTK_PARAMETERS;

  //Set the velocities u^i/u^0 to zero:
  CCTK_REAL VrL = -betarL;
  CCTK_REAL Vth = -betathL;
  CCTK_REAL Vph = -betaphL;

  //Coordinate transformation from Cartesian to spherical:
  GiRaFFEfood_MagnW_Sph_to_Cart_XformVector(VxL,VyL,VzL, VrL,Vth,Vph, xL,yL,zL);

}

