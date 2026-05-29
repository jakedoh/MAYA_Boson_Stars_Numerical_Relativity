// Please consult GiRaFFEfood/README for the list of authors and current maintainers //
/*
 * This program implements the initial configuration called Split Monopole,
 * which is a 3D test, derived from the BZ force free monopole solution,
 * by inverting it in the lower hemisphere. This solution is exact but first order accurate.
 * A stable equatorial current sheet is required to maintain this configuration.
 * The vector potential is calculated with formulas (134, 135), arXiv:1310.3274v2.
 * The velocity is calculated with formula (75), arXiv:1310.3274v2, then transformed
 * to the Valencia formalism, arXiv:1501.07276v2, eq.(11)
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

#include "gsl/gsl_sf_dilog.h"
#ifndef M_PI
#define M_PI 3.141592653589793238463
#endif


using namespace std;
//header for the defined functions
void GiRaFFEfood_splitA_mu(CCTK_REAL &xL, CCTK_REAL &yL, CCTK_REAL &zL,
                           CCTK_REAL &AtL, CCTK_REAL &AxL, CCTK_REAL &AyL, CCTK_REAL &AzL);
void GiRaFFEfood_driftV_iB_i(CCTK_REAL &xL, CCTK_REAL &yL, CCTK_REAL &zL,
                             CCTK_REAL &VxL, CCTK_REAL &VyL, CCTK_REAL &VzL,
                             CCTK_REAL &BxL, CCTK_REAL &ByL, CCTK_REAL &BzL,
                             CCTK_REAL &alpL,
                             CCTK_REAL &betarL, CCTK_REAL &betathL, CCTK_REAL &betaphL,
                             CCTK_REAL &grrL, CCTK_REAL &grthL, CCTK_REAL &grphL,
                             CCTK_REAL &gththL, CCTK_REAL &gthphL, CCTK_REAL &gphphL);

void GiRaFFEfood_SplitM_spherical_to_Cartesian_XformVector(CCTK_REAL &VxL, CCTK_REAL &VyL, CCTK_REAL &VzL,
                                                           CCTK_REAL &Vrh, CCTK_REAL &Vth, CCTK_REAL &Vph,
                                                           CCTK_REAL xL, CCTK_REAL yL, CCTK_REAL zL);

void GiRaFFEfood_SplitM_spherical_to_Cartesian_XformOneForm(CCTK_REAL &VxL, CCTK_REAL &VyL, CCTK_REAL &VzL,
                                                            CCTK_REAL &Vrh, CCTK_REAL &Vth, CCTK_REAL &Vph,
                                                            CCTK_REAL xL, CCTK_REAL yL, CCTK_REAL zL);


extern "C" void GiRaFFEfood_SplitMonopole(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT imin = 0, jmin = 0, kmin = 0;
  CCTK_INT imax = cctk_lsh[0], jmax = cctk_lsh[1], kmax = cctk_lsh[2];

#pragma omp parallel for
  for(int k=kmin;k<kmax;k++)
    for(int j=jmin;j<jmax;j++)
      for(int i=imin;i<imax;i++) {
        int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        CCTK_REAL dx = CCTK_DELTA_SPACE(0);
        CCTK_REAL dy = CCTK_DELTA_SPACE(1);
        CCTK_REAL dz = CCTK_DELTA_SPACE(2);

        //Initialize the potential
        CCTK_REAL AtL = 0.0;
        CCTK_REAL AxL = 0.0;
        CCTK_REAL AyL = 0.0;
        CCTK_REAL AzL = 0.0;

        //Initialize the Valencia formulation velocity
        CCTK_REAL VxL = 0.0;
        CCTK_REAL VyL = 0.0;
        CCTK_REAL VzL = 0.0;

        //Cartesian coordinates
        CCTK_REAL xL = x[index];
        CCTK_REAL yL = y[index];
        CCTK_REAL zL = z[index];

        //Avoid divisions by zero (later we divide by sqrt(x^2 + y^2 + z^2):
        CCTK_REAL fudge_factor = 1e-5;
        if(fabs(xL) < fudge_factor*dx) {xL = fudge_factor*dx;}
        if(fabs(yL) < fudge_factor*dy) {yL = fudge_factor*dy;}
        if(fabs(zL) < fudge_factor*dz) {zL = fudge_factor*dz;}

        //The lapse and shift of the background metric -- Shifted Kerr-Schild
        CCTK_REAL alpL   = alp[index];
        CCTK_REAL betaxL = betax[index];
        CCTK_REAL betayL = betay[index];
        CCTK_REAL betazL = betaz[index];
        CCTK_REAL gxxL = gxx[index];
        CCTK_REAL gxyL = gxy[index];
        CCTK_REAL gxzL = gxz[index];
        CCTK_REAL gyyL = gyy[index];
        CCTK_REAL gyzL = gyz[index];
        CCTK_REAL gzzL = gzz[index];

        CCTK_REAL betarL = SKSbetar[index];
        CCTK_REAL betathL = SKSbetath[index];
        CCTK_REAL betaphL = SKSbetaph[index];
        CCTK_REAL grrL = SKSgrr[index];
        CCTK_REAL grthL = SKSgrth[index];
        CCTK_REAL grphL = SKSgrph[index];
        CCTK_REAL gththL = SKSgthth[index];
        CCTK_REAL gthphL = SKSgthph[index];
        CCTK_REAL gphphL = SKSgphph[index];

        //Call the function that computes B-field and IllinoisGRMHD velocity (i.e., velocity found in induction equation).
        CCTK_REAL BxLdumb,ByLdumb,BzLdumb;
        GiRaFFEfood_driftV_iB_i(xL, yL, zL,
                                VxL, VyL, VzL,
                                BxLdumb,ByLdumb,BzLdumb,
                                alpL,
                                betarL, betathL, betaphL,
                                grrL, grthL, grphL,
                                gththL, gthphL, gphphL);

        /** START: CHECK FOR SUPERLUMINAL VELOCITIES **/
#define SQR(x) ((x) * (x))
        //Fill in the interface variables:
        CCTK_REAL sqrtg = sqrt(gxxL * gyyL * gzzL + gxyL * gyzL * gxzL + gxzL * gxyL * gyzL -
                               gxzL * gyyL * gxzL - gxyL * gxyL * gzzL - gxxL * gyzL * gyzL);

        //CCTK_REAL one_minus_one_over_alpha_u0_squared = pow(sqrtg,4.0/6.0)*(gxxL* SQR(VxL + betaxL) +
        CCTK_REAL one_minus_one_over_alpha_u0_squared = (gxxL* SQR(VxL + betaxL) +
                                                         2.0*gxyL*(VxL + betaxL)*(VyL + betayL) +
                                                         2.0*gxzL*(VxL + betaxL)*(VzL + betazL) +
                                                         gyyL* SQR(VyL + betayL) +
                                                         2.0*gyzL*(VyL + betayL)*(VzL + betazL) +
                                                         gzzL* SQR(VzL + betazL) )/(alpL*alpL);
        CCTK_REAL ONE_MINUS_ONE_OVER_GAMMA_SPEED_LIMIT_SQUARED = 1.0-1.0/SQR(GAMMA_SPEED_LIMIT);

        if(one_minus_one_over_alpha_u0_squared > ONE_MINUS_ONE_OVER_GAMMA_SPEED_LIMIT_SQUARED) {
          CCTK_VWarn(CCTK_WARN_ALERT,__LINE__, __FILE__, CCTK_THORNSTRING,
                     "%e %e %e | %e %e %e RUH ROH. EXCEEDED SPEED LIMIT! %e > %e\n",xL,yL,zL,  VxL,VyL,VzL, one_minus_one_over_alpha_u0_squared, ONE_MINUS_ONE_OVER_GAMMA_SPEED_LIMIT_SQUARED);
        }
        /** END: CHECK FOR SUPERLUMINAL VELOCITIES **/


        //Calculates the Valencia formulation velocities.
        CCTK_REAL ETVxL = (VxL + betaxL)/alpL;
        CCTK_REAL ETVyL = (VyL + betayL)/alpL;
        CCTK_REAL ETVzL = (VzL + betazL)/alpL;

        //Fill in the Valencia formulation velocities.
        vel[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = ETVxL;
        vel[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = ETVyL;
        vel[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] = ETVzL;

        // Next compute the 4-vector potential in staggered coordinates.
        CCTK_REAL xh = xL + 0.5*dx;
        CCTK_REAL yh = yL + 0.5*dy;
        CCTK_REAL zh = zL + 0.5*dz;

        //Call the function that computes the vector potentials in staggered coordinates.
        CCTK_REAL dummy0,dummy1,dummy2,dummy3;
        //At(i+0.5,j+0.5,k+0.5)
        GiRaFFEfood_splitA_mu(xh, yh, zh, AtL, dummy1, dummy2, dummy3);
        //Ax(i,j+0.5,k+0.5)
        GiRaFFEfood_splitA_mu(xL, yh, zh, dummy0, AxL, dummy2, dummy3);
        //Ay(i+0.5,j,k+0.5)
        GiRaFFEfood_splitA_mu(xh, yL, zh, dummy0, dummy1, AyL, dummy3);
        //Az(i+0.5,j+0.5,k)
        GiRaFFEfood_splitA_mu(xh, yh, zL, dummy0, dummy1, dummy2, AzL);

        Aphi[index] = alpL*AtL*sqrtg;
        Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = AxL;
        Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = AyL;
        Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] = AzL;

        psi6phi[index] = alpL*AtL*sqrtg;
        Ax[index] = AxL;
        Ay[index] = AyL;
        Az[index] = AzL;
      }
}

void GiRaFFEfood_splitA_mu(CCTK_REAL &xL, CCTK_REAL &yL, CCTK_REAL &zL,
                           CCTK_REAL &AtL, CCTK_REAL &AxL, CCTK_REAL &AyL, CCTK_REAL &AzL) {
  DECLARE_CCTK_PARAMETERS;

  //Transformation from Cartesian to spherical coordinates.
  CCTK_REAL r_KS = sqrt(xL*xL + yL*yL + zL*zL) + KerrSchild_radial_shift;
  CCTK_REAL th = atan2(sqrt(xL*xL + yL*yL),zL);

  CCTK_REAL sinth  = sin(th);
  CCTK_REAL costh  = cos(th);

  //Parameters
  CCTK_REAL mL = BH_mass;
  CCTK_REAL roM = r_KS/mL;
  CCTK_REAL twoMor = 2.0/roM;
  CCTK_REAL Li2x = gsl_sf_dilog(twoMor);
  CCTK_REAL Lx = 0.0;
  CCTK_REAL fr;
  CCTK_REAL cL = split_C;
  CCTK_REAL aL = BH_spin;

  //The dilogarithmic function:
  //$L(x)=Li_2(x)+\frac{1}{2}\ln x\ln(1-x)$ for $0 < x < 1$
  //where $x=\frac{2M}{r}$.
  if (( twoMor > 0.0) and (twoMor < 1.0)) {
    Lx = Li2x + 0.5*log(twoMor)*log(1.0 - twoMor); }
  //    Lx = (Li2x - log(1.0/twoMor*log(1.0-twoMor); }

  //The radial function:
  //$f(r)=\frac{r^2(2r-3M)}{8M^3}L(\frac{2M}{r}) + \frac{M^2+3Mr-6r^2}{12M^2}\ln\frac{r}{2M}
  //      +\frac{11}{72}+\frac{M}{3r}+\frac{r}{2M}-\frac{r^2}{2M^2}$
  //fr is replaced with equation(41) from astro-ph/0404512

  //switch for 1st order solution
  if(drop_fr_SplitM) {
    fr = 0.0;
  } else {

    fr = roM*roM*(2.0*roM - 3.0)/8.0 * (Lx)
      + (1.0 + 3.0*roM - 6.0*roM*roM)/12.0*log(1.0/twoMor)
      + 11.0/72.0 + 1.0/(3.0*roM) + 1.0/twoMor - roM*roM/2.0;
  }

  /** HERE WE USE EQS 134 AND 135 FROM https://arxiv.org/pdf/1310.3274.pdf . IGNORE THE f(r) AND f'(r) TERMS; THEY MUST BE SET TO ZERO ANYWAY **/
  //The 4-vector potential in spherical coordinates:
  //Calculate A_r.
  //$A_r=-\frac{C a}{8}|\cos\theta|(1+\frac{4M}{r})\sqrt{1+\frac{2M}{r}}
  CCTK_REAL Arh = - aL*cL/8.0*fabs(costh)*(1.0+2.0*twoMor)*sqrt(1.0+twoMor);
  //$A_\theta=0$
  CCTK_REAL Ath = 0.0;

  //Calculate A_{\phi}.
  //$A_{\phi}=M^2 C[1-|\cos\theta|+a^2 f(r)\cos\theta\sin^2\theta]+O(a_{*}^4$
  //CCTK_REAL Aph = mL*mL*cL*(1.0 - costh + aL*aL*fr*costh*sin(th)*sin(th));
  //Aph is replaced with equation(35) from astro-ph/0404512
  CCTK_REAL Aph0 =  cL - cL*fabs(costh);
  CCTK_REAL Aph2 =   cL*fr*costh*sinth*sinth;
  CCTK_REAL Aph  =   mL*mL* (Aph0 + aL*aL*Aph2);

  GiRaFFEfood_SplitM_spherical_to_Cartesian_XformOneForm(AxL,AyL,AzL, Arh,Ath,Aph,xL,yL,zL);

  /* Finally, set the t component: A_t */
  AtL = -aL/(8.0*mL*mL)*Aph;
}

void GiRaFFEfood_SplitM_spherical_to_Cartesian_XformVector(CCTK_REAL &VxL, CCTK_REAL &VyL, CCTK_REAL &VzL,
                                                           CCTK_REAL &Vrh, CCTK_REAL &Vth, CCTK_REAL &Vph,
                                                           CCTK_REAL xL, CCTK_REAL yL, CCTK_REAL zL) {

  //Coordinate transformation from Cartesian to spherical:
  CCTK_REAL rh = sqrt(xL*xL + yL*yL + zL*zL);
  CCTK_REAL th = atan2(sqrt(xL*xL + yL*yL),zL);
  CCTK_REAL ph = atan2(yL,xL);

  CCTK_REAL sinth  = sin(th);
  CCTK_REAL costh  = cos(th);
  CCTK_REAL sinph  = sin(ph);
  CCTK_REAL cosph  = cos(ph);

  //Transform. Note that we must be careful since we define v^i (the vector):
  /*
    v^x = v^r dx/dr + v^ph dx/dph + v^th dx/dth
    v^y = v^r dy/dr + v^ph dy/dph + v^th dy/dth
    v^z = v^r dz/dr + v^ph dz/dph + v^th dz/dth
  */

  // E.g., can find transformation matrix here: https://en.wikipedia.org/w/index.php?title=List_of_common_coordinate_transformations&oldid=845649164#From_spherical_coordinates

  CCTK_REAL dx_dr = sinth*cosph;
  CCTK_REAL dy_dr = sinth*sinph;
  CCTK_REAL dz_dr = costh;

  CCTK_REAL dx_dth = rh*costh*cosph;
  CCTK_REAL dy_dth = rh*costh*sinph;
  CCTK_REAL dz_dth =-rh*sinth;

  CCTK_REAL dx_dph =-rh*sinth*sinph;
  CCTK_REAL dy_dph = rh*sinth*cosph;
  CCTK_REAL dz_dph = 0.0;

  //Transform.
  VxL = Vrh*dx_dr + Vth*dx_dth + Vph*dx_dph;
  VyL = Vrh*dy_dr + Vth*dy_dth + Vph*dy_dph;
  VzL = Vrh*dz_dr + Vth*dz_dth + Vph*dz_dph;
}

void GiRaFFEfood_SplitM_spherical_to_Cartesian_XformOneForm(CCTK_REAL &VxL, CCTK_REAL &VyL, CCTK_REAL &VzL,
                                                            CCTK_REAL &Vrh, CCTK_REAL &Vth, CCTK_REAL &Vph,
                                                            CCTK_REAL xL, CCTK_REAL yL, CCTK_REAL zL) {

  //Coordinate transformation from Cartesian to spherical:
  CCTK_REAL rho = sqrt(xL*xL + yL*yL + zL*zL);

  //Transform. Note that we must be careful since we define v^i (the vector):
  /*
    v_x = v_r dr/dx + v_ph dph/dx + v_th dth/dx
    v_y = v_r dr/dy + v_ph dph/dy + v_th dth/dy
    v_z = v_r dr/dz + v_ph dph/dz + v_th dth/dz
  */

  // E.g., can find transformation matrix here: https://en.wikipedia.org/w/index.php?title=List_of_common_coordinate_transformations&oldid=845649164#From_Cartesian_coordinates_2

  CCTK_REAL dr_dx = xL/rho;
  CCTK_REAL dr_dy = yL/rho;
  CCTK_REAL dr_dz = zL/rho;

  CCTK_REAL dth_dx = xL*zL / (rho*rho*sqrt(xL*xL + yL*yL));
  CCTK_REAL dth_dy = yL*zL / (rho*rho*sqrt(xL*xL + yL*yL));
  CCTK_REAL dth_dz = - sqrt(xL*xL + yL*yL) / (rho*rho);

  CCTK_REAL dph_dx = -yL/(xL*xL + yL*yL);
  CCTK_REAL dph_dy =  xL/(xL*xL + yL*yL);
  CCTK_REAL dph_dz = 0.0;

  //Transform.
  VxL = Vrh*dr_dx + Vth*dth_dx + Vph*dph_dx;
  VyL = Vrh*dr_dy + Vth*dth_dy + Vph*dph_dy;
  VzL = Vrh*dr_dz + Vth*dth_dz + Vph*dph_dz;
}

void GiRaFFEfood_driftV_iB_i(CCTK_REAL &xL, CCTK_REAL &yL, CCTK_REAL &zL,
                             CCTK_REAL &VxL, CCTK_REAL &VyL, CCTK_REAL &VzL,
                             CCTK_REAL &BxL, CCTK_REAL &ByL, CCTK_REAL &BzL,
                             CCTK_REAL &alpL,
                             CCTK_REAL &betarL, CCTK_REAL &betathL, CCTK_REAL &betaphL,
                             CCTK_REAL &grrL, CCTK_REAL &grthL, CCTK_REAL &grphL,
                             CCTK_REAL &gththL, CCTK_REAL &gthphL, CCTK_REAL &gphphL) {
  DECLARE_CCTK_PARAMETERS;

  //Coordinate transformation from Cartesian to spherical:
  CCTK_REAL r_KS = sqrt(xL*xL + yL*yL + zL*zL) + KerrSchild_radial_shift;
  CCTK_REAL th = atan2(sqrt(xL*xL + yL*yL),zL);

  CCTK_REAL sinth  = sin(th);
  CCTK_REAL costh  = cos(th);
  CCTK_REAL cos2th = cos(2.0*th);

  //Parameters
  CCTK_REAL mL = BH_mass;
  CCTK_REAL roM = r_KS/mL;
  CCTK_REAL twoMor = 2.0/roM;
  CCTK_REAL Li2x = gsl_sf_dilog(twoMor);
  CCTK_REAL Li1x = -log(1.0 - twoMor)/twoMor;
  CCTK_REAL Lx = 0.0, Lxp = 0.0;
  CCTK_REAL fr, fpr;
  CCTK_REAL cL = split_C;
  CCTK_REAL aL = BH_spin;

  CCTK_REAL sqgm = sqrt(grrL * gththL * gphphL + grthL * gthphL * grphL + grphL * grthL * gthphL -
                        grphL * gththL * grphL - grthL * grthL * gphphL - grrL * gthphL * gthphL);

  //The dilogarithmic function:
  //$L(x)=Li_2(x)+\frac{1}{2}\ln x\ln(1-x)$ for $0 < x < 1$
  //where $x=\frac{2M}{r}$.
  if (( twoMor > 0.0) and (twoMor < 1.0)) {
    //    Lx = Li2x + 0.5*log(twoMor)*log(1.0 - twoMor); }
    Lx = (Li2x - log(1.0-twoMor)*log(1.0/twoMor));
    Lxp = -1.0/r_KS*(Li1x*twoMor + log(1.0-twoMor) + log(1.0/twoMor)*twoMor/(1.0-twoMor));
  }

  //The radial function:
  //$f(r)=\frac{r^2(2r-3M)}{8M^3}L(\frac{2M}{r}) + \frac{M^2+3Mr-6r^2}{12M^2}\ln\frac{r}{2M}
  //      +\frac{11}{72}+\frac{M}{3r}+\frac{r}{2M}-\frac{r^2}{2M^2}$
  //fr is replaced with equation(41) from astro-ph/0404512

  //switch for 1st order solution
  if(drop_fr_SplitM) {
    fr = 0.0, fpr = 0.0;
  } else {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING ,"drop_fr_SplitM SHOULD BE SET TO ONE. SHOULDN'T BE HERE.\n");

    fr = roM*roM*(2.0*roM - 3.0)/8.0 * (Lx)
      + (1.0 + 3.0*roM - 6.0*roM*roM)/12.0*log(1.0/twoMor)
      + 11.0/72.0 + 1.0/(3.0*roM) + 1.0/twoMor - roM*roM/2.0;
    //The derivarive of the radial function
    fpr = (roM*(2.0*roM - 3.0) + roM*roM ) * (Lx)/(4.0*mL)
      + roM*roM*(2.0*roM - 3.0)/8.0 * (Lxp)
      + (3.0 - 12.0*roM)/(12.0*mL)*log(1.0/twoMor)
      + (1.0 + 3.0*roM - 6.0*roM*roM)/(12.0*r_KS)
      - 1.0/mL*(1.0/(3.0*roM*roM) + 1.0/2.0 - 1.0*roM);
  }

  //The split monopole condition C<0 for z<0:
  if(zL<0.0) cL*=-1.0;

  // Eq. (128) arXiv:1310.3274v2 was used to calculate the magnetic field
  //Calculate B^r: eq.128
  CCTK_REAL Burh = cL*alpL/(roM*roM)*(1.0 + aL*aL/(2.0*r_KS*r_KS)
                                      *(-2.0*costh + roM*roM*(1.0 + 3.0*cos2th)*fr));

  //Calculate B^\theta$
  CCTK_REAL Buth = -cL*alpL*aL*aL/(r_KS*r_KS)*sinth*costh*fpr;

  //Calculate B^\phi.
  CCTK_REAL Buph = - cL*alpL*aL/(8.0*roM*r_KS)*(1.0 + 2.0*twoMor);

  CCTK_REAL Bdrh = grrL*Burh  + grthL*Buth  + grphL*Buph;
  CCTK_REAL Bdth = grthL*Burh + gththL*Buth + gthphL*Buph;
  CCTK_REAL Bdph = grphL*Burh + gthphL*Buth + gphphL*Buph;

  // Eq. (131, 132, 133) arXiv:1310.3274v2 were used to calculate the electric field
  //Calculate E_r.
  CCTK_REAL Edrh = -cL*aL*aL*aL/(8.0*alpL*mL*mL*mL)*costh*sinth*sinth*fpr;

  //Calculate E_\theta$
  CCTK_REAL Edth = -cL*aL/(8.0*alpL)*(sinth + aL*aL*fr*sinth
                                      *(2.0*costh*costh - sinth*sinth))
    - betarL*sqgm*aL*cL/(8.0*r_KS*r_KS)*(1.0 + 2.0*twoMor);

  //Calculate E_\phi.
  CCTK_REAL Edph = betarL/(alpL*mL)*cL*aL*aL*costh*sinth*sinth*fpr;

  /* START: CHECK B2 < E2 CONDITION. NOTE THAT IN UNSHIFTED KERR-SCHILD, THIS IS VIOLATED DEEP IN THE HORIZON. */
  CCTK_REAL guprr = (-2*r_KS)/(r_KS*(2 + r_KS) + pow(aL,2)*pow(cos(th),2)) + 1/(1 - (pow(aL,2)*pow(sin(th),2))/(pow(aL,2) + pow(r_KS,2)));
  CCTK_REAL guprth = 0.0;
  CCTK_REAL guprph = (2*aL)/(pow(aL,2) + 2*pow(r_KS,2) + pow(aL,2)*cos(2*th));
  CCTK_REAL gupthth = 1/(pow(r_KS,2) + pow(aL,2)*pow(cos(th),2));
  CCTK_REAL gupthph = 0.0;
  CCTK_REAL gupphph = pow(sin(th),-2)/(pow(aL,2) + pow(r_KS,2) - pow(aL,2)*pow(sin(th),2));


  CCTK_REAL Eurh = guprr*Edrh + guprth*Edth + guprph*Edph;
  CCTK_REAL Euth = guprth*Edrh + gupthth*Edth + gupthph*Edph;
  CCTK_REAL Euph = guprph*Edrh + gupthph*Edth + gupphph*Edph;

  CCTK_REAL B2 = Bdrh*Burh + Bdth*Buth + Bdph*Buph;
  CCTK_REAL E2 = Edrh*Eurh + Edth*Euth + Edph*Euph;

  if(B2 < E2) {
    CCTK_REAL ph = atan2(yL,xL);
    CCTK_VWarn(CCTK_WARN_ALERT,__LINE__, __FILE__, CCTK_THORNSTRING,
               "WARNING: FOUND B2 >= E2: %e %e %e | %e < %e | B2-E2 = %e . DID YOU SET KerrSchild_radial_shift=1 ???\n",r_KS,th,ph,B2,E2,B2-E2);
  }
  /* END: CHECK B2 < E2 CONDITION. NOTE THAT IN UNSHIFTED KERR-SCHILD, THIS IS VIOLATED DEEP IN THE HORIZON. */

  // Only holds in the Newtonian (flat space) limit.
  //if(sqrt(B2/(B2-E2)) > 2000.0) printf("BADGAMMA: %e %e %e | %e < %e | B2-E2 = %e\n",r_KS,th,ph,B2,E2,B2-E2);


  CCTK_REAL Vurh  = alpL * (Edth*Bdph - Edph*Bdth)/B2/sqgm - betarL;
  CCTK_REAL Vuth  = alpL * (Edph*Bdrh - Edrh*Bdph)/B2/sqgm - betathL;
  CCTK_REAL Vuph  = alpL * (Edrh*Bdth - Edth*Bdrh)/B2/sqgm - betaphL;

  /* SET VELOCITY, B-FIELDS */
  GiRaFFEfood_SplitM_spherical_to_Cartesian_XformVector(VxL,VyL,VzL, Vurh,Vuth,Vuph,xL,yL,zL);
  GiRaFFEfood_SplitM_spherical_to_Cartesian_XformVector(BxL,ByL,BzL, Burh,Buth,Buph,xL,yL,zL);
  CCTK_REAL ExL,EyL,EzL;
  GiRaFFEfood_SplitM_spherical_to_Cartesian_XformVector(ExL,EyL,EzL, Eurh,Euth,Euph,xL,yL,zL);


  /* START: CHECK B2 < E2 CONDITION AFTER THE COORDINATE TRANSFORMATION FROM SPHERICAL TO CARTESIAN. NOTE THAT IN UNSHIFTED KERR-SCHILD, THIS IS VIOLATED DEEP IN THE HORIZON. */
  CCTK_REAL B_xL,B_yL,B_zL;
  CCTK_REAL E_xL,E_yL,E_zL;
  GiRaFFEfood_SplitM_spherical_to_Cartesian_XformOneForm(B_xL,B_yL,B_zL, Bdrh,Bdth,Bdph,xL,yL,zL);
  GiRaFFEfood_SplitM_spherical_to_Cartesian_XformOneForm(E_xL,E_yL,E_zL, Edrh,Edth,Edph,xL,yL,zL);

  CCTK_REAL B2C = BxL*B_xL + ByL*B_yL + BzL*B_zL;
  CCTK_REAL E2C = ExL*E_xL + EyL*E_yL + EzL*E_zL;

  /* NOW CHECK B^2 < E^2 AFTER THE COORDINATE TRANSFORMATION FROM SPHERICAL TO CARTESIAN */
  if(B2C < E2C) {
    CCTK_REAL ph = atan2(yL,xL);
    CCTK_VWarn(CCTK_WARN_ALERT,__LINE__, __FILE__, CCTK_THORNSTRING,
               "WARNING: FOUND B2 >= E2: %e %e %e | %e < %e | B2-E2 = %e . DID YOU SET KerrSchild_radial_shift=1 ???\n",r_KS,th,ph,B2C,E2C,B2C-E2C);
  }
  /* END: CHECK B2 < E2 CONDITION AFTER THE COORDINATE TRANSFORMATION FROM SPHERICAL TO CARTESIAN. NOTE THAT IN UNSHIFTED KERR-SCHILD, THIS IS VIOLATED DEEP IN THE HORIZON. */
}

/*
 * This program calculates the error between the exact and evolved variables - Bi, Vi and Ai - for the Split Monopole.
 * The magnetic field is calculated with formulas (140-145), arXiv:1310.3274v2
 * */
extern "C" void GiRaFFEfood_ErrorSplitM(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if( !(Compute_Exact_Every > 0 && (cctk_iteration%Compute_Exact_Every) == 0) ) return;

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

        //Cartesian Coordinates
        CCTK_REAL xL = x[index];
        CCTK_REAL yL = y[index];
        CCTK_REAL zL = z[index];

        CCTK_REAL dx = CCTK_DELTA_SPACE(0);
        CCTK_REAL dy = CCTK_DELTA_SPACE(1);
        CCTK_REAL dz = CCTK_DELTA_SPACE(2);

        //Staggered coordinates
        CCTK_REAL xh = xL + 0.5*dx;
        CCTK_REAL yh = yL + 0.5*dy;
        CCTK_REAL zh = zL + 0.5*dz;

        //Calculate the magnetic and velocity fields
        CCTK_REAL alpL   = alp[index];

        CCTK_REAL betarL = SKSbetar[index];
        CCTK_REAL betathL = SKSbetath[index];
        CCTK_REAL betaphL = SKSbetaph[index];
        CCTK_REAL grrL = SKSgrr[index];
        CCTK_REAL grthL = SKSgrth[index];
        CCTK_REAL grphL = SKSgrph[index];
        CCTK_REAL gththL = SKSgthth[index];
        CCTK_REAL gthphL = SKSgthph[index];
        CCTK_REAL gphphL = SKSgphph[index];

        CCTK_REAL VxL,VyL,VzL;
        CCTK_REAL BxL, ByL, BzL;
        // Compute B-field and IllinoisGRMHD velocity (i.e., velocity found in induction equation).
        GiRaFFEfood_driftV_iB_i(xL, yL, zL,
                                VxL, VyL, VzL,
                                BxL,ByL,BzL,
                                alpL,
                                betarL, betathL, betaphL,
                                grrL, grthL, grphL,
                                gththL, gthphL, gphphL);

        exactBx[index] = BxL;
        exactBy[index] = ByL;
        exactBz[index] = BzL;

        exactVx[index] = VxL;
        exactVy[index] = VyL;
        exactVz[index] = VzL;

        //Calculate the potentials
        CCTK_REAL AtE, AxE, AyE, AzE;
        CCTK_REAL dummy0,dummy1,dummy2,dummy3;
        //At(i+0.5,j+0.5,k+0.5)
        GiRaFFEfood_splitA_mu(xh, yh, zh, AtE, dummy1, dummy2, dummy3);
        //Ax(i,j+0.5,k+0.5)
        GiRaFFEfood_splitA_mu(xL, yh, zh, dummy0, AxE, dummy2, dummy3);
        //Ay(i+0.5,j,k+0.5)
        GiRaFFEfood_splitA_mu(xh, yL, zh, dummy0, dummy1, AyE, dummy3);
        //Az(i+0.5,j+0.5,k)
        GiRaFFEfood_splitA_mu(xh, yh, zL, dummy0, dummy1, dummy2, AzE);

        //Fill in the interface variables:
        delBx[index] = Bx[index] - BxL;
        delBy[index] = By[index] - ByL;
        delBz[index] = Bz[index] - BzL;

        delvx[index] = vx[index] - exactVx[index];
        delvy[index] = vy[index] - exactVy[index];
        delvz[index] = vz[index] - exactVz[index];

        delAx[index] = Ax[index] - AxE;
        delAy[index] = Ay[index] - AyE;
        delAz[index] = Az[index] - AzE;

        errorsumx += (delBx[index])*(delBy[index]);//Bx[index] - BxL);
        errorsumy += (delBy[index])*(delBy[index]);//By[index] - ByL);
        errorsumz += (delBz[index])*(delBz[index]);//Bz[index] - BzL);
      }

  /* Useful diagnostic:
     CCTK_REAL dx = CCTK_DELTA_SPACE(0);
     CCTK_REAL dy = CCTK_DELTA_SPACE(1);
     CCTK_REAL dz = CCTK_DELTA_SPACE(2);

     CCTK_REAL norm2Bx = sqrt(errorsumx*dx*dy*dz);
     CCTK_REAL norm2By = sqrt(errorsumy*dx*dy*dz);
     CCTK_REAL norm2Bz = sqrt(errorsumz*dx*dy*dz);
  */
  CCTK_VInfo(CCTK_THORNSTRING,
             "Lev:\t %d\t norm2BxByBz: %e %e %e | %e\n",GetRefinementLevel(cctkGH),errorsumx,errorsumy,errorsumz,delBx[CCTK_GFINDEX3D(cctkGH,3,3,3)]);
}

