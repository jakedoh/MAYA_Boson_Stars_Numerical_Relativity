// Please consult GiRaFFEfood/README for the list of authors and current maintainers //
/*
 * This program implements the initial configuration called Exact Wald,
 * an exact electrovacuum solution to Maxwell's equations, for a black hole
 * immersed in uniform magnetic field aligned with the axis of rotation.
 * In the non-spinning case this is also a force-free solution, but only outside the hole.
 * The vector potential is calculated with equation (138) arXiv:1310.3274v2.
 * The magnetic field is set with the formulas (140, 141, 142).
 * The components of the velocity are given by the equations (153, 154, 155),
 * arXiv:1310.3274v2 then transformed to the Valencia formalism, arXiv:1501.07276v2, eq.(11)
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

// header for the defined functions
void GiRaFFEfood_ExactWald_exactV_i(CCTK_REAL &VxL, CCTK_REAL &VyL, CCTK_REAL &VzL,
                                    CCTK_REAL xL, CCTK_REAL yL, CCTK_REAL zL);

void GiRaFFEfood_ExactWald_exactA_i(CCTK_REAL &AxL, CCTK_REAL &AyL, CCTK_REAL &AzL,
                                    CCTK_REAL xL, CCTK_REAL yL, CCTK_REAL zL);

void GiRaFFEfood_ExactWald_exactB_i(CCTK_REAL &BxL, CCTK_REAL &ByL, CCTK_REAL &BzL,
                                    CCTK_REAL xL, CCTK_REAL yL, CCTK_REAL zL);

void GiRaFFEfood_ExactWald_spherical_to_Cartesian_XformVector(CCTK_REAL &VxL, CCTK_REAL &VyL, CCTK_REAL &VzL,
                                                              CCTK_REAL &Vrh, CCTK_REAL &Vth, CCTK_REAL &Vph,
                                                              CCTK_REAL xL, CCTK_REAL yL, CCTK_REAL zL);

extern "C" void GiRaFFEfood_ExactWald(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT imin = 0, jmin = 0, kmin = 0;
  CCTK_INT imax = cctk_lsh[0], jmax = cctk_lsh[1], kmax = cctk_lsh[2];

  CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  CCTK_REAL dz = CCTK_DELTA_SPACE(2);

#pragma omp parallel for
  for(int k=kmin;k<kmax;k++)
    for(int j=jmin;j<jmax;j++)
      for(int i=imin;i<imax;i++) {
        int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        //Initialize the potential
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

        //Avoid divisions by zero (later we divide by sqrt(x^2 + y^2 + z^2):
        CCTK_REAL fudge_factor = 1e-15;
        if(fabs(xL) < fudge_factor*dx) {xL = fudge_factor*dx;}
        if(fabs(yL) < fudge_factor*dy) {yL = fudge_factor*dy;}
        if(fabs(zL) < fudge_factor*dz) {zL = fudge_factor*dz;}

        //The vector potential is stagered so we need staggered coordinates:
        CCTK_REAL xh = xL + 0.5*dx;
        CCTK_REAL yh = yL + 0.5*dy;
        CCTK_REAL zh = zL + 0.5*dz;

        //Call the function that calculates the vector potential. The vector potential is staggered.
        CCTK_REAL dummy0,dummy1,dummy2;
        //Ax(i,j+0.5,k+0.5)
        GiRaFFEfood_ExactWald_exactA_i(AxL,dummy1,dummy2, xL,yh,zh);
        //Ay(i+0.5,j,k+0.5)
        GiRaFFEfood_ExactWald_exactA_i(dummy0,AyL,dummy2, xh,yL,zh);
        //Az(i+0.5,j+0.5,k)
        GiRaFFEfood_ExactWald_exactA_i(dummy0,dummy1,AzL, xh,yh,zL);

        //Fill in the interface variables:
        Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = AxL;
        Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = AyL;
        Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] = AzL;

        Ax[index] = AxL;
        Ay[index] = AyL;
        Az[index] = AzL;

        //Call the function that calculates the velocity. Note that the velocity is not staggered.
        GiRaFFEfood_ExactWald_exactV_i(VxL,VyL,VzL,xL,yL,zL);

        CCTK_REAL alpL = alp[index];
        CCTK_REAL betaxL = betax[index];
        CCTK_REAL betayL = betay[index];
        CCTK_REAL betazL = betaz[index];
        /*
        // Confirm that beta^i A_i = 0:
        //Ax(i,j+0.5,k+0.5)
        GiRaFFEfood_ExactWald_exactA_i(AxL,dummy1,dummy2, xL,yL,zL);
        //Ay(i+0.5,j,k+0.5)
        GiRaFFEfood_ExactWald_exactA_i(dummy0,AyL,dummy2, xL,yL,zL);
        //Az(i+0.5,j+0.5,k)
        GiRaFFEfood_ExactWald_exactA_i(dummy0,dummy1,AzL, xL,yL,zL);

        CCTK_VInfo(CCTK_THORNSTRING,"iii %e | %e ?= %e %e %e ?= %e | %e %e %e\n",AxL*betaxL + AyL*betayL + AzL*betazL,2*sqrt(1-zL*zL/(rh*rh))/((2+rh)*sqrt(1+yL*yL/(xL*xL))),betaxL,betayL,betazL ,2*zL/(2*rh+rh*rh) , AxL,AyL,AzL);
        */

        // GiRaFFEfood_ExactWald_exactV_i will yield velocities consistent with those found in induction equation.
        //  HydroBase's vel[] gridfunctions require VALENCIA FORMULATION version of velocities.
        // Here we convert velocities to be consistent with Valencia formulation: eq.(11) arXiv:1501.07276
        CCTK_REAL ETVx = (VxL + betaxL)/alpL;
        CCTK_REAL ETVy = (VyL + betayL)/alpL;
        CCTK_REAL ETVz = (VzL + betazL)/alpL;

        //Fill in the interface variables:
        vel[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = ETVx;
        vel[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = ETVy;
        vel[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] = ETVz;
      }
}


//Convergence routine that implements the difference between the numerical and exact
//vector potential and magnetic field for the exact Wald solution, outside the hole.
extern "C" void GiRaFFEfood_ConvergenceWald(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if( !(Compute_Exact_Every > 0 && (cctk_iteration%Compute_Exact_Every) == 0) ) return;

  CCTK_INT imin = 0, jmin = 0, kmin = 0;
  CCTK_INT imax = cctk_lsh[0], jmax = cctk_lsh[1], kmax = cctk_lsh[2];

  CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  CCTK_REAL dz = CCTK_DELTA_SPACE(2);

  CCTK_REAL errorsumt = 0.0;
  CCTK_REAL errorsumx = 0.0;
  CCTK_REAL errorsumy = 0.0;
  CCTK_REAL errorsumz = 0.0;

#pragma omp parallel for reduction(+:errorsumt,errorsumx,errorsumy,errorsumz)
  for(int k=kmin;k<kmax;k++)
    for(int j=jmin;j<jmax;j++)
      for(int i=imin;i<imax;i++) {
        int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        CCTK_REAL rL = r[index];

        //Cartesian coordinates
        CCTK_REAL xL = x[index];
        CCTK_REAL yL = y[index];
        CCTK_REAL zL = z[index];

        //Avoid divisions by zero (later we divide by sqrt(x^2 + y^2 + z^2):
        CCTK_REAL fudge_factor = 1e-15;
        if(fabs(xL) < fudge_factor*dx) {xL = fudge_factor*dx;}
        if(fabs(yL) < fudge_factor*dy) {yL = fudge_factor*dy;}
        if(fabs(zL) < fudge_factor*dz) {zL = fudge_factor*dz;}

        //Staggered coordinates
        CCTK_REAL xh = xL + 0.5*dx;
        CCTK_REAL yh = yL + 0.5*dy;
        CCTK_REAL zh = zL + 0.5*dz;

        //Call the function that calculates the vector potential. The vector potential is staggered.
        CCTK_REAL AxE,  AyE, AzE;
        CCTK_REAL dummy0,dummy1,dummy2;
        //Ax(i,j+0.5,k+0.5)
        GiRaFFEfood_ExactWald_exactA_i(AxE,dummy1,dummy2, xL,yh,zh);
        //Ay(i+0.5,j,k+0.5)
        GiRaFFEfood_ExactWald_exactA_i(dummy0,AyE,dummy2, xh,yL,zh);
        //Az(i+0.5,j+0.5,k)
        GiRaFFEfood_ExactWald_exactA_i(dummy0,dummy1,AzE, xh,yh,zL);

        //Fill in the interface variables:
        delpsi6phi[index] = psi6phi[index];
        delAx[index] = Ax[index] - AxE;
        delAy[index] = Ay[index] - AyE;
        delAz[index] = Az[index] - AzE;

        //The magnetic field
        CCTK_REAL BxE,  ByE, BzE;
        GiRaFFEfood_ExactWald_exactB_i(BxE,ByE,BzE, xL,yL,zL);
        exactBx[index] = BxE;
        exactBy[index] = ByE;
        exactBz[index] = BzE;

        delBx[index] = Bx[index] - BxE;
        delBy[index] = By[index] - ByE;
        delBz[index] = Bz[index] - BzE;

        //Calculating errors for velocities
        CCTK_REAL VxL,VyL,VzL;
        GiRaFFEfood_ExactWald_exactV_i(VxL,VyL,VzL, xL,yL,zL);
        exactVx[index] = VxL;
        exactVy[index] = VyL;
        exactVz[index] = VzL;

        delvx[index] = vx[index] - VxL;
        delvy[index] = vy[index] - VyL;
        delvz[index] = vz[index] - VzL;

        //Here is applied the excision of the domain that is too close to the black hole.
        /*
          Our goal is to integrate these "del" quantities, but note that FFE
          B^2 > E^2 condition is *violated* inside horizon.
          So to make the integral meaningful, we excise the BH from integration domain:
        */
        if(rL < 8.0*(2.0-KerrSchild_radial_shift) || rL > 90.0) {
          delpsi6phi[index] = 0.0;
          delAx[index] = 0.0;
          delAy[index] = 0.0;
          delAz[index] = 0.0;

          delBx[index] = 0.0;
          delBy[index] = 0.0;
          delBz[index] = 0.0;

          delvx[index] = 0.0;
          delvy[index] = 0.0;
          delvz[index] = 0.0;
        }

        if(rL > 8.0*(2.0-KerrSchild_radial_shift) && rL < 90.0) {
          errorsumt += (delpsi6phi[index]*delpsi6phi[index]);
          errorsumx += (delAx[index])*(delAx[index]);
          errorsumy += (delAy[index])*(delAy[index]);
          errorsumz += (delAz[index])*(delAz[index]);
        }

      }

  CCTK_REAL norm2At = sqrt(errorsumt);
  CCTK_REAL norm2Ax = sqrt(errorsumx);
  CCTK_REAL norm2Ay = sqrt(errorsumy);
  CCTK_REAL norm2Az = sqrt(errorsumz);
  CCTK_VInfo(CCTK_THORNSTRING,"norm2AtAxAyAz: %e %e %e %e | %e\n",norm2At, norm2Ax, norm2Ay, norm2Az, delAx[CCTK_GFINDEX3D(cctkGH,3,3,3)]);
}


//Transform vectors from spherical to Cartesian coordinates
void GiRaFFEfood_ExactWald_spherical_to_Cartesian_XformVector(CCTK_REAL &VxL, CCTK_REAL &VyL, CCTK_REAL &VzL,
                                                              CCTK_REAL &Vrh, CCTK_REAL &Vth, CCTK_REAL &Vph,
                                                              CCTK_REAL xL, CCTK_REAL yL, CCTK_REAL zL) {

  //Coordinate transformation from Cartesian to spherical:
  CCTK_REAL rh = sqrt(xL*xL + yL*yL + zL*zL);
  CCTK_REAL th = atan2(sqrt(xL*xL + yL*yL),zL);
  CCTK_REAL ph = atan2(yL,xL);

  CCTK_REAL sinth = sin(th);
  CCTK_REAL costh = cos(th);
  CCTK_REAL sinph = sin(ph);
  CCTK_REAL cosph = cos(ph);

  //Transform. Note that we must be careful since we define v^i (the vector):
  /*
    v^x = v^r dx/dr + v^ph dx/dph + v^th dx/dth
    v^y = v^r dy/dr + v^ph dy/dph + v^th dy/dth
    v^z = v^r dz/dr + v^ph dz/dph + v^th dz/dth
  */

  // E.g., can find transformation matrix here: https://en.wikipedia.org/wiki/List_of_common_coordinate_transformations#From_spherical_coordinates

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

//Equation (138), arXiv:1310.3274v2.
//Calculate the vector potential in spherical coordinates then transform to cartesian
void GiRaFFEfood_ExactWald_exactA_i(CCTK_REAL &AxL, CCTK_REAL &AyL, CCTK_REAL &AzL,
                                    CCTK_REAL xL, CCTK_REAL yL, CCTK_REAL zL) {
  DECLARE_CCTK_PARAMETERS;

  //Coordinate transformation form Cartesian to spherical
  CCTK_REAL rL = sqrt(xL*xL + yL*yL + zL*zL);
  CCTK_REAL costh = zL/rL;
  CCTK_REAL sinth2 = 1.0 - costh*costh;

  // K-S coordinate r. Note that after the radial shift, theta and phi values at a given point REMAIN THE SAME.
  CCTK_REAL r_KS = rL + KerrSchild_radial_shift;

  //Equation (138). Calculate the vector potential in spherical coordinates.
  CCTK_REAL Aph = 0.5*Wald_B0*r_KS*r_KS*sinth2;

  /*
    A_x = A_ph * dph/dx = A_ph * [-y / (x^2 + y^2) ]
    A_y = A_ph * dph/dy = A_ph * [+x / (x^2 + y^2) ]
    A_z = A_ph * dph/dz = A_ph * [ 0 ] = 0
  */

  /* You will notice that when KerrSchild_radial_shift=0,
     the below expressions for Ax,Ay,Az reduce to:
     Ax = - y * Wald_B0/2;
     Ay = + x * Wald_B0/2;
     Az = 0;
     Thus prolongation & restriction operations will result in exactly the same A-field!
     THE SAME IS NOT TRUE FOR KerrSchild_radial_shift!=0,
     so expect larger initial errors in this case when performing AMR runs.
  */

  //Coordinate transformation from spherical to cartesian:
  //Note that we must be careful since we define A_i (the one-form) and not A^i (the vector).
  AxL = -yL*Aph/(xL*xL + yL*yL);
  AyL =  xL*Aph/(xL*xL + yL*yL);
  AzL =  0.0;
}

//Equations (153), (154), (155), arXiv:1310.3274v2.
//Calculate the velocity in spherical coordinates then transform to cartesian
void GiRaFFEfood_ExactWald_exactV_i(CCTK_REAL &VxL, CCTK_REAL &VyL, CCTK_REAL &VzL,
                                    CCTK_REAL xL, CCTK_REAL yL, CCTK_REAL zL) {
  DECLARE_CCTK_PARAMETERS;

  //Coordinate transformation from Cartesian to spherical:
  CCTK_REAL rL = sqrt(xL*xL + yL*yL + zL*zL);
  CCTK_REAL costh = zL/rL;
  CCTK_REAL sinth = sqrt(1.0 - costh*costh);
  CCTK_REAL sintwoth = 2.0*costh*sinth;

  // K-S coordinate r. Note that after the radial shift, theta and phi values at a given point REMAIN THE SAME.
  CCTK_REAL r_KS = rL + KerrSchild_radial_shift;

  //Set the velocities, this is Illinois formalism -- formulas (153, 154, 155), arXiv:1310.3274v2.
  //$v^r=-\frac{2M\cos^2\theta}{r+2M\cos^2\theta}$
  CCTK_REAL Vrh = - 2.0*BH_mass*costh*costh/(r_KS + 2.0*BH_mass*costh*costh);

  //$v^{\theta}=\frac{M\sin 2\theta}{r(r+2M\cos^2\theta}$
  CCTK_REAL Vth = BH_mass*sintwoth/(r_KS*r_KS + 2.0*r_KS*BH_mass*costh*costh);

  CCTK_REAL Vph = 0.0;

  //Now transform to Cartesian coordinates
  GiRaFFEfood_ExactWald_spherical_to_Cartesian_XformVector(VxL,VyL,VzL, Vrh,Vth,Vph,xL,yL,zL);
}

//Equation (140), (141), (142), arXiv:1310.3274v2.
//Calculate the magnetic field in spherical coordinates then transform to cartesian coordinates
void GiRaFFEfood_ExactWald_exactB_i(CCTK_REAL &BxL, CCTK_REAL &ByL, CCTK_REAL &BzL,
                                    CCTK_REAL xL, CCTK_REAL yL, CCTK_REAL zL) {
  DECLARE_CCTK_PARAMETERS;

  //Coordinate transformation from Cartesian to spherical:
  CCTK_REAL rL = sqrt(xL*xL + yL*yL + zL*zL);
  CCTK_REAL th = atan2(sqrt(xL*xL + yL*yL),zL);

  CCTK_REAL sinth = sin(th);
  CCTK_REAL costh = cos(th);

  // K-S coordinate r. Note that after the radial shift, theta and phi values at a given point REMAIN THE SAME.
  CCTK_REAL r_KS = rL + KerrSchild_radial_shift;

  //Set the magnetic field -- formulas (140, 141, 142)
  //$v^r=-\frac{2M\cos^2\theta}{r+2M\cos^2\theta}$
  CCTK_REAL Brh = + Wald_B0/sqrt(1.0 + 2.0*BH_mass/r_KS)*costh;

  //$v^{\theta}=\frac{M\sin 2\theta}{r(r+2M\cos^2\theta}$
  CCTK_REAL Bth = - Wald_B0/sqrt(1.0 + 2.0*BH_mass/r_KS)*sinth/r_KS;

  CCTK_REAL Bph = 0.0;

  //Now transform vector B^i from spherical to Cartesian coordinates
  GiRaFFEfood_ExactWald_spherical_to_Cartesian_XformVector(BxL,ByL,BzL, Brh,Bth,Bph,xL,yL,zL);
}
