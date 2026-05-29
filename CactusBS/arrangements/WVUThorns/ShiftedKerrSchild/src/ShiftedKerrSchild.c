/*
 * This function implements shifted Kerr-Schild coordinates,
 * which are the same as regular K-S coordinates, except with
 *    r_{shifted} = r_{KS} + KerrSchild_radial_shift*BH_mass
 * We draw extensively from the K-S coordinate prescription
 * provided in Baumgarte & Shapiro, "Numerical Relativity:
 * Solving Einstein's Equations on the Computer"
 * */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void ShiftedKS_spherical_to_Cartesian_XformVector(CCTK_REAL *VxL, CCTK_REAL *VyL, CCTK_REAL *VzL,
                                                  const CCTK_REAL Vrh,const CCTK_REAL Vth, const CCTK_REAL Vph,
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
  *VxL = Vrh*dx_dr + Vth*dx_dth + Vph*dx_dph;
  *VyL = Vrh*dy_dr + Vth*dy_dth + Vph*dy_dph;
  *VzL = Vrh*dz_dr + Vth*dz_dth + Vph*dz_dph;
}

void ShiftedKS_spherical_to_Cartesian_Xform_gij(CCTK_REAL *gxx,CCTK_REAL *gxy,CCTK_REAL *gxz,CCTK_REAL *gyy,CCTK_REAL *gyz,CCTK_REAL *gzz,
                                                const CCTK_REAL grr,const CCTK_REAL grth,const CCTK_REAL grph,const CCTK_REAL gthth,const CCTK_REAL gthph,const CCTK_REAL gphph,
                                                CCTK_REAL xL, CCTK_REAL yL, CCTK_REAL zL) {

  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL gijSph[4][4];
  gijSph[1][1]             =grr;
  gijSph[1][2]=gijSph[2][1]=grth;
  gijSph[1][3]=gijSph[3][1]=grph;
  gijSph[2][2]             =gthth;
  gijSph[2][3]=gijSph[3][2]=gthph;
  gijSph[3][3]             =gphph;

  //Coordinate transformation from Cartesian to spherical:
  CCTK_REAL rho = sqrt(xL*xL + yL*yL + zL*zL);

  //Transform. Note that we must be careful since we define v^i (the vector):
  /*
    v_x = v_r dr/dx + v_ph dph/dx + v_th dth/dx
    v_y = v_r dr/dy + v_ph dph/dy + v_th dth/dy
    v_z = v_r dr/dz + v_ph dph/dz + v_th dth/dz
  */

  // E.g., can find transformation matrix here: https://en.wikipedia.org/wiki/List_of_common_coordinate_transformations#From_spherical_coordinates

  CCTK_REAL Conv[4][4];
  Conv[1][1] = xL/rho; // dr_dx
  Conv[1][2] = yL/rho; // dr_dy
  Conv[1][3] = zL/rho; // dr_dz
  Conv[2][1] = xL*zL / (rho*rho*sqrt(xL*xL + yL*yL)); // dth_dx
  Conv[2][2] = yL*zL / (rho*rho*sqrt(xL*xL + yL*yL)); // dth_dy
  Conv[2][3] = - sqrt(xL*xL + yL*yL) / (rho*rho);     // dth_dz
  Conv[3][1] = -yL/(xL*xL + yL*yL); // dph_dx
  Conv[3][2] =  xL/(xL*xL + yL*yL); // dph_dy
  Conv[3][3] =  0.0;                // dph_dz

  //Transform.
  CCTK_REAL gijCart[4][4];
  for(int i=1;i<=3;i++) for(int j=1;j<=3;j++) gijCart[i][j] = 0.0;

  //gijCart_{ij} = Conv^k_i Conv^l_j gijSph_{kl}
  for(int i=1;i<=3;i++) for(int j=1;j<=3;j++) {
      for(int k=1;k<=3;k++) for(int l=1;l<=3;l++) {
          gijCart[i][j] += Conv[k][i]*Conv[l][j]*gijSph[k][l];
        }
    }

  *gxx = gijCart[1][1];
  *gxy = gijCart[1][2];
  *gxz = gijCart[1][3];
  *gyy = gijCart[2][2];
  *gyz = gijCart[2][3];
  *gzz = gijCart[3][3];
}

void ShiftedKS_ID_onept_all_but_Kij(const CCTK_REAL radial_shift,
                                    const CCTK_REAL xL, const CCTK_REAL yL, const CCTK_REAL zL,
                                    const CCTK_REAL a,  const CCTK_REAL BH_mass,
                                    CCTK_REAL *SKSgrrL, CCTK_REAL *SKSgrthL, CCTK_REAL *SKSgrphL,
                                    CCTK_REAL *SKSgththL, CCTK_REAL *SKSgthphL, CCTK_REAL *SKSgphphL,
                                    CCTK_REAL *SKSbetarL, CCTK_REAL *SKSbetathL, CCTK_REAL *SKSbetaphL,
                                    CCTK_REAL *gxxL, CCTK_REAL *gxyL, CCTK_REAL *gxzL,
                                    CCTK_REAL *gyyL, CCTK_REAL *gyzL, CCTK_REAL *gzzL,
                                    CCTK_REAL *alpL, CCTK_REAL *betaxL, CCTK_REAL *betayL, CCTK_REAL *betazL) {

  CCTK_REAL r_shiftedL = sqrt(xL*xL + yL*yL + zL*zL);
  CCTK_REAL r_KS = r_shiftedL + radial_shift;

  CCTK_REAL r_KS2 = r_KS*r_KS;
  // Angles are unchanged after a radial shift:
  CCTK_REAL costh = zL/r_shiftedL;
  CCTK_REAL sinth2 = 1.0 - costh*costh;

  // Kerr-Schild Coordinates, from:
  // http://relativity.livingreviews.org/open?pubNo=lrr-2000-5&page=articlesu8.html
  CCTK_REAL rho2 = r_KS2 + a*a*costh*costh;

  *alpL   = 1.0 / sqrt( 1.0 + 2.0*BH_mass*r_KS/rho2 );
  CCTK_REAL betar  = (*alpL)*(*alpL)*2.0*BH_mass*r_KS/rho2;
  CCTK_REAL betath = 0.0;
  CCTK_REAL betaph = 0.0;
  CCTK_REAL grr    = 1.0 + 2.0*BH_mass*r_KS/rho2;
  CCTK_REAL grth   = 0.0;
  CCTK_REAL grph   =-grr*a*sinth2;
  CCTK_REAL gthth  = rho2;
  CCTK_REAL gthph  = 0.0;
  CCTK_REAL gphph  = ( r_KS2 + a*a + 2.0*BH_mass*r_KS/rho2 * a*a*sinth2 ) * sinth2;

  *SKSbetarL  = betar;
  *SKSbetathL = betath;
  *SKSbetaphL = betaph;

  *SKSgrrL   = grr;
  *SKSgrthL  = grth;
  *SKSgrphL  = grph;
  *SKSgththL = gthth;
  *SKSgthphL = gthph;
  *SKSgphphL = gphph;

  ShiftedKS_spherical_to_Cartesian_XformVector(betaxL, betayL, betazL,
                                               betar,betath,betaph,
                                               xL,yL,zL);


  ShiftedKS_spherical_to_Cartesian_Xform_gij(gxxL,gxyL,gxzL,gyyL,gyzL,gzzL,
                                             grr,grth,grph,gthth,gthph,gphph,
                                             xL,yL,zL);
}

void ShiftedKS_ID_onept_Kij_only(const CCTK_REAL radial_shift,
                                 const CCTK_REAL xL, const CCTK_REAL yL, const CCTK_REAL zL,
                                 const CCTK_REAL a,  const CCTK_REAL BH_mass,
                                 CCTK_REAL *kxxL, CCTK_REAL *kxyL, CCTK_REAL *kxzL,
                                 CCTK_REAL *kyyL, CCTK_REAL *kyzL, CCTK_REAL *kzzL) {

  /*
    Eq 13 in http://link.springer.com/content/pdf/10.12942%2Flrr-2000-5.pdf ; Cook, G.B. Living Rev. Relativ. (2000) 3: 5. doi:10.12942/lrr-2000-5:
    *     Extrinsic curvature is given by K_{ij} = (Del_i N_j + Del_j N_i)/(2 alpha)*)
    "BoxRules" -> System`Convert`TeXFormDump`$GreekWords;
    a := a;
    BHmass := BHmass;
    rho2 = r^2 + a*a*Cos[theta]^2;
    alpL = 1/Sqrt[1 + 2*BHmass*r/rho2];
    beta[1] = alpL*alpL*2*BHmass*r/rho2;
    beta[2] = 0;
    beta[3] = 0;
    g[1][1] = 1 + 2*BHmass*r/rho2;
    g[1][2] = g[2][1] = 0;
    g[1][3] = g[3][1] = -g[1][1]*a*Sin[theta]^2;
    g[2][2] = rho2;
    g[2][3] = g[3][2] = 0;
    g[3][3] = (r^2 + a*a + 2*BHmass*r/rho2*a*a*Sin[theta]^2)*
    Sin[theta]^2;
    x[1] = r;
    x[2] = theta;
    x[3] = phi;
    gDN = Table[Subscript[m, i, j], {i, 3}, {j, 3}];
    Do[gDN[[i, j]] = g[i][j], {i, 1, 3}, {j, 1, 3}]
    gUP = FullSimplify[Inverse[gDN]];
    gUP // MatrixForm
    Do[GamD[m][k][l] =
    D[g[m][k], x[l]] + D[g[m][l], x[k]] - D[g[k][l], x[m]], {m, 1,
    3}, {k, 1, 3}, {l, 1, 3}]
    Do[Gam[i][k][l] = 0, {i, 1, 3}, {k, 1, 3}, {l, 1, 3}]
    Do[Gam[i][k][l] =
    Gam[i][k][l] +
    1/2*gUP[[i, m]]*
    GamD[m][k][l];(*Print[i,k,l,":  ",Gam[i][k][l]];*), {i, 1,
    3}, {k, 1, 3}, {l, 1, 3}, {m, 1, 3}]
    (*FullSimplify[Gam[3][1][2]-Gam[3][2][1]]*)
    Do[NN[i] = 0, {i, 1, 3}]
    Do[NN[i] = NN[i] + beta[j]*g[j][i], {i, 1, 3}, {j, 1, 3}]
    Do[DelN[i][j] = D[NN[i], x[j]], {i, 1, 3}, {j, 1, 3}]
    Do[DelN[aa][c] = DelN[aa][c] - Gam[b][c][aa]*NN[b], {aa, 1, 3}, {b, 1,
    3}, {c, 1, 3}]
    Do[K[i][j] = FullSimplify[(DelN[i][j] + DelN[j][i])/(2*alpL)], {i, 1,
    3}, {j, 1, 3}]
    (*Simplify[DelN[3][3]+DelN[3][3]]*)
    (*FullSimplify[K[1][1]]
    FullSimplify[K[1][2]]
    FullSimplify[K[1][3]]
    FullSimplify[K[2][2]]
    FullSimplify[K[2][3]]
    FullSimplify[K[3][3]]*)

    r := rKS;
    Do[Print["CCTK_REAL ", K, i, j, "=" CForm[FullSimplify[K[i][j]]],
    ";"], {i, 1, 3}, {j, i, 3}]
    (*Do[Print["\\fl K_{",i,j,"}&=&" TeXForm[FullSimplify[K[i][j]]]," \
    \\\\"],{i,1,3},{j,i,3}]*)
    r := r;
    BHmass := M;
    AA = (a^2 Cos[2 theta] + a^2 + 2 r^2);
    BB = AA + 4 BHmass*r;
    DD = Sqrt[(2 M r)/(a^2 Cos[theta]^2 + r^2) + 1];

    Print["\\fl K_{rr}&=&\\frac{D (A+2 M r)}{A^2 B}\\left[" TeXForm[
    FullSimplify[
    AA^2*BB/(DD*(AA + 2*BHmass*r))*K[1][1]]], "\\right] \\\\"]
    Print["\\fl K_{r\\theta}&=&\\frac{D}{AB}\\left[" TeXForm[
    FullSimplify[AA*BB/DD*K[1][2]]], "\\right] \\\\"]
    Print["\\fl K_{r\\phi}&=&\\frac{D}{A^2}\\left[" TeXForm[
    FullSimplify[AA^2/DD*K[1][3]]], "\\right] \\\\"]
    Print["\\fl K_{\\theta\\theta}&=&\\frac{D}{B}\\left[" TeXForm[
    FullSimplify[BB/DD*K[2][2]]], "\\right] \\\\"]
    Print["\\fl K_{\\theta\\phi}&=&\\frac{D}{AB}\\left[" TeXForm[
    FullSimplify[AA*BB/DD*K[2][3]]], "\\right] \\\\"]
    Print["\\fl K_{\\phi\\phi}&=&\\frac{D}{A^2 B}\\left[" TeXForm[
    FullSimplify[AA^2*BB/DD*K[3][3]]], "\\right] \\\\"]
  */

  CCTK_REAL theta = atan2(sqrt(xL*xL + yL*yL),zL);
  CCTK_REAL r_shiftedL = sqrt(xL*xL + yL*yL + zL*zL);
  CCTK_REAL rKS = r_shiftedL + radial_shift;

  CCTK_REAL Krr,Krth,Krph,Kthth,Kthph,Kphph;
  {
    // BEGIN MATHEMATICA CODE CHUNK. See above (large) code comment for full Mathematica notebook.
    /* Mathematica dislikes underscores in variable names.
     * Also it prefers to call the nonexistent Power(), Sqrt(), Sin(), and Cos() functions.
     * These #defines fix that.
     */
#define BHmass BH_mass
#define Power pow
#define Sqrt sqrt
#define Sin sin
#define Cos cos
    CCTK_REAL K11= (-4*BHmass*Sqrt(1 + (2*BHmass*rKS)/(Power(rKS,2) + Power(a,2)*Power(Cos(theta),2)))*(-Power(a,2) + 2*Power(rKS,2) - Power(a,2)*Cos(2*theta))*(Power(a,2) + 2*rKS*(BHmass + rKS) + Power(a,2)*Cos(2*theta)))/(Power(Power(a,2) + 2*Power(rKS,2) + Power(a,2)*Cos(2*theta),2)*(Power(a,2) + 2*rKS*(2*BHmass + rKS) + Power(a,2)*Cos(2*theta)));
    CCTK_REAL K12= (4*Power(a,2)*BHmass*rKS*Sqrt(1 + (2*BHmass*rKS)/(Power(rKS,2) + Power(a,2)*Power(Cos(theta),2)))*Sin(2*theta))/((Power(a,2) + 2*Power(rKS,2) + Power(a,2)*Cos(2*theta))*(Power(a,2) + 2*rKS*(2*BHmass + rKS) + Power(a,2)*Cos(2*theta)));
    CCTK_REAL K13= (-2*a*BHmass*Sqrt(1 + (2*BHmass*rKS)/(Power(rKS,2) + Power(a,2)*Power(Cos(theta),2)))*(Power(a,2) - 2*Power(rKS,2) + Power(a,2)*Cos(2*theta))*Power(Sin(theta),2))/Power(Power(a,2) + 2*Power(rKS,2) + Power(a,2)*Cos(2*theta),2);
    CCTK_REAL K22= (4*BHmass*Power(rKS,2)*Sqrt(1 + (2*BHmass*rKS)/(Power(rKS,2) + Power(a,2)*Power(Cos(theta),2))))/(Power(a,2) + 2*rKS*(2*BHmass + rKS) + Power(a,2)*Cos(2*theta));
    CCTK_REAL K23= (-8*Power(a,3)*BHmass*rKS*Cos(theta)*Sqrt(1 + (2*BHmass*rKS)/(Power(rKS,2) + Power(a,2)*Power(Cos(theta),2)))*Power(Sin(theta),3))/((Power(a,2) + 2*Power(rKS,2) + Power(a,2)*Cos(2*theta))*(Power(a,2) + 2*rKS*(2*BHmass + rKS) + Power(a,2)*Cos(2*theta)));
    CCTK_REAL K33= (2*BHmass*rKS*Sqrt(1 + (2*BHmass*rKS)/(Power(rKS,2) + Power(a,2)*Power(Cos(theta),2)))*(8*Power(rKS,5) + 4*Power(a,2)*Power(rKS,2)*(-BHmass + 2*rKS) + Power(a,4)*(BHmass + 3*rKS) + 4*Power(a,2)*rKS*(Power(a,2) + rKS*(BHmass + 2*rKS))*Cos(2*theta) + Power(a,4)*(-BHmass + rKS)*Cos(4*theta))*Power(Sin(theta),2))/(Power(Power(a,2) + 2*Power(rKS,2) + Power(a,2)*Cos(2*theta),2)*(Power(a,2) + 2*rKS*(2*BHmass + rKS) + Power(a,2)*Cos(2*theta)));
    // Undefine the Mathematica #defines, so they do not affect other areas of the code.
#undef BHmass
#undef Power
#undef Sqrt
#undef Sin
#undef Cos
    Krr=K11;Krth=K12;Krph=K13;Kthth=K22;Kthph=K23;Kphph=K33;
    // END MATHEMATICA CODE CHUNK
  }
  ShiftedKS_spherical_to_Cartesian_Xform_gij(kxxL,kxyL,kxzL,kyyL,kyzL,kzzL,
                                             Krr,Krth,Krph,Kthth,Kthph,Kphph,
                                             xL,yL,zL);

}

void ShiftedKS_ID(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const int do_metric = CCTK_EQUALS(initial_data , "ShiftedKerrSchild");
  const int do_lapse = CCTK_EQUALS(initial_lapse, "ShiftedKerrSchild");
  const int do_shift = CCTK_EQUALS(initial_shift, "ShiftedKerrSchild");

#pragma omp parallel for
  for(int idx=0;idx<cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; idx++) {
    /* CCTK_REAL a  = BH_spin; */

    CCTK_REAL xL = x[idx];
    CCTK_REAL yL = y[idx];
    CCTK_REAL zL = z[idx];

    CCTK_REAL alpL,betaxL,betayL,betazL,gxxL,gxyL,gxzL,gyyL,gyzL,gzzL;
    CCTK_REAL grrL,grthL,grphL,gththL,gthphL,gphphL,betarL,betathL,betaphL;

    ShiftedKS_ID_onept_all_but_Kij(KerrSchild_radial_shift,
                                   xL, yL, zL,
                                   BH_spin,  BH_mass,
                                   &grrL, &grthL, &grphL,
                                   &gththL, &gthphL, &gphphL,
                                   &betarL, &betathL, &betaphL,
                                   &gxxL, &gxyL, &gxzL,
                                   &gyyL, &gyzL, &gzzL,
                                   &alpL, &betaxL, &betayL, &betazL);

    SKSgrr[idx]   = grrL;
    SKSgrth[idx]  = grthL;
    SKSgrph[idx]  = grphL;
    SKSgthth[idx] = gththL;
    SKSgthph[idx] = gthphL;
    SKSgphph[idx] = gphphL;

    SKSbetar[idx]  = betarL;
    SKSbetath[idx] = betathL;
    SKSbetaph[idx] = betaphL;

    if(do_lapse) { alp[idx] = alpL; }

    //#define DEBUGSHIFTEDKS

#ifdef DEBUGSHIFTEDKS
    // GREAT FOR DEBUGGING; just enable KerrSchild or Exact/KerrSchild & verify we're in agreement with zero radial shift.
    if(fabs(alp[idx] - alpL)/alp[idx] > 1e-8) { printf("%e %e %e : %e | %.15e %.15e BADxx \n",xL,yL,zL, fabs(alp[idx] - alpL)/alp[idx], alp[idx], alpL);  }
    if(fabs(betax[idx] - betaxL)/betax[idx] > 1e-8) { printf("%e %e %e : %e | %.15e %.15e BADxx \n",xL,yL,zL, fabs(betax[idx] - betaxL)/betax[idx], betax[idx], betaxL);  }
    if(fabs(betay[idx] - betayL)/betay[idx] > 1e-8) { printf("%e %e %e : %e | %.15e %.15e BADxy \n",xL,yL,zL, fabs(betay[idx] - betayL)/betay[idx], betay[idx], betayL);  }
    if(fabs(betaz[idx] - betazL)/betaz[idx] > 1e-8) { printf("%e %e %e : %e | %.15e %.15e BADxz \n",xL,yL,zL, fabs(betaz[idx] - betazL)/betaz[idx], betaz[idx], betazL);  }
#endif

#if 0
#define BHmass BH_mass
#define Power pow
#define Sqrt sqrt
#define Sin sin
#define Cos cos
    // MATHEMATICA WORKSHEET
/* BHmass = 1; */
/* r0 = 0.4; */

/* rh = Sqrt[xL*xL + yL*yL + zL*zL]; */
/* rs = rh + r0; */
/* rho2 = rs*rs; */
/* th = ArcTan[zL, Sqrt[xL*xL + yL*yL]]; */
/* (*th=ArcTan[Sqrt[xL*xL+yL*yL]/zL]*) */
/* ph = ArcTan[xL, yL]; */
/* (*ph=ArcTan[yL/xL];*) */

/* alpL = 1/Sqrt[1 + 2*BHmass*rs/rho2]; */
/* betar = alpL*alpL*2*BHmass*rs/rho2; */
/* betath = 0; */
/* betaph = 0; */

/* sinth = Sin[th]; */
/* costh = Cos[th]; */
/* sinph = Sin[ph]; */
/* cosph = Cos[ph]; */

/* dxdr = sinth*cosph; */
/* dydr = sinth*sinph; */
/* dzdr = costh; */

/* dxdth = rh*costh*cosph; */
/* dydth = rh*costh*sinph; */
/* dzdth = -rh*sinth; */

/* dxdph = -rh*sinth*sinph; */
/* dydph = rh*sinth*cosph; */
/* dzdph = 0; */

/*   VxL = betar*dxdr + betath*dxdth + betaph*dxdph; */
/*   VyL = betar*dydr + betath*dydth + betaph*dydph; */
/*   VzL = betar*dzdr + betath*dzdth + betaph*dzdph; */

/* Print["CCTK_REAL betaxLtest=", CForm[FullSimplify[VxL]], ";\n", */
/*   "CCTK_REAL betayLtest=", CForm[FullSimplify[VyL]], ";\n", */
/*   "CCTK_REAL betazLtest=", CForm[FullSimplify[VzL]], ";"]; */
CCTK_REAL betaxLtest=(2.*xL)/(Power(xL,2) + Power(yL,2) + Power(zL,2) + 2.4*Sqrt(Power(xL,2) + Power(yL,2) + Power(zL,2)));
CCTK_REAL betayLtest=(2.*yL)/(Power(xL,2) + Power(yL,2) + Power(zL,2) + 2.4*Sqrt(Power(xL,2) + Power(yL,2) + Power(zL,2)));
CCTK_REAL betazLtest=(2.*zL)/(Power(xL,2) + Power(yL,2) + Power(zL,2) + 2.4*Sqrt(Power(xL,2) + Power(yL,2) + Power(zL,2)));
#undef BHmass
#undef Power
#undef Sqrt
#undef Sin
#undef Cos
#endif
    if(do_shift) {
      betax[idx] = betaxL;
      betay[idx] = betayL;
      betaz[idx] = betazL;
    }

    // Eq. 2.122 in Baumgarte & Shapiro "Numerical Relativity: Solving Einstein's Equations on the Computer"
#ifdef DEBUGSHIFTEDKS
    // GREAT FOR DEBUGGING; just enable KerrSchild or Exact/KerrSchild & verify we're in agreement with zero radial shift.
    if(fabs(gxx[idx] - gxxL)/gxx[idx] > 1e-8) { printf("%e %e %e : %e | %.15e %.15e BADxx \n",xL,yL,zL, fabs(gxx[idx] - gxxL)/gxx[idx], gxx[idx], gxxL);  }
    if(fabs(gxy[idx] - gxyL)/gxy[idx] > 1e-8) { printf("%e %e %e : %e | %.15e %.15e BADxy \n",xL,yL,zL, fabs(gxy[idx] - gxyL)/gxy[idx], gxy[idx], gxyL);  }
    if(fabs(gxz[idx] - gxzL)/gxz[idx] > 1e-8) { printf("%e %e %e : %e | %.15e %.15e BADxz \n",xL,yL,zL, fabs(gxz[idx] - gxzL)/gxz[idx], gxz[idx], gxzL);  }
    if(fabs(gyy[idx] - gyyL)/gyy[idx] > 1e-8) { printf("%e %e %e : %e | %.15e %.15e BADyy \n",xL,yL,zL, fabs(gyy[idx] - gyyL)/gyy[idx], gyy[idx], gyyL);  }
    if(fabs(gyz[idx] - gyzL)/gyz[idx] > 1e-8) { printf("%e %e %e : %e | %.15e %.15e BADyz \n",xL,yL,zL, fabs(gyz[idx] - gyzL)/gyz[idx], gyz[idx], gyzL);  }
    if(fabs(gzz[idx] - gzzL)/gzz[idx] > 1e-8) { printf("%e %e %e : %e | %.15e %.15e BADzz \n",xL,yL,zL, fabs(gzz[idx] - gzzL)/gzz[idx], gzz[idx], gzzL);  }
#endif

    if(do_metric) {
      gxx[idx] = gxxL;
      gxy[idx] = gxyL;
      gxz[idx] = gxzL;
      gyy[idx] = gyyL;
      gyz[idx] = gyzL;
      gzz[idx] = gzzL;
    }

    CCTK_REAL kxxL,kxyL,kxzL,kyyL,kyzL,kzzL;
    ShiftedKS_ID_onept_Kij_only(KerrSchild_radial_shift,
                                xL, yL, zL,
                                BH_spin,  BH_mass,
                                &kxxL, &kxyL, &kxzL,
                                &kyyL, &kyzL, &kzzL);

#ifdef DEBUGSHIFTEDKS
      // GREAT FOR DEBUGGING; just enable KerrSchild or Exact/KerrSchild & verify we're in agreement with zero radial shift.
      if(fabs(kxx[idx] - kxxL)/kxx[idx] > 1e-8) { printf("%e %e %e : %e | %.15e %.15e BADkxx \n",xL,yL,zL, fabs(kxx[idx] - kxxL)/kxx[idx], kxx[idx], kxxL);  }
      if(fabs(kxy[idx] - kxyL)/kxy[idx] > 1e-8) { printf("%e %e %e : %e | %.15e %.15e BADkxy \n",xL,yL,zL, fabs(kxy[idx] - kxyL)/kxy[idx], kxy[idx], kxyL);  }
      if(fabs(kxz[idx] - kxzL)/kxz[idx] > 1e-8) { printf("%e %e %e : %e | %.15e %.15e BADkxz \n",xL,yL,zL, fabs(kxz[idx] - kxzL)/kxz[idx], kxz[idx], kxzL);  }
      if(fabs(kyy[idx] - kyyL)/kyy[idx] > 1e-8) { printf("%e %e %e : %e | %.15e %.15e BADkyy \n",xL,yL,zL, fabs(kyy[idx] - kyyL)/kyy[idx], kyy[idx], kyyL);  }
      if(fabs(kyz[idx] - kyzL)/kyz[idx] > 1e-8) { printf("%e %e %e : %e | %.15e %.15e BADkyz \n",xL,yL,zL, fabs(kyz[idx] - kyzL)/kyz[idx], kyz[idx], kyzL);  }
      if(fabs(kzz[idx] - kzzL)/kzz[idx] > 1e-8) { printf("%e %e %e : %e | %.15e %.15e BADkzz \n",xL,yL,zL, fabs(kzz[idx] - kzzL)/kzz[idx], kzz[idx], kzzL);  }
#endif

    // Set these to NaN, as they're not needed for our tests.
    if(do_metric) {
      kxx[idx] = kxxL;
      kxy[idx] = kxyL;
      kxz[idx] = kxzL;
      kyy[idx] = kyyL;
      kyz[idx] = kyzL;
      kzz[idx] = kzzL;
    }
  }
}

void ShiftedKerrSchild_ParamCheck(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(!CCTK_EQUALS(metric_type, "physical")) {
    CCTK_VPARAMWARN("Only the 'physical' metric_type is supported, not '%s'",
                    metric_type);
  }
}
