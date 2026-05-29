/* We evolve forward in time a set of functions called the
 * "conservative variables" (magnetic field and Poynting vector),
 * and any time the conserv's are updated, we must recover the
 * primitive variables (velocities), before reconstructing & evaluating
 * the RHSs of the MHD equations again.
 *
 * This file contains the routine for this algebraic calculation.
 * The velocity is calculated with formula (85), arXiv:1310.3274v2
 * $v^i = 4 \pi \alpha \gamma^{ij} {\tilde S}_j \gamma{-1/2} B^{-2} - \beta^i$
 * The force-free condition: $B^2>E^2$ is checked before computing the velocity.
 * and after imposing the constraint ${\tilde B}^i {\tilde S}_i = 0$

 * The procedure is as described in arXiv:1310.3274v2:
 * 1. ${\tilde S}_i ->{\tilde S}_i - ({\tilde S}_j {\tilde B}^j) {\tilde B}^i/{\tilde B}^2$
 * 2. $f = \sqrt{(1-\gamma_{max}^{-2}){\tilde B}^4/(16 \pi^2 \gamma {\tilde S}^2)}$
 * 3. ${\tilde S}_i -> {\tilde S}_i min(1,f)
 * 4. $v^i = 4 \pi \alpha \gamma^{ij} {\tilde S}_j \gamma{-1/2} B^{-2} - \beta^i$
 * 5. ${\tilde n}_i v^i = 0$
 *
 * All equations are from: http://arxiv.org/pdf/1310.3274.pdf (v2)
 * */

#include "cctk.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sys/time.h>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

#ifndef M_PI
#define M_PI 3.141592653589793238463
#endif

#include "GiRaFFE_headers.h"
#include "inlined_functions.C"

extern "C" void GiRaFFE_conserv_to_prims_FFE(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // We use proper C++ here, for file I/O later.
  using namespace std;

  // Here we convert ADM variables (from ADMBase) to the BSSN-based variables expected by this routine.
  GiRaFFE_convert_ADM_to_BSSN__enforce_detgtij_eq_1__and_compute_gtupij(cctkGH,cctk_lsh,  gxx,gxy,gxz,gyy,gyz,gzz,alp,
                                                                        gtxx,gtxy,gtxz,gtyy,gtyz,gtzz,
                                                                        gtupxx,gtupxy,gtupxz,gtupyy,gtupyz,gtupzz,
                                                                        phi_bssn,psi_bssn,lapm1);
  const int imin=0,jmin=0,kmin=0;
  const int imax=cctk_lsh[0],jmax=cctk_lsh[1],kmax=cctk_lsh[2];

  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);

  CCTK_REAL error_int_numer=0,error_int_denom=0;

  CCTK_INT num_vel_limits=0,num_vel_nulls_current_sheet=0;
#pragma omp parallel for reduction(+:error_int_numer,error_int_denom,num_vel_limits,num_vel_nulls_current_sheet) schedule(static)
  for(int k=kmin;k<kmax;k++)
    for(int j=jmin;j<jmax;j++)
      for(int i=imin;i<imax;i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        const CCTK_REAL rL = r[index];
        if(rL>min_radius_inside_of_which_conserv_to_prims_FFE_and_FFE_evolution_is_DISABLED) {

          const CCTK_REAL psi_bssnL = psi_bssn[index];
          const CCTK_REAL psi2 = psi_bssnL*psi_bssnL;
          const CCTK_REAL psi4 = psi2*psi2;
          const CCTK_REAL psim4 = 1.0/psi4;
          const CCTK_REAL sqrtg = psi4*psi2; // Determinant of 3-metric

          // \gamma_{ij}, computed from \tilde{\gamma}_{ij}
          const CCTK_REAL gxxL = psi4*gtxx[index];
          const CCTK_REAL gxyL = psi4*gtxy[index];
          const CCTK_REAL gxzL = psi4*gtxz[index];
          const CCTK_REAL gyyL = psi4*gtyy[index];
          const CCTK_REAL gyzL = psi4*gtyz[index];
          const CCTK_REAL gzzL = psi4*gtzz[index];

          // \gamma^{ij} = psim4 * \tilde{\gamma}^{ij}
          const CCTK_REAL gupxxL = psim4*gtupxx[index];
          const CCTK_REAL gupxyL = psim4*gtupxy[index];
          const CCTK_REAL gupxzL = psim4*gtupxz[index];
          const CCTK_REAL gupyyL = psim4*gtupyy[index];
          const CCTK_REAL gupyzL = psim4*gtupyz[index];
          const CCTK_REAL gupzzL = psim4*gtupzz[index];

          // Read in magnetic field and momentum variables once from memory, since memory access is expensive:
          const CCTK_REAL BxL = Bx[index];
          const CCTK_REAL ByL = By[index];
          const CCTK_REAL BzL = Bz[index];

          // End of page 7 on http://arxiv.org/pdf/1310.3274.pdf
          const CCTK_REAL BtildexL = BxL*sqrtg;
          const CCTK_REAL BtildeyL = ByL*sqrtg;
          const CCTK_REAL BtildezL = BzL*sqrtg;

          const CCTK_REAL Btilde_xL = gxxL*BtildexL + gxyL*BtildeyL + gxzL*BtildezL;
          const CCTK_REAL Btilde_yL = gxyL*BtildexL + gyyL*BtildeyL + gyzL*BtildezL;
          const CCTK_REAL Btilde_zL = gxzL*BtildexL + gyzL*BtildeyL + gzzL*BtildezL;

          CCTK_REAL mhd_st_xL = mhd_st_x[index];
          CCTK_REAL mhd_st_yL = mhd_st_y[index];
          CCTK_REAL mhd_st_zL = mhd_st_z[index];

          const CCTK_REAL mhd_st_x_orig = mhd_st_xL;
          const CCTK_REAL mhd_st_y_orig = mhd_st_yL;
          const CCTK_REAL mhd_st_z_orig = mhd_st_zL;

          const CCTK_REAL alpL = alp[index];
          const CCTK_REAL fourpialpha = 4.0*M_PI*alpL;

          const CCTK_REAL betaxL = betax[index];
          const CCTK_REAL betayL = betay[index];
          const CCTK_REAL betazL = betaz[index];

          //* 1. Just below Eq 90: Enforce orthogonality of B^i & S^i, so that B^i S_i = 0
          //*    Correction ${\tilde S}_i ->{\tilde S}_i - ({\tilde S}_j {\tilde B}^j) {\tilde B}_i/{\tilde B}^2$
          //*    NOTICE THAT THE {\tilde B}_i IS LOWERED, AS IT SHOULD BE. THIS IS A TYPO IN PASCHALIDIS ET AL.

          // First compute Btilde^i Stilde_i:
          const CCTK_REAL BtildeiSt_i = mhd_st_xL*BtildexL + mhd_st_yL*BtildeyL + mhd_st_zL*BtildezL;

          // Then compute (Btilde)^2
          const CCTK_REAL Btilde2 = gxxL*BtildexL*BtildexL + gyyL*BtildeyL*BtildeyL + gzzL*BtildezL*BtildezL
            + 2.0*(gxyL*BtildexL*BtildeyL + gxzL*BtildexL*BtildezL + gyzL*BtildeyL*BtildezL);

#define APPLY_GRFFE_FIXES

          // Now apply constraint: Stilde_i = Stilde_i - (Btilde^i Stilde_i) / (Btilde)^2
#ifdef APPLY_GRFFE_FIXES
          mhd_st_xL -= BtildeiSt_i*Btilde_xL/Btilde2;
          mhd_st_yL -= BtildeiSt_i*Btilde_yL/Btilde2;
          mhd_st_zL -= BtildeiSt_i*Btilde_zL/Btilde2;
#endif
          // Now that tildeS_i has been fixed, let's compute tildeS^i:
          CCTK_REAL mhd_st_upx = gupxxL*mhd_st_xL + gupxyL*mhd_st_yL + gupxzL*mhd_st_zL;
          CCTK_REAL mhd_st_upy = gupxyL*mhd_st_xL + gupyyL*mhd_st_yL + gupyzL*mhd_st_zL;
          CCTK_REAL mhd_st_upz = gupxzL*mhd_st_xL + gupyzL*mhd_st_yL + gupzzL*mhd_st_zL;

          // Just below Eq. 86 in http://arxiv.org/pdf/1310.3274.pdf:
          CCTK_REAL St2 = mhd_st_xL*mhd_st_upx + mhd_st_yL*mhd_st_upy + mhd_st_zL*mhd_st_upz;

          //* 2. Eq. 92: Factor $f = \sqrt{(1-\gamma_{max}^{-2}){\tilde B}^4/(16 \pi^2 \gamma {\tilde S}^2)}$

#ifdef APPLY_GRFFE_FIXES
          const CCTK_REAL gmax = GAMMA_SPEED_LIMIT;
          if(St2 > (1.0 - 1.0/(gmax*gmax))*Btilde2*Btilde2/ (16.0*M_PI*M_PI*sqrtg*sqrtg)) {
            const CCTK_REAL fact = sqrt((1.0 - 1.0/(gmax*gmax))/St2)*Btilde2/(4.0*M_PI*sqrtg);

            //* 3. ${\tilde S}_i -> {\tilde S}_i min(1,f)
            mhd_st_xL *= MIN(1.0,fact);
            mhd_st_yL *= MIN(1.0,fact);
            mhd_st_zL *= MIN(1.0,fact);

            // Recompute S^i
            mhd_st_upx = gupxxL*mhd_st_xL + gupxyL*mhd_st_yL + gupxzL*mhd_st_zL;
            mhd_st_upy = gupxyL*mhd_st_xL + gupyyL*mhd_st_yL + gupyzL*mhd_st_zL;
            mhd_st_upz = gupxzL*mhd_st_xL + gupyzL*mhd_st_yL + gupzzL*mhd_st_zL;
            /*
            printf("%e %e %e | %e %e %e | %e %e %e | oldgamma: %e %e should be > %e vfix\n",x[index],y[index],z[index],
                   BxL,ByL,BzL,
                   St2,(1.0 - 1.0/(gmax*gmax))*Btilde2*Btilde2/ (16.0*M_PI*M_PI*sqrtg*sqrtg),gmax,
                   sqrt(Btilde2 / (Btilde2 - 16*M_PI*M_PI*sqrtg*sqrtg * St2 / Btilde2) ) , Btilde2,16*M_PI*M_PI*sqrtg*sqrtg * St2 / Btilde2  );
            //exit(1);
            */
            // Recompute Stilde^2:
            St2 = mhd_st_xL*mhd_st_upx + mhd_st_yL*mhd_st_upy + mhd_st_zL*mhd_st_upz;

            if( St2 >= Btilde2*Btilde2/ (16.0*M_PI*M_PI*sqrtg*sqrtg) ) {
              printf("ERROR: Velocity cap fix wasn't effective; still have B^2 > E^2\n"); exit(1);
            }
            num_vel_limits++;
          }
#endif

          //* 4. Eq. 85: $v^i = 4 pi \alpha \gamma^{ij} {\tilde S}_j \gamma{-1/2} B^{-2} - \beta^i$:

          // See, e.g., Eq 71 in http://arxiv.org/pdf/1310.3274.pdf
          // ... or end of page 7 on http://arxiv.org/pdf/1310.3274.pdf:
          const CCTK_REAL B2 = Btilde2/(sqrtg*sqrtg);
          /*
             Eq. 75:
             v^i = \alpha \gamma^{ij} S_j / \mathcal{B}^2 - \beta^i
             Eq. 7: \mathcal{B}^{\mu} = B^{\mu}/\sqrt{4 \pi}
             -> v^i = 4 \pi \alpha \gamma^{ij} S_j / B^2 - \beta^i
             Eq. 79: \tilde{S_i} = \sqrt{\gamma} S_i
             -> v^i = 4 \pi \alpha \gamma^{ij} \tilde{S}_j / (\sqrt{\gamma} B^2) - \beta^i
          */
          const CCTK_REAL vxL = fourpialpha*mhd_st_upx/(sqrtg*B2) - betaxL;
          const CCTK_REAL vyL = fourpialpha*mhd_st_upy/(sqrtg*B2) - betayL;
          /* vzL not necessarily const! See below. */
          CCTK_REAL vzL = fourpialpha*mhd_st_upz/(sqrtg*B2) - betazL;

          //* 5. Eq. 94: ${\tilde n}_i v^i = 0$ in the current sheet region
          //     n^i is defined as the normal from the current sheet, which lies in the
          //     xy-plane (z=0). So n = (0,0,1)
#ifdef APPLY_GRFFE_FIXES
          if(current_sheet_null_v) {
            CCTK_REAL zL = z[index];
            if (fabs(zL) <= (4.0 + 1.0e-2)*dz ) {
              //vzL = 0.0;
              vzL = - (vxL*gxzL + vyL*gyzL) / gzzL;

              // vzL reset: TYPICALLY WOULD RESET CONSERVATIVES TO BE CONSISTENT. LET'S NOT DO THAT, TO AVOID MESSING UP B-FIELDS

              if(1==1) {
                CCTK_REAL PRIMS[MAXNUMVARS];
                int ww=0;
                PRIMS[ww] = vxL;    ww++;
                PRIMS[ww] = vyL;    ww++;
                PRIMS[ww] = vzL;    ww++;
                PRIMS[ww] = BxL;    ww++;
                PRIMS[ww] = ByL;    ww++;
                PRIMS[ww] = BzL;    ww++;

                CCTK_REAL METRIC[NUMVARS_FOR_METRIC],dummy=0;
                ww=0;
                // FIXME: NECESSARY?
                //psi_bssn[index] = exp(phi[index]);
                METRIC[ww] = phi_bssn[index];ww++;
                METRIC[ww] = dummy;          ww++; // Don't need to set psi.
                METRIC[ww] = gtxx[index];    ww++;
                METRIC[ww] = gtxy[index];    ww++;
                METRIC[ww] = gtxz[index];    ww++;
                METRIC[ww] = gtyy[index];    ww++;
                METRIC[ww] = gtyz[index];    ww++;
                METRIC[ww] = gtzz[index];    ww++;
                METRIC[ww] = lapm1[index];   ww++;
                METRIC[ww] = betax[index];   ww++;
                METRIC[ww] = betay[index];   ww++;
                METRIC[ww] = betaz[index];   ww++;
                METRIC[ww] = gtupxx[index];  ww++;
                METRIC[ww] = gtupyy[index];  ww++;
                METRIC[ww] = gtupzz[index];  ww++;
                METRIC[ww] = gtupxy[index];  ww++;
                METRIC[ww] = gtupxz[index];  ww++;
                METRIC[ww] = gtupyz[index];  ww++;

                CCTK_REAL CONSERVS[NUM_CONSERVS] = {0.0, 0.0, 0.0}; // 3 conservative variables: Stilde_x, Stilde_y, Stilde_z
                GiRaFFE_compute_conservatives(PRIMS,METRIC, CONSERVS);

                mhd_st_xL = CONSERVS[STILDEX];
                mhd_st_yL = CONSERVS[STILDEY];
                mhd_st_zL = CONSERVS[STILDEZ];
              }
              num_vel_nulls_current_sheet++;
            }
          }
#endif

          vx[index] = vxL;
          vy[index] = vyL;
          vz[index] = vzL;

          //Now we compute the difference between original & new conservatives, for diagnostic purposes:
          error_int_numer += fabs(mhd_st_xL - mhd_st_x_orig) + fabs(mhd_st_yL - mhd_st_y_orig) + fabs(mhd_st_zL - mhd_st_z_orig);
          error_int_denom += fabs(mhd_st_x_orig) + fabs(mhd_st_y_orig) + fabs(mhd_st_z_orig);
          /*
            if(fabs(vx_orig) > 1e-13 && fabs(vxL-vx_orig)/vx_orig > 1e-2) printf("BAD vx: %e %e | %e %e %e\n",vxL,vx_orig,x[index],y[index],z[index]);
            if(fabs(vy_orig) > 1e-13 && fabs(vyL-vy_orig)/vy_orig > 1e-2) printf("BAD vy: %e %e | %e %e %e\n",vyL,vy_orig,x[index],y[index],z[index]);
            if(fabs(vz_orig) > 1e-13 && fabs(vzL-vz_orig)/vz_orig > 1e-2) printf("BAD vz: %e %e | %e %e %e\n",vzL,vz_orig,x[index],y[index],z[index]);
            error_int_numer += fabs(vxL - vx_orig) + fabs(vyL - vy_orig) + fabs(vzL - vz_orig);
            error_int_denom += fabs(vx_orig) + fabs(vy_orig) + fabs(vz_orig);
          */


          mhd_st_x[index] = mhd_st_xL;
          mhd_st_y[index] = mhd_st_yL;
          mhd_st_z[index] = mhd_st_zL;
        }
      }

  CCTK_VInfo(CCTK_THORNSTRING,"FFEC2P: Lev: %d NumPts= %d | Error: %.3e, ErrDenom: %.3e, v_limits: %d / %d = %.3e, v_nulls: %d / %d = %.3e",
             (int)GetRefinementLevel(cctkGH),
             cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2],
             error_int_numer/(error_int_denom+1e-300),error_int_denom,
             /**/       num_vel_limits,            cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2],
             (CCTK_REAL)num_vel_limits/((CCTK_REAL)cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]),
             /**/       num_vel_nulls_current_sheet,            cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2],
             (CCTK_REAL)num_vel_nulls_current_sheet/((CCTK_REAL)cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]));
}

