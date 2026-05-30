#include <assert.h>
#include <math.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "hll.h"

// Equ. (41) of Del Zanna, actually used Equ (18) of Giacomazzo arXiv:gr-qc/0701109v2
// I use flux equation of Giacomazzo (differs by a factor of the lapse from that of Del Zanna)
void hll_numerical_flux(CCTK_REAL dens, 
    CCTK_REAL sx, CCTK_REAL sy, CCTK_REAL sz, CCTK_REAL tau,
    CCTK_REAL velx, CCTK_REAL vely, CCTK_REAL velz, 
    CCTK_REAL press, 
    CCTK_REAL alp, CCTK_REAL betax, CCTK_REAL betay, CCTK_REAL betaz,
    CCTK_REAL gxx, CCTK_REAL gxy, CCTK_REAL gxz, CCTK_REAL gyy, CCTK_REAL gyz, CCTK_REAL gzz, 
    CCTK_REAL sqrt_detg, 
    CCTK_REAL Bnx, CCTK_REAL Bny, CCTK_REAL Bnz, CCTK_REAL psi,
    CCTK_REAL *densflux, CCTK_REAL *sxflux, CCTK_REAL *syflux, CCTK_REAL *szflux,
    CCTK_REAL *tauflux, CCTK_REAL *Bnxflux, CCTK_REAL *Bnyflux, CCTK_REAL *Bnzflux, CCTK_REAL *psiflux)
{
  DECLARE_CCTK_PARAMETERS;
  CCTK_REAL bx, by, bz;               // b_i = F^*_{munu} u^mu, undensititzed
  CCTK_REAL Bnlowx, Bnlowy, Bnlowz;   // densitized
  CCTK_REAL vlowx, vlowy, vlowz;
  CCTK_REAL vdotB;                    // densitized
  CCTK_REAL udotB;                    // alpha b^0 = W B^i v_i = u_i B^i, densitized
  CCTK_REAL b2, sqrt_detg_press_star; // b_mu b^mu, P + b^2/2, densitized
  CCTK_REAL w_lorentz;

  vlowx = gxx*velx + gxy*vely + gxz*velz;
  vlowy = gxy*velx + gyy*vely + gyz*velz;
  vlowz = gxz*velx + gyz*vely + gzz*velz;

  w_lorentz = 1/sqrt(1-(vlowx*velx+vlowy*vely+vlowz*velz));
  
  Bnlowx = gxx*Bnx + gxy*Bny + gxz*Bnz;
  Bnlowy = gxy*Bnx + gyy*Bny + gyz*Bnz;
  Bnlowz = gxz*Bnx + gyz*Bny + gzz*Bnz;
  
  vdotB = velx*Bnlowx + vely*Bnlowy + velz*Bnlowz;
  udotB = w_lorentz * vdotB;

  // I use an undensitized b_i similar to the primitives which are undensitized
  bx = (Bnlowx/w_lorentz + vlowx*udotB) / sqrt_detg;
  by = (Bnlowy/w_lorentz + vlowy*udotB) / sqrt_detg;
  bz = (Bnlowz/w_lorentz + vlowz*udotB) / sqrt_detg;

  b2 = (Bnlowx*Bnx + Bnlowy*Bny + Bnlowz*Bnz) / SQR(w_lorentz) + SQR(vdotB);
  // very messy: P is undensitized, b2 contains two factors of sqrt_detg
  sqrt_detg_press_star = sqrt_detg * press + 0.5 * b2 / sqrt_detg; 

  *densflux = dens * ( velx - betax/alp );
  *sxflux = sx * ( velx - betax/alp ) - bx*Bnx/w_lorentz + sqrt_detg_press_star;
  *syflux = sy * ( velx - betax/alp ) - by*Bnx/w_lorentz;
  *szflux = sz * ( velx - betax/alp ) - bz*Bnx/w_lorentz;
  *tauflux = tau * ( velx - betax/alp ) + sqrt_detg_press_star*velx - vdotB / sqrt_detg * Bnx;

  // Whisky splits off a factor of alpha from the fluxes
  if (clean_divergence) {

        // This is the flux as derived from Penner's thesis and in GRHydro (remember factor alpha) 
        CCTK_REAL uxx = (-SQR(gyz) + gyy*gzz)/SQR(sqrt_detg);
        CCTK_REAL uxy = (gxz*gyz - gxy*gzz)/ SQR(sqrt_detg); 
        CCTK_REAL uxz = (gxy*gyz - gxz*gyy)/ SQR(sqrt_detg);
        *psiflux = Bnx - psi*betax/alp;
        *Bnxflux = - Bnx * betax/alp + sqrt_detg*uxx*psi;
        *Bnyflux = Bny * ( velx - betax/alp ) - vely*Bnx + sqrt_detg*uxy*psi;
        *Bnzflux = Bnz * ( velx - betax/alp ) - velz*Bnx + sqrt_detg*uxz*psi;
        // This is the flux from Lehner arXiv:1001.0575
        // *psiflux = Bnx / sqrt_detg;

  } else {
        *Bnxflux = 0.;
        *Bnyflux = Bny * ( velx - betax/alp ) - Bnx * ( vely - betay/alp );
        *Bnzflux = Bnz * ( velx - betax/alp ) - Bnx * ( velz - betaz/alp );
  }
}
