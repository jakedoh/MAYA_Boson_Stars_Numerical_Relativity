#include <assert.h>
#include <math.h>

#include "cctk.h"
#include "EOS_Base.h"

#include "hll.h"

// Eq (76) of Del Zanna, Langrangian signal speed
void hll_lambda_prime(CCTK_REAL rho, CCTK_REAL eps, CCTK_REAL press, 
  CCTK_REAL velx, CCTK_REAL vely, CCTK_REAL velz, CCTK_REAL Bx, CCTK_REAL By, CCTK_REAL Bz,
  CCTK_REAL gxx, CCTK_REAL gxy, CCTK_REAL gxz, CCTK_REAL gyy, CCTK_REAL gyz, CCTK_REAL gzz,
  CCTK_REAL uxx, CCTK_INT eos_handle,
  CCTK_REAL *lambda_minus, CCTK_REAL *lambda_plus)
{
  CCTK_REAL a2, b2;   // DelZanna's a^2 from Eq 76, B^2
  CCTK_REAL cs2, ca2; // speed of sound, Alfv\'en speed
  CCTK_REAL v2, inv_w2, vB; // v^2, 1/w_lorentz^2, v^i*g_ij*B^j
  CCTK_REAL dpdrho, dpdeps, rho_enthalpy;
  CCTK_REAL sqrt_term;

  // speed of sound
  dpdrho = EOS_DPressByDRho(eos_handle,rho,eps);
  dpdeps = EOS_DPressByDEps(eos_handle,rho,eps);
  rho_enthalpy = rho*(1 + eps) + press;

  // checked that it agrees with Eq 73 for ideal gas
  cs2 = (rho*dpdrho + press*dpdeps/rho) / rho_enthalpy;
  
  v2 = gxx*SQR(velx) + gyy*SQR(vely) + gzz*SQR(velz)
        + 2.*gxy*velx*vely + 2.*gxz*velx*velz + 2.*gyz*vely*velz;
  inv_w2 = 1 - v2;

  // Alfv\'en speed
  vB = velx*(gxx*Bx + gxy*By + gxz*Bz) + 
       vely*(gxy*Bx + gyy*By + gyz*Bz) +
       velz*(gxz*Bx + gyz*By + gzz*Bz);
  b2 = (gxx*SQR(Bx) + gyy*SQR(By) + gzz*SQR(Bz)
        + 2.*gxy*Bx*By + 2.*gxz*Bx*Bz + 2.*gyz*By*Bz) * inv_w2 + SQR(vB);

  ca2 = b2/(rho_enthalpy + b2);

  a2 = cs2 + ca2 - cs2*ca2;

  // compute signal speeds Equ (76)
  sqrt_term = sqrt( a2*(1-v2) * ( (1-v2*a2)*uxx - (1-a2)*SQR(velx) ) );
  *lambda_minus = ( (1-a2)*velx - sqrt_term ) / (1-v2*a2);
  *lambda_plus  = ( (1-a2)*velx + sqrt_term ) / (1-v2*a2);
}
