#ifndef _HLL_H
#define _HLL_H_
#include <assert.h>
#include <stdio.h>

#include "cctk.h"
#include "cctk_Arguments.h"

// access one of Cactus' vector grid functions
#define vector_component(vec, comp) (&vec[comp*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])

// in hll.c
void RiemannSolverHLL_Solve(CCTK_ARGUMENTS);

// in flux.c
void hll_numerical_flux(CCTK_REAL dens, 
    CCTK_REAL sx, CCTK_REAL sy, CCTK_REAL sz, CCTK_REAL tau,
    CCTK_REAL velx, CCTK_REAL vely, CCTK_REAL velz, 
    CCTK_REAL press, 
    CCTK_REAL alp, CCTK_REAL betax, CCTK_REAL betay, CCTK_REAL betaz,
    CCTK_REAL gxx, CCTK_REAL gxy, CCTK_REAL gxz, CCTK_REAL gyy, CCTK_REAL gyz, CCTK_REAL gzz, 
    CCTK_REAL sqrt_detg, 
    CCTK_REAL Bx, CCTK_REAL By, CCTK_REAL Bz, CCTK_REAL psi,
    CCTK_REAL *densflux, CCTK_REAL *sxflux, CCTK_REAL *syflux, CCTK_REAL *szflux,
    CCTK_REAL *tauflux, CCTK_REAL *Bxflux, CCTK_REAL *Byflux, CCTK_REAL *Bzflux, CCTK_REAL *psiflux);

// in signal_speed.c
void hll_lambda_prime(CCTK_REAL rho, CCTK_REAL eps, CCTK_REAL press, 
  CCTK_REAL velx, CCTK_REAL vely, CCTK_REAL velz, CCTK_REAL Bx, CCTK_REAL By, CCTK_REAL Bz,
  CCTK_REAL gxx, CCTK_REAL gxy, CCTK_REAL gxz, CCTK_REAL gyy, CCTK_REAL gyz, CCTK_REAL gzz,
  CCTK_REAL uXX, CCTK_INT eoshandle,
  CCTK_REAL *lambda_minus, CCTK_REAL *lambda_plus);


// some math
static inline CCTK_REAL SQR(const CCTK_REAL x)
{
  return x*x;
}
static inline CCTK_REAL max3(CCTK_REAL a, CCTK_REAL b, CCTK_REAL c)
{
  CCTK_REAL ans=a;
  if ( b > ans )
	ans=b;
  if ( c > ans )
	ans=c;

  return ans;
}
static inline CCTK_REAL min3(CCTK_REAL a, CCTK_REAL b, CCTK_REAL c)
{

  CCTK_REAL ans=a;
  if ( b < ans )
	ans = b;
  if ( c < ans )
	ans = c;
  return ans;
}
#endif /* _HLL_H_ */
