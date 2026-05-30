/*  file produced by user hinder, 15/11/2006 */
/*  Produced with Mathematica Version 5.2 for Linux (June 20, 2005) */

/*  Mathematica script written by Ian Hinder and Sascha Husa */

/*  $Id$ */

#define KRANC_C

#include <math.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "GenericFD.h"
#include "Differencing.h"

extern void WeylScal4_levelinfo(const cGH * cctkGH,
                              CCTK_INT * my_level,
                              CCTK_INT * n_levels);
extern CCTK_INT WeylScal4_rnd_iteration(CCTK_INT it);

/* Define macros used in calculations */
#define INITVALUE  (42)
#define INV(x) ((1.0) / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))
#define QAD(x) ((x) * (x) * (x) * (x))

void psis_calc_2nd_Body(cGH *cctkGH, CCTK_INT dir, CCTK_INT face, CCTK_REAL normal[3], CCTK_REAL tangent1[3], CCTK_REAL tangent2[3], CCTK_INT min[3], CCTK_INT max[3], CCTK_INT n_subblock_gfs, CCTK_REAL *subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  /* Refinement level we are on */
  CCTK_INT reflevel, nlevels;
  
  /* Declare the variables used for looping over grid points */
  CCTK_INT i = INITVALUE, j = INITVALUE, k = INITVALUE;
  CCTK_INT index = INITVALUE;
  CCTK_INT subblock_index = INITVALUE;
  
  /* Declare finite differencing variables */
  CCTK_REAL dx = INITVALUE, dy = INITVALUE, dz = INITVALUE;
  CCTK_REAL dxi = INITVALUE, dyi = INITVALUE, dzi = INITVALUE;
  CCTK_REAL khalf = INITVALUE, kthird = INITVALUE, ktwothird = INITVALUE, kfourthird = INITVALUE, keightthird = INITVALUE;
  CCTK_REAL hdxi = INITVALUE, hdyi = INITVALUE, hdzi = INITVALUE;
  
  /* Declare shorthands */
  CCTK_REAL detg = INITVALUE;
  CCTK_REAL gamma111 = INITVALUE, gamma121 = INITVALUE, gamma122 = INITVALUE, gamma131 = INITVALUE, gamma132 = INITVALUE, gamma133 = INITVALUE;
  CCTK_REAL gamma211 = INITVALUE, gamma221 = INITVALUE, gamma222 = INITVALUE, gamma231 = INITVALUE, gamma232 = INITVALUE, gamma233 = INITVALUE;
  CCTK_REAL gamma311 = INITVALUE, gamma321 = INITVALUE, gamma322 = INITVALUE, gamma331 = INITVALUE, gamma332 = INITVALUE, gamma333 = INITVALUE;
  CCTK_REAL gInv11 = INITVALUE, gInv12 = INITVALUE, gInv13 = INITVALUE, gInv21 = INITVALUE, gInv22 = INITVALUE, gInv23 = INITVALUE;
  CCTK_REAL gInv31 = INITVALUE, gInv32 = INITVALUE, gInv33 = INITVALUE;
  CCTK_REAL im1 = INITVALUE, im2 = INITVALUE, im3 = INITVALUE;
  CCTK_REAL imbar1 = INITVALUE, imbar2 = INITVALUE, imbar3 = INITVALUE;
  CCTK_REAL invdetg = INITVALUE;
  CCTK_REAL isqrt2 = INITVALUE;
  CCTK_REAL n0 = INITVALUE, n1 = INITVALUE, n2 = INITVALUE, n3 = INITVALUE;
  CCTK_REAL omega11 = INITVALUE, omega12 = INITVALUE, omega13 = INITVALUE, omega22 = INITVALUE, omega23 = INITVALUE, omega33 = INITVALUE;
  CCTK_REAL R1111 = INITVALUE, R1112 = INITVALUE, R1113 = INITVALUE, R1121 = INITVALUE, R1122 = INITVALUE, R1123 = INITVALUE;
  CCTK_REAL R1131 = INITVALUE, R1132 = INITVALUE, R1133 = INITVALUE, R1211 = INITVALUE, R1212 = INITVALUE, R1213 = INITVALUE;
  CCTK_REAL R1221 = INITVALUE, R1222 = INITVALUE, R1223 = INITVALUE, R1231 = INITVALUE, R1232 = INITVALUE, R1233 = INITVALUE;
  CCTK_REAL R1311 = INITVALUE, R1312 = INITVALUE, R1313 = INITVALUE, R1321 = INITVALUE, R1322 = INITVALUE, R1323 = INITVALUE;
  CCTK_REAL R1331 = INITVALUE, R1332 = INITVALUE, R1333 = INITVALUE, R2111 = INITVALUE, R2112 = INITVALUE, R2113 = INITVALUE;
  CCTK_REAL R2121 = INITVALUE, R2122 = INITVALUE, R2123 = INITVALUE, R2131 = INITVALUE, R2132 = INITVALUE, R2133 = INITVALUE;
  CCTK_REAL R2211 = INITVALUE, R2212 = INITVALUE, R2213 = INITVALUE, R2221 = INITVALUE, R2222 = INITVALUE, R2223 = INITVALUE;
  CCTK_REAL R2231 = INITVALUE, R2232 = INITVALUE, R2233 = INITVALUE, R2311 = INITVALUE, R2312 = INITVALUE, R2313 = INITVALUE;
  CCTK_REAL R2321 = INITVALUE, R2322 = INITVALUE, R2323 = INITVALUE, R2331 = INITVALUE, R2332 = INITVALUE, R2333 = INITVALUE;
  CCTK_REAL R3111 = INITVALUE, R3112 = INITVALUE, R3113 = INITVALUE, R3121 = INITVALUE, R3122 = INITVALUE, R3123 = INITVALUE;
  CCTK_REAL R3131 = INITVALUE, R3132 = INITVALUE, R3133 = INITVALUE, R3211 = INITVALUE, R3212 = INITVALUE, R3213 = INITVALUE;
  CCTK_REAL R3221 = INITVALUE, R3222 = INITVALUE, R3223 = INITVALUE, R3231 = INITVALUE, R3232 = INITVALUE, R3233 = INITVALUE;
  CCTK_REAL R3311 = INITVALUE, R3312 = INITVALUE, R3313 = INITVALUE, R3321 = INITVALUE, R3322 = INITVALUE, R3323 = INITVALUE;
  CCTK_REAL R3331 = INITVALUE, R3332 = INITVALUE, R3333 = INITVALUE;
  CCTK_REAL R4p1111 = INITVALUE, R4p1112 = INITVALUE, R4p1113 = INITVALUE, R4p1121 = INITVALUE, R4p1122 = INITVALUE, R4p1123 = INITVALUE;
  CCTK_REAL R4p1131 = INITVALUE, R4p1132 = INITVALUE, R4p1133 = INITVALUE, R4p1211 = INITVALUE, R4p1212 = INITVALUE, R4p1213 = INITVALUE;
  CCTK_REAL R4p1221 = INITVALUE, R4p1222 = INITVALUE, R4p1223 = INITVALUE, R4p1231 = INITVALUE, R4p1232 = INITVALUE, R4p1233 = INITVALUE;
  CCTK_REAL R4p1311 = INITVALUE, R4p1312 = INITVALUE, R4p1313 = INITVALUE, R4p1321 = INITVALUE, R4p1322 = INITVALUE, R4p1323 = INITVALUE;
  CCTK_REAL R4p1331 = INITVALUE, R4p1332 = INITVALUE, R4p1333 = INITVALUE, R4p2111 = INITVALUE, R4p2112 = INITVALUE, R4p2113 = INITVALUE;
  CCTK_REAL R4p2121 = INITVALUE, R4p2122 = INITVALUE, R4p2123 = INITVALUE, R4p2131 = INITVALUE, R4p2132 = INITVALUE, R4p2133 = INITVALUE;
  CCTK_REAL R4p2211 = INITVALUE, R4p2212 = INITVALUE, R4p2213 = INITVALUE, R4p2221 = INITVALUE, R4p2222 = INITVALUE, R4p2223 = INITVALUE;
  CCTK_REAL R4p2231 = INITVALUE, R4p2232 = INITVALUE, R4p2233 = INITVALUE, R4p2311 = INITVALUE, R4p2312 = INITVALUE, R4p2313 = INITVALUE;
  CCTK_REAL R4p2321 = INITVALUE, R4p2322 = INITVALUE, R4p2323 = INITVALUE, R4p2331 = INITVALUE, R4p2332 = INITVALUE, R4p2333 = INITVALUE;
  CCTK_REAL R4p3111 = INITVALUE, R4p3112 = INITVALUE, R4p3113 = INITVALUE, R4p3121 = INITVALUE, R4p3122 = INITVALUE, R4p3123 = INITVALUE;
  CCTK_REAL R4p3131 = INITVALUE, R4p3132 = INITVALUE, R4p3133 = INITVALUE, R4p3211 = INITVALUE, R4p3212 = INITVALUE, R4p3213 = INITVALUE;
  CCTK_REAL R4p3221 = INITVALUE, R4p3222 = INITVALUE, R4p3223 = INITVALUE, R4p3231 = INITVALUE, R4p3232 = INITVALUE, R4p3233 = INITVALUE;
  CCTK_REAL R4p3311 = INITVALUE, R4p3312 = INITVALUE, R4p3313 = INITVALUE, R4p3321 = INITVALUE, R4p3322 = INITVALUE, R4p3323 = INITVALUE;
  CCTK_REAL R4p3331 = INITVALUE, R4p3332 = INITVALUE, R4p3333 = INITVALUE;
  CCTK_REAL rm1 = INITVALUE, rm2 = INITVALUE, rm3 = INITVALUE;
  CCTK_REAL rmbar1 = INITVALUE, rmbar2 = INITVALUE, rmbar3 = INITVALUE;
  CCTK_REAL Ro111 = INITVALUE, Ro112 = INITVALUE, Ro113 = INITVALUE, Ro121 = INITVALUE, Ro122 = INITVALUE, Ro123 = INITVALUE;
  CCTK_REAL Ro131 = INITVALUE, Ro132 = INITVALUE, Ro133 = INITVALUE, Ro211 = INITVALUE, Ro212 = INITVALUE, Ro213 = INITVALUE;
  CCTK_REAL Ro221 = INITVALUE, Ro222 = INITVALUE, Ro223 = INITVALUE, Ro231 = INITVALUE, Ro232 = INITVALUE, Ro233 = INITVALUE;
  CCTK_REAL Ro311 = INITVALUE, Ro312 = INITVALUE, Ro313 = INITVALUE, Ro321 = INITVALUE, Ro322 = INITVALUE, Ro323 = INITVALUE;
  CCTK_REAL Ro331 = INITVALUE, Ro332 = INITVALUE, Ro333 = INITVALUE;
  CCTK_REAL Rojo11 = INITVALUE, Rojo12 = INITVALUE, Rojo13 = INITVALUE, Rojo21 = INITVALUE, Rojo22 = INITVALUE, Rojo23 = INITVALUE;
  CCTK_REAL Rojo31 = INITVALUE, Rojo32 = INITVALUE, Rojo33 = INITVALUE;
  CCTK_REAL va1 = INITVALUE, va2 = INITVALUE, va3 = INITVALUE;
  CCTK_REAL vb1 = INITVALUE, vb2 = INITVALUE, vb3 = INITVALUE;
  CCTK_REAL vc1 = INITVALUE, vc2 = INITVALUE, vc3 = INITVALUE;
  CCTK_REAL xmoved = INITVALUE, ymoved = INITVALUE, zmoved = INITVALUE;
  
  /* Declare local copies of grid functions */
  CCTK_REAL gxxL = INITVALUE;
  CCTK_REAL gxyL = INITVALUE;
  CCTK_REAL gxzL = INITVALUE;
  CCTK_REAL gyyL = INITVALUE;
  CCTK_REAL gyzL = INITVALUE;
  CCTK_REAL gzzL = INITVALUE;
  CCTK_REAL kxxL = INITVALUE;
  CCTK_REAL kxyL = INITVALUE;
  CCTK_REAL kxzL = INITVALUE;
  CCTK_REAL kyyL = INITVALUE;
  CCTK_REAL kyzL = INITVALUE;
  CCTK_REAL kzzL = INITVALUE;
  CCTK_REAL Psi4iL = INITVALUE;
  CCTK_REAL Psi4rL = INITVALUE;
  CCTK_REAL xL = INITVALUE;
  CCTK_REAL yL = INITVALUE;
  CCTK_REAL zL = INITVALUE;
  
  
  CCTK_REAL PDstandard2nd1gxx = INITVALUE;
  CCTK_REAL PDstandard2nd2gxx = INITVALUE;
  CCTK_REAL PDstandard2nd3gxx = INITVALUE;
  CCTK_REAL PDstandard2nd22gxx = INITVALUE;
  CCTK_REAL PDstandard2nd33gxx = INITVALUE;
  CCTK_REAL PDstandard2nd12gxx = INITVALUE;
  CCTK_REAL PDstandard2nd13gxx = INITVALUE;
  CCTK_REAL PDstandard2nd21gxx = INITVALUE;
  CCTK_REAL PDstandard2nd23gxx = INITVALUE;
  CCTK_REAL PDstandard2nd31gxx = INITVALUE;
  CCTK_REAL PDstandard2nd32gxx = INITVALUE;
  CCTK_REAL PDstandard2nd1gxy = INITVALUE;
  CCTK_REAL PDstandard2nd2gxy = INITVALUE;
  CCTK_REAL PDstandard2nd3gxy = INITVALUE;
  CCTK_REAL PDstandard2nd33gxy = INITVALUE;
  CCTK_REAL PDstandard2nd12gxy = INITVALUE;
  CCTK_REAL PDstandard2nd13gxy = INITVALUE;
  CCTK_REAL PDstandard2nd21gxy = INITVALUE;
  CCTK_REAL PDstandard2nd23gxy = INITVALUE;
  CCTK_REAL PDstandard2nd31gxy = INITVALUE;
  CCTK_REAL PDstandard2nd32gxy = INITVALUE;
  CCTK_REAL PDstandard2nd1gxz = INITVALUE;
  CCTK_REAL PDstandard2nd2gxz = INITVALUE;
  CCTK_REAL PDstandard2nd3gxz = INITVALUE;
  CCTK_REAL PDstandard2nd22gxz = INITVALUE;
  CCTK_REAL PDstandard2nd12gxz = INITVALUE;
  CCTK_REAL PDstandard2nd13gxz = INITVALUE;
  CCTK_REAL PDstandard2nd21gxz = INITVALUE;
  CCTK_REAL PDstandard2nd23gxz = INITVALUE;
  CCTK_REAL PDstandard2nd31gxz = INITVALUE;
  CCTK_REAL PDstandard2nd32gxz = INITVALUE;
  CCTK_REAL PDstandard2nd1gyy = INITVALUE;
  CCTK_REAL PDstandard2nd2gyy = INITVALUE;
  CCTK_REAL PDstandard2nd3gyy = INITVALUE;
  CCTK_REAL PDstandard2nd11gyy = INITVALUE;
  CCTK_REAL PDstandard2nd33gyy = INITVALUE;
  CCTK_REAL PDstandard2nd12gyy = INITVALUE;
  CCTK_REAL PDstandard2nd13gyy = INITVALUE;
  CCTK_REAL PDstandard2nd21gyy = INITVALUE;
  CCTK_REAL PDstandard2nd23gyy = INITVALUE;
  CCTK_REAL PDstandard2nd31gyy = INITVALUE;
  CCTK_REAL PDstandard2nd32gyy = INITVALUE;
  CCTK_REAL PDstandard2nd1gyz = INITVALUE;
  CCTK_REAL PDstandard2nd2gyz = INITVALUE;
  CCTK_REAL PDstandard2nd3gyz = INITVALUE;
  CCTK_REAL PDstandard2nd11gyz = INITVALUE;
  CCTK_REAL PDstandard2nd12gyz = INITVALUE;
  CCTK_REAL PDstandard2nd13gyz = INITVALUE;
  CCTK_REAL PDstandard2nd21gyz = INITVALUE;
  CCTK_REAL PDstandard2nd23gyz = INITVALUE;
  CCTK_REAL PDstandard2nd31gyz = INITVALUE;
  CCTK_REAL PDstandard2nd32gyz = INITVALUE;
  CCTK_REAL PDstandard2nd1gzz = INITVALUE;
  CCTK_REAL PDstandard2nd2gzz = INITVALUE;
  CCTK_REAL PDstandard2nd3gzz = INITVALUE;
  CCTK_REAL PDstandard2nd11gzz = INITVALUE;
  CCTK_REAL PDstandard2nd22gzz = INITVALUE;
  CCTK_REAL PDstandard2nd12gzz = INITVALUE;
  CCTK_REAL PDstandard2nd13gzz = INITVALUE;
  CCTK_REAL PDstandard2nd21gzz = INITVALUE;
  CCTK_REAL PDstandard2nd23gzz = INITVALUE;
  CCTK_REAL PDstandard2nd31gzz = INITVALUE;
  CCTK_REAL PDstandard2nd32gzz = INITVALUE;
  CCTK_REAL PDstandard2nd2kxx = INITVALUE;
  CCTK_REAL PDstandard2nd3kxx = INITVALUE;
  CCTK_REAL PDstandard2nd1kxy = INITVALUE;
  CCTK_REAL PDstandard2nd2kxy = INITVALUE;
  CCTK_REAL PDstandard2nd3kxy = INITVALUE;
  CCTK_REAL PDstandard2nd1kxz = INITVALUE;
  CCTK_REAL PDstandard2nd2kxz = INITVALUE;
  CCTK_REAL PDstandard2nd3kxz = INITVALUE;
  CCTK_REAL PDstandard2nd1kyy = INITVALUE;
  CCTK_REAL PDstandard2nd3kyy = INITVALUE;
  CCTK_REAL PDstandard2nd1kyz = INITVALUE;
  CCTK_REAL PDstandard2nd2kyz = INITVALUE;
  CCTK_REAL PDstandard2nd3kyz = INITVALUE;
  CCTK_REAL PDstandard2nd1kzz = INITVALUE;
  CCTK_REAL PDstandard2nd2kzz = INITVALUE;
  
  /* Declare predefined quantities */
  CCTK_REAL p1o12dx = INITVALUE;
  CCTK_REAL p1o12dy = INITVALUE;
  CCTK_REAL p1o12dz = INITVALUE;
  CCTK_REAL p1o144dxdy = INITVALUE;
  CCTK_REAL p1o144dxdz = INITVALUE;
  CCTK_REAL p1o144dydz = INITVALUE;
  CCTK_REAL p1o2dx = INITVALUE;
  CCTK_REAL p1o2dy = INITVALUE;
  CCTK_REAL p1o2dz = INITVALUE;
  CCTK_REAL p1o4dxdy = INITVALUE;
  CCTK_REAL p1o4dxdz = INITVALUE;
  CCTK_REAL p1o4dydz = INITVALUE;
  CCTK_REAL p1odx = INITVALUE;
  CCTK_REAL p1odx2 = INITVALUE;
  CCTK_REAL p1ody = INITVALUE;
  CCTK_REAL p1ody2 = INITVALUE;
  CCTK_REAL p1odz = INITVALUE;
  CCTK_REAL p1odz2 = INITVALUE;
  CCTK_REAL pm1o12dx2 = INITVALUE;
  CCTK_REAL pm1o12dy2 = INITVALUE;
  CCTK_REAL pm1o12dz2 = INITVALUE;
  CCTK_REAL pm1o2dx = INITVALUE;
  CCTK_REAL pm1o2dy = INITVALUE;
  CCTK_REAL pm1o2dz = INITVALUE;
  CCTK_REAL pm1odxdy = INITVALUE;
  CCTK_REAL pm1odxdz = INITVALUE;
  CCTK_REAL pm1odydz = INITVALUE;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering psis_calc_2nd_Body");
  }
  
  /* Carpet lies about the iteration number during EVOL, so we have to round to
   * the nearest "valid" iteration number 
   */
  if (WeylScal4_rnd_iteration(cctk_iteration) % psis_calc_2nd_calc_every != psis_calc_2nd_calc_offset)
  {
    return;
  }

  WeylScal4_levelinfo(cctkGH, &reflevel, &nlevels);
  if (calc_max_reflevel != -1 && reflevel > calc_max_reflevel)
  {
    return;
  }
  
  /* Include user-supplied include files */
  
  /* Initialise finite differencing variables */
  dx = CCTK_DELTA_SPACE(0);
  dy = CCTK_DELTA_SPACE(1);
  dz = CCTK_DELTA_SPACE(2);
  dxi = 1.0 / dx;
  dyi = 1.0 / dy;
  dzi = 1.0 / dz;
  khalf = 0.5;
  kthird = 1/3.0;
  ktwothird = 2.0/3.0;
  kfourthird = 4.0/3.0;
  keightthird = 8.0/3.0;
  hdxi = 0.5 * dxi;
  hdyi = 0.5 * dyi;
  hdzi = 0.5 * dzi;
  
  /* Initialize predefined quantities */
  p1o12dx = INV(dx)/12.;
  p1o12dy = INV(dy)/12.;
  p1o12dz = INV(dz)/12.;
  p1o144dxdy = (INV(dx)*INV(dy))/144.;
  p1o144dxdz = (INV(dx)*INV(dz))/144.;
  p1o144dydz = (INV(dy)*INV(dz))/144.;
  p1o2dx = khalf*INV(dx);
  p1o2dy = khalf*INV(dy);
  p1o2dz = khalf*INV(dz);
  p1o4dxdy = (INV(dx)*INV(dy))/4.;
  p1o4dxdz = (INV(dx)*INV(dz))/4.;
  p1o4dydz = (INV(dy)*INV(dz))/4.;
  p1odx = INV(dx);
  p1odx2 = pow(dx,-2);
  p1ody = INV(dy);
  p1ody2 = pow(dy,-2);
  p1odz = INV(dz);
  p1odz2 = pow(dz,-2);
  pm1o12dx2 = -pow(dx,-2)/12.;
  pm1o12dy2 = -pow(dy,-2)/12.;
  pm1o12dz2 = -pow(dz,-2)/12.;
  pm1o2dx = -(khalf*INV(dx));
  pm1o2dy = -(khalf*INV(dy));
  pm1o2dz = -(khalf*INV(dz));
  pm1odxdy = -(INV(dx)*INV(dy));
  pm1odxdz = -(INV(dx)*INV(dz));
  pm1odydz = -(INV(dy)*INV(dz));
  
  /* Loop over the grid points */
  for (k = min[2]; k < max[2]; k++)
  {
    for (j = min[1]; j < max[1]; j++)
    {
      for (i = min[0]; i < max[0]; i++)
      {
         index  =  CCTK_GFINDEX3D(cctkGH,i,j,k) ;
         subblock_index  =  i - min[0] + (max[0] - min[0]) * (j - min[1] + (max[1]-min[1]) * (k - min[2])) ;
        
        /* Assign local copies of grid functions */
        gxxL = gxx[index];
        gxyL = gxy[index];
        gxzL = gxz[index];
        gyyL = gyy[index];
        gyzL = gyz[index];
        gzzL = gzz[index];
        kxxL = kxx[index];
        kxyL = kxy[index];
        kxzL = kxz[index];
        kyyL = kyy[index];
        kyzL = kyz[index];
        kzzL = kzz[index];
        xL = x[index];
        yL = y[index];
        zL = z[index];
        
        /* Assign local copies of subblock grid functions */
        
        /* Include user supplied include files */
        
        /* Precompute derivatives (new style) */
        PDstandard2nd1gxx = PDstandard2nd1(gxx, i, j, k);
        PDstandard2nd2gxx = PDstandard2nd2(gxx, i, j, k);
        PDstandard2nd3gxx = PDstandard2nd3(gxx, i, j, k);
        PDstandard2nd22gxx = PDstandard2nd22(gxx, i, j, k);
        PDstandard2nd33gxx = PDstandard2nd33(gxx, i, j, k);
        PDstandard2nd23gxx = PDstandard2nd23(gxx, i, j, k);
        PDstandard2nd1gxy = PDstandard2nd1(gxy, i, j, k);
        PDstandard2nd2gxy = PDstandard2nd2(gxy, i, j, k);
        PDstandard2nd3gxy = PDstandard2nd3(gxy, i, j, k);
        PDstandard2nd33gxy = PDstandard2nd33(gxy, i, j, k);
        PDstandard2nd12gxy = PDstandard2nd12(gxy, i, j, k);
        PDstandard2nd13gxy = PDstandard2nd13(gxy, i, j, k);
        PDstandard2nd23gxy = PDstandard2nd23(gxy, i, j, k);
        PDstandard2nd1gxz = PDstandard2nd1(gxz, i, j, k);
        PDstandard2nd2gxz = PDstandard2nd2(gxz, i, j, k);
        PDstandard2nd3gxz = PDstandard2nd3(gxz, i, j, k);
        PDstandard2nd22gxz = PDstandard2nd22(gxz, i, j, k);
        PDstandard2nd12gxz = PDstandard2nd12(gxz, i, j, k);
        PDstandard2nd13gxz = PDstandard2nd13(gxz, i, j, k);
        PDstandard2nd23gxz = PDstandard2nd23(gxz, i, j, k);
        PDstandard2nd1gyy = PDstandard2nd1(gyy, i, j, k);
        PDstandard2nd2gyy = PDstandard2nd2(gyy, i, j, k);
        PDstandard2nd3gyy = PDstandard2nd3(gyy, i, j, k);
        PDstandard2nd11gyy = PDstandard2nd11(gyy, i, j, k);
        PDstandard2nd33gyy = PDstandard2nd33(gyy, i, j, k);
        PDstandard2nd13gyy = PDstandard2nd13(gyy, i, j, k);
        PDstandard2nd1gyz = PDstandard2nd1(gyz, i, j, k);
        PDstandard2nd2gyz = PDstandard2nd2(gyz, i, j, k);
        PDstandard2nd3gyz = PDstandard2nd3(gyz, i, j, k);
        PDstandard2nd11gyz = PDstandard2nd11(gyz, i, j, k);
        PDstandard2nd12gyz = PDstandard2nd12(gyz, i, j, k);
        PDstandard2nd13gyz = PDstandard2nd13(gyz, i, j, k);
        PDstandard2nd23gyz = PDstandard2nd23(gyz, i, j, k);
        PDstandard2nd1gzz = PDstandard2nd1(gzz, i, j, k);
        PDstandard2nd2gzz = PDstandard2nd2(gzz, i, j, k);
        PDstandard2nd3gzz = PDstandard2nd3(gzz, i, j, k);
        PDstandard2nd11gzz = PDstandard2nd11(gzz, i, j, k);
        PDstandard2nd22gzz = PDstandard2nd22(gzz, i, j, k);
        PDstandard2nd12gzz = PDstandard2nd12(gzz, i, j, k);
        PDstandard2nd2kxx = PDstandard2nd2(kxx, i, j, k);
        PDstandard2nd3kxx = PDstandard2nd3(kxx, i, j, k);
        PDstandard2nd1kxy = PDstandard2nd1(kxy, i, j, k);
        PDstandard2nd2kxy = PDstandard2nd2(kxy, i, j, k);
        PDstandard2nd3kxy = PDstandard2nd3(kxy, i, j, k);
        PDstandard2nd1kxz = PDstandard2nd1(kxz, i, j, k);
        PDstandard2nd2kxz = PDstandard2nd2(kxz, i, j, k);
        PDstandard2nd3kxz = PDstandard2nd3(kxz, i, j, k);
        PDstandard2nd1kyy = PDstandard2nd1(kyy, i, j, k);
        PDstandard2nd3kyy = PDstandard2nd3(kyy, i, j, k);
        PDstandard2nd1kyz = PDstandard2nd1(kyz, i, j, k);
        PDstandard2nd2kyz = PDstandard2nd2(kyz, i, j, k);
        PDstandard2nd3kyz = PDstandard2nd3(kyz, i, j, k);
        PDstandard2nd1kzz = PDstandard2nd1(kzz, i, j, k);
        PDstandard2nd2kzz = PDstandard2nd2(kzz, i, j, k);
        
        /* Precompute derivatives (old style) */
        
        /* Calculate temporaries and grid functions */
        detg  =  2*gxyL*gxzL*gyzL + gzzL*(gxxL*gyyL - SQR(gxyL)) - gyyL*SQR(gxzL) - gxxL*SQR(gyzL);
        
        invdetg  =  INV(detg);
        
        gInv11  =  invdetg*(gyyL*gzzL - SQR(gyzL));
        
        gInv12  =  (gxzL*gyzL - gxyL*gzzL)*invdetg;
        
        gInv13  =  (-(gxzL*gyyL) + gxyL*gyzL)*invdetg;
        
        gInv21  =  (gxzL*gyzL - gxyL*gzzL)*invdetg;
        
        gInv22  =  invdetg*(gxxL*gzzL - SQR(gxzL));
        
        gInv23  =  (gxyL*gxzL - gxxL*gyzL)*invdetg;
        
        gInv31  =  (-(gxzL*gyyL) + gxyL*gyzL)*invdetg;
        
        gInv32  =  (gxyL*gxzL - gxxL*gyzL)*invdetg;
        
        gInv33  =  invdetg*(gxxL*gyyL - SQR(gxyL));
        
        gamma111  =  khalf*(gInv11*PDstandard2nd1gxx + 2*(gInv12*PDstandard2nd1gxy + gInv13*PDstandard2nd1gxz) - 
              gInv12*PDstandard2nd2gxx - gInv13*PDstandard2nd3gxx);
        
        gamma211  =  khalf*(gInv21*PDstandard2nd1gxx + 2*(gInv22*PDstandard2nd1gxy + gInv23*PDstandard2nd1gxz) - 
              gInv22*PDstandard2nd2gxx - gInv23*PDstandard2nd3gxx);
        
        gamma311  =  khalf*(gInv31*PDstandard2nd1gxx + 2*(gInv32*PDstandard2nd1gxy + gInv33*PDstandard2nd1gxz) - 
              gInv32*PDstandard2nd2gxx - gInv33*PDstandard2nd3gxx);
        
        gamma121  =  khalf*(gInv12*PDstandard2nd1gyy + gInv11*PDstandard2nd2gxx + 
              gInv13*(PDstandard2nd1gyz + PDstandard2nd2gxz - PDstandard2nd3gxy));
        
        gamma221  =  khalf*(gInv22*PDstandard2nd1gyy + gInv21*PDstandard2nd2gxx + 
              gInv23*(PDstandard2nd1gyz + PDstandard2nd2gxz - PDstandard2nd3gxy));
        
        gamma321  =  khalf*(gInv32*PDstandard2nd1gyy + gInv31*PDstandard2nd2gxx + 
              gInv33*(PDstandard2nd1gyz + PDstandard2nd2gxz - PDstandard2nd3gxy));
        
        gamma131  =  khalf*(gInv13*PDstandard2nd1gzz + gInv11*PDstandard2nd3gxx + 
              gInv12*(PDstandard2nd1gyz - PDstandard2nd2gxz + PDstandard2nd3gxy));
        
        gamma231  =  khalf*(gInv23*PDstandard2nd1gzz + gInv21*PDstandard2nd3gxx + 
              gInv22*(PDstandard2nd1gyz - PDstandard2nd2gxz + PDstandard2nd3gxy));
        
        gamma331  =  khalf*(gInv33*PDstandard2nd1gzz + gInv31*PDstandard2nd3gxx + 
              gInv32*(PDstandard2nd1gyz - PDstandard2nd2gxz + PDstandard2nd3gxy));
        
        gamma122  =  khalf*(gInv11*(-PDstandard2nd1gyy + 2*PDstandard2nd2gxy) + gInv12*PDstandard2nd2gyy + 
              gInv13*(2*PDstandard2nd2gyz - PDstandard2nd3gyy));
        
        gamma222  =  khalf*(gInv21*(-PDstandard2nd1gyy + 2*PDstandard2nd2gxy) + gInv22*PDstandard2nd2gyy + 
              gInv23*(2*PDstandard2nd2gyz - PDstandard2nd3gyy));
        
        gamma322  =  khalf*(gInv31*(-PDstandard2nd1gyy + 2*PDstandard2nd2gxy) + gInv32*PDstandard2nd2gyy + 
              gInv33*(2*PDstandard2nd2gyz - PDstandard2nd3gyy));
        
        gamma132  =  khalf*(gInv13*PDstandard2nd2gzz + gInv11*(-PDstandard2nd1gyz + PDstandard2nd2gxz + PDstandard2nd3gxy) + 
              gInv12*PDstandard2nd3gyy);
        
        gamma232  =  khalf*(gInv23*PDstandard2nd2gzz + gInv21*(-PDstandard2nd1gyz + PDstandard2nd2gxz + PDstandard2nd3gxy) + 
              gInv22*PDstandard2nd3gyy);
        
        gamma332  =  khalf*(gInv33*PDstandard2nd2gzz + gInv31*(-PDstandard2nd1gyz + PDstandard2nd2gxz + PDstandard2nd3gxy) + 
              gInv32*PDstandard2nd3gyy);
        
        gamma133  =  khalf*(-(gInv11*PDstandard2nd1gzz) - gInv12*PDstandard2nd2gzz + 2*gInv11*PDstandard2nd3gxz + 
              2*gInv12*PDstandard2nd3gyz + gInv13*PDstandard2nd3gzz);
        
        gamma233  =  khalf*(-(gInv21*PDstandard2nd1gzz) - gInv22*PDstandard2nd2gzz + 2*gInv21*PDstandard2nd3gxz + 
              2*gInv22*PDstandard2nd3gyz + gInv23*PDstandard2nd3gzz);
        
        gamma333  =  khalf*(-(gInv31*PDstandard2nd1gzz) - gInv32*PDstandard2nd2gzz + 2*gInv31*PDstandard2nd3gxz + 
              2*gInv32*PDstandard2nd3gyz + gInv33*PDstandard2nd3gzz);
        
        xmoved  =  xL - xorig;
        
        ymoved  =  yL - yorig;
        
        zmoved  =  zL - zorig;
        
        va1  =  -ymoved;
        
        va2  =  offset + xmoved;
        
        va3  =  0;
        
        vb1  =  offset + xmoved;
        
        vb2  =  ymoved;
        
        vb3  =  zmoved;
        
        vc1  =  ((-(gInv13*va2) + gInv12*va3)*vb1 + (gInv13*va1 - gInv11*va3)*vb2 + (-(gInv12*va1) + gInv11*va2)*vb3)*
            pow(detg,0.5);
        
        vc2  =  ((-(gInv23*va2) + gInv22*va3)*vb1 + (gInv23*va1 - gInv21*va3)*vb2 + (-(gInv22*va1) + gInv21*va2)*vb3)*
            pow(detg,0.5);
        
        vc3  =  ((-(gInv33*va2) + gInv32*va3)*vb1 + (gInv33*va1 - gInv31*va3)*vb2 + (-(gInv32*va1) + gInv31*va2)*vb3)*
            pow(detg,0.5);
        
        omega11  =  2*(gyzL*va2*va3 + va1*(gxyL*va2 + gxzL*va3)) + gxxL*SQR(va1) + gyyL*SQR(va2) + gzzL*SQR(va3);
        
        va1  =  va1*pow(omega11,-khalf);
        
        va2  =  va2*pow(omega11,-khalf);
        
        va3  =  va3*pow(omega11,-khalf);
        
        omega12  =  (gxxL*va1 + gxyL*va2 + gxzL*va3)*vb1 + (gxyL*va1 + gyyL*va2 + gyzL*va3)*vb2 + 
            (gxzL*va1 + gyzL*va2 + gzzL*va3)*vb3;
        
        vb1  =  -(omega12*va1) + vb1;
        
        vb2  =  -(omega12*va2) + vb2;
        
        vb3  =  -(omega12*va3) + vb3;
        
        omega22  =  2*(gyzL*vb2*vb3 + vb1*(gxyL*vb2 + gxzL*vb3)) + gxxL*SQR(vb1) + gyyL*SQR(vb2) + gzzL*SQR(vb3);
        
        vb1  =  vb1*pow(omega22,-khalf);
        
        vb2  =  vb2*pow(omega22,-khalf);
        
        vb3  =  vb3*pow(omega22,-khalf);
        
        omega13  =  (gxxL*va1 + gxyL*va2 + gxzL*va3)*vc1 + (gxyL*va1 + gyyL*va2 + gyzL*va3)*vc2 + 
            (gxzL*va1 + gyzL*va2 + gzzL*va3)*vc3;
        
        omega23  =  (gxxL*vb1 + gxyL*vb2 + gxzL*vb3)*vc1 + (gxyL*vb1 + gyyL*vb2 + gyzL*vb3)*vc2 + 
            (gxzL*vb1 + gyzL*vb2 + gzzL*vb3)*vc3;
        
        vc1  =  -(omega13*va1) - omega23*vb1 + vc1;
        
        vc2  =  -(omega13*va2) - omega23*vb2 + vc2;
        
        vc3  =  -(omega13*va3) - omega23*vb3 + vc3;
        
        omega33  =  2*(gyzL*vc2*vc3 + vc1*(gxyL*vc2 + gxzL*vc3)) + gxxL*SQR(vc1) + gyyL*SQR(vc2) + gzzL*SQR(vc3);
        
        vc1  =  vc1*pow(omega33,-khalf);
        
        vc2  =  vc2*pow(omega33,-khalf);
        
        vc3  =  vc3*pow(omega33,-khalf);
        
        isqrt2  =  0.7071067811865475244;
        
        n1  =  -(isqrt2*vb1);
        
        n2  =  -(isqrt2*vb2);
        
        n3  =  -(isqrt2*vb3);
        
        rm1  =  isqrt2*vc1;
        
        rm2  =  isqrt2*vc2;
        
        rm3  =  isqrt2*vc3;
        
        im1  =  isqrt2*va1;
        
        im2  =  isqrt2*va2;
        
        im3  =  isqrt2*va3;
        
        rmbar1  =  isqrt2*vc1;
        
        rmbar2  =  isqrt2*vc2;
        
        rmbar3  =  isqrt2*vc3;
        
        imbar1  =  -(isqrt2*va1);
        
        imbar2  =  -(isqrt2*va2);
        
        imbar3  =  -(isqrt2*va3);
        
        n0  =  isqrt2;
        
        R1111  =  0;
        
        R1112  =  0;
        
        R1113  =  0;
        
        R1121  =  0;
        
        R1122  =  0;
        
        R1123  =  0;
        
        R1131  =  0;
        
        R1132  =  0;
        
        R1133  =  0;
        
        R1211  =  0;
        
        R1212  =  khalf*((4*gamma121*gamma221 - 2*gamma111*gamma222)*gxyL + (4*gamma121*gamma321 - 2*gamma111*gamma322)*gxzL + 
              (4*gamma221*gamma321 - 2*gamma211*gamma322)*gyzL - 
              2*(gamma111*gamma122*gxxL + gamma122*gamma211*gxyL + gamma122*gamma311*gxzL + gamma211*gamma222*gyyL + 
                 gamma222*gamma311*gyzL + gamma311*gamma322*gzzL) - PDstandard2nd11gyy + 2*PDstandard2nd12gxy - 
              PDstandard2nd22gxx + 2*gxxL*SQR(gamma121) + 2*gyyL*SQR(gamma221) + 2*gzzL*SQR(gamma321));
        
        R1213  =  (gamma121*gamma131 - gamma111*gamma132)*gxxL + 
            (-(gamma132*gamma211) + gamma131*gamma221 + gamma121*gamma231 - gamma111*gamma232)*gxyL + 
            (-(gamma132*gamma311) + gamma131*gamma321 + gamma121*gamma331 - gamma111*gamma332)*gxzL + 
            (gamma221*gamma231 - gamma211*gamma232)*gyyL + 
            (-(gamma232*gamma311) + gamma231*gamma321 + gamma221*gamma331 - gamma211*gamma332)*gyzL + 
            (gamma321*gamma331 - gamma311*gamma332)*gzzL + 
            khalf*(-PDstandard2nd11gyz + PDstandard2nd12gxz + PDstandard2nd13gxy - PDstandard2nd23gxx);
        
        R1221  =  (-2*gamma121*gamma221 + gamma111*gamma222)*gxyL + (-2*gamma121*gamma321 + gamma111*gamma322)*gxzL + 
            gamma122*(gamma111*gxxL + gamma211*gxyL + gamma311*gxzL) + (-2*gamma221*gamma321 + gamma211*gamma322)*gyzL + 
            gamma222*(gamma211*gyyL + gamma311*gyzL) - PDstandard2nd12gxy + khalf*(PDstandard2nd11gyy + PDstandard2nd22gxx) - 
            gxxL*SQR(gamma121) - gyyL*SQR(gamma221) + gzzL*(gamma311*gamma322 - SQR(gamma321));
        
        R1222  =  0;
        
        R1223  =  (gamma122*gamma131 - gamma121*gamma132)*gxxL + 
            (-(gamma132*gamma221) + gamma131*gamma222 + gamma122*gamma231 - gamma121*gamma232)*gxyL + 
            (-(gamma132*gamma321) + gamma131*gamma322 + gamma122*gamma331 - gamma121*gamma332)*gxzL + 
            (gamma222*gamma231 - gamma221*gamma232)*gyyL + 
            (-(gamma232*gamma321) + gamma231*gamma322 + gamma222*gamma331 - gamma221*gamma332)*gyzL + 
            (gamma322*gamma331 - gamma321*gamma332)*gzzL + 
            khalf*(-PDstandard2nd12gyz + PDstandard2nd13gyy + PDstandard2nd22gxz - PDstandard2nd23gxy);
        
        R1231  =  (-(gamma121*gamma131) + gamma111*gamma132)*gxxL + 
            (gamma132*gamma211 - gamma131*gamma221 - gamma121*gamma231 + gamma111*gamma232)*gxyL + 
            (gamma132*gamma311 - gamma131*gamma321 - gamma121*gamma331 + gamma111*gamma332)*gxzL + 
            (-(gamma221*gamma231) + gamma211*gamma232)*gyyL + 
            (gamma232*gamma311 - gamma231*gamma321 - gamma221*gamma331 + gamma211*gamma332)*gyzL + 
            (-(gamma321*gamma331) + gamma311*gamma332)*gzzL + 
            khalf*(PDstandard2nd11gyz - PDstandard2nd12gxz - PDstandard2nd13gxy + PDstandard2nd23gxx);
        
        R1232  =  (-(gamma122*gamma131) + gamma121*gamma132)*gxxL + 
            (gamma132*gamma221 - gamma131*gamma222 - gamma122*gamma231 + gamma121*gamma232)*gxyL + 
            (gamma132*gamma321 - gamma131*gamma322 - gamma122*gamma331 + gamma121*gamma332)*gxzL + 
            (-(gamma222*gamma231) + gamma221*gamma232)*gyyL + 
            (gamma232*gamma321 - gamma231*gamma322 - gamma222*gamma331 + gamma221*gamma332)*gyzL + 
            (-(gamma322*gamma331) + gamma321*gamma332)*gzzL + 
            khalf*(PDstandard2nd12gyz - PDstandard2nd13gyy - PDstandard2nd22gxz + PDstandard2nd23gxy);
        
        R1233  =  0;
        
        R1311  =  0;
        
        R1312  =  (gamma121*gamma131 - gamma111*gamma132)*gxxL + 
            (-(gamma132*gamma211) + gamma131*gamma221 + gamma121*gamma231 - gamma111*gamma232)*gxyL + 
            (-(gamma132*gamma311) + gamma131*gamma321 + gamma121*gamma331 - gamma111*gamma332)*gxzL + 
            (gamma221*gamma231 - gamma211*gamma232)*gyyL + 
            (-(gamma232*gamma311) + gamma231*gamma321 + gamma221*gamma331 - gamma211*gamma332)*gyzL + 
            (gamma321*gamma331 - gamma311*gamma332)*gzzL + 
            khalf*(-PDstandard2nd11gyz + PDstandard2nd12gxz + PDstandard2nd13gxy - PDstandard2nd23gxx);
        
        R1313  =  khalf*((4*gamma131*gamma231 - 2*gamma111*gamma233)*gxyL + (4*gamma131*gamma331 - 2*gamma111*gamma333)*gxzL + 
              (4*gamma231*gamma331 - 2*gamma211*gamma333)*gyzL - 
              2*(gamma111*gamma133*gxxL + gamma133*gamma211*gxyL + gamma133*gamma311*gxzL + gamma211*gamma233*gyyL + 
                 gamma233*gamma311*gyzL + gamma311*gamma333*gzzL) - PDstandard2nd11gzz + 2*PDstandard2nd13gxz - 
              PDstandard2nd33gxx + 2*gxxL*SQR(gamma131) + 2*gyyL*SQR(gamma231) + 2*gzzL*SQR(gamma331));
        
        R1321  =  (-(gamma121*gamma131) + gamma111*gamma132)*gxxL + 
            (gamma132*gamma211 - gamma131*gamma221 - gamma121*gamma231 + gamma111*gamma232)*gxyL + 
            (gamma132*gamma311 - gamma131*gamma321 - gamma121*gamma331 + gamma111*gamma332)*gxzL + 
            (-(gamma221*gamma231) + gamma211*gamma232)*gyyL + 
            (gamma232*gamma311 - gamma231*gamma321 - gamma221*gamma331 + gamma211*gamma332)*gyzL + 
            (-(gamma321*gamma331) + gamma311*gamma332)*gzzL + 
            khalf*(PDstandard2nd11gyz - PDstandard2nd12gxz - PDstandard2nd13gxy + PDstandard2nd23gxx);
        
        R1322  =  0;
        
        R1323  =  (gamma131*gamma132 - gamma121*gamma133)*gxxL + 
            (-(gamma133*gamma221) + gamma132*gamma231 + gamma131*gamma232 - gamma121*gamma233)*gxyL + 
            (-(gamma133*gamma321) + gamma132*gamma331 + gamma131*gamma332 - gamma121*gamma333)*gxzL + 
            (gamma231*gamma232 - gamma221*gamma233)*gyyL + 
            (-(gamma233*gamma321) + gamma232*gamma331 + gamma231*gamma332 - gamma221*gamma333)*gyzL + 
            (gamma331*gamma332 - gamma321*gamma333)*gzzL + 
            khalf*(-PDstandard2nd12gzz + PDstandard2nd13gyz + PDstandard2nd23gxz - PDstandard2nd33gxy);
        
        R1331  =  (-2*gamma131*gamma231 + gamma111*gamma233)*gxyL + (-2*gamma131*gamma331 + gamma111*gamma333)*gxzL + 
            gamma133*(gamma111*gxxL + gamma211*gxyL + gamma311*gxzL) + (-2*gamma231*gamma331 + gamma211*gamma333)*gyzL + 
            gamma233*(gamma211*gyyL + gamma311*gyzL) - PDstandard2nd13gxz + khalf*(PDstandard2nd11gzz + PDstandard2nd33gxx) - 
            gxxL*SQR(gamma131) - gyyL*SQR(gamma231) + gzzL*(gamma311*gamma333 - SQR(gamma331));
        
        R1332  =  (-(gamma131*gamma132) + gamma121*gamma133)*gxxL + 
            (gamma133*gamma221 - gamma132*gamma231 - gamma131*gamma232 + gamma121*gamma233)*gxyL + 
            (gamma133*gamma321 - gamma132*gamma331 - gamma131*gamma332 + gamma121*gamma333)*gxzL + 
            (-(gamma231*gamma232) + gamma221*gamma233)*gyyL + 
            (gamma233*gamma321 - gamma232*gamma331 - gamma231*gamma332 + gamma221*gamma333)*gyzL + 
            (-(gamma331*gamma332) + gamma321*gamma333)*gzzL + 
            khalf*(PDstandard2nd12gzz - PDstandard2nd13gyz - PDstandard2nd23gxz + PDstandard2nd33gxy);
        
        R1333  =  0;
        
        R2111  =  0;
        
        R2112  =  (-2*gamma121*gamma221 + gamma111*gamma222)*gxyL + (-2*gamma121*gamma321 + gamma111*gamma322)*gxzL + 
            gamma122*(gamma111*gxxL + gamma211*gxyL + gamma311*gxzL) + (-2*gamma221*gamma321 + gamma211*gamma322)*gyzL + 
            gamma222*(gamma211*gyyL + gamma311*gyzL) - PDstandard2nd12gxy + khalf*(PDstandard2nd11gyy + PDstandard2nd22gxx) - 
            gxxL*SQR(gamma121) - gyyL*SQR(gamma221) + gzzL*(gamma311*gamma322 - SQR(gamma321));
        
        R2113  =  (-(gamma121*gamma131) + gamma111*gamma132)*gxxL + 
            (gamma132*gamma211 - gamma131*gamma221 - gamma121*gamma231 + gamma111*gamma232)*gxyL + 
            (gamma132*gamma311 - gamma131*gamma321 - gamma121*gamma331 + gamma111*gamma332)*gxzL + 
            (-(gamma221*gamma231) + gamma211*gamma232)*gyyL + 
            (gamma232*gamma311 - gamma231*gamma321 - gamma221*gamma331 + gamma211*gamma332)*gyzL + 
            (-(gamma321*gamma331) + gamma311*gamma332)*gzzL + 
            khalf*(PDstandard2nd11gyz - PDstandard2nd12gxz - PDstandard2nd13gxy + PDstandard2nd23gxx);
        
        R2121  =  khalf*((4*gamma121*gamma221 - 2*gamma111*gamma222)*gxyL + (4*gamma121*gamma321 - 2*gamma111*gamma322)*gxzL + 
              (4*gamma221*gamma321 - 2*gamma211*gamma322)*gyzL - 
              2*(gamma111*gamma122*gxxL + gamma122*gamma211*gxyL + gamma122*gamma311*gxzL + gamma211*gamma222*gyyL + 
                 gamma222*gamma311*gyzL + gamma311*gamma322*gzzL) - PDstandard2nd11gyy + 2*PDstandard2nd12gxy - 
              PDstandard2nd22gxx + 2*gxxL*SQR(gamma121) + 2*gyyL*SQR(gamma221) + 2*gzzL*SQR(gamma321));
        
        R2122  =  0;
        
        R2123  =  (-(gamma122*gamma131) + gamma121*gamma132)*gxxL + 
            (gamma132*gamma221 - gamma131*gamma222 - gamma122*gamma231 + gamma121*gamma232)*gxyL + 
            (gamma132*gamma321 - gamma131*gamma322 - gamma122*gamma331 + gamma121*gamma332)*gxzL + 
            (-(gamma222*gamma231) + gamma221*gamma232)*gyyL + 
            (gamma232*gamma321 - gamma231*gamma322 - gamma222*gamma331 + gamma221*gamma332)*gyzL + 
            (-(gamma322*gamma331) + gamma321*gamma332)*gzzL + 
            khalf*(PDstandard2nd12gyz - PDstandard2nd13gyy - PDstandard2nd22gxz + PDstandard2nd23gxy);
        
        R2131  =  (gamma121*gamma131 - gamma111*gamma132)*gxxL + 
            (-(gamma132*gamma211) + gamma131*gamma221 + gamma121*gamma231 - gamma111*gamma232)*gxyL + 
            (-(gamma132*gamma311) + gamma131*gamma321 + gamma121*gamma331 - gamma111*gamma332)*gxzL + 
            (gamma221*gamma231 - gamma211*gamma232)*gyyL + 
            (-(gamma232*gamma311) + gamma231*gamma321 + gamma221*gamma331 - gamma211*gamma332)*gyzL + 
            (gamma321*gamma331 - gamma311*gamma332)*gzzL + 
            khalf*(-PDstandard2nd11gyz + PDstandard2nd12gxz + PDstandard2nd13gxy - PDstandard2nd23gxx);
        
        R2132  =  (gamma122*gamma131 - gamma121*gamma132)*gxxL + 
            (-(gamma132*gamma221) + gamma131*gamma222 + gamma122*gamma231 - gamma121*gamma232)*gxyL + 
            (-(gamma132*gamma321) + gamma131*gamma322 + gamma122*gamma331 - gamma121*gamma332)*gxzL + 
            (gamma222*gamma231 - gamma221*gamma232)*gyyL + 
            (-(gamma232*gamma321) + gamma231*gamma322 + gamma222*gamma331 - gamma221*gamma332)*gyzL + 
            (gamma322*gamma331 - gamma321*gamma332)*gzzL + 
            khalf*(-PDstandard2nd12gyz + PDstandard2nd13gyy + PDstandard2nd22gxz - PDstandard2nd23gxy);
        
        R2133  =  0;
        
        R2211  =  0;
        
        R2212  =  0;
        
        R2213  =  0;
        
        R2221  =  0;
        
        R2222  =  0;
        
        R2223  =  0;
        
        R2231  =  0;
        
        R2232  =  0;
        
        R2233  =  0;
        
        R2311  =  0;
        
        R2312  =  (gamma122*gamma131 - gamma121*gamma132)*gxxL + 
            (-(gamma132*gamma221) + gamma131*gamma222 + gamma122*gamma231 - gamma121*gamma232)*gxyL + 
            (-(gamma132*gamma321) + gamma131*gamma322 + gamma122*gamma331 - gamma121*gamma332)*gxzL + 
            (gamma222*gamma231 - gamma221*gamma232)*gyyL + 
            (-(gamma232*gamma321) + gamma231*gamma322 + gamma222*gamma331 - gamma221*gamma332)*gyzL + 
            (gamma322*gamma331 - gamma321*gamma332)*gzzL + 
            khalf*(-PDstandard2nd12gyz + PDstandard2nd13gyy + PDstandard2nd22gxz - PDstandard2nd23gxy);
        
        R2313  =  (gamma131*gamma132 - gamma121*gamma133)*gxxL + 
            (-(gamma133*gamma221) + gamma132*gamma231 + gamma131*gamma232 - gamma121*gamma233)*gxyL + 
            (-(gamma133*gamma321) + gamma132*gamma331 + gamma131*gamma332 - gamma121*gamma333)*gxzL + 
            (gamma231*gamma232 - gamma221*gamma233)*gyyL + 
            (-(gamma233*gamma321) + gamma232*gamma331 + gamma231*gamma332 - gamma221*gamma333)*gyzL + 
            (gamma331*gamma332 - gamma321*gamma333)*gzzL + 
            khalf*(-PDstandard2nd12gzz + PDstandard2nd13gyz + PDstandard2nd23gxz - PDstandard2nd33gxy);
        
        R2321  =  (-(gamma122*gamma131) + gamma121*gamma132)*gxxL + 
            (gamma132*gamma221 - gamma131*gamma222 - gamma122*gamma231 + gamma121*gamma232)*gxyL + 
            (gamma132*gamma321 - gamma131*gamma322 - gamma122*gamma331 + gamma121*gamma332)*gxzL + 
            (-(gamma222*gamma231) + gamma221*gamma232)*gyyL + 
            (gamma232*gamma321 - gamma231*gamma322 - gamma222*gamma331 + gamma221*gamma332)*gyzL + 
            (-(gamma322*gamma331) + gamma321*gamma332)*gzzL + 
            khalf*(PDstandard2nd12gyz - PDstandard2nd13gyy - PDstandard2nd22gxz + PDstandard2nd23gxy);
        
        R2322  =  0;
        
        R2323  =  khalf*((4*gamma132*gamma232 - 2*gamma122*gamma233)*gxyL + (4*gamma132*gamma332 - 2*gamma122*gamma333)*gxzL + 
              (4*gamma232*gamma332 - 2*gamma222*gamma333)*gyzL - 
              2*(gamma122*gamma133*gxxL + gamma133*gamma222*gxyL + gamma133*gamma322*gxzL + gamma222*gamma233*gyyL + 
                 gamma233*gamma322*gyzL + gamma322*gamma333*gzzL) - PDstandard2nd22gzz + 2*PDstandard2nd23gyz - 
              PDstandard2nd33gyy + 2*gxxL*SQR(gamma132) + 2*gyyL*SQR(gamma232) + 2*gzzL*SQR(gamma332));
        
        R2331  =  (-(gamma131*gamma132) + gamma121*gamma133)*gxxL + 
            (gamma133*gamma221 - gamma132*gamma231 - gamma131*gamma232 + gamma121*gamma233)*gxyL + 
            (gamma133*gamma321 - gamma132*gamma331 - gamma131*gamma332 + gamma121*gamma333)*gxzL + 
            (-(gamma231*gamma232) + gamma221*gamma233)*gyyL + 
            (gamma233*gamma321 - gamma232*gamma331 - gamma231*gamma332 + gamma221*gamma333)*gyzL + 
            (-(gamma331*gamma332) + gamma321*gamma333)*gzzL + 
            khalf*(PDstandard2nd12gzz - PDstandard2nd13gyz - PDstandard2nd23gxz + PDstandard2nd33gxy);
        
        R2332  =  (-2*gamma132*gamma232 + gamma122*gamma233)*gxyL + (-2*gamma132*gamma332 + gamma122*gamma333)*gxzL + 
            gamma133*(gamma122*gxxL + gamma222*gxyL + gamma322*gxzL) + (-2*gamma232*gamma332 + gamma222*gamma333)*gyzL + 
            gamma233*(gamma222*gyyL + gamma322*gyzL) - PDstandard2nd23gyz + khalf*(PDstandard2nd22gzz + PDstandard2nd33gyy) - 
            gxxL*SQR(gamma132) - gyyL*SQR(gamma232) + gzzL*(gamma322*gamma333 - SQR(gamma332));
        
        R2333  =  0;
        
        R3111  =  0;
        
        R3112  =  (-(gamma121*gamma131) + gamma111*gamma132)*gxxL + 
            (gamma132*gamma211 - gamma131*gamma221 - gamma121*gamma231 + gamma111*gamma232)*gxyL + 
            (gamma132*gamma311 - gamma131*gamma321 - gamma121*gamma331 + gamma111*gamma332)*gxzL + 
            (-(gamma221*gamma231) + gamma211*gamma232)*gyyL + 
            (gamma232*gamma311 - gamma231*gamma321 - gamma221*gamma331 + gamma211*gamma332)*gyzL + 
            (-(gamma321*gamma331) + gamma311*gamma332)*gzzL + 
            khalf*(PDstandard2nd11gyz - PDstandard2nd12gxz - PDstandard2nd13gxy + PDstandard2nd23gxx);
        
        R3113  =  (-2*gamma131*gamma231 + gamma111*gamma233)*gxyL + (-2*gamma131*gamma331 + gamma111*gamma333)*gxzL + 
            gamma133*(gamma111*gxxL + gamma211*gxyL + gamma311*gxzL) + (-2*gamma231*gamma331 + gamma211*gamma333)*gyzL + 
            gamma233*(gamma211*gyyL + gamma311*gyzL) - PDstandard2nd13gxz + khalf*(PDstandard2nd11gzz + PDstandard2nd33gxx) - 
            gxxL*SQR(gamma131) - gyyL*SQR(gamma231) + gzzL*(gamma311*gamma333 - SQR(gamma331));
        
        R3121  =  (gamma121*gamma131 - gamma111*gamma132)*gxxL + 
            (-(gamma132*gamma211) + gamma131*gamma221 + gamma121*gamma231 - gamma111*gamma232)*gxyL + 
            (-(gamma132*gamma311) + gamma131*gamma321 + gamma121*gamma331 - gamma111*gamma332)*gxzL + 
            (gamma221*gamma231 - gamma211*gamma232)*gyyL + 
            (-(gamma232*gamma311) + gamma231*gamma321 + gamma221*gamma331 - gamma211*gamma332)*gyzL + 
            (gamma321*gamma331 - gamma311*gamma332)*gzzL + 
            khalf*(-PDstandard2nd11gyz + PDstandard2nd12gxz + PDstandard2nd13gxy - PDstandard2nd23gxx);
        
        R3122  =  0;
        
        R3123  =  (-(gamma131*gamma132) + gamma121*gamma133)*gxxL + 
            (gamma133*gamma221 - gamma132*gamma231 - gamma131*gamma232 + gamma121*gamma233)*gxyL + 
            (gamma133*gamma321 - gamma132*gamma331 - gamma131*gamma332 + gamma121*gamma333)*gxzL + 
            (-(gamma231*gamma232) + gamma221*gamma233)*gyyL + 
            (gamma233*gamma321 - gamma232*gamma331 - gamma231*gamma332 + gamma221*gamma333)*gyzL + 
            (-(gamma331*gamma332) + gamma321*gamma333)*gzzL + 
            khalf*(PDstandard2nd12gzz - PDstandard2nd13gyz - PDstandard2nd23gxz + PDstandard2nd33gxy);
        
        R3131  =  khalf*((4*gamma131*gamma231 - 2*gamma111*gamma233)*gxyL + (4*gamma131*gamma331 - 2*gamma111*gamma333)*gxzL + 
              (4*gamma231*gamma331 - 2*gamma211*gamma333)*gyzL - 
              2*(gamma111*gamma133*gxxL + gamma133*gamma211*gxyL + gamma133*gamma311*gxzL + gamma211*gamma233*gyyL + 
                 gamma233*gamma311*gyzL + gamma311*gamma333*gzzL) - PDstandard2nd11gzz + 2*PDstandard2nd13gxz - 
              PDstandard2nd33gxx + 2*gxxL*SQR(gamma131) + 2*gyyL*SQR(gamma231) + 2*gzzL*SQR(gamma331));
        
        R3132  =  (gamma131*gamma132 - gamma121*gamma133)*gxxL + 
            (-(gamma133*gamma221) + gamma132*gamma231 + gamma131*gamma232 - gamma121*gamma233)*gxyL + 
            (-(gamma133*gamma321) + gamma132*gamma331 + gamma131*gamma332 - gamma121*gamma333)*gxzL + 
            (gamma231*gamma232 - gamma221*gamma233)*gyyL + 
            (-(gamma233*gamma321) + gamma232*gamma331 + gamma231*gamma332 - gamma221*gamma333)*gyzL + 
            (gamma331*gamma332 - gamma321*gamma333)*gzzL + 
            khalf*(-PDstandard2nd12gzz + PDstandard2nd13gyz + PDstandard2nd23gxz - PDstandard2nd33gxy);
        
        R3133  =  0;
        
        R3211  =  0;
        
        R3212  =  (-(gamma122*gamma131) + gamma121*gamma132)*gxxL + 
            (gamma132*gamma221 - gamma131*gamma222 - gamma122*gamma231 + gamma121*gamma232)*gxyL + 
            (gamma132*gamma321 - gamma131*gamma322 - gamma122*gamma331 + gamma121*gamma332)*gxzL + 
            (-(gamma222*gamma231) + gamma221*gamma232)*gyyL + 
            (gamma232*gamma321 - gamma231*gamma322 - gamma222*gamma331 + gamma221*gamma332)*gyzL + 
            (-(gamma322*gamma331) + gamma321*gamma332)*gzzL + 
            khalf*(PDstandard2nd12gyz - PDstandard2nd13gyy - PDstandard2nd22gxz + PDstandard2nd23gxy);
        
        R3213  =  (-(gamma131*gamma132) + gamma121*gamma133)*gxxL + 
            (gamma133*gamma221 - gamma132*gamma231 - gamma131*gamma232 + gamma121*gamma233)*gxyL + 
            (gamma133*gamma321 - gamma132*gamma331 - gamma131*gamma332 + gamma121*gamma333)*gxzL + 
            (-(gamma231*gamma232) + gamma221*gamma233)*gyyL + 
            (gamma233*gamma321 - gamma232*gamma331 - gamma231*gamma332 + gamma221*gamma333)*gyzL + 
            (-(gamma331*gamma332) + gamma321*gamma333)*gzzL + 
            khalf*(PDstandard2nd12gzz - PDstandard2nd13gyz - PDstandard2nd23gxz + PDstandard2nd33gxy);
        
        R3221  =  (gamma122*gamma131 - gamma121*gamma132)*gxxL + 
            (-(gamma132*gamma221) + gamma131*gamma222 + gamma122*gamma231 - gamma121*gamma232)*gxyL + 
            (-(gamma132*gamma321) + gamma131*gamma322 + gamma122*gamma331 - gamma121*gamma332)*gxzL + 
            (gamma222*gamma231 - gamma221*gamma232)*gyyL + 
            (-(gamma232*gamma321) + gamma231*gamma322 + gamma222*gamma331 - gamma221*gamma332)*gyzL + 
            (gamma322*gamma331 - gamma321*gamma332)*gzzL + 
            khalf*(-PDstandard2nd12gyz + PDstandard2nd13gyy + PDstandard2nd22gxz - PDstandard2nd23gxy);
        
        R3222  =  0;
        
        R3223  =  (-2*gamma132*gamma232 + gamma122*gamma233)*gxyL + (-2*gamma132*gamma332 + gamma122*gamma333)*gxzL + 
            gamma133*(gamma122*gxxL + gamma222*gxyL + gamma322*gxzL) + (-2*gamma232*gamma332 + gamma222*gamma333)*gyzL + 
            gamma233*(gamma222*gyyL + gamma322*gyzL) - PDstandard2nd23gyz + khalf*(PDstandard2nd22gzz + PDstandard2nd33gyy) - 
            gxxL*SQR(gamma132) - gyyL*SQR(gamma232) + gzzL*(gamma322*gamma333 - SQR(gamma332));
        
        R3231  =  (gamma131*gamma132 - gamma121*gamma133)*gxxL + 
            (-(gamma133*gamma221) + gamma132*gamma231 + gamma131*gamma232 - gamma121*gamma233)*gxyL + 
            (-(gamma133*gamma321) + gamma132*gamma331 + gamma131*gamma332 - gamma121*gamma333)*gxzL + 
            (gamma231*gamma232 - gamma221*gamma233)*gyyL + 
            (-(gamma233*gamma321) + gamma232*gamma331 + gamma231*gamma332 - gamma221*gamma333)*gyzL + 
            (gamma331*gamma332 - gamma321*gamma333)*gzzL + 
            khalf*(-PDstandard2nd12gzz + PDstandard2nd13gyz + PDstandard2nd23gxz - PDstandard2nd33gxy);
        
        R3232  =  khalf*((4*gamma132*gamma232 - 2*gamma122*gamma233)*gxyL + (4*gamma132*gamma332 - 2*gamma122*gamma333)*gxzL + 
              (4*gamma232*gamma332 - 2*gamma222*gamma333)*gyzL - 
              2*(gamma122*gamma133*gxxL + gamma133*gamma222*gxyL + gamma133*gamma322*gxzL + gamma222*gamma233*gyyL + 
                 gamma233*gamma322*gyzL + gamma322*gamma333*gzzL) - PDstandard2nd22gzz + 2*PDstandard2nd23gyz - 
              PDstandard2nd33gyy + 2*gxxL*SQR(gamma132) + 2*gyyL*SQR(gamma232) + 2*gzzL*SQR(gamma332));
        
        R3233  =  0;
        
        R3311  =  0;
        
        R3312  =  0;
        
        R3313  =  0;
        
        R3321  =  0;
        
        R3322  =  0;
        
        R3323  =  0;
        
        R3331  =  0;
        
        R3332  =  0;
        
        R3333  =  0;
        
        R4p1111  =  R1111;
        
        R4p1112  =  R1112;
        
        R4p1113  =  R1113;
        
        R4p1121  =  R1121;
        
        R4p1122  =  R1122;
        
        R4p1123  =  R1123;
        
        R4p1131  =  R1131;
        
        R4p1132  =  R1132;
        
        R4p1133  =  R1133;
        
        R4p1211  =  R1211;
        
        R4p1212  =  kxxL*kyyL + R1212 - SQR(kxyL);
        
        R4p1213  =  -(kxyL*kxzL) + kxxL*kyzL + R1213;
        
        R4p1221  =  -(kxxL*kyyL) + R1221 + SQR(kxyL);
        
        R4p1222  =  R1222;
        
        R4p1223  =  -(kxzL*kyyL) + kxyL*kyzL + R1223;
        
        R4p1231  =  kxyL*kxzL - kxxL*kyzL + R1231;
        
        R4p1232  =  kxzL*kyyL - kxyL*kyzL + R1232;
        
        R4p1233  =  R1233;
        
        R4p1311  =  R1311;
        
        R4p1312  =  -(kxyL*kxzL) + kxxL*kyzL + R1312;
        
        R4p1313  =  kxxL*kzzL + R1313 - SQR(kxzL);
        
        R4p1321  =  kxyL*kxzL - kxxL*kyzL + R1321;
        
        R4p1322  =  R1322;
        
        R4p1323  =  -(kxzL*kyzL) + kxyL*kzzL + R1323;
        
        R4p1331  =  -(kxxL*kzzL) + R1331 + SQR(kxzL);
        
        R4p1332  =  kxzL*kyzL - kxyL*kzzL + R1332;
        
        R4p1333  =  R1333;
        
        R4p2111  =  R2111;
        
        R4p2112  =  -(kxxL*kyyL) + R2112 + SQR(kxyL);
        
        R4p2113  =  kxyL*kxzL - kxxL*kyzL + R2113;
        
        R4p2121  =  kxxL*kyyL + R2121 - SQR(kxyL);
        
        R4p2122  =  R2122;
        
        R4p2123  =  kxzL*kyyL - kxyL*kyzL + R2123;
        
        R4p2131  =  -(kxyL*kxzL) + kxxL*kyzL + R2131;
        
        R4p2132  =  -(kxzL*kyyL) + kxyL*kyzL + R2132;
        
        R4p2133  =  R2133;
        
        R4p2211  =  R2211;
        
        R4p2212  =  R2212;
        
        R4p2213  =  R2213;
        
        R4p2221  =  R2221;
        
        R4p2222  =  R2222;
        
        R4p2223  =  R2223;
        
        R4p2231  =  R2231;
        
        R4p2232  =  R2232;
        
        R4p2233  =  R2233;
        
        R4p2311  =  R2311;
        
        R4p2312  =  -(kxzL*kyyL) + kxyL*kyzL + R2312;
        
        R4p2313  =  -(kxzL*kyzL) + kxyL*kzzL + R2313;
        
        R4p2321  =  kxzL*kyyL - kxyL*kyzL + R2321;
        
        R4p2322  =  R2322;
        
        R4p2323  =  kyyL*kzzL + R2323 - SQR(kyzL);
        
        R4p2331  =  kxzL*kyzL - kxyL*kzzL + R2331;
        
        R4p2332  =  -(kyyL*kzzL) + R2332 + SQR(kyzL);
        
        R4p2333  =  R2333;
        
        R4p3111  =  R3111;
        
        R4p3112  =  kxyL*kxzL - kxxL*kyzL + R3112;
        
        R4p3113  =  -(kxxL*kzzL) + R3113 + SQR(kxzL);
        
        R4p3121  =  -(kxyL*kxzL) + kxxL*kyzL + R3121;
        
        R4p3122  =  R3122;
        
        R4p3123  =  kxzL*kyzL - kxyL*kzzL + R3123;
        
        R4p3131  =  kxxL*kzzL + R3131 - SQR(kxzL);
        
        R4p3132  =  -(kxzL*kyzL) + kxyL*kzzL + R3132;
        
        R4p3133  =  R3133;
        
        R4p3211  =  R3211;
        
        R4p3212  =  kxzL*kyyL - kxyL*kyzL + R3212;
        
        R4p3213  =  kxzL*kyzL - kxyL*kzzL + R3213;
        
        R4p3221  =  -(kxzL*kyyL) + kxyL*kyzL + R3221;
        
        R4p3222  =  R3222;
        
        R4p3223  =  -(kyyL*kzzL) + R3223 + SQR(kyzL);
        
        R4p3231  =  -(kxzL*kyzL) + kxyL*kzzL + R3231;
        
        R4p3232  =  kyyL*kzzL + R3232 - SQR(kyzL);
        
        R4p3233  =  R3233;
        
        R4p3311  =  R3311;
        
        R4p3312  =  R3312;
        
        R4p3313  =  R3313;
        
        R4p3321  =  R3321;
        
        R4p3322  =  R3322;
        
        R4p3323  =  R3323;
        
        R4p3331  =  R3331;
        
        R4p3332  =  R3332;
        
        R4p3333  =  R3333;
        
        Ro111  =  0;
        
        Ro112  =  gamma121*kxxL + (-gamma111 + gamma221)*kxyL + gamma321*kxzL - gamma211*kyyL - gamma311*kyzL + 
            PDstandard2nd1kxy - PDstandard2nd2kxx;
        
        Ro113  =  gamma131*kxxL + gamma231*kxyL + (-gamma111 + gamma331)*kxzL - gamma211*kyzL - gamma311*kzzL + 
            PDstandard2nd1kxz - PDstandard2nd3kxx;
        
        Ro121  =  -(gamma121*kxxL) + gamma111*kxyL - gamma221*kxyL - gamma321*kxzL + gamma211*kyyL + gamma311*kyzL - 
            PDstandard2nd1kxy + PDstandard2nd2kxx;
        
        Ro122  =  0;
        
        Ro123  =  gamma131*kxyL - gamma121*kxzL + gamma231*kyyL - gamma221*kyzL + gamma331*kyzL - gamma321*kzzL + 
            PDstandard2nd2kxz - PDstandard2nd3kxy;
        
        Ro131  =  -(gamma131*kxxL) - gamma231*kxyL + gamma111*kxzL - gamma331*kxzL + gamma211*kyzL + gamma311*kzzL - 
            PDstandard2nd1kxz + PDstandard2nd3kxx;
        
        Ro132  =  -(gamma131*kxyL) + gamma121*kxzL - gamma231*kyyL + gamma221*kyzL - gamma331*kyzL + gamma321*kzzL - 
            PDstandard2nd2kxz + PDstandard2nd3kxy;
        
        Ro133  =  0;
        
        Ro211  =  0;
        
        Ro212  =  gamma122*kxxL + (-gamma121 + gamma222)*kxyL + gamma322*kxzL - gamma221*kyyL - gamma321*kyzL + 
            PDstandard2nd1kyy - PDstandard2nd2kxy;
        
        Ro213  =  gamma132*kxxL + gamma232*kxyL + (-gamma121 + gamma332)*kxzL - gamma221*kyzL - gamma321*kzzL + 
            PDstandard2nd1kyz - PDstandard2nd3kxy;
        
        Ro221  =  -(gamma122*kxxL) + gamma121*kxyL - gamma222*kxyL - gamma322*kxzL + gamma221*kyyL + gamma321*kyzL - 
            PDstandard2nd1kyy + PDstandard2nd2kxy;
        
        Ro222  =  0;
        
        Ro223  =  gamma132*kxyL - gamma122*kxzL + gamma232*kyyL - gamma222*kyzL + gamma332*kyzL - gamma322*kzzL + 
            PDstandard2nd2kyz - PDstandard2nd3kyy;
        
        Ro231  =  -(gamma132*kxxL) - gamma232*kxyL + gamma121*kxzL - gamma332*kxzL + gamma221*kyzL + gamma321*kzzL - 
            PDstandard2nd1kyz + PDstandard2nd3kxy;
        
        Ro232  =  -(gamma132*kxyL) + gamma122*kxzL - gamma232*kyyL + gamma222*kyzL - gamma332*kyzL + gamma322*kzzL - 
            PDstandard2nd2kyz + PDstandard2nd3kyy;
        
        Ro233  =  0;
        
        Ro311  =  0;
        
        Ro312  =  gamma132*kxxL + (-gamma131 + gamma232)*kxyL + gamma332*kxzL - gamma231*kyyL - gamma331*kyzL + 
            PDstandard2nd1kyz - PDstandard2nd2kxz;
        
        Ro313  =  gamma133*kxxL + gamma233*kxyL + (-gamma131 + gamma333)*kxzL - gamma231*kyzL - gamma331*kzzL + 
            PDstandard2nd1kzz - PDstandard2nd3kxz;
        
        Ro321  =  -(gamma132*kxxL) + gamma131*kxyL - gamma232*kxyL - gamma332*kxzL + gamma231*kyyL + gamma331*kyzL - 
            PDstandard2nd1kyz + PDstandard2nd2kxz;
        
        Ro322  =  0;
        
        Ro323  =  gamma133*kxyL - gamma132*kxzL + gamma233*kyyL - gamma232*kyzL + gamma333*kyzL - gamma332*kzzL + 
            PDstandard2nd2kzz - PDstandard2nd3kyz;
        
        Ro331  =  -(gamma133*kxxL) - gamma233*kxyL + gamma131*kxzL - gamma333*kxzL + gamma231*kyzL + gamma331*kzzL - 
            PDstandard2nd1kzz + PDstandard2nd3kxz;
        
        Ro332  =  -(gamma133*kxyL) + gamma132*kxzL - gamma233*kyyL + gamma232*kyzL - gamma333*kyzL + gamma332*kzzL - 
            PDstandard2nd2kzz + PDstandard2nd3kyz;
        
        Ro333  =  0;
        
        Rojo11  =  gInv11*R1111 + gInv12*R1112 + gInv13*R1113 + gInv21*R1211 + gInv23*(-(kxyL*kxzL) + kxxL*kyzL + R1213) + 
            gInv31*R1311 + gInv32*(-(kxyL*kxzL) + kxxL*kyzL + R1312) + gInv22*(kxxL*kyyL + R1212 - SQR(kxyL)) + 
            gInv33*(kxxL*kzzL + R1313 - SQR(kxzL));
        
        Rojo12  =  gInv11*R1121 + gInv13*(kxyL*kxzL - kxxL*kyzL + R1123) + gInv21*R1221 + gInv22*R1222 + gInv23*R1223 + 
            gInv31*R1321 + gInv32*(-(kxzL*kyyL) + kxyL*kyzL + R1322) + gInv33*(-(kxzL*kyzL) + kxyL*kzzL + R1323) + 
            gInv12*(-(kxxL*kyyL) + R1122 + SQR(kxyL));
        
        Rojo13  =  gInv11*R1131 + gInv12*(kxyL*kxzL - kxxL*kyzL + R1132) + gInv21*R1231 + 
            gInv22*(kxzL*kyyL - kxyL*kyzL + R1232) + gInv23*(kxzL*kyzL - kxyL*kzzL + R1233) + gInv31*R1331 + gInv32*R1332 + 
            gInv33*R1333 + gInv13*(-(kxxL*kzzL) + R1133 + SQR(kxzL));
        
        Rojo21  =  gInv11*R2111 + gInv12*R2112 + gInv13*R2113 + gInv22*R2212 + gInv23*(-(kxzL*kyyL) + kxyL*kyzL + R2213) + 
            gInv31*(kxyL*kxzL - kxxL*kyzL + R2311) + gInv32*R2312 + gInv33*(-(kxzL*kyzL) + kxyL*kzzL + R2313) + 
            gInv21*(-(kxxL*kyyL) + R2211 + SQR(kxyL));
        
        Rojo22  =  gInv12*R2122 + gInv13*(kxzL*kyyL - kxyL*kyzL + R2123) + gInv21*R2221 + gInv22*R2222 + gInv23*R2223 + 
            gInv31*(kxzL*kyyL - kxyL*kyzL + R2321) + gInv32*R2322 + gInv11*(kxxL*kyyL + R2121 - SQR(kxyL)) + 
            gInv33*(kyyL*kzzL + R2323 - SQR(kyzL));
        
        Rojo23  =  gInv11*(-(kxyL*kxzL) + kxxL*kyzL + R2131) + gInv12*R2132 + gInv13*(kxzL*kyzL - kxyL*kzzL + R2133) + 
            gInv21*(-(kxzL*kyyL) + kxyL*kyzL + R2231) + gInv22*R2232 + gInv31*R2331 + gInv32*R2332 + gInv33*R2333 + 
            gInv23*(-(kyyL*kzzL) + R2233 + SQR(kyzL));
        
        Rojo31  =  gInv11*R3111 + gInv12*R3112 + gInv13*R3113 + gInv21*(kxyL*kxzL - kxxL*kyzL + R3211) + 
            gInv22*(kxzL*kyyL - kxyL*kyzL + R3212) + gInv23*R3213 + gInv32*(kxzL*kyzL - kxyL*kzzL + R3312) + gInv33*R3313 + 
            gInv31*(-(kxxL*kzzL) + R3311 + SQR(kxzL));
        
        Rojo32  =  gInv11*(-(kxyL*kxzL) + kxxL*kyzL + R3121) + gInv12*(-(kxzL*kyyL) + kxyL*kyzL + R3122) + gInv13*R3123 + 
            gInv21*R3221 + gInv22*R3222 + gInv23*R3223 + gInv31*(kxzL*kyzL - kxyL*kzzL + R3321) + gInv33*R3323 + 
            gInv32*(-(kyyL*kzzL) + R3322 + SQR(kyzL));
        
        Rojo33  =  gInv12*(-(kxzL*kyzL) + kxyL*kzzL + R3132) + gInv13*R3133 + gInv21*(-(kxzL*kyzL) + kxyL*kzzL + R3231) + 
            gInv23*R3233 + gInv31*R3331 + gInv32*R3332 + gInv33*R3333 + gInv11*(kxxL*kzzL + R3131 - SQR(kxzL)) + 
            gInv22*(kyyL*kzzL + R3232 - SQR(kyzL));
        
        Psi4rL  =  n1*(-(imbar1*(imbar2*(n1*(R4p1112 + R4p1211) + n2*(R4p1122 + R4p1221) + n3*(R4p1132 + R4p1231)) + 
                    imbar3*(n1*(R4p1113 + R4p1311) + n2*(R4p1123 + R4p1321) + n3*(R4p1133 + R4p1331)))) - 
               imbar2*imbar3*(n1*(R4p1213 + R4p1312) + n2*(R4p1223 + R4p1322) + n3*(R4p1233 + R4p1332)) + 
               n1*R4p1112*rmbar1*rmbar2 + n2*R4p1122*rmbar1*rmbar2 + n3*R4p1132*rmbar1*rmbar2 + n1*R4p1211*rmbar1*rmbar2 + 
               n2*R4p1221*rmbar1*rmbar2 + n3*R4p1231*rmbar1*rmbar2 + n1*R4p1113*rmbar1*rmbar3 + n2*R4p1123*rmbar1*rmbar3 + 
               n3*R4p1133*rmbar1*rmbar3 + n1*R4p1311*rmbar1*rmbar3 + n2*R4p1321*rmbar1*rmbar3 + n3*R4p1331*rmbar1*rmbar3 + 
               n1*R4p1213*rmbar2*rmbar3 + n2*R4p1223*rmbar2*rmbar3 + n3*R4p1233*rmbar2*rmbar3 + n1*R4p1312*rmbar2*rmbar3 + 
               n2*R4p1322*rmbar2*rmbar3 + n3*R4p1332*rmbar2*rmbar3 - (n1*R4p1111 + n2*R4p1121 + n3*R4p1131)*SQR(imbar1) - 
               (n1*R4p1212 + n2*R4p1222 + n3*R4p1232)*SQR(imbar2) - n1*R4p1313*SQR(imbar3) - n2*R4p1323*SQR(imbar3) - 
               n3*R4p1333*SQR(imbar3) + n1*R4p1111*SQR(rmbar1) + n2*R4p1121*SQR(rmbar1) + n3*R4p1131*SQR(rmbar1) + 
               n1*R4p1212*SQR(rmbar2) + n2*R4p1222*SQR(rmbar2) + n3*R4p1232*SQR(rmbar2) + n1*R4p1313*SQR(rmbar3) + 
               n2*R4p1323*SQR(rmbar3) + n3*R4p1333*SQR(rmbar3)) + 
            n2*(-(imbar1*(imbar2*(n1*(R4p2112 + R4p2211) + n2*(R4p2122 + R4p2221) + n3*(R4p2132 + R4p2231)) + 
                    imbar3*(n1*(R4p2113 + R4p2311) + n2*(R4p2123 + R4p2321) + n3*(R4p2133 + R4p2331)))) - 
               imbar2*imbar3*(n1*(R4p2213 + R4p2312) + n2*(R4p2223 + R4p2322) + n3*(R4p2233 + R4p2332)) + 
               n1*R4p2112*rmbar1*rmbar2 + n2*R4p2122*rmbar1*rmbar2 + n3*R4p2132*rmbar1*rmbar2 + n1*R4p2211*rmbar1*rmbar2 + 
               n2*R4p2221*rmbar1*rmbar2 + n3*R4p2231*rmbar1*rmbar2 + n1*R4p2113*rmbar1*rmbar3 + n2*R4p2123*rmbar1*rmbar3 + 
               n3*R4p2133*rmbar1*rmbar3 + n1*R4p2311*rmbar1*rmbar3 + n2*R4p2321*rmbar1*rmbar3 + n3*R4p2331*rmbar1*rmbar3 + 
               n1*R4p2213*rmbar2*rmbar3 + n2*R4p2223*rmbar2*rmbar3 + n3*R4p2233*rmbar2*rmbar3 + n1*R4p2312*rmbar2*rmbar3 + 
               n2*R4p2322*rmbar2*rmbar3 + n3*R4p2332*rmbar2*rmbar3 - (n1*R4p2111 + n2*R4p2121 + n3*R4p2131)*SQR(imbar1) - 
               (n1*R4p2212 + n2*R4p2222 + n3*R4p2232)*SQR(imbar2) - n1*R4p2313*SQR(imbar3) - n2*R4p2323*SQR(imbar3) - 
               n3*R4p2333*SQR(imbar3) + n1*R4p2111*SQR(rmbar1) + n2*R4p2121*SQR(rmbar1) + n3*R4p2131*SQR(rmbar1) + 
               n1*R4p2212*SQR(rmbar2) + n2*R4p2222*SQR(rmbar2) + n3*R4p2232*SQR(rmbar2) + n1*R4p2313*SQR(rmbar3) + 
               n2*R4p2323*SQR(rmbar3) + n3*R4p2333*SQR(rmbar3)) + 
            n3*(-(imbar1*(imbar2*(n1*(R4p3112 + R4p3211) + n2*(R4p3122 + R4p3221) + n3*(R4p3132 + R4p3231)) + 
                    imbar3*(n1*(R4p3113 + R4p3311) + n2*(R4p3123 + R4p3321) + n3*(R4p3133 + R4p3331)))) - 
               imbar2*imbar3*(n1*(R4p3213 + R4p3312) + n2*(R4p3223 + R4p3322) + n3*(R4p3233 + R4p3332)) + 
               n1*R4p3112*rmbar1*rmbar2 + n2*R4p3122*rmbar1*rmbar2 + n3*R4p3132*rmbar1*rmbar2 + n1*R4p3211*rmbar1*rmbar2 + 
               n2*R4p3221*rmbar1*rmbar2 + n3*R4p3231*rmbar1*rmbar2 + n1*R4p3113*rmbar1*rmbar3 + n2*R4p3123*rmbar1*rmbar3 + 
               n3*R4p3133*rmbar1*rmbar3 + n1*R4p3311*rmbar1*rmbar3 + n2*R4p3321*rmbar1*rmbar3 + n3*R4p3331*rmbar1*rmbar3 + 
               n1*R4p3213*rmbar2*rmbar3 + n2*R4p3223*rmbar2*rmbar3 + n3*R4p3233*rmbar2*rmbar3 + n1*R4p3312*rmbar2*rmbar3 + 
               n2*R4p3322*rmbar2*rmbar3 + n3*R4p3332*rmbar2*rmbar3 - (n1*R4p3111 + n2*R4p3121 + n3*R4p3131)*SQR(imbar1) - 
               (n1*R4p3212 + n2*R4p3222 + n3*R4p3232)*SQR(imbar2) - n1*R4p3313*SQR(imbar3) - n2*R4p3323*SQR(imbar3) - 
               n3*R4p3333*SQR(imbar3) + n1*R4p3111*SQR(rmbar1) + n2*R4p3121*SQR(rmbar1) + n3*R4p3131*SQR(rmbar1) + 
               n1*R4p3212*SQR(rmbar2) + n2*R4p3222*SQR(rmbar2) + n3*R4p3232*SQR(rmbar2) + n1*R4p3313*SQR(rmbar3) + 
               n2*R4p3323*SQR(rmbar3) + n3*R4p3333*SQR(rmbar3)) - 
            SQR(n0)*((imbar2*imbar3 - rmbar2*rmbar3)*Rojo23 + imbar1*(imbar2*(Rojo12 + Rojo21) + imbar3*(Rojo13 + Rojo31)) - 
               rmbar1*(rmbar2*(Rojo12 + Rojo21) + rmbar3*(Rojo13 + Rojo31)) + imbar2*imbar3*Rojo32 - rmbar2*rmbar3*Rojo32 + 
               Rojo11*SQR(imbar1) + Rojo22*SQR(imbar2) + Rojo33*SQR(imbar3) - Rojo11*SQR(rmbar1) - Rojo22*SQR(rmbar2) - 
               Rojo33*SQR(rmbar3)) + 2*n0*(-(n1*(imbar2*imbar3*Ro213 - rmbar2*rmbar3*Ro213 + 
                    imbar1*(imbar2*(Ro112 + Ro211) + imbar3*(Ro113 + Ro311)) - 
                    rmbar1*(rmbar2*(Ro112 + Ro211) + rmbar3*(Ro113 + Ro311)) + imbar2*imbar3*Ro312 - rmbar2*rmbar3*Ro312 + 
                    Ro111*SQR(imbar1) + Ro212*SQR(imbar2) + Ro313*SQR(imbar3) - Ro111*SQR(rmbar1) - Ro212*SQR(rmbar2) - 
                    Ro313*SQR(rmbar3))) - n2*(imbar2*imbar3*Ro223 - rmbar2*rmbar3*Ro223 + 
                  imbar1*(imbar2*(Ro122 + Ro221) + imbar3*(Ro123 + Ro321)) - 
                  rmbar1*(rmbar2*(Ro122 + Ro221) + rmbar3*(Ro123 + Ro321)) + imbar2*imbar3*Ro322 - rmbar2*rmbar3*Ro322 + 
                  Ro121*SQR(imbar1) + Ro222*SQR(imbar2) + Ro323*SQR(imbar3) - Ro121*SQR(rmbar1) - Ro222*SQR(rmbar2) - 
                  Ro323*SQR(rmbar3)) - n3*(imbar2*imbar3*Ro233 - rmbar2*rmbar3*Ro233 + 
                  imbar1*(imbar2*(Ro132 + Ro231) + imbar3*(Ro133 + Ro331)) - 
                  rmbar1*(rmbar2*(Ro132 + Ro231) + rmbar3*(Ro133 + Ro331)) + imbar2*imbar3*Ro332 - rmbar2*rmbar3*Ro332 + 
                  Ro131*SQR(imbar1) + Ro232*SQR(imbar2) + Ro333*SQR(imbar3) - Ro131*SQR(rmbar1) - Ro232*SQR(rmbar2) - 
                  Ro333*SQR(rmbar3)));
        
        Psi4iL  =  -(n1*(im1*(n1*(2*R4p1111*rm1 + R4p1112*rm2 + R4p1211*rm2 + R4p1113*rm3 + R4p1311*rm3) + 
                    n2*(2*R4p1121*rm1 + R4p1122*rm2 + R4p1221*rm2 + R4p1123*rm3 + R4p1321*rm3) + 
                    n3*(2*R4p1131*rm1 + R4p1132*rm2 + R4p1231*rm2 + R4p1133*rm3 + R4p1331*rm3)) + 
                 im2*(n1*(R4p1112*rm1 + R4p1211*rm1 + 2*R4p1212*rm2 + R4p1213*rm3 + R4p1312*rm3) + 
                    n2*(R4p1122*rm1 + R4p1221*rm1 + 2*R4p1222*rm2 + R4p1223*rm3 + R4p1322*rm3) + 
                    n3*(R4p1132*rm1 + R4p1231*rm1 + 2*R4p1232*rm2 + R4p1233*rm3 + R4p1332*rm3)) + 
                 im3*(n1*(R4p1113*rm1 + R4p1311*rm1 + R4p1213*rm2 + R4p1312*rm2 + 2*R4p1313*rm3) + 
                    n2*(R4p1123*rm1 + R4p1321*rm1 + R4p1223*rm2 + R4p1322*rm2 + 2*R4p1323*rm3) + 
                    n3*(R4p1133*rm1 + R4p1331*rm1 + R4p1233*rm2 + R4p1332*rm2 + 2*R4p1333*rm3)))) - 
            n2*(im1*(n1*(2*R4p2111*rm1 + R4p2112*rm2 + R4p2211*rm2 + R4p2113*rm3 + R4p2311*rm3) + 
                  n2*(2*R4p2121*rm1 + R4p2122*rm2 + R4p2221*rm2 + R4p2123*rm3 + R4p2321*rm3) + 
                  n3*(2*R4p2131*rm1 + R4p2132*rm2 + R4p2231*rm2 + R4p2133*rm3 + R4p2331*rm3)) + 
               im2*(n1*(R4p2112*rm1 + R4p2211*rm1 + 2*R4p2212*rm2 + R4p2213*rm3 + R4p2312*rm3) + 
                  n2*(R4p2122*rm1 + R4p2221*rm1 + 2*R4p2222*rm2 + R4p2223*rm3 + R4p2322*rm3) + 
                  n3*(R4p2132*rm1 + R4p2231*rm1 + 2*R4p2232*rm2 + R4p2233*rm3 + R4p2332*rm3)) + 
               im3*(n1*(R4p2113*rm1 + R4p2311*rm1 + R4p2213*rm2 + R4p2312*rm2 + 2*R4p2313*rm3) + 
                  n2*(R4p2123*rm1 + R4p2321*rm1 + R4p2223*rm2 + R4p2322*rm2 + 2*R4p2323*rm3) + 
                  n3*(R4p2133*rm1 + R4p2331*rm1 + R4p2233*rm2 + R4p2332*rm2 + 2*R4p2333*rm3))) - 
            n3*(im1*(n1*(2*R4p3111*rm1 + R4p3112*rm2 + R4p3211*rm2 + R4p3113*rm3 + R4p3311*rm3) + 
                  n2*(2*R4p3121*rm1 + R4p3122*rm2 + R4p3221*rm2 + R4p3123*rm3 + R4p3321*rm3) + 
                  n3*(2*R4p3131*rm1 + R4p3132*rm2 + R4p3231*rm2 + R4p3133*rm3 + R4p3331*rm3)) + 
               im2*(n1*(R4p3112*rm1 + R4p3211*rm1 + 2*R4p3212*rm2 + R4p3213*rm3 + R4p3312*rm3) + 
                  n2*(R4p3122*rm1 + R4p3221*rm1 + 2*R4p3222*rm2 + R4p3223*rm3 + R4p3322*rm3) + 
                  n3*(R4p3132*rm1 + R4p3231*rm1 + 2*R4p3232*rm2 + R4p3233*rm3 + R4p3332*rm3)) + 
               im3*(n1*(R4p3113*rm1 + R4p3311*rm1 + R4p3213*rm2 + R4p3312*rm2 + 2*R4p3313*rm3) + 
                  n2*(R4p3123*rm1 + R4p3321*rm1 + R4p3223*rm2 + R4p3322*rm2 + 2*R4p3323*rm3) + 
                  n3*(R4p3133*rm1 + R4p3331*rm1 + R4p3233*rm2 + R4p3332*rm2 + 2*R4p3333*rm3))) - 
            2*n0*(im1*(n1*(2*rm1*Ro111 + rm2*(Ro112 + Ro211) + rm3*(Ro113 + Ro311)) + 
                  n2*(2*rm1*Ro121 + rm2*(Ro122 + Ro221) + rm3*(Ro123 + Ro321)) + 
                  n3*(2*rm1*Ro131 + rm2*(Ro132 + Ro231) + rm3*(Ro133 + Ro331))) + 
               im2*(n1*(rm1*(Ro112 + Ro211) + 2*rm2*Ro212 + rm3*(Ro213 + Ro312)) + 
                  n2*(rm1*(Ro122 + Ro221) + 2*rm2*Ro222 + rm3*(Ro223 + Ro322)) + 
                  n3*(rm1*(Ro132 + Ro231) + 2*rm2*Ro232 + rm3*(Ro233 + Ro332))) + 
               im3*(n1*(rm1*(Ro113 + Ro311) + rm2*(Ro213 + Ro312) + 2*rm3*Ro313) + 
                  n2*(rm1*(Ro123 + Ro321) + rm2*(Ro223 + Ro322) + 2*rm3*Ro323) + 
                  n3*(rm1*(Ro133 + Ro331) + rm2*(Ro233 + Ro332) + 2*rm3*Ro333))) - 
            (im1*(2*rm1*Rojo11 + rm2*(Rojo12 + Rojo21) + rm3*(Rojo13 + Rojo31)) + 
               im2*(rm1*(Rojo12 + Rojo21) + 2*rm2*Rojo22 + rm3*(Rojo23 + Rojo32)) + 
               im3*(rm1*(Rojo13 + Rojo31) + rm2*(Rojo23 + Rojo32) + 2*rm3*Rojo33))*SQR(n0);
        
        
        /* Copy local copies back to grid functions */
        Psi4i[index] = Psi4iL;
        Psi4r[index] = Psi4rL;
        x[index] = xL;
        y[index] = yL;
        z[index] = zL;
        
        /* Copy local copies back to subblock grid functions */
      }
    }
  }
}

void psis_calc_2nd(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  GenericFD_LoopOverInterior(cctkGH, &psis_calc_2nd_Body);
}
