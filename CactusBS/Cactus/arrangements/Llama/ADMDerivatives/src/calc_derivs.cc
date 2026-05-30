
/* Copyright 2013 Peter Diener, Nils Dorband, Roland Haas, Ian Hinder,
Christian Ott, Denis Pollney, Thomas Radke, Christian Reisswig, Erik
Schnetter, Barry Wardell and Burkhard Zink

This file is part of Llama.

Llama is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 2 of the License, or (at your
option) any later version.

Llama is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with Llama.  If not, see <http://www.gnu.org/licenses/>. */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include <assert.h>
#include <iostream>

#include "loopcontrol.h"

#include "GlobalDerivative.h"
#include "grids.hh"

using namespace std;


#define invert_metric(g11, g12, g13, g22, g23, g33, detg, \
                              ig11, ig12, ig13, ig22, ig23, ig33) \
  do                                                              \
    {                                                             \
      CCTK_REAL detg11 = g22*g33 - g23*g23;                       \
      CCTK_REAL detg12 = g13*g23 - g12*g33;                       \
      CCTK_REAL detg13 = g12*g23 - g13*g22;                       \
      CCTK_REAL detg22 = g11*g33 - g13*g13;                       \
      CCTK_REAL detg23 = g12*g13 - g11*g23;                       \
      CCTK_REAL detg33 = g11*g22 - g12*g12;                       \
                                                                  \
      detg = detg11*g11 + detg12*g12 + detg13*g13;                \
                                                                  \
      ig11 = detg11 / detg;                                       \
      ig12 = detg12 / detg;                                       \
      ig13 = detg13 / detg;                                       \
      ig22 = detg22 / detg;                                       \
      ig23 = detg23 / detg;                                       \
      ig33 = detg33 / detg;                                       \
    }                                                             \
  while (0)



extern "C" void ADMDerivatives_CalcDerivatives(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS
   DECLARE_CCTK_PARAMETERS
   
   bool radial_coordinate = false;
   
   // calculate radial derivative
   /*if (CCTK_IsFunctionAliased("MultiPatch_ProvidesThornburg04"))
   {
      if (MultiPatch_ProvidesThornburg04())
      {
         int map = MultiPatch_GetMap(cctkGH);
         // don't do anything on Cartesian map
         //if (map != 0)  // we would also have to check for sphere_inner_radius...for now we assume the user knows that the extraction slice has to be on spherical grid
            radial_coordinate = true;
         
      }
   }*/
   
   
   CCTK_INT* imin[3];
   CCTK_INT* imax[3];
   CCTK_REAL* q[3];
   int ni = cctk_lsh[0]; 
   int nj = cctk_lsh[1];
   int nk = cctk_lsh[2];
   
   CCTK_REAL ihx = 1.0 / CCTK_DELTA_SPACE(0);
   CCTK_REAL ihy = 1.0 / CCTK_DELTA_SPACE(1);
   CCTK_REAL ihz = 1.0 / CCTK_DELTA_SPACE(2);
   
   int orig_fd_order =
    *(const CCTK_INT*)CCTK_ParameterGet("order", "SummationByParts", NULL);
   
   if (spatial_deriv_order > 0)
   {
      ostringstream fd_order_str;
      fd_order_str << spatial_deriv_order;
      CCTK_ParameterSet("order", "SummationByParts", fd_order_str.str().c_str());
   }
   
   for (int i=0; i<3; ++i)
   {
      int n = cctk_lsh[i];
      imin[i] = new CCTK_INT[n];
      imax[i] = new CCTK_INT[n];
      q[i] = new CCTK_REAL[n*n];

      if ((imin[i]==NULL) || (imax[i]==NULL) || (q[i]==NULL))
      CCTK_WARN(0, "Could not allocated derivative coefficient arrays.");

      Diff_coeff(cctkGH, i, n, imin[i], imax[i], q[i], -1);

      for (int j=0; j<n; ++j)
      {
         imin[i][j] -= 1;
         imax[i][j] -= 1;
      }
   }
   
   if (spatial_deriv_order > 0)
   {
      ostringstream fd_order_str;
      fd_order_str << orig_fd_order;
      CCTK_ParameterSet("order", "SummationByParts", fd_order_str.str().c_str());
   }
   
   CCTK_INT istart[3], iend[3];

   ADMDerivatives_GetGridRanges(cctkGH, istart, iend);
   
   if (store_radial_derivatives || store_cartesian_derivatives || store_time_derivatives)
   {
      if (!radial_coordinate)
      {     
         // we would have to make use of Jacobains to transform to a radial coordinate system
#pragma omp parallel
         LC_LOOP3 (ADMDerivatives_CalcDerivatives,
               i, j, k,
               istart[0], istart[1], istart[2],
               iend[0], iend[1], iend[2],
               cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
         {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
         
            // calculate derivatives \partial_i \beta^j
            CCTK_REAL dxBetax(0), dxBetay(0), dyBetax(0), dxBetaz(0), dzBetax(0);
            CCTK_REAL dyBetay(0), dyBetaz(0), dzBetay(0), dzBetaz(0);
            CCTK_REAL dxalp(0), dyalp(0), dzalp(0);
   
            CCTK_REAL dadx=1, dbdx=0, dcdx=0, dady=0, dbdy=1, dcdy=0, dadz=0,
                      dbdz=0, dcdz=1;

            CCTK_REAL dg[3][3][3], detg;

#include "derivatives.h"
            
            if (store_radial_derivatives)
            {
               const CCTK_REAL irad = 1.0/sqrt(x[ijk]*x[ijk] + y[ijk]*y[ijk] + z[ijk]*z[ijk]);
               
               alp_dr[ijk] = irad * (x[ijk]*dxalp + y[ijk]*dyalp + z[ijk]*dzalp);
               
               betax_dr[ijk] = irad * (x[ijk]*dxBetax + y[ijk]*dyBetax + z[ijk]*dzBetax);
               betay_dr[ijk] = irad * (x[ijk]*dxBetay + y[ijk]*dyBetay + z[ijk]*dzBetay);
               betaz_dr[ijk] = irad * (x[ijk]*dxBetaz + y[ijk]*dyBetaz + z[ijk]*dzBetaz);
                           
               gxx_dr[ijk] = irad * (x[ijk]*dg[0][0][0] + y[ijk]*dg[0][0][1] + z[ijk]*dg[0][0][2]);
               gxy_dr[ijk] = irad * (x[ijk]*dg[0][1][0] + y[ijk]*dg[0][1][1] + z[ijk]*dg[0][1][2]);
               gxz_dr[ijk] = irad * (x[ijk]*dg[0][2][0] + y[ijk]*dg[0][2][1] + z[ijk]*dg[0][2][2]);
               gyy_dr[ijk] = irad * (x[ijk]*dg[1][1][0] + y[ijk]*dg[1][1][1] + z[ijk]*dg[1][1][2]);
               gyz_dr[ijk] = irad * (x[ijk]*dg[1][2][0] + y[ijk]*dg[1][2][1] + z[ijk]*dg[1][2][2]);
               gzz_dr[ijk] = irad * (x[ijk]*dg[2][2][0] + y[ijk]*dg[2][2][1] + z[ijk]*dg[2][2][2]);
            }
            
            if (store_cartesian_derivatives)
            {
               alp_dx[ijk] = dxalp;
               
               betax_dx[ijk] = dxBetax;
               betay_dx[ijk] = dxBetay;
               betaz_dx[ijk] = dxBetaz;
                           
               gxx_dx[ijk] = dg[0][0][0];
               gxy_dx[ijk] = dg[0][1][0];
               gxz_dx[ijk] = dg[0][2][0];
               gyy_dx[ijk] = dg[1][1][0];
               gyz_dx[ijk] = dg[1][2][0];
               gzz_dx[ijk] = dg[2][2][0];
               
               alp_dy[ijk] = dyalp;
               
               betax_dy[ijk] = dyBetax;
               betay_dy[ijk] = dyBetay;
               betaz_dy[ijk] = dyBetaz;
                           
               gxx_dy[ijk] = dg[0][0][1];
               gxy_dy[ijk] = dg[0][1][1];
               gxz_dy[ijk] = dg[0][2][1];
               gyy_dy[ijk] = dg[1][1][1];
               gyz_dy[ijk] = dg[1][2][1];
               gzz_dy[ijk] = dg[2][2][1];
               
               alp_dz[ijk] = dzalp;
               
               betax_dz[ijk] = dzBetax;
               betay_dz[ijk] = dzBetay;
               betaz_dz[ijk] = dzBetaz;
                           
               gxx_dz[ijk] = dg[0][0][2];
               gxy_dz[ijk] = dg[0][1][2];
               gxz_dz[ijk] = dg[0][2][2];
               gyy_dz[ijk] = dg[1][1][2];
               gyz_dz[ijk] = dg[1][2][2];
               gzz_dz[ijk] = dg[2][2][2];
            }
            
            if (store_time_derivatives)
            {
               for (int ii=0; ii < 3; ++ii)
               {
                  dg[1][0][ii] = dg[0][1][ii];
                  dg[2][0][ii] = dg[0][2][ii];
                  dg[2][1][ii] = dg[1][2][ii];
               } 
               
               // calculate inverse metric compoenents
               CCTK_REAL ig[3][3];
               
               invert_metric(gxx[ijk], gxy[ijk], gxz[ijk], gyy[ijk], gyz[ijk], gzz[ijk], detg,
                           ig[0][0], ig[0][1], ig[0][2], ig[1][1], ig[1][2], ig[2][2]);
               
               ig[1][0] = ig[0][1];
               ig[2][0] = ig[0][2];
               ig[2][1] = ig[1][2];
               
               // calculate Christoffels
               
               CCTK_REAL C[3][3][3]; // last index is the upper index of Christoffel symbol
                     
               for (int ii=0; ii < 3; ++ii)
                  for (int jj=ii; jj < 3; ++jj)
                     for (int kk=0; kk < 3; ++kk)
                     {
                        C[ii][jj][kk] = 0.5* ( ig[kk][0] * (dg[0][ii][jj] + dg[0][jj][ii] - dg[ii][jj][0])
                                             + ig[kk][1] * (dg[1][ii][jj] + dg[1][jj][ii] - dg[ii][jj][1])
                                             + ig[kk][2] * (dg[2][ii][jj] + dg[2][jj][ii] - dg[ii][jj][2]) );
                        if (ii != jj)
                           C[jj][ii][kk] = C[ii][jj][kk];
                     }
               
               
               // calculate covariant derivative D_i \beta^j
               CCTK_REAL DxBetax(0), DxBetay(0), DyBetax(0), DxBetaz(0), DzBetax(0);
               CCTK_REAL DyBetay(0), DyBetaz(0), DzBetay(0), DzBetaz(0);
               
               DxBetax = dxBetax + (C[0][0][0]*betax[ijk] + C[0][1][0]*betay[ijk] + C[0][2][0]*betaz[ijk]);
               DxBetay = dxBetay + (C[0][0][1]*betax[ijk] + C[0][1][1]*betay[ijk] + C[0][2][1]*betaz[ijk]);
               DxBetaz = dxBetaz + (C[0][0][2]*betax[ijk] + C[0][1][2]*betay[ijk] + C[0][2][2]*betaz[ijk]);
         
               DyBetax = dyBetax + (C[1][0][0]*betax[ijk] + C[1][1][0]*betay[ijk] + C[1][2][0]*betaz[ijk]);
               DyBetay = dyBetay + (C[1][0][1]*betax[ijk] + C[1][1][1]*betay[ijk] + C[1][2][1]*betaz[ijk]);
               DyBetaz = dyBetaz + (C[1][0][2]*betax[ijk] + C[1][1][2]*betay[ijk] + C[1][2][2]*betaz[ijk]);
         
               DzBetax = dzBetax + (C[2][0][0]*betax[ijk] + C[2][1][0]*betay[ijk] + C[2][2][0]*betaz[ijk]);
               DzBetay = dzBetay + (C[2][0][1]*betax[ijk] + C[2][1][1]*betay[ijk] + C[2][2][1]*betaz[ijk]);
               DzBetaz = dzBetaz + (C[2][0][2]*betax[ijk] + C[2][1][2]*betay[ijk] + C[2][2][2]*betaz[ijk]);
         
               // get index down: D_i \beta_j 
               CCTK_REAL Dxbetax(0), Dxbetay(0), Dybetax(0), Dxbetaz(0), Dzbetax(0);
               CCTK_REAL Dybetay(0), Dybetaz(0), Dzbetay(0), Dzbetaz(0);
               
               Dxbetax = gxx[ijk]*DxBetax + gxy[ijk]*DxBetay + gxz[ijk]*DxBetaz;
               Dxbetay = gxy[ijk]*DxBetax + gyy[ijk]*DxBetay + gyz[ijk]*DxBetaz;
               Dxbetaz = gxz[ijk]*DxBetax + gyz[ijk]*DxBetay + gzz[ijk]*DxBetaz;
               
               Dybetax = gxx[ijk]*DyBetax + gxy[ijk]*DyBetay + gxz[ijk]*DyBetaz;
               Dybetay = gxy[ijk]*DyBetax + gyy[ijk]*DyBetay + gyz[ijk]*DyBetaz;
               Dybetaz = gxz[ijk]*DyBetax + gyz[ijk]*DyBetay + gzz[ijk]*DyBetaz;
               
               Dzbetax = gxx[ijk]*DzBetax + gxy[ijk]*DzBetay + gxz[ijk]*DzBetaz;
               Dzbetay = gxy[ijk]*DzBetax + gyy[ijk]*DzBetay + gyz[ijk]*DzBetaz;
               Dzbetaz = gxz[ijk]*DzBetax + gyz[ijk]*DzBetay + gzz[ijk]*DzBetaz;
               
               // get time derivative of metric
               gxx_dt[ijk] = -2.0*alp[ijk]*kxx[ijk] + Dxbetax + Dxbetax;
               gxy_dt[ijk] = -2.0*alp[ijk]*kxy[ijk] + Dxbetay + Dybetax;
               gxz_dt[ijk] = -2.0*alp[ijk]*kxz[ijk] + Dxbetaz + Dzbetax;
               gyy_dt[ijk] = -2.0*alp[ijk]*kyy[ijk] + Dybetay + Dybetay;
               gyz_dt[ijk] = -2.0*alp[ijk]*kyz[ijk] + Dybetaz + Dzbetay;
               gzz_dt[ijk] = -2.0*alp[ijk]*kzz[ijk] + Dzbetaz + Dzbetaz;
            }
 
 
         }
         LC_ENDLOOP3 (ADMDerivatives_CalcDerivatives);
         
      }
/*      else
      {
         // we need to do only radial Jacobian trafo since we have one radial coordinate direction which might be stretched!
   
         assert(*general_coordinates);
   
#pragma omp parallel
         LC_LOOP3 (ADMDerivatives_CalcRadialDerivatives,
               i, j, k,
               istart[0], istart[1], istart[2],
               iend[0], iend[1], iend[2],
               cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
         {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
   
            // FIXME: Get Jacobian right!
            const CCTK_REAL Jac = sqrt(J31[ijk]*J31[ijk] + J32[ijk]*J32[ijk] + J33[ijk]*J33[ijk]);
   
            alp_dr[ijk] = Jac * calc_dc(cctkGH, alp, i, j, k, nk, imin[2], imax[2], q[2], ihz);
            
            betax_dr[ijk] = Jac * calc_dc(cctkGH, betax, i, j, k, nk, imin[2], imax[2], q[2], ihz);
            betay_dr[ijk] = Jac * calc_dc(cctkGH, betay, i, j, k, nk, imin[2], imax[2], q[2], ihz);
            betaz_dr[ijk] = Jac * calc_dc(cctkGH, betaz, i, j, k, nk, imin[2], imax[2], q[2], ihz);
            
            gxx_dr[ijk] = Jac * calc_dc(cctkGH, gxx, i, j, k, nk, imin[2], imax[2], q[2], ihz);
            gxy_dr[ijk] = Jac * calc_dc(cctkGH, gxy, i, j, k, nk, imin[2], imax[2], q[2], ihz);
            gxz_dr[ijk] = Jac * calc_dc(cctkGH, gxz, i, j, k, nk, imin[2], imax[2], q[2], ihz);
            gyy_dr[ijk] = Jac * calc_dc(cctkGH, gyy, i, j, k, nk, imin[2], imax[2], q[2], ihz);
            gyz_dr[ijk] = Jac * calc_dc(cctkGH, gyz, i, j, k, nk, imin[2], imax[2], q[2], ihz);
            gzz_dr[ijk] = Jac * calc_dc(cctkGH, gzz, i, j, k, nk, imin[2], imax[2], q[2], ihz);
         }
         LC_ENDLOOP3 (ADMDerivatives_CalcRadialDerivatives);
      }*/
   }  



   delete [] imin[0];
   delete [] imin[1];
   delete [] imin[2];
   delete [] imax[0];
   delete [] imax[1];
   delete [] imax[2];
   delete [] q[0];
   delete [] q[1];
   delete [] q[2];
   
}
