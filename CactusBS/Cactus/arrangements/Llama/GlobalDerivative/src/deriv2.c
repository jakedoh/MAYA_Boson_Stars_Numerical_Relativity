
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

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Functions.h>
#include <cctk_Parameters.h>
#include <stdlib.h>

#include "GlobalDerivative.h"

/*
 * Applies a global second derivative to an entire gridvariable
 */
void globalDiff2Gv(const CCTK_POINTER_TO_CONST cctkGH_, CCTK_INT dir1,
		   CCTK_INT dir2, const CCTK_REAL* restrict const var,
		   CCTK_REAL* dvar,
                   const CCTK_REAL* restrict const J_dadx,
		   const CCTK_REAL* restrict const J_dbdx,
		   const CCTK_REAL* restrict const J_dcdx,
		   const CCTK_REAL* restrict const J_dady,
		   const CCTK_REAL* restrict const J_dbdy,
		   const CCTK_REAL* restrict const J_dcdy,
                   const CCTK_REAL* restrict const J_dadz,
		   const CCTK_REAL* restrict const J_dbdz,
		   const CCTK_REAL* restrict const J_dcdz,
		   const CCTK_REAL* restrict const dJ_dadxdx,
		   const CCTK_REAL* restrict const dJ_dbdxdx,
		   const CCTK_REAL* restrict const dJ_dcdxdx,
                   const CCTK_REAL* restrict const dJ_dadxdy,
		   const CCTK_REAL* restrict const dJ_dbdxdy,
		   const CCTK_REAL* restrict const dJ_dcdxdy,
		   const CCTK_REAL* restrict const dJ_dadxdz,
		   const CCTK_REAL* restrict const dJ_dbdxdz,
		   const CCTK_REAL* restrict const dJ_dcdxdz,
                   const CCTK_REAL* restrict const dJ_dadydy,
		   const CCTK_REAL* restrict const dJ_dbdydy,
		   const CCTK_REAL* restrict const dJ_dcdydy,
		   const CCTK_REAL* restrict const dJ_dadydz,
		   const CCTK_REAL* restrict const dJ_dbdydz,
		   const CCTK_REAL* restrict const dJ_dcdydz,
                   const CCTK_REAL* restrict const dJ_dadzdz,
		   const CCTK_REAL* restrict const dJ_dbdzdz,
		   const CCTK_REAL* restrict const dJ_dcdzdz,
		   const CCTK_INT table_handle)
{
  cGH const * restrict const cctkGH = cctkGH_;
  DECLARE_CCTK_PARAMETERS
    DECLARE_CCTK_ARGUMENTS
  int ni, nj, nk;
   
  CCTK_INT *imin, *jmin, *kmin;
  CCTK_INT *imax, *jmax, *kmax;
  CCTK_REAL *qa, *qb, *qc;
  CCTK_INT *imin2, *jmin2, *kmin2;
  CCTK_INT *imax2, *jmax2, *kmax2;
  CCTK_REAL *qa2, *qb2, *qc2;
  CCTK_REAL iha, ihb, ihc, ihaa, ihab, ihac, ihbb, ihbc, ihcc;
   
   
  ni = cctk_lsh[0];
  nj = cctk_lsh[1];
  nk = cctk_lsh[2];
   
  /*
   * Grid spacings required by the finite difference operators.
   */
  iha = 1.0 / CCTK_DELTA_SPACE(0);
  ihb = 1.0 / CCTK_DELTA_SPACE(1);
  ihc = 1.0 / CCTK_DELTA_SPACE(2);
   
  ihaa = iha*iha;
  ihab = iha*ihb;
  ihac = iha*ihc;
  ihbb = ihb*ihb;
  ihbc = ihb*ihc;
  ihcc = ihc*ihc;
   
  /*
   * Call SummationByParts for finite-difference operator coefficients.
   */
  // first derivatives
   
  imin = (CCTK_INT *) malloc(ni*sizeof(CCTK_INT));
  jmin = (CCTK_INT *) malloc(nj*sizeof(CCTK_INT));
  kmin = (CCTK_INT *) malloc(nk*sizeof(CCTK_INT));
   
  imax = (CCTK_INT *) malloc(ni*sizeof(CCTK_INT));
  jmax = (CCTK_INT *) malloc(nj*sizeof(CCTK_INT));
  kmax = (CCTK_INT *) malloc(nk*sizeof(CCTK_INT));
   
  qa = (CCTK_REAL *) malloc(ni*ni*sizeof(CCTK_REAL));
  qb = (CCTK_REAL *) malloc(nj*nj*sizeof(CCTK_REAL));
  qc = (CCTK_REAL *) malloc(nk*nk*sizeof(CCTK_REAL));
   
  Diff_coeff(cctkGH, 0, ni, imin, imax, qa, -1);
  Diff_coeff(cctkGH, 1, nj, jmin, jmax, qb, -1);
  Diff_coeff(cctkGH, 2, nk, kmin, kmax, qc, -1);
   
  // second derivatives
   
  imin2 = (CCTK_INT *) malloc(ni*sizeof(CCTK_INT));
  jmin2 = (CCTK_INT *) malloc(nj*sizeof(CCTK_INT));
  kmin2 = (CCTK_INT *) malloc(nk*sizeof(CCTK_INT));
   
  imax2 = (CCTK_INT *) malloc(ni*sizeof(CCTK_INT));
  jmax2 = (CCTK_INT *) malloc(nj*sizeof(CCTK_INT));
  kmax2 = (CCTK_INT *) malloc(nk*sizeof(CCTK_INT));
   
  qa2 = (CCTK_REAL *) malloc(ni*ni*sizeof(CCTK_REAL));
  qb2 = (CCTK_REAL *) malloc(nj*nj*sizeof(CCTK_REAL));
  qc2 = (CCTK_REAL *) malloc(nk*nk*sizeof(CCTK_REAL));
   
  Diff2_coeff(cctkGH, 0, ni, imin2, imax2, qa2, -1);
  Diff2_coeff(cctkGH, 1, nj, jmin2, jmax2, qb2, -1);
  Diff2_coeff(cctkGH, 2, nk, kmin2, kmax2, qc2, -1);
   
  /*
   * Convert Fortran->C indexing in the returned Diff_coeff arrays.
   */
#pragma omp parallel
  {
#pragma omp for nowait
  for (int i=0; i<ni; ++i)
    {
      imin[i] -= 1;
      imax[i] -= 1;
    }
#pragma omp for nowait
  for (int j=0; j<nj; ++j)
    {
      jmin[j] -= 1;
      jmax[j] -= 1;
    }
#pragma omp for nowait
  for (int k=0; k<nk; ++k)
    {
      kmin[k] -= 1;
      kmax[k] -= 1;
    }
#pragma omp for nowait
  for (int i=0; i<ni; ++i)
    {
      imin2[i] -= 1;
      imax2[i] -= 1;
    }
#pragma omp for nowait
  for (int j=0; j<nj; ++j)
    {
      jmin2[j] -= 1;
      jmax2[j] -= 1;
    }
#pragma omp for
  for (int k=0; k<nk; ++k)
    {
      kmin2[k] -= 1;
      kmax2[k] -= 1;
    }
  } 
   
  if (dir2 < dir1)
    {
      // swap directions due to symmetry 
   
      int d = dir2;
      dir2 = dir1;
      dir1 = d;
    }
   
   
  switch (dir1)
    {
    case 0:
      {
      
	switch (dir2)
	  {
	  case 0:
            {
	      // loop over all points and create derivative
#pragma omp parallel for collapse(3)
	      for (int k=0; k<nk; ++k)
		for (int j=0; j<nj; ++j)
		  for (int i=0; i<ni; ++i)
		    {
		      int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
                        
		      dvar[ijk] = g_diff_dxdx(cctkGH, var,
					      *general_coordinates,
					      J_dadx[ijk], J_dbdx[ijk],
					      J_dcdx[ijk], dJ_dadxdx[ijk],
					      dJ_dbdxdx[ijk], dJ_dcdxdx[ijk],
					      i, j, k, ni, nj, nk, 
					      imin, imax, jmin, jmax, kmin,
					      kmax, qa, qb, qc, imin2, imax2,
					      jmin2, jmax2, kmin2, kmax2,
					      qa2, qb2, qc2, iha, ihb, ihc,
					      ihaa, ihab, ihac, ihbb, ihbc,
					      ihcc);
		    }
	      break;
            }
            
	  case 1:
            {
	      // loop over all points and create derivative
#pragma omp parallel for collapse(3)
	      for (int k=0; k<nk; ++k)
		for (int j=0; j<nj; ++j)
		  for (int i=0; i<ni; ++i)
		    {
		      int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
                        
		      dvar[ijk] = g_diff_dxdy(cctkGH, var,
					      *general_coordinates, 
					      J_dadx[ijk], J_dbdx[ijk],
					      J_dcdx[ijk], J_dady[ijk],
					      J_dbdy[ijk], J_dcdy[ijk], 
					      dJ_dadxdy[ijk], dJ_dbdxdy[ijk],
					      dJ_dcdxdy[ijk],
					      i, j, k, ni, nj, nk, 
					      imin, imax, jmin, jmax,
					      kmin, kmax, qa, qb, qc,
					      imin2, imax2, jmin2, jmax2,
					      kmin2, kmax2, qa2, qb2, qc2,
					      iha, ihb, ihc, ihaa, ihab, ihac,
					      ihbb, ihbc, ihcc);
		    }
	      break;
            }
            
	  case 2:
            {
	      // loop over all points and create derivative
#pragma omp parallel for collapse(3)
	      for (int k=0; k<nk; ++k)
		for (int j=0; j<nj; ++j)
		  for (int i=0; i<ni; ++i)
		    {
		      int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
                        
		      dvar[ijk] = g_diff_dxdz(cctkGH, var,
					      *general_coordinates,
					      J_dadx[ijk], J_dbdx[ijk],
					      J_dcdx[ijk], J_dadz[ijk],
					      J_dbdz[ijk], J_dcdz[ijk], 
					      dJ_dadxdz[ijk], dJ_dbdxdz[ijk],
					      dJ_dcdxdz[ijk],
					      i, j, k, ni, nj, nk, 
					      imin, imax, jmin, jmax,
					      kmin, kmax, qa, qb, qc,
					      imin2, imax2, jmin2, jmax2,
					      kmin2, kmax2, qa2, qb2, qc2,
					      iha, ihb, ihc, ihaa, ihab, ihac,
					      ihbb, ihbc, ihcc);
		    }
	      break;
            }
	  }
	break;
      }
      
    case 1:
      {
      
	switch (dir2)
	  {
	  case 1:
            {
	      // loop over all points and create derivative
#pragma omp parallel for collapse(3)
	      for (int k=0; k<nk; ++k)
		for (int j=0; j<nj; ++j)
		  for (int i=0; i<ni; ++i)
		    {
		      int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
                        
		      dvar[ijk] = g_diff_dydy(cctkGH, var,
					      *general_coordinates,
					      J_dady[ijk], J_dbdy[ijk],
					      J_dcdy[ijk], dJ_dadydy[ijk],
					      dJ_dbdydy[ijk], dJ_dcdydy[ijk],
					      i, j, k, ni, nj, nk, 
					      imin, imax, jmin, jmax,
					      kmin, kmax, qa, qb, qc,
					      imin2, imax2, jmin2, jmax2,
					      kmin2, kmax2, qa2, qb2, qc2,
					      iha, ihb, ihc, ihaa, ihab, ihac,
					      ihbb, ihbc, ihcc);
		    }
	      break;
            }
            
	  case 2:
            {
	      // loop over all points and create derivative
#pragma omp parallel for collapse(3)
	      for (int k=0; k<nk; ++k)
		for (int j=0; j<nj; ++j)
		  for (int i=0; i<ni; ++i)
		    {
		      int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
                        
		      dvar[ijk] = g_diff_dydz(cctkGH, var,
					      *general_coordinates,
					      J_dady[ijk], J_dbdy[ijk],
					      J_dcdy[ijk], J_dadz[ijk],
					      J_dbdz[ijk], J_dcdz[ijk], 
					      dJ_dadydz[ijk], dJ_dbdydz[ijk],
					      dJ_dcdydz[ijk],
					      i, j, k, ni, nj, nk, 
					      imin, imax, jmin, jmax,
					      kmin, kmax, qa, qb, qc,
					      imin2, imax2, jmin2, jmax2,
					      kmin2, kmax2, qa2, qb2, qc2,
					      iha, ihb, ihc, ihaa, ihab, ihac,
					      ihbb, ihbc, ihcc);
		    }
	      break;
            }
            
            
	  }
	break;
      }
      
    case 2:
      {
	switch (dir2)
	  {
	  case 2:
            {
	      // loop over all points and create derivative
#pragma omp parallel for collapse(3)
	      for (int k=0; k<nk; ++k)
		for (int j=0; j<nj; ++j)
		  for (int i=0; i<ni; ++i)
		    {
		      int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
                        
		      dvar[ijk] = g_diff_dzdz(cctkGH, var,
					      *general_coordinates,
					      J_dadz[ijk], J_dbdz[ijk],
					      J_dcdz[ijk], dJ_dadzdz[ijk],
					      dJ_dbdzdz[ijk], dJ_dcdzdz[ijk],
					      i, j, k, ni, nj, nk, 
					      imin, imax, jmin, jmax,
					      kmin, kmax, qa, qb, qc,
					      imin2, imax2, jmin2, jmax2,
					      kmin2, kmax2, qa2, qb2, qc2,
					      iha, ihb, ihc, ihaa, ihab, ihac,
					      ihbb, ihbc, ihcc);
		    }
	      break;
            }
	  }
	break;
      }
    }
   
  free(imin);
  free(jmin);
  free(kmin);
   
  free(imax);
  free(jmax);
  free(kmax);
   
  free(qa);
  free(qb);
  free(qc);
   
  free(imin2);
  free(jmin2);
  free(kmin2);
   
  free(imax2);
  free(jmax2);
  free(kmax2);
   
  free(qa2);
  free(qb2);
  free(qc2);
}









