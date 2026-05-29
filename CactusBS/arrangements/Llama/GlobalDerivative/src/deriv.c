
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
 *  Applies a global first derivative to an entire gridvariable
 */
void globalDiffGv(const CCTK_POINTER_TO_CONST cctkGH_, const CCTK_INT dir,
                  const CCTK_REAL *var, CCTK_REAL *dvar,
                  const CCTK_REAL* restrict const J_dadx,
		  const CCTK_REAL* restrict const J_dbdx,
		  const CCTK_REAL* restrict const J_dcdx,
		  const CCTK_REAL* restrict const J_dady,
		  const CCTK_REAL* restrict const J_dbdy,
		  const CCTK_REAL* restrict const J_dcdy,
                  const CCTK_REAL* restrict const J_dadz,
		  const CCTK_REAL* restrict const J_dbdz,
		  const CCTK_REAL* restrict const J_dcdz,
		  const CCTK_INT table_handle)
{
  cGH const * restrict const cctkGH = cctkGH_;
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;
//  int i, j, k, ijk;
  int ni, nj, nk;
   
  CCTK_INT *imin, *jmin, *kmin;
  CCTK_INT *imax, *jmax, *kmax;
  CCTK_REAL *qa, *qb, *qc;
  CCTK_REAL iha, ihb, ihc;
   
  ni = cctk_lsh[0];
  nj = cctk_lsh[1];
  nk = cctk_lsh[2];
   
  /*
   * Grid spacings required by the finite difference operators.
   */
  iha = 1.0 / CCTK_DELTA_SPACE(0);
  ihb = 1.0 / CCTK_DELTA_SPACE(1);
  ihc = 1.0 / CCTK_DELTA_SPACE(2);
   
  /*
   * Call SummationByParts for finite-difference operator coefficients.
   */
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
#pragma omp for
  for (int k=0; k<nk; ++k)
    {
      kmin[k] -= 1;
      kmax[k] -= 1;
    }
  }
   
   
  switch (dir)
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
                  
		dvar[ijk] = g_diff_dx(cctkGH, var, *general_coordinates,
				      J_dadx[ijk], J_dbdx[ijk], J_dcdx[ijk],
				      i, j, k, ni, nj, nk, 
				      imin, imax, jmin, jmax, kmin, kmax,
				      qa, qb, qc, iha, ihb, ihc);
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
                  
		dvar[ijk] = g_diff_dy(cctkGH, var, *general_coordinates,
				      J_dady[ijk], J_dbdy[ijk], J_dcdy[ijk],
				      i, j, k, ni, nj, nk, 
				      imin, imax, jmin, jmax, kmin, kmax,
				      qa, qb, qc, iha, ihb, ihc);
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
                  
		dvar[ijk] = g_diff_dz(cctkGH, var, *general_coordinates,
				      J_dadz[ijk], J_dbdz[ijk], J_dcdz[ijk],
				      i, j, k, ni, nj, nk, 
				      imin, imax, jmin, jmax, kmin, kmax,
				      qa, qb, qc, iha, ihb, ihc);
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
}
