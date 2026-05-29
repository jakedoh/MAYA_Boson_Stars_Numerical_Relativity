
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

#ifndef ALLDERIV_H
#define ALLDERIV_H


#include "cctk.h"

static inline void all_diff (
          const cGH* restrict const cctkGH,
          const CCTK_INT nvars,
          const CCTK_REAL* const vars[],
          CCTK_REAL dvars[][3], 
          const CCTK_INT i, const CCTK_INT j, const CCTK_INT k,
          const CCTK_INT ni, const CCTK_INT nj, const CCTK_INT nk,
          const CCTK_INT* restrict const imin,
          const CCTK_INT* restrict const imax,
          const CCTK_INT* restrict const jmin,
          const CCTK_INT* restrict const jmax,
          const CCTK_INT* restrict const kmin,
          const CCTK_INT* restrict const kmax,
          const CCTK_REAL* restrict const qa,
          const CCTK_REAL* restrict const qb,
          const CCTK_REAL* restrict const qc,
          const CCTK_REAL ha,
          const CCTK_REAL hb,
          const CCTK_REAL hc)
{
  for (int v=0; v<nvars; v++)
  {
    for (int l=0; l<3; l++)
    {
      dvars[v][l] = 0.0;
    }
  }

  for (int a=imin[i]; a<=imax[i]; ++a)
  {
    int index = CCTK_GFINDEX3D(cctkGH, a, j, k);
    for (int v=0; v<nvars; v++)
    {
      dvars[v][0] += qa[a+i*ni] * (vars[v])[index];
    }
  }

  for (int b=jmin[j]; b<=jmax[j]; ++b)
  {
    int index = CCTK_GFINDEX3D(cctkGH, i, b, k);
    for (int v=0; v<nvars; v++)
    {
      dvars[v][1] += qb[b+j*nj] * (vars[v])[index];
    }
  }

  for (int c=kmin[k]; c<=kmax[k]; ++c)
  {
    int index = CCTK_GFINDEX3D(cctkGH, i, j, c);
    for (int v=0; v<nvars; v++)
    {
      dvars[v][2] += qc[c+k*nk] * (vars[v])[index];
    }
  }               

  for (int v=0; v<nvars; v++)
  {
    dvars[v][0] = ha * dvars[v][0];
    dvars[v][1] = hb * dvars[v][1];
    dvars[v][2] = hc * dvars[v][2];
  }

}

static inline void all_diff2 (
          const cGH* restrict const cctkGH,
          const CCTK_INT nvars,
          const CCTK_REAL* const vars[],
          CCTK_REAL dvars[][6], 
          const CCTK_INT i, const CCTK_INT j, const CCTK_INT k,
          const CCTK_INT ni, const CCTK_INT nj, const CCTK_INT nk,
          const CCTK_INT* restrict const imin,
          const CCTK_INT* restrict const imax,
          const CCTK_INT* restrict const jmin,
          const CCTK_INT* restrict const jmax,
          const CCTK_INT* restrict const kmin,
          const CCTK_INT* restrict const kmax,
          const CCTK_REAL* restrict const qa,
          const CCTK_REAL* restrict const qb,
          const CCTK_REAL* restrict const qc,
          const CCTK_INT* restrict const imin2,
          const CCTK_INT* restrict const imax2,
          const CCTK_INT* restrict const jmin2,
          const CCTK_INT* restrict const jmax2,
          const CCTK_INT* restrict const kmin2,
          const CCTK_INT* restrict const kmax2,
          const CCTK_REAL* restrict const qa2,
          const CCTK_REAL* restrict const qb2,
          const CCTK_REAL* restrict const qc2,
          const CCTK_REAL hhaa,
          const CCTK_REAL hhab,
          const CCTK_REAL hhac,
          const CCTK_REAL hhbb,
          const CCTK_REAL hhbc,
          const CCTK_REAL hhcc)
{
  for (int v=0; v<nvars; v++)
  {
    for (int l=0; l<6; l++)
    {
      dvars[v][l] = 0.0;
    }
  }

  for (int a=imin2[i]; a<=imax2[i]; ++a)
  {
    int index = CCTK_GFINDEX3D(cctkGH, a, j, k);
    for (int v=0; v<nvars; v++)
    {
      dvars[v][0] += qa2[a+i*ni] * (vars[v])[index];
    }
  }

  for (int b=jmin[j]; b<=jmax[j]; ++b)
  {
    CCTK_REAL tmp[nvars];
    for (int v=0; v<nvars; v++)
    {
      tmp[v] = 0.0;
    }
    for (int a=imin[i]; a<=imax[i]; ++a)
    {
      int index = CCTK_GFINDEX3D(cctkGH, a, b, k);
      for (int v=0; v<nvars; v++)
      {
        tmp[v] += qa[a+i*ni] * (vars[v])[index];
      }
    }
    for (int v=0; v<nvars; v++)
    {
      dvars[v][1] += qb[b+j*nj] * tmp[v];
    }
  }

  for (int c=kmin[k]; c<=kmax[k]; ++c)
  {
    CCTK_REAL tmp[nvars];
    for (int v=0; v<nvars; v++)
    {
      tmp[v] = 0.0;
    }
    for (int a=imin[i]; a<=imax[i]; ++a)
    {
      int index = CCTK_GFINDEX3D(cctkGH, a, j, c);
      for (int v=0; v<nvars; v++)
      {
        tmp[v] += qa[a+i*ni] * (vars[v])[index];
      }
    }
    for (int v=0; v<nvars; v++)
    {
      dvars[v][2] += qc[c+k*nk] * tmp[v];
    }
  }

  for (int b=jmin2[j]; b<=jmax2[j]; ++b)
  {
    int index = CCTK_GFINDEX3D(cctkGH, i, b, k);
    for (int v=0; v<nvars; v++)
    {
      dvars[v][3] += qb2[b+j*nj] * (vars[v])[index];
    }
  }

  for (int c=kmin[k]; c<=kmax[k]; ++c)
  {
    CCTK_REAL tmp[nvars];
    for (int v=0; v<nvars; v++)
    {
      tmp[v] = 0.0;
    }
    for (int b=jmin[j]; b<=jmax[j]; ++b)
    {
      int index = CCTK_GFINDEX3D(cctkGH, i, b, c);
      for (int v=0; v<nvars; v++)
      {
        tmp[v] += qb[b+j*nj] * (vars[v])[index];
      }
    }
    for (int v=0; v<nvars; v++)
    {
      dvars[v][4] += qc[c+k*nk] * tmp[v];
    }
  }

  for (int c=kmin2[k]; c<=kmax2[k]; ++c)
  {
    int index = CCTK_GFINDEX3D(cctkGH, i, j, c);
    for (int v=0; v<nvars; v++)
    {
      dvars[v][5] += qc2[c+k*nk] * (vars[v])[index];
    }
  }

  for (int v=0; v<nvars; v++)
  {
    dvars[v][0] = hhaa * dvars[v][0];
    dvars[v][1] = hhab * dvars[v][1];
    dvars[v][2] = hhac * dvars[v][2];
    dvars[v][3] = hhbb * dvars[v][3];
    dvars[v][4] = hhbc * dvars[v][4];
    dvars[v][5] = hhcc * dvars[v][5];
  }
}

#endif
