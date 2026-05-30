
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

#ifndef ALLDERIV_8_H
#define ALLDERIV_8_H


#include "cctk.h"

#define UPWIND_SMALL 1e-10

static inline void all_diff_8 (
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
          const CCTK_INT* restrict const bsize,
          const CCTK_REAL ha,
          const CCTK_REAL hb,
          const CCTK_REAL hc)
{
  if (i >= bsize[0] && i < cctkGH->cctk_lsh[0]-bsize[1])
  {
    CCTK_REAL hafac = ha*0.0011904761904761904762;  // 1/840
    CCTK_INT indip4 = CCTK_GFINDEX3D(cctkGH, i+4, j, k);
    CCTK_INT indip3 = CCTK_GFINDEX3D(cctkGH, i+3, j, k);
    CCTK_INT indip2 = CCTK_GFINDEX3D(cctkGH, i+2, j, k);
    CCTK_INT indip1 = CCTK_GFINDEX3D(cctkGH, i+1, j, k);
    CCTK_INT indim1 = CCTK_GFINDEX3D(cctkGH, i-1, j, k);
    CCTK_INT indim2 = CCTK_GFINDEX3D(cctkGH, i-2, j, k);
    CCTK_INT indim3 = CCTK_GFINDEX3D(cctkGH, i-3, j, k);
    CCTK_INT indim4 = CCTK_GFINDEX3D(cctkGH, i-4, j, k);
    for (int v=0; v<nvars; v++)
    {
      dvars[v][0] = hafac*( -3.0 *((vars[v])[indip4]
                                 - (vars[v])[indim4])
                           +32.0 *((vars[v])[indip3]
                                 - (vars[v])[indim3])
                          -168.0 *((vars[v])[indip2]
                                 - (vars[v])[indim2])
                          +672.0 *((vars[v])[indip1]
                                 - (vars[v])[indim1])); 
    }
  } else {

    for (int v=0; v<nvars; v++)
    {
      dvars[v][0] = 0.0;
    }

    for (int a=imin[i]; a<=imax[i]; ++a)
    {
      int index = CCTK_GFINDEX3D(cctkGH, a, j, k);
      for (int v=0; v<nvars; v++)
      {
        dvars[v][0] += qa[a+i*ni] * (vars[v])[index];
      }
    }
    for (int v=0; v<nvars; v++)
    {
      dvars[v][0] = ha * dvars[v][0];
    }

  }

  if (j >= bsize[2] && j < cctkGH->cctk_lsh[1]-bsize[3])
  {
    CCTK_REAL hbfac = hb*0.0011904761904761904762;  // 1/840
    CCTK_INT indjp4 = CCTK_GFINDEX3D(cctkGH, i, j+4, k);
    CCTK_INT indjp3 = CCTK_GFINDEX3D(cctkGH, i, j+3, k);
    CCTK_INT indjp2 = CCTK_GFINDEX3D(cctkGH, i, j+2, k);
    CCTK_INT indjp1 = CCTK_GFINDEX3D(cctkGH, i, j+1, k);
    CCTK_INT indjm1 = CCTK_GFINDEX3D(cctkGH, i, j-1, k);
    CCTK_INT indjm2 = CCTK_GFINDEX3D(cctkGH, i, j-2, k);
    CCTK_INT indjm3 = CCTK_GFINDEX3D(cctkGH, i, j-3, k);
    CCTK_INT indjm4 = CCTK_GFINDEX3D(cctkGH, i, j-4, k);
    for (int v=0; v<nvars; v++)
    {
      dvars[v][1] = hbfac*( -3.0 *((vars[v])[indjp4]
                                 - (vars[v])[indjm4])
                           +32.0 *((vars[v])[indjp3]
                                 - (vars[v])[indjm3])
                          -168.0 *((vars[v])[indjp2]
                                 - (vars[v])[indjm2])
                          +672.0 *((vars[v])[indjp1]
                                 - (vars[v])[indjm1]));
    }
  } else {

    for (int v=0; v<nvars; v++)
    {
      dvars[v][1] = 0.0;
    }
    for (int b=jmin[j]; b<=jmax[j]; ++b)
    {
      int index = CCTK_GFINDEX3D(cctkGH, i, b, k);
      for (int v=0; v<nvars; v++)
      {
        dvars[v][1] += qb[b+j*nj] * (vars[v])[index];
      }
    }
    for (int v=0; v<nvars; v++)
    {
      dvars[v][1] = hb * dvars[v][1];
    }

  }

  if (k >= bsize[4] && k < cctkGH->cctk_lsh[2]-bsize[5])
  {
    CCTK_REAL hcfac = hc*0.0011904761904761904762;  // 1/840
    CCTK_INT indkp4 = CCTK_GFINDEX3D(cctkGH, i, j, k+4);
    CCTK_INT indkp3 = CCTK_GFINDEX3D(cctkGH, i, j, k+3);
    CCTK_INT indkp2 = CCTK_GFINDEX3D(cctkGH, i, j, k+2);
    CCTK_INT indkp1 = CCTK_GFINDEX3D(cctkGH, i, j, k+1);
    CCTK_INT indkm1 = CCTK_GFINDEX3D(cctkGH, i, j, k-1);
    CCTK_INT indkm2 = CCTK_GFINDEX3D(cctkGH, i, j, k-2);
    CCTK_INT indkm3 = CCTK_GFINDEX3D(cctkGH, i, j, k-3);
    CCTK_INT indkm4 = CCTK_GFINDEX3D(cctkGH, i, j, k-4);
    for (int v=0; v<nvars; v++)
    {
      dvars[v][2] = hcfac*( -3.0 *((vars[v])[indkp4]
                                 - (vars[v])[indkm4])
                           +32.0 *((vars[v])[indkp3]
                                 - (vars[v])[indkm3])
                          -168.0 *((vars[v])[indkp2]
                                 - (vars[v])[indkm2])
                          +672.0 *((vars[v])[indkp1]
                                 - (vars[v])[indkm1]));
    }
  } else {

    for (int v=0; v<nvars; v++)
    {
      dvars[v][2] = 0.0;
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
      dvars[v][2] = hc * dvars[v][2];
    }
  }

}

static inline void all_diff2_8 (
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
          const CCTK_INT* restrict const bsize,
          const CCTK_REAL hhaa,
          const CCTK_REAL hhab,
          const CCTK_REAL hhac,
          const CCTK_REAL hhbb,
          const CCTK_REAL hhbc,
          const CCTK_REAL hhcc)
{

  if (i >= 4 && i < cctkGH->cctk_lsh[0]-4)
  {
    CCTK_REAL hhaafac = hhaa*0.00019841269841269841270; // 1/5040
    CCTK_INT indip4 = CCTK_GFINDEX3D(cctkGH, i+4, j, k);
    CCTK_INT indip3 = CCTK_GFINDEX3D(cctkGH, i+3, j, k);
    CCTK_INT indip2 = CCTK_GFINDEX3D(cctkGH, i+2, j, k);
    CCTK_INT indip1 = CCTK_GFINDEX3D(cctkGH, i+1, j, k);
    CCTK_INT ind    = CCTK_GFINDEX3D(cctkGH, i, j, k);
    CCTK_INT indim1 = CCTK_GFINDEX3D(cctkGH, i-1, j, k);
    CCTK_INT indim2 = CCTK_GFINDEX3D(cctkGH, i-2, j, k);
    CCTK_INT indim3 = CCTK_GFINDEX3D(cctkGH, i-3, j, k);
    CCTK_INT indim4 = CCTK_GFINDEX3D(cctkGH, i-4, j, k);
    for (int v=0; v<nvars; v++)
    {
      dvars[v][0] = hhaafac*(  - 9.0 * ((vars[v])[indip4]
                                      + (vars[v])[indim4])
                            +  128.0 * ((vars[v])[indip3]
                                      + (vars[v])[indim3])
                            - 1008.0 * ((vars[v])[indip2]
                                      + (vars[v])[indim2])
                            + 8064.0 * ((vars[v])[indip1]
                                      + (vars[v])[indim1])
                            - 14350.0 * (vars[v])[ind]
                            );
    }
  } else {

    for (int v=0; v<nvars; v++)
    {
      dvars[v][0] = 0.0;
    }

    for (int a=imin2[i]; a<=imax2[i]; ++a)
    {
      int index = CCTK_GFINDEX3D(cctkGH, a, j, k);
      for (int v=0; v<nvars; v++)
      {
        dvars[v][0] += qa2[a+i*ni] * (vars[v])[index];
      }
    }

    for (int v=0; v<nvars; v++)
    {
      dvars[v][0] = hhaa * dvars[v][0];
    }

  }

  if (i >= bsize[0] && i < cctkGH->cctk_lsh[0]-bsize[1] &&
      j >= bsize[2] && j < cctkGH->cctk_lsh[1]-bsize[3])
  {
    CCTK_REAL hhabfac = hhab*0.0000014172335600907029478; // 1/705600
    CCTK_INT indip4jp4 = CCTK_GFINDEX3D(cctkGH, i+4, j+4, k);
    CCTK_INT indip4jp3 = CCTK_GFINDEX3D(cctkGH, i+4, j+3, k);
    CCTK_INT indip4jp2 = CCTK_GFINDEX3D(cctkGH, i+4, j+2, k);
    CCTK_INT indip4jp1 = CCTK_GFINDEX3D(cctkGH, i+4, j+1, k);
    CCTK_INT indip4jm1 = CCTK_GFINDEX3D(cctkGH, i+4, j-1, k);
    CCTK_INT indip4jm2 = CCTK_GFINDEX3D(cctkGH, i+4, j-2, k);
    CCTK_INT indip4jm3 = CCTK_GFINDEX3D(cctkGH, i+4, j-3, k);
    CCTK_INT indip4jm4 = CCTK_GFINDEX3D(cctkGH, i+4, j-4, k);
    CCTK_INT indip3jp4 = CCTK_GFINDEX3D(cctkGH, i+3, j+4, k);
    CCTK_INT indip3jp3 = CCTK_GFINDEX3D(cctkGH, i+3, j+3, k);
    CCTK_INT indip3jp2 = CCTK_GFINDEX3D(cctkGH, i+3, j+2, k);
    CCTK_INT indip3jp1 = CCTK_GFINDEX3D(cctkGH, i+3, j+1, k);
    CCTK_INT indip3jm1 = CCTK_GFINDEX3D(cctkGH, i+3, j-1, k);
    CCTK_INT indip3jm2 = CCTK_GFINDEX3D(cctkGH, i+3, j-2, k);
    CCTK_INT indip3jm3 = CCTK_GFINDEX3D(cctkGH, i+3, j-3, k);
    CCTK_INT indip3jm4 = CCTK_GFINDEX3D(cctkGH, i+3, j-4, k);
    CCTK_INT indip2jp4 = CCTK_GFINDEX3D(cctkGH, i+2, j+4, k);
    CCTK_INT indip2jp3 = CCTK_GFINDEX3D(cctkGH, i+2, j+3, k);
    CCTK_INT indip2jp2 = CCTK_GFINDEX3D(cctkGH, i+2, j+2, k);
    CCTK_INT indip2jp1 = CCTK_GFINDEX3D(cctkGH, i+2, j+1, k);
    CCTK_INT indip2jm1 = CCTK_GFINDEX3D(cctkGH, i+2, j-1, k);
    CCTK_INT indip2jm2 = CCTK_GFINDEX3D(cctkGH, i+2, j-2, k);
    CCTK_INT indip2jm3 = CCTK_GFINDEX3D(cctkGH, i+2, j-3, k);
    CCTK_INT indip2jm4 = CCTK_GFINDEX3D(cctkGH, i+2, j-4, k);
    CCTK_INT indip1jp4 = CCTK_GFINDEX3D(cctkGH, i+1, j+4, k);
    CCTK_INT indip1jp3 = CCTK_GFINDEX3D(cctkGH, i+1, j+3, k);
    CCTK_INT indip1jp2 = CCTK_GFINDEX3D(cctkGH, i+1, j+2, k);
    CCTK_INT indip1jp1 = CCTK_GFINDEX3D(cctkGH, i+1, j+1, k);
    CCTK_INT indip1jm1 = CCTK_GFINDEX3D(cctkGH, i+1, j-1, k);
    CCTK_INT indip1jm2 = CCTK_GFINDEX3D(cctkGH, i+1, j-2, k);
    CCTK_INT indip1jm3 = CCTK_GFINDEX3D(cctkGH, i+1, j-3, k);
    CCTK_INT indip1jm4 = CCTK_GFINDEX3D(cctkGH, i+1, j-4, k);
    CCTK_INT indim1jp4 = CCTK_GFINDEX3D(cctkGH, i-1, j+4, k);
    CCTK_INT indim1jp3 = CCTK_GFINDEX3D(cctkGH, i-1, j+3, k);
    CCTK_INT indim1jp2 = CCTK_GFINDEX3D(cctkGH, i-1, j+2, k);
    CCTK_INT indim1jp1 = CCTK_GFINDEX3D(cctkGH, i-1, j+1, k);
    CCTK_INT indim1jm1 = CCTK_GFINDEX3D(cctkGH, i-1, j-1, k);
    CCTK_INT indim1jm2 = CCTK_GFINDEX3D(cctkGH, i-1, j-2, k);
    CCTK_INT indim1jm3 = CCTK_GFINDEX3D(cctkGH, i-1, j-3, k);
    CCTK_INT indim1jm4 = CCTK_GFINDEX3D(cctkGH, i-1, j-4, k);
    CCTK_INT indim2jp4 = CCTK_GFINDEX3D(cctkGH, i-2, j+4, k);
    CCTK_INT indim2jp3 = CCTK_GFINDEX3D(cctkGH, i-2, j+3, k);
    CCTK_INT indim2jp2 = CCTK_GFINDEX3D(cctkGH, i-2, j+2, k);
    CCTK_INT indim2jp1 = CCTK_GFINDEX3D(cctkGH, i-2, j+1, k);
    CCTK_INT indim2jm1 = CCTK_GFINDEX3D(cctkGH, i-2, j-1, k);
    CCTK_INT indim2jm2 = CCTK_GFINDEX3D(cctkGH, i-2, j-2, k);
    CCTK_INT indim2jm3 = CCTK_GFINDEX3D(cctkGH, i-2, j-3, k);
    CCTK_INT indim2jm4 = CCTK_GFINDEX3D(cctkGH, i-2, j-4, k);
    CCTK_INT indim3jp4 = CCTK_GFINDEX3D(cctkGH, i-3, j+4, k);
    CCTK_INT indim3jp3 = CCTK_GFINDEX3D(cctkGH, i-3, j+3, k);
    CCTK_INT indim3jp2 = CCTK_GFINDEX3D(cctkGH, i-3, j+2, k);
    CCTK_INT indim3jp1 = CCTK_GFINDEX3D(cctkGH, i-3, j+1, k);
    CCTK_INT indim3jm1 = CCTK_GFINDEX3D(cctkGH, i-3, j-1, k);
    CCTK_INT indim3jm2 = CCTK_GFINDEX3D(cctkGH, i-3, j-2, k);
    CCTK_INT indim3jm3 = CCTK_GFINDEX3D(cctkGH, i-3, j-3, k);
    CCTK_INT indim3jm4 = CCTK_GFINDEX3D(cctkGH, i-3, j-4, k);
    CCTK_INT indim4jp4 = CCTK_GFINDEX3D(cctkGH, i-4, j+4, k);
    CCTK_INT indim4jp3 = CCTK_GFINDEX3D(cctkGH, i-4, j+3, k);
    CCTK_INT indim4jp2 = CCTK_GFINDEX3D(cctkGH, i-4, j+2, k);
    CCTK_INT indim4jp1 = CCTK_GFINDEX3D(cctkGH, i-4, j+1, k);
    CCTK_INT indim4jm1 = CCTK_GFINDEX3D(cctkGH, i-4, j-1, k);
    CCTK_INT indim4jm2 = CCTK_GFINDEX3D(cctkGH, i-4, j-2, k);
    CCTK_INT indim4jm3 = CCTK_GFINDEX3D(cctkGH, i-4, j-3, k);
    CCTK_INT indim4jm4 = CCTK_GFINDEX3D(cctkGH, i-4, j-4, k);
    for (int v=0; v<nvars; v++)
    {
      dvars[v][1] = hhabfac*(    9.0 * ((vars[v])[indip4jp4]
                                      + (vars[v])[indim4jm4]
                                      - (vars[v])[indip4jm4]
                                      - (vars[v])[indim4jp4])
                              - 96.0 * ((vars[v])[indip4jp3]
                                      + (vars[v])[indip3jp4]
                                      - (vars[v])[indip4jm3]
                                      - (vars[v])[indip3jm4]
                                      - (vars[v])[indim3jp4]
                                      - (vars[v])[indim4jp3]
                                      + (vars[v])[indim3jm4]
                                      + (vars[v])[indim4jm3])
                            + 1024.0 * ((vars[v])[indip3jp3]
                                      + (vars[v])[indim3jm3]
                                      - (vars[v])[indip3jm3]
                                      - (vars[v])[indim3jp3])
                             + 504.0 * ((vars[v])[indip4jp2]
                                      + (vars[v])[indip2jp4]
                                      - (vars[v])[indip4jm2]
                                      - (vars[v])[indip2jm4]
                                      - (vars[v])[indim2jp4]
                                      - (vars[v])[indim4jp2]
                                      + (vars[v])[indim2jm4]
                                      + (vars[v])[indim4jm2])
                            - 5376.0 * ((vars[v])[indip3jp2]
                                      + (vars[v])[indip2jp3]
                                      - (vars[v])[indip3jm2]
                                      - (vars[v])[indip2jm3]
                                      - (vars[v])[indim2jp3]
                                      - (vars[v])[indim3jp2]
                                      + (vars[v])[indim2jm3]
                                      + (vars[v])[indim3jm2])
                           + 28224.0 * ((vars[v])[indip2jp2]
                                      + (vars[v])[indim2jm2]
                                      - (vars[v])[indip2jm2]
                                      - (vars[v])[indim2jp2])
                            - 2016.0 * ((vars[v])[indip4jp1]
                                      + (vars[v])[indip1jp4]
                                      - (vars[v])[indip4jm1]
                                      - (vars[v])[indip1jm4]
                                      - (vars[v])[indim1jp4]
                                      - (vars[v])[indim4jp1]
                                      + (vars[v])[indim1jm4]
                                      + (vars[v])[indim4jm1])
                           + 21504.0 * ((vars[v])[indip3jp1]
                                      + (vars[v])[indip1jp3]
                                      - (vars[v])[indip3jm1]
                                      - (vars[v])[indip1jm3]
                                      - (vars[v])[indim1jp3]
                                      - (vars[v])[indim3jp1]
                                      + (vars[v])[indim1jm3]
                                      + (vars[v])[indim3jm1])
                          - 112896.0 * ((vars[v])[indip2jp1]
                                      + (vars[v])[indip1jp2]
                                      - (vars[v])[indip2jm1]
                                      - (vars[v])[indip1jm2]
                                      - (vars[v])[indim1jp2]
                                      - (vars[v])[indim2jp1]
                                      + (vars[v])[indim1jm2]
                                      + (vars[v])[indim2jm1])
                          + 451584.0 * ((vars[v])[indip1jp1]
                                      + (vars[v])[indim1jm1]
                                      - (vars[v])[indip1jm1]
                                      - (vars[v])[indim1jp1])
                            );
    }
  } else {

    for (int v=0; v<nvars; v++)
    {
      dvars[v][1] = 0.0;
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

    for (int v=0; v<nvars; v++)
    {
      dvars[v][1] = hhab * dvars[v][1];
    }

  }

  if (i >= bsize[0] && i < cctkGH->cctk_lsh[0]-bsize[1] &&
      k >= bsize[4] && k < cctkGH->cctk_lsh[2]-bsize[5])
  {
    CCTK_REAL hhacfac = hhac*0.0000014172335600907029478; // 1/705600
    CCTK_INT indip4kp4 = CCTK_GFINDEX3D(cctkGH, i+4, j, k+4);
    CCTK_INT indip4kp3 = CCTK_GFINDEX3D(cctkGH, i+4, j, k+3);
    CCTK_INT indip4kp2 = CCTK_GFINDEX3D(cctkGH, i+4, j, k+2);
    CCTK_INT indip4kp1 = CCTK_GFINDEX3D(cctkGH, i+4, j, k+1);
    CCTK_INT indip4km1 = CCTK_GFINDEX3D(cctkGH, i+4, j, k-1);
    CCTK_INT indip4km2 = CCTK_GFINDEX3D(cctkGH, i+4, j, k-2);
    CCTK_INT indip4km3 = CCTK_GFINDEX3D(cctkGH, i+4, j, k-3);
    CCTK_INT indip4km4 = CCTK_GFINDEX3D(cctkGH, i+4, j, k-4);
    CCTK_INT indip3kp4 = CCTK_GFINDEX3D(cctkGH, i+3, j, k+4);
    CCTK_INT indip3kp3 = CCTK_GFINDEX3D(cctkGH, i+3, j, k+3);
    CCTK_INT indip3kp2 = CCTK_GFINDEX3D(cctkGH, i+3, j, k+2);
    CCTK_INT indip3kp1 = CCTK_GFINDEX3D(cctkGH, i+3, j, k+1);
    CCTK_INT indip3km1 = CCTK_GFINDEX3D(cctkGH, i+3, j, k-1);
    CCTK_INT indip3km2 = CCTK_GFINDEX3D(cctkGH, i+3, j, k-2);
    CCTK_INT indip3km3 = CCTK_GFINDEX3D(cctkGH, i+3, j, k-3);
    CCTK_INT indip3km4 = CCTK_GFINDEX3D(cctkGH, i+3, j, k-4);
    CCTK_INT indip2kp4 = CCTK_GFINDEX3D(cctkGH, i+2, j, k+4);
    CCTK_INT indip2kp3 = CCTK_GFINDEX3D(cctkGH, i+2, j, k+3);
    CCTK_INT indip2kp2 = CCTK_GFINDEX3D(cctkGH, i+2, j, k+2);
    CCTK_INT indip2kp1 = CCTK_GFINDEX3D(cctkGH, i+2, j, k+1);
    CCTK_INT indip2km1 = CCTK_GFINDEX3D(cctkGH, i+2, j, k-1);
    CCTK_INT indip2km2 = CCTK_GFINDEX3D(cctkGH, i+2, j, k-2);
    CCTK_INT indip2km3 = CCTK_GFINDEX3D(cctkGH, i+2, j, k-3);
    CCTK_INT indip2km4 = CCTK_GFINDEX3D(cctkGH, i+2, j, k-4);
    CCTK_INT indip1kp4 = CCTK_GFINDEX3D(cctkGH, i+1, j, k+4);
    CCTK_INT indip1kp3 = CCTK_GFINDEX3D(cctkGH, i+1, j, k+3);
    CCTK_INT indip1kp2 = CCTK_GFINDEX3D(cctkGH, i+1, j, k+2);
    CCTK_INT indip1kp1 = CCTK_GFINDEX3D(cctkGH, i+1, j, k+1);
    CCTK_INT indip1km1 = CCTK_GFINDEX3D(cctkGH, i+1, j, k-1);
    CCTK_INT indip1km2 = CCTK_GFINDEX3D(cctkGH, i+1, j, k-2);
    CCTK_INT indip1km3 = CCTK_GFINDEX3D(cctkGH, i+1, j, k-3);
    CCTK_INT indip1km4 = CCTK_GFINDEX3D(cctkGH, i+1, j, k-4);
    CCTK_INT indim1kp4 = CCTK_GFINDEX3D(cctkGH, i-1, j, k+4);
    CCTK_INT indim1kp3 = CCTK_GFINDEX3D(cctkGH, i-1, j, k+3);
    CCTK_INT indim1kp2 = CCTK_GFINDEX3D(cctkGH, i-1, j, k+2);
    CCTK_INT indim1kp1 = CCTK_GFINDEX3D(cctkGH, i-1, j, k+1);
    CCTK_INT indim1km1 = CCTK_GFINDEX3D(cctkGH, i-1, j, k-1);
    CCTK_INT indim1km2 = CCTK_GFINDEX3D(cctkGH, i-1, j, k-2);
    CCTK_INT indim1km3 = CCTK_GFINDEX3D(cctkGH, i-1, j, k-3);
    CCTK_INT indim1km4 = CCTK_GFINDEX3D(cctkGH, i-1, j, k-4);
    CCTK_INT indim2kp4 = CCTK_GFINDEX3D(cctkGH, i-2, j, k+4);
    CCTK_INT indim2kp3 = CCTK_GFINDEX3D(cctkGH, i-2, j, k+3);
    CCTK_INT indim2kp2 = CCTK_GFINDEX3D(cctkGH, i-2, j, k+2);
    CCTK_INT indim2kp1 = CCTK_GFINDEX3D(cctkGH, i-2, j, k+1);
    CCTK_INT indim2km1 = CCTK_GFINDEX3D(cctkGH, i-2, j, k-1);
    CCTK_INT indim2km2 = CCTK_GFINDEX3D(cctkGH, i-2, j, k-2);
    CCTK_INT indim2km3 = CCTK_GFINDEX3D(cctkGH, i-2, j, k-3);
    CCTK_INT indim2km4 = CCTK_GFINDEX3D(cctkGH, i-2, j, k-4);
    CCTK_INT indim3kp4 = CCTK_GFINDEX3D(cctkGH, i-3, j, k+4);
    CCTK_INT indim3kp3 = CCTK_GFINDEX3D(cctkGH, i-3, j, k+3);
    CCTK_INT indim3kp2 = CCTK_GFINDEX3D(cctkGH, i-3, j, k+2);
    CCTK_INT indim3kp1 = CCTK_GFINDEX3D(cctkGH, i-3, j, k+1);
    CCTK_INT indim3km1 = CCTK_GFINDEX3D(cctkGH, i-3, j, k-1);
    CCTK_INT indim3km2 = CCTK_GFINDEX3D(cctkGH, i-3, j, k-2);
    CCTK_INT indim3km3 = CCTK_GFINDEX3D(cctkGH, i-3, j, k-3);
    CCTK_INT indim3km4 = CCTK_GFINDEX3D(cctkGH, i-3, j, k-4);
    CCTK_INT indim4kp4 = CCTK_GFINDEX3D(cctkGH, i-4, j, k+4);
    CCTK_INT indim4kp3 = CCTK_GFINDEX3D(cctkGH, i-4, j, k+3);
    CCTK_INT indim4kp2 = CCTK_GFINDEX3D(cctkGH, i-4, j, k+2);
    CCTK_INT indim4kp1 = CCTK_GFINDEX3D(cctkGH, i-4, j, k+1);
    CCTK_INT indim4km1 = CCTK_GFINDEX3D(cctkGH, i-4, j, k-1);
    CCTK_INT indim4km2 = CCTK_GFINDEX3D(cctkGH, i-4, j, k-2);
    CCTK_INT indim4km3 = CCTK_GFINDEX3D(cctkGH, i-4, j, k-3);
    CCTK_INT indim4km4 = CCTK_GFINDEX3D(cctkGH, i-4, j, k-4);
    for (int v=0; v<nvars; v++)
    {
      dvars[v][2] = hhacfac*(    9.0 * ((vars[v])[indip4kp4]
                                      + (vars[v])[indim4km4]
                                      - (vars[v])[indip4km4]
                                      - (vars[v])[indim4kp4])
                              - 96.0 * ((vars[v])[indip4kp3]
                                      + (vars[v])[indip3kp4]
                                      - (vars[v])[indip4km3]
                                      - (vars[v])[indip3km4]
                                      - (vars[v])[indim3kp4]
                                      - (vars[v])[indim4kp3]
                                      + (vars[v])[indim3km4]
                                      + (vars[v])[indim4km3])
                            + 1024.0 * ((vars[v])[indip3kp3]
                                      + (vars[v])[indim3km3]
                                      - (vars[v])[indip3km3]
                                      - (vars[v])[indim3kp3])
                             + 504.0 * ((vars[v])[indip4kp2]
                                      + (vars[v])[indip2kp4]
                                      - (vars[v])[indip4km2]
                                      - (vars[v])[indip2km4]
                                      - (vars[v])[indim2kp4]
                                      - (vars[v])[indim4kp2]
                                      + (vars[v])[indim2km4]
                                      + (vars[v])[indim4km2])
                            - 5376.0 * ((vars[v])[indip3kp2]
                                      + (vars[v])[indip2kp3]
                                      - (vars[v])[indip3km2]
                                      - (vars[v])[indip2km3]
                                      - (vars[v])[indim2kp3]
                                      - (vars[v])[indim3kp2]
                                      + (vars[v])[indim2km3]
                                      + (vars[v])[indim3km2])
                           + 28224.0 * ((vars[v])[indip2kp2]
                                      + (vars[v])[indim2km2]
                                      - (vars[v])[indip2km2]
                                      - (vars[v])[indim2kp2])
                            - 2016.0 * ((vars[v])[indip4kp1]
                                      + (vars[v])[indip1kp4]
                                      - (vars[v])[indip4km1]
                                      - (vars[v])[indip1km4]
                                      - (vars[v])[indim1kp4]
                                      - (vars[v])[indim4kp1]
                                      + (vars[v])[indim1km4]
                                      + (vars[v])[indim4km1])
                           + 21504.0 * ((vars[v])[indip3kp1]
                                      + (vars[v])[indip1kp3]
                                      - (vars[v])[indip3km1]
                                      - (vars[v])[indip1km3]
                                      - (vars[v])[indim1kp3]
                                      - (vars[v])[indim3kp1]
                                      + (vars[v])[indim1km3]
                                      + (vars[v])[indim3km1])
                          - 112896.0 * ((vars[v])[indip2kp1]
                                      + (vars[v])[indip1kp2]
                                      - (vars[v])[indip2km1]
                                      - (vars[v])[indip1km2]
                                      - (vars[v])[indim1kp2]
                                      - (vars[v])[indim2kp1]
                                      + (vars[v])[indim1km2]
                                      + (vars[v])[indim2km1])
                          + 451584.0 * ((vars[v])[indip1kp1]
                                      + (vars[v])[indim1km1]
                                      - (vars[v])[indip1km1]
                                      - (vars[v])[indim1kp1])
                            );
    }
  } else {

    for (int v=0; v<nvars; v++)
    {
      dvars[v][2] = 0.0;
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

    for (int v=0; v<nvars; v++)
    {
      dvars[v][2] = hhac * dvars[v][2];
    }

  }

  if (j >= 4 && j < cctkGH->cctk_lsh[1]-4)
  {
    CCTK_REAL hhbbfac = hhbb*0.00019841269841269841270; // 1/5040;
    CCTK_INT indjp4 = CCTK_GFINDEX3D(cctkGH, i, j+4, k);
    CCTK_INT indjp3 = CCTK_GFINDEX3D(cctkGH, i, j+3, k);
    CCTK_INT indjp2 = CCTK_GFINDEX3D(cctkGH, i, j+2, k);
    CCTK_INT indjp1 = CCTK_GFINDEX3D(cctkGH, i, j+1, k);
    CCTK_INT ind    = CCTK_GFINDEX3D(cctkGH, i, j, k);
    CCTK_INT indjm1 = CCTK_GFINDEX3D(cctkGH, i, j-1, k);
    CCTK_INT indjm2 = CCTK_GFINDEX3D(cctkGH, i, j-2, k);
    CCTK_INT indjm3 = CCTK_GFINDEX3D(cctkGH, i, j-3, k);
    CCTK_INT indjm4 = CCTK_GFINDEX3D(cctkGH, i, j-4, k);
    for (int v=0; v<nvars; v++)
    {
      dvars[v][3] = hhbbfac*(  - 9.0 * ((vars[v])[indjp4]
                                      + (vars[v])[indjm4])
                            +  128.0 * ((vars[v])[indjp3]
                                      + (vars[v])[indjm3])
                            - 1008.0 * ((vars[v])[indjp2]
                                      + (vars[v])[indjm2])
                            + 8064.0 * ((vars[v])[indjp1]
                                      + (vars[v])[indjm1])
                            - 14350.0 * (vars[v])[ind]
                            );
    }
  } else {

    for (int v=0; v<nvars; v++)
    {
      dvars[v][3] = 0.0;
    }

    for (int b=jmin2[j]; b<=jmax2[j]; ++b)
    {
      int index = CCTK_GFINDEX3D(cctkGH, i, b, k);
      for (int v=0; v<nvars; v++)
      {
        dvars[v][3] += qb2[b+j*nj] * (vars[v])[index];
      }
    }

    for (int v=0; v<nvars; v++)
    {
      dvars[v][3] = hhbb * dvars[v][3];
    }

  }

  if (j >= bsize[2] && j < cctkGH->cctk_lsh[1]-bsize[3] &&
      k >= bsize[4] && k < cctkGH->cctk_lsh[2]-bsize[5])
  {
    CCTK_REAL hhbcfac = hhbc*0.0000014172335600907029478; // 1/705600;
    CCTK_INT indjp4kp4 = CCTK_GFINDEX3D(cctkGH, i, j+4, k+4);
    CCTK_INT indjp4kp3 = CCTK_GFINDEX3D(cctkGH, i, j+4, k+3);
    CCTK_INT indjp4kp2 = CCTK_GFINDEX3D(cctkGH, i, j+4, k+2);
    CCTK_INT indjp4kp1 = CCTK_GFINDEX3D(cctkGH, i, j+4, k+1);
    CCTK_INT indjp4km1 = CCTK_GFINDEX3D(cctkGH, i, j+4, k-1);
    CCTK_INT indjp4km2 = CCTK_GFINDEX3D(cctkGH, i, j+4, k-2);
    CCTK_INT indjp4km3 = CCTK_GFINDEX3D(cctkGH, i, j+4, k-3);
    CCTK_INT indjp4km4 = CCTK_GFINDEX3D(cctkGH, i, j+4, k-4);
    CCTK_INT indjp3kp4 = CCTK_GFINDEX3D(cctkGH, i, j+3, k+4);
    CCTK_INT indjp3kp3 = CCTK_GFINDEX3D(cctkGH, i, j+3, k+3);
    CCTK_INT indjp3kp2 = CCTK_GFINDEX3D(cctkGH, i, j+3, k+2);
    CCTK_INT indjp3kp1 = CCTK_GFINDEX3D(cctkGH, i, j+3, k+1);
    CCTK_INT indjp3km1 = CCTK_GFINDEX3D(cctkGH, i, j+3, k-1);
    CCTK_INT indjp3km2 = CCTK_GFINDEX3D(cctkGH, i, j+3, k-2);
    CCTK_INT indjp3km3 = CCTK_GFINDEX3D(cctkGH, i, j+3, k-3);
    CCTK_INT indjp3km4 = CCTK_GFINDEX3D(cctkGH, i, j+3, k-4);
    CCTK_INT indjp2kp4 = CCTK_GFINDEX3D(cctkGH, i, j+2, k+4);
    CCTK_INT indjp2kp3 = CCTK_GFINDEX3D(cctkGH, i, j+2, k+3);
    CCTK_INT indjp2kp2 = CCTK_GFINDEX3D(cctkGH, i, j+2, k+2);
    CCTK_INT indjp2kp1 = CCTK_GFINDEX3D(cctkGH, i, j+2, k+1);
    CCTK_INT indjp2km1 = CCTK_GFINDEX3D(cctkGH, i, j+2, k-1);
    CCTK_INT indjp2km2 = CCTK_GFINDEX3D(cctkGH, i, j+2, k-2);
    CCTK_INT indjp2km3 = CCTK_GFINDEX3D(cctkGH, i, j+2, k-3);
    CCTK_INT indjp2km4 = CCTK_GFINDEX3D(cctkGH, i, j+2, k-4);
    CCTK_INT indjp1kp4 = CCTK_GFINDEX3D(cctkGH, i, j+1, k+4);
    CCTK_INT indjp1kp3 = CCTK_GFINDEX3D(cctkGH, i, j+1, k+3);
    CCTK_INT indjp1kp2 = CCTK_GFINDEX3D(cctkGH, i, j+1, k+2);
    CCTK_INT indjp1kp1 = CCTK_GFINDEX3D(cctkGH, i, j+1, k+1);
    CCTK_INT indjp1km1 = CCTK_GFINDEX3D(cctkGH, i, j+1, k-1);
    CCTK_INT indjp1km2 = CCTK_GFINDEX3D(cctkGH, i, j+1, k-2);
    CCTK_INT indjp1km3 = CCTK_GFINDEX3D(cctkGH, i, j+1, k-3);
    CCTK_INT indjp1km4 = CCTK_GFINDEX3D(cctkGH, i, j+1, k-4);
    CCTK_INT indjm1kp4 = CCTK_GFINDEX3D(cctkGH, i, j-1, k+4);
    CCTK_INT indjm1kp3 = CCTK_GFINDEX3D(cctkGH, i, j-1, k+3);
    CCTK_INT indjm1kp2 = CCTK_GFINDEX3D(cctkGH, i, j-1, k+2);
    CCTK_INT indjm1kp1 = CCTK_GFINDEX3D(cctkGH, i, j-1, k+1);
    CCTK_INT indjm1km1 = CCTK_GFINDEX3D(cctkGH, i, j-1, k-1);
    CCTK_INT indjm1km2 = CCTK_GFINDEX3D(cctkGH, i, j-1, k-2);
    CCTK_INT indjm1km3 = CCTK_GFINDEX3D(cctkGH, i, j-1, k-3);
    CCTK_INT indjm1km4 = CCTK_GFINDEX3D(cctkGH, i, j-1, k-4);
    CCTK_INT indjm2kp4 = CCTK_GFINDEX3D(cctkGH, i, j-2, k+4);
    CCTK_INT indjm2kp3 = CCTK_GFINDEX3D(cctkGH, i, j-2, k+3);
    CCTK_INT indjm2kp2 = CCTK_GFINDEX3D(cctkGH, i, j-2, k+2);
    CCTK_INT indjm2kp1 = CCTK_GFINDEX3D(cctkGH, i, j-2, k+1);
    CCTK_INT indjm2km1 = CCTK_GFINDEX3D(cctkGH, i, j-2, k-1);
    CCTK_INT indjm2km2 = CCTK_GFINDEX3D(cctkGH, i, j-2, k-2);
    CCTK_INT indjm2km3 = CCTK_GFINDEX3D(cctkGH, i, j-2, k-3);
    CCTK_INT indjm2km4 = CCTK_GFINDEX3D(cctkGH, i, j-2, k-4);
    CCTK_INT indjm3kp4 = CCTK_GFINDEX3D(cctkGH, i, j-3, k+4);
    CCTK_INT indjm3kp3 = CCTK_GFINDEX3D(cctkGH, i, j-3, k+3);
    CCTK_INT indjm3kp2 = CCTK_GFINDEX3D(cctkGH, i, j-3, k+2);
    CCTK_INT indjm3kp1 = CCTK_GFINDEX3D(cctkGH, i, j-3, k+1);
    CCTK_INT indjm3km1 = CCTK_GFINDEX3D(cctkGH, i, j-3, k-1);
    CCTK_INT indjm3km2 = CCTK_GFINDEX3D(cctkGH, i, j-3, k-2);
    CCTK_INT indjm3km3 = CCTK_GFINDEX3D(cctkGH, i, j-3, k-3);
    CCTK_INT indjm3km4 = CCTK_GFINDEX3D(cctkGH, i, j-3, k-4);
    CCTK_INT indjm4kp4 = CCTK_GFINDEX3D(cctkGH, i, j-4, k+4);
    CCTK_INT indjm4kp3 = CCTK_GFINDEX3D(cctkGH, i, j-4, k+3);
    CCTK_INT indjm4kp2 = CCTK_GFINDEX3D(cctkGH, i, j-4, k+2);
    CCTK_INT indjm4kp1 = CCTK_GFINDEX3D(cctkGH, i, j-4, k+1);
    CCTK_INT indjm4km1 = CCTK_GFINDEX3D(cctkGH, i, j-4, k-1);
    CCTK_INT indjm4km2 = CCTK_GFINDEX3D(cctkGH, i, j-4, k-2);
    CCTK_INT indjm4km3 = CCTK_GFINDEX3D(cctkGH, i, j-4, k-3);
    CCTK_INT indjm4km4 = CCTK_GFINDEX3D(cctkGH, i, j-4, k-4);
    for (int v=0; v<nvars; v++)
    {
      dvars[v][4] = hhbcfac*(    9.0 * ((vars[v])[indjp4kp4]
                                      + (vars[v])[indjm4km4]
                                      - (vars[v])[indjp4km4]
                                      - (vars[v])[indjm4kp4])
                              - 96.0 * ((vars[v])[indjp4kp3]
                                      + (vars[v])[indjp3kp4]
                                      - (vars[v])[indjp4km3]
                                      - (vars[v])[indjp3km4]
                                      - (vars[v])[indjm3kp4]
                                      - (vars[v])[indjm4kp3]
                                      + (vars[v])[indjm3km4]
                                      + (vars[v])[indjm4km3])
                            + 1024.0 * ((vars[v])[indjp3kp3]
                                      + (vars[v])[indjm3km3]
                                      - (vars[v])[indjp3km3]
                                      - (vars[v])[indjm3kp3])
                             + 504.0 * ((vars[v])[indjp4kp2]
                                      + (vars[v])[indjp2kp4]
                                      - (vars[v])[indjp4km2]
                                      - (vars[v])[indjp2km4]
                                      - (vars[v])[indjm2kp4]
                                      - (vars[v])[indjm4kp2]
                                      + (vars[v])[indjm2km4]
                                      + (vars[v])[indjm4km2])
                            - 5376.0 * ((vars[v])[indjp3kp2]
                                      + (vars[v])[indjp2kp3]
                                      - (vars[v])[indjp3km2]
                                      - (vars[v])[indjp2km3]
                                      - (vars[v])[indjm2kp3]
                                      - (vars[v])[indjm3kp2]
                                      + (vars[v])[indjm2km3]
                                      + (vars[v])[indjm3km2])
                           + 28224.0 * ((vars[v])[indjp2kp2]
                                      + (vars[v])[indjm2km2]
                                      - (vars[v])[indjp2km2]
                                      - (vars[v])[indjm2kp2])
                            - 2016.0 * ((vars[v])[indjp4kp1]
                                      + (vars[v])[indjp1kp4]
                                      - (vars[v])[indjp4km1]
                                      - (vars[v])[indjp1km4]
                                      - (vars[v])[indjm1kp4]
                                      - (vars[v])[indjm4kp1]
                                      + (vars[v])[indjm1km4]
                                      + (vars[v])[indjm4km1])
                           + 21504.0 * ((vars[v])[indjp3kp1]
                                      + (vars[v])[indjp1kp3]
                                      - (vars[v])[indjp3km1]
                                      - (vars[v])[indjp1km3]
                                      - (vars[v])[indjm1kp3]
                                      - (vars[v])[indjm3kp1]
                                      + (vars[v])[indjm1km3]
                                      + (vars[v])[indjm3km1])
                          - 112896.0 * ((vars[v])[indjp2kp1]
                                      + (vars[v])[indjp1kp2]
                                      - (vars[v])[indjp2km1]
                                      - (vars[v])[indjp1km2]
                                      - (vars[v])[indjm1kp2]
                                      - (vars[v])[indjm2kp1]
                                      + (vars[v])[indjm1km2]
                                      + (vars[v])[indjm2km1])
                          + 451584.0 * ((vars[v])[indjp1kp1]
                                      + (vars[v])[indjm1km1]
                                      - (vars[v])[indjp1km1]
                                      - (vars[v])[indjm1kp1])
                            );
    }
  } else {

    for (int v=0; v<nvars; v++)
    {
      dvars[v][4] = 0.0;
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

    for (int v=0; v<nvars; v++)
    {
      dvars[v][4] = hhbc * dvars[v][4];
    }

  }

  if (k >= 4 && k < cctkGH->cctk_lsh[2]-4)
  {
    CCTK_REAL hhccfac = hhcc*0.00019841269841269841270; // 1/5040 
    CCTK_INT indkp4 = CCTK_GFINDEX3D(cctkGH, i, j, k+4);
    CCTK_INT indkp3 = CCTK_GFINDEX3D(cctkGH, i, j, k+3);
    CCTK_INT indkp2 = CCTK_GFINDEX3D(cctkGH, i, j, k+2);
    CCTK_INT indkp1 = CCTK_GFINDEX3D(cctkGH, i, j, k+1);
    CCTK_INT ind    = CCTK_GFINDEX3D(cctkGH, i, j, k);
    CCTK_INT indkm1 = CCTK_GFINDEX3D(cctkGH, i, j, k-1);
    CCTK_INT indkm2 = CCTK_GFINDEX3D(cctkGH, i, j, k-2);
    CCTK_INT indkm3 = CCTK_GFINDEX3D(cctkGH, i, j, k-3);
    CCTK_INT indkm4 = CCTK_GFINDEX3D(cctkGH, i, j, k-4);
    for (int v=0; v<nvars; v++)
    {
      dvars[v][5] = hhccfac*(  - 9.0 * ((vars[v])[indkp4]
                                      + (vars[v])[indkm4])
                            +  128.0 * ((vars[v])[indkp3]
                                      + (vars[v])[indkm3])
                            - 1008.0 * ((vars[v])[indkp2]
                                      + (vars[v])[indkm2])
                            + 8064.0 * ((vars[v])[indkp1]
                                      + (vars[v])[indkm1])
                            - 14350.0 * (vars[v])[ind]);
    }
  } else {

    for (int v=0; v<nvars; v++)
    {
      dvars[v][5] = 0.0;
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
      dvars[v][5] = hhcc * dvars[v][5];
    }

  }

}

static inline void adv_diff_8 (
          const cGH* restrict const cctkGH,
          const CCTK_INT nvars,
          const CCTK_REAL* const vars[],
          CCTK_REAL dvars[][3],
          const CCTK_REAL upwind_dir[3],
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
          const CCTK_INT* restrict const bsize,
          const CCTK_REAL ha,
          const CCTK_REAL hb,
          const CCTK_REAL hc)
{

  if (upwind_dir[0] > UPWIND_SMALL)
  {

    if (i>=3 && i < cctkGH->cctk_lsh[0]-5)
    {
      CCTK_REAL hafac = ha*0.0011904761904761904762; // 1/840
      CCTK_INT indip5 = CCTK_GFINDEX3D(cctkGH, i+5, j, k);
      CCTK_INT indip4 = CCTK_GFINDEX3D(cctkGH, i+4, j, k);
      CCTK_INT indip3 = CCTK_GFINDEX3D(cctkGH, i+3, j, k);
      CCTK_INT indip2 = CCTK_GFINDEX3D(cctkGH, i+2, j, k);
      CCTK_INT indip1 = CCTK_GFINDEX3D(cctkGH, i+1, j, k);
      CCTK_INT ind    = CCTK_GFINDEX3D(cctkGH, i, j, k);
      CCTK_INT indim1 = CCTK_GFINDEX3D(cctkGH, i-1, j, k);
      CCTK_INT indim2 = CCTK_GFINDEX3D(cctkGH, i-2, j, k);
      CCTK_INT indim3 = CCTK_GFINDEX3D(cctkGH, i-3, j, k);
      for (int v=0; v<nvars; v++)
      {
        dvars[v][0] = hafac*(    3.0 * (vars[v])[indip5]
                            -   30.0 * (vars[v])[indip4]
                            +  140.0 * (vars[v])[indip3]
                            -  420.0 * (vars[v])[indip2]
                            + 1050.0 * (vars[v])[indip1]
                            -  378.0 * (vars[v])[ind]
                            -  420.0 * (vars[v])[indim1]
                            +   60.0 * (vars[v])[indim2]
                            -    5.0 * (vars[v])[indim3]);
      }
    } else {

      for (int v=0; v<nvars; v++)
      {
        dvars[v][0] = 0.0;
      }

      for (int a=imin[i]; a<=imax[i]; ++a)
      {
        int index = CCTK_GFINDEX3D(cctkGH, a, j, k);
        for (int v=0; v<nvars; v++)
        {
          dvars[v][0] += qa[a+i*ni] * (vars[v])[index];
        }
      }
      for (int v=0; v<nvars; v++)
      {
        dvars[v][0] = ha * dvars[v][0];
      }

    }

  } 
  else if (upwind_dir[0] < -UPWIND_SMALL)
  {

    if (i>=5 && i < cctkGH->cctk_lsh[0]-3)
    {
      CCTK_REAL hafac = ha*0.0011904761904761904762; // 1/840
      CCTK_INT indip3 = CCTK_GFINDEX3D(cctkGH, i+3, j, k);
      CCTK_INT indip2 = CCTK_GFINDEX3D(cctkGH, i+2, j, k);
      CCTK_INT indip1 = CCTK_GFINDEX3D(cctkGH, i+1, j, k);
      CCTK_INT ind    = CCTK_GFINDEX3D(cctkGH, i, j, k);
      CCTK_INT indim1 = CCTK_GFINDEX3D(cctkGH, i-1, j, k);
      CCTK_INT indim2 = CCTK_GFINDEX3D(cctkGH, i-2, j, k);
      CCTK_INT indim3 = CCTK_GFINDEX3D(cctkGH, i-3, j, k);
      CCTK_INT indim4 = CCTK_GFINDEX3D(cctkGH, i-4, j, k);
      CCTK_INT indim5 = CCTK_GFINDEX3D(cctkGH, i-5, j, k);
      for (int v=0; v<nvars; v++)
      {
        dvars[v][0] = hafac*(    5.0 * (vars[v])[indip3]
                            -   60.0 * (vars[v])[indip2]
                            +  420.0 * (vars[v])[indip1]
                            +  378.0 * (vars[v])[ind]
                            - 1050.0 * (vars[v])[indim1]
                            +  420.0 * (vars[v])[indim2]
                            -  140.0 * (vars[v])[indim3]
                            +   30.0 * (vars[v])[indim4]
                            -    3.0 * (vars[v])[indim5]);
      }
    } else {

      for (int v=0; v<nvars; v++)
      {
        dvars[v][0] = 0.0;
      }

      for (int a=imin[i]; a<=imax[i]; ++a)
      {
        int index = CCTK_GFINDEX3D(cctkGH, a, j, k);
        for (int v=0; v<nvars; v++)
        {
          dvars[v][0] += qa[a+i*ni] * (vars[v])[index];
        }
      }
      for (int v=0; v<nvars; v++)
      {
        dvars[v][0] = ha * dvars[v][0];
      }

    }

  }
  else
  {

    if (i >= bsize[0] && i < cctkGH->cctk_lsh[0]-bsize[1])
    {
      CCTK_REAL hafac = ha*0.0011904761904761904762;  // 1/840
      CCTK_INT indip4 = CCTK_GFINDEX3D(cctkGH, i+4, j, k);
      CCTK_INT indip3 = CCTK_GFINDEX3D(cctkGH, i+3, j, k);
      CCTK_INT indip2 = CCTK_GFINDEX3D(cctkGH, i+2, j, k);
      CCTK_INT indip1 = CCTK_GFINDEX3D(cctkGH, i+1, j, k);
      CCTK_INT indim1 = CCTK_GFINDEX3D(cctkGH, i-1, j, k);
      CCTK_INT indim2 = CCTK_GFINDEX3D(cctkGH, i-2, j, k);
      CCTK_INT indim3 = CCTK_GFINDEX3D(cctkGH, i-3, j, k);
      CCTK_INT indim4 = CCTK_GFINDEX3D(cctkGH, i-4, j, k);
      for (int v=0; v<nvars; v++)
      {
        dvars[v][0] = hafac*( -3.0 *((vars[v])[indip4]
                                   - (vars[v])[indim4])
                             +32.0 *((vars[v])[indip3]
                                   - (vars[v])[indim3])
                            -168.0 *((vars[v])[indip2]
                                   - (vars[v])[indim2])
                            +672.0 *((vars[v])[indip1]
                                   - (vars[v])[indim1]));
      }
    } else {

      for (int v=0; v<nvars; v++)
      {
        dvars[v][0] = 0.0;
      }

      for (int a=imin[i]; a<=imax[i]; ++a)
      {
        int index = CCTK_GFINDEX3D(cctkGH, a, j, k);
        for (int v=0; v<nvars; v++)
        {
          dvars[v][0] += qa[a+i*ni] * (vars[v])[index];
        }
      }
      for (int v=0; v<nvars; v++)
      {
        dvars[v][0] = ha * dvars[v][0];
      }

    }

  }

  if (upwind_dir[1] > UPWIND_SMALL)
  {

    if (j>=3 && j < cctkGH->cctk_lsh[1]-5)
    {
      CCTK_REAL hbfac = hb*0.0011904761904761904762; // 1/840
      CCTK_INT indjp5 = CCTK_GFINDEX3D(cctkGH, i, j+5, k);
      CCTK_INT indjp4 = CCTK_GFINDEX3D(cctkGH, i, j+4, k);
      CCTK_INT indjp3 = CCTK_GFINDEX3D(cctkGH, i, j+3, k);
      CCTK_INT indjp2 = CCTK_GFINDEX3D(cctkGH, i, j+2, k);
      CCTK_INT indjp1 = CCTK_GFINDEX3D(cctkGH, i, j+1, k);
      CCTK_INT ind    = CCTK_GFINDEX3D(cctkGH, i, j, k);
      CCTK_INT indjm1 = CCTK_GFINDEX3D(cctkGH, i, j-1, k);
      CCTK_INT indjm2 = CCTK_GFINDEX3D(cctkGH, i, j-2, k);
      CCTK_INT indjm3 = CCTK_GFINDEX3D(cctkGH, i, j-3, k);
      for (int v=0; v<nvars; v++)
      {
        dvars[v][1] = hbfac*(    3.0 * (vars[v])[indjp5]
                            -   30.0 * (vars[v])[indjp4]
                            +  140.0 * (vars[v])[indjp3]
                            -  420.0 * (vars[v])[indjp2]
                            + 1050.0 * (vars[v])[indjp1]
                            -  378.0 * (vars[v])[ind]
                            -  420.0 * (vars[v])[indjm1]
                            +   60.0 * (vars[v])[indjm2]
                            -    5.0 * (vars[v])[indjm3]);
      }
    } else {

      for (int v=0; v<nvars; v++)
      {
        dvars[v][1] = 0.0;
      }

      for (int b=jmin[j]; b<=jmax[j]; ++b)
      {
        int index = CCTK_GFINDEX3D(cctkGH, i, b, k);
        for (int v=0; v<nvars; v++)
        {
          dvars[v][1] += qb[b+j*nj] * (vars[v])[index];
        }
      }
      for (int v=0; v<nvars; v++)
      {
        dvars[v][1] = hb * dvars[v][1];
      }

    }

  }
  else if (upwind_dir[1] < -UPWIND_SMALL)
  {

    if (j>=5 && j < cctkGH->cctk_lsh[1]-3)
    {
      CCTK_REAL hbfac = hb*0.0011904761904761904762; // 1/840
      CCTK_INT indjp3 = CCTK_GFINDEX3D(cctkGH, i, j+3, k);
      CCTK_INT indjp2 = CCTK_GFINDEX3D(cctkGH, i, j+2, k);
      CCTK_INT indjp1 = CCTK_GFINDEX3D(cctkGH, i, j+1, k);
      CCTK_INT ind    = CCTK_GFINDEX3D(cctkGH, i, j, k);
      CCTK_INT indjm1 = CCTK_GFINDEX3D(cctkGH, i, j-1, k);
      CCTK_INT indjm2 = CCTK_GFINDEX3D(cctkGH, i, j-2, k);
      CCTK_INT indjm3 = CCTK_GFINDEX3D(cctkGH, i, j-3, k);
      CCTK_INT indjm4 = CCTK_GFINDEX3D(cctkGH, i, j-4, k);
      CCTK_INT indjm5 = CCTK_GFINDEX3D(cctkGH, i, j-5, k);
      for (int v=0; v<nvars; v++)
      {
        dvars[v][1] = hbfac*(    5.0 * (vars[v])[indjp3]
                            -   60.0 * (vars[v])[indjp2]
                            +  420.0 * (vars[v])[indjp1]
                            +  378.0 * (vars[v])[ind]
                            - 1050.0 * (vars[v])[indjm1]
                            +  420.0 * (vars[v])[indjm2]
                            -  140.0 * (vars[v])[indjm3]
                            +   30.0 * (vars[v])[indjm4]
                            -    3.0 * (vars[v])[indjm5]);
      }
    } else {

      for (int v=0; v<nvars; v++)
      {
        dvars[v][1] = 0.0;
      }

      for (int b=jmin[j]; b<=jmax[j]; ++b)
      {
        int index = CCTK_GFINDEX3D(cctkGH, i, b, k);
        for (int v=0; v<nvars; v++)
        {
          dvars[v][1] += qb[b+j*nj] * (vars[v])[index];
        }
      }
      for (int v=0; v<nvars; v++)
      {
        dvars[v][1] = hb * dvars[v][1];
      }

    }

  }
  else
  {

    if (j >= bsize[2] && j < cctkGH->cctk_lsh[1]-bsize[3])
    {
      CCTK_REAL hbfac = hb*0.0011904761904761904762;  // 1/840
      CCTK_INT indjp4 = CCTK_GFINDEX3D(cctkGH, i, j+4, k);
      CCTK_INT indjp3 = CCTK_GFINDEX3D(cctkGH, i, j+3, k);
      CCTK_INT indjp2 = CCTK_GFINDEX3D(cctkGH, i, j+2, k);
      CCTK_INT indjp1 = CCTK_GFINDEX3D(cctkGH, i, j+1, k);
      CCTK_INT indjm1 = CCTK_GFINDEX3D(cctkGH, i, j-1, k);
      CCTK_INT indjm2 = CCTK_GFINDEX3D(cctkGH, i, j-2, k);
      CCTK_INT indjm3 = CCTK_GFINDEX3D(cctkGH, i, j-3, k);
      CCTK_INT indjm4 = CCTK_GFINDEX3D(cctkGH, i, j-4, k);
      for (int v=0; v<nvars; v++)
      {
        dvars[v][1] = hbfac*( -3.0 *((vars[v])[indjp4]
                                   - (vars[v])[indjm4])
                             +32.0 *((vars[v])[indjp3]
                                   - (vars[v])[indjm3])
                            -168.0 *((vars[v])[indjp2]
                                   - (vars[v])[indjm2])
                            +672.0 *((vars[v])[indjp1]
                                   - (vars[v])[indjm1]));
      }
    } else {

      for (int v=0; v<nvars; v++)
      {
        dvars[v][1] = 0.0;
      }
      for (int b=jmin[j]; b<=jmax[j]; ++b)
      {
        int index = CCTK_GFINDEX3D(cctkGH, i, b, k);
        for (int v=0; v<nvars; v++)
        {
          dvars[v][1] += qb[b+j*nj] * (vars[v])[index];
        }
      }
      for (int v=0; v<nvars; v++)
      {
        dvars[v][1] = hb * dvars[v][1];
      }

    }

  }

  if (upwind_dir[2] > UPWIND_SMALL)
  {

    if (k>=3 && k < cctkGH->cctk_lsh[2]-5)
    {
      CCTK_REAL hcfac = hc*0.0011904761904761904762; // 1/840
      CCTK_INT indkp5 = CCTK_GFINDEX3D(cctkGH, i, j, k+5);
      CCTK_INT indkp4 = CCTK_GFINDEX3D(cctkGH, i, j, k+4);
      CCTK_INT indkp3 = CCTK_GFINDEX3D(cctkGH, i, j, k+3);
      CCTK_INT indkp2 = CCTK_GFINDEX3D(cctkGH, i, j, k+2);
      CCTK_INT indkp1 = CCTK_GFINDEX3D(cctkGH, i, j, k+1);
      CCTK_INT ind    = CCTK_GFINDEX3D(cctkGH, i, j, k);
      CCTK_INT indkm1 = CCTK_GFINDEX3D(cctkGH, i, j, k-1);
      CCTK_INT indkm2 = CCTK_GFINDEX3D(cctkGH, i, j, k-2);
      CCTK_INT indkm3 = CCTK_GFINDEX3D(cctkGH, i, j, k-3);
      for (int v=0; v<nvars; v++)
      {
        dvars[v][2] = hcfac*(    3.0 * (vars[v])[indkp5]
                            -   30.0 * (vars[v])[indkp4]
                            +  140.0 * (vars[v])[indkp3]
                            -  420.0 * (vars[v])[indkp2]
                            + 1050.0 * (vars[v])[indkp1]
                            -  378.0 * (vars[v])[ind]
                            -  420.0 * (vars[v])[indkm1]
                            +   60.0 * (vars[v])[indkm2]
                            -    5.0 * (vars[v])[indkm3]);
      }
    } else {

      for (int v=0; v<nvars; v++)
      {
        dvars[v][2] = 0.0;
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
        dvars[v][2] = hc * dvars[v][2];
      }

    }

  }
  else if (upwind_dir[2] < -UPWIND_SMALL)
  {

    if (k>=5 && k < cctkGH->cctk_lsh[2]-3)
    {
      CCTK_REAL hcfac = hc*0.0011904761904761904762; // 1/840
      CCTK_INT indkp3 = CCTK_GFINDEX3D(cctkGH, i, j, k+3);
      CCTK_INT indkp2 = CCTK_GFINDEX3D(cctkGH, i, j, k+2);
      CCTK_INT indkp1 = CCTK_GFINDEX3D(cctkGH, i, j, k+1);
      CCTK_INT ind    = CCTK_GFINDEX3D(cctkGH, i, j, k);
      CCTK_INT indkm1 = CCTK_GFINDEX3D(cctkGH, i, j, k-1);
      CCTK_INT indkm2 = CCTK_GFINDEX3D(cctkGH, i, j, k-2);
      CCTK_INT indkm3 = CCTK_GFINDEX3D(cctkGH, i, j, k-3);
      CCTK_INT indkm4 = CCTK_GFINDEX3D(cctkGH, i, j, k-4);
      CCTK_INT indkm5 = CCTK_GFINDEX3D(cctkGH, i, j, k-5);
      for (int v=0; v<nvars; v++)
      {
        dvars[v][2] = hcfac*(    5.0 * (vars[v])[indkp3]
                            -   60.0 * (vars[v])[indkp2]
                            +  420.0 * (vars[v])[indkp1]
                            +  378.0 * (vars[v])[ind]
                            - 1050.0 * (vars[v])[indkm1]
                            +  420.0 * (vars[v])[indkm2]
                            -  140.0 * (vars[v])[indkm3]
                            +   30.0 * (vars[v])[indkm4]
                            -    3.0 * (vars[v])[indkm5]);
      }
    } else {

      for (int v=0; v<nvars; v++)
      {
        dvars[v][2] = 0.0;
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
        dvars[v][2] = hc * dvars[v][2];
      }
    }

  }
  else
  {

    if (k >= bsize[4] && k < cctkGH->cctk_lsh[2]-bsize[5])
    {
      CCTK_REAL hcfac = hc*0.0011904761904761904762;  // 1/840
      CCTK_INT indkp4 = CCTK_GFINDEX3D(cctkGH, i, j, k+4);
      CCTK_INT indkp3 = CCTK_GFINDEX3D(cctkGH, i, j, k+3);
      CCTK_INT indkp2 = CCTK_GFINDEX3D(cctkGH, i, j, k+2);
      CCTK_INT indkp1 = CCTK_GFINDEX3D(cctkGH, i, j, k+1);
      CCTK_INT indkm1 = CCTK_GFINDEX3D(cctkGH, i, j, k-1);
      CCTK_INT indkm2 = CCTK_GFINDEX3D(cctkGH, i, j, k-2);
      CCTK_INT indkm3 = CCTK_GFINDEX3D(cctkGH, i, j, k-3);
      CCTK_INT indkm4 = CCTK_GFINDEX3D(cctkGH, i, j, k-4);
      for (int v=0; v<nvars; v++)
      {
        dvars[v][2] = hcfac*( -3.0 *((vars[v])[indkp4]
                                   - (vars[v])[indkm4])
                             +32.0 *((vars[v])[indkp3]
                                   - (vars[v])[indkm3])
                            -168.0 *((vars[v])[indkp2]
                                   - (vars[v])[indkm2])
                            +672.0 *((vars[v])[indkp1]
                                   - (vars[v])[indkm1]));
      }
    } else {

      for (int v=0; v<nvars; v++)
      {
        dvars[v][2] = 0.0;
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
        dvars[v][2] = hc * dvars[v][2];
      }
    }

  }
}

#endif
