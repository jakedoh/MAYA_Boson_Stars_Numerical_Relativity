
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

#ifndef GLOBALDERIV_H
#define GLOBALDERIV_H

#include "cctk.h"

#undef HARDCODED_STENCILS_2
#undef HARDCODED_STENCILS_4
#define SBP_STENCILS

/*
 * Poison checking within the derivative operators can be implementing
 * by setting
 *  #define CHECK_FOR_POISON
 * This is somewhat inefficient, so it is disabled by default. The poison
 * value is the standard one used by Carpet.
 */
#undef CHECK_FOR_POISON

#define POISON_VAL -424242.0
#define Check_For_Poison(i, j, k, iiimin, iiimax, stride, dir)          \
  do                                                                    \
    {                                                                   \
      for (int iii=iiimin; iii<=iiimax; ++iii)                           \
        {                                                               \
          int ijk;                                                      \
          if (dir==0)                                                   \
            ijk = CCTK_GFINDEX3D(cctkGH, iii, j, k);                    \
          else if (dir==1)                                              \
            ijk = CCTK_GFINDEX3D(cctkGH, i, iii, k);                    \
          else                                                          \
            ijk = CCTK_GFINDEX3D(cctkGH, i, j, iii);                    \
          if (fn[ijk] == POISON_VAL)                                    \
            {                                                           \
              CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,       \
                         "Derivative stencil attempted to use poisoned data at point (%d %d %d), dir=%d", i, j, k, dir); \
            }                                                           \
        }                                                               \
    }                                                                   \
  while (0)

#define Check_For_Poison2(i, j, k, iiimin, iiimax, jjjmin, jjjmax,      \
                          istride, jstride, kdir)                       \
  do                                                                    \
    {                                                                   \
      for (int iii=iiimin; iii<=iiimax; ++iii)                           \
        {                                                               \
          for (int jjj=jjjmin; jjj<=jjjmax; ++jjj)                       \
            {                                                           \
              int ijk;                                                  \
              if (kdir==0)                                              \
                ijk = CCTK_GFINDEX3D(cctkGH, i, iii, jjj);              \
              else if (kdir==1)                                         \
                ijk = CCTK_GFINDEX3D(cctkGH, iii, j, jjj);              \
              else                                                      \
                ijk = CCTK_GFINDEX3D(cctkGH, iii, jjj, k);              \
              if (fn[ijk] == POISON_VAL)                                \
                {                                                       \
                  CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,   \
                             "Derivative stencil attempted to use poisoned data at point (%d %d %d), dir=%d", i, j, k, kdir); \
                }                                                       \
            }                                                           \
        }                                                               \
    }                                                                   \
  while (0)



#ifdef SBP_STENCILS

static inline __attribute__((always_inline)) CCTK_REAL
calc_da(const cGH* restrict const cctkGH, const CCTK_REAL* const fn,
        const CCTK_INT i, const CCTK_INT j, const CCTK_INT k,
        const CCTK_INT stride, const CCTK_INT* restrict const imin,
        const CCTK_INT* restrict const imax, const CCTK_REAL* restrict const q,
        const CCTK_REAL h)
{
  CCTK_REAL dval = 0.0;

#ifdef CHECK_FOR_POISON
    Check_For_Poison(i, j, k, imin[i], imax[i], stride, 0);
#endif

  for (int a=imin[i]; a<=imax[i]; ++a)
    dval += q[a+i*stride] * fn[CCTK_GFINDEX3D(cctkGH, a, j, k)];
  return h*dval;
}



static inline __attribute__((always_inline)) CCTK_REAL
calc_db(const cGH* restrict const cctkGH, const CCTK_REAL* const fn,
        const CCTK_INT i, const CCTK_INT j, const CCTK_INT k,
        const CCTK_INT stride, const CCTK_INT* restrict const jmin,
        const CCTK_INT* restrict const jmax, const CCTK_REAL* restrict const q,
        const CCTK_REAL h)
{
  CCTK_REAL dval = 0.0;

#ifdef CHECK_FOR_POISON
    Check_For_Poison(i, j, k, jmin[j], jmax[j], stride, 1);
#endif

  for (int a=jmin[j]; a<=jmax[j]; ++a)
    dval += q[a+j*stride] * fn[CCTK_GFINDEX3D(cctkGH, i, a, k)];
  return h*dval;
}


static inline __attribute__((always_inline)) CCTK_REAL
calc_dc(const cGH* restrict const cctkGH, const CCTK_REAL* const fn,
        const CCTK_INT i, const CCTK_INT j, const CCTK_INT k,
        const CCTK_INT stride, const CCTK_INT* restrict const kmin,
        const CCTK_INT* restrict const kmax, const CCTK_REAL* restrict const q,
        const CCTK_REAL h)
{
  CCTK_REAL dval = 0.0;

#ifdef CHECK_FOR_POISON
    Check_For_Poison(i, j, k, kmin[k], kmax[k], stride, 2);
#endif

  for (int a=kmin[k]; a<=kmax[k]; ++a)
    dval += q[a+k*stride] * fn[CCTK_GFINDEX3D(cctkGH, i, j, a)];
  return h*dval;
}


static inline __attribute__((always_inline)) CCTK_REAL
calc_dada(const cGH* restrict const cctkGH, const CCTK_REAL* const fn, 
          const CCTK_INT i, const CCTK_INT j, const CCTK_INT k,
          const CCTK_INT stride, const CCTK_INT* restrict const imin,
          const CCTK_INT* restrict const imax,
          const CCTK_REAL* restrict const q, const CCTK_REAL hh)
{
  CCTK_REAL ddval = 0.0;

#ifdef CHECK_FOR_POISON
    Check_For_Poison(i, j, k, imin[i], imax[i], stride, 0);
#endif

  for (int a=imin[i]; a<=imax[i]; ++a)
    ddval += q[a+i*stride] * fn[CCTK_GFINDEX3D(cctkGH, a, j, k)];
  return hh*ddval;
}


static inline __attribute__((always_inline)) CCTK_REAL
calc_dadb(const cGH* restrict const cctkGH, const CCTK_REAL* const fn, 
          const CCTK_INT i, const CCTK_INT j, const CCTK_INT k,
          const CCTK_INT istride, const CCTK_INT jstride, 
          const CCTK_INT* restrict const imin,
          const CCTK_INT* restrict const imax, 
          const CCTK_INT* restrict const jmin,
          const CCTK_INT* restrict const jmax, 
          const CCTK_REAL* restrict const iq,
          const CCTK_REAL* restrict const jq, const CCTK_REAL hh)
{
  CCTK_REAL ddval=0.0;

#ifdef CHECK_FOR_POISON
  Check_For_Poison2(i, j, k, imin[i], imax[i], jmin[j], jmax[j], istride,
		    jstride, 2);
#endif

  for (int b=jmin[j]; b<=jmax[j]; ++b)
    {
      CCTK_REAL dval = 0.0;
      for (int a=imin[i]; a<=imax[i]; ++a)
        dval += iq[a+i*istride] * fn[CCTK_GFINDEX3D(cctkGH, a, b, k)];
      ddval += jq[b+j*jstride] * dval;
    }
  return hh*ddval;
}

static inline __attribute__((always_inline)) CCTK_REAL
calc_dadc(const cGH* restrict const cctkGH, const CCTK_REAL* const fn, 
          const CCTK_INT i, const CCTK_INT j, const CCTK_INT k,
          const CCTK_INT istride, const CCTK_INT kstride, 
          const CCTK_INT* restrict const imin,
          const CCTK_INT* restrict const imax,
          const CCTK_INT* restrict const kmin,
          const CCTK_INT* restrict const kmax,
          const CCTK_REAL* restrict const iq,
          const CCTK_REAL* restrict const kq, const CCTK_REAL hh)
{
  CCTK_REAL ddval=0.0;

#ifdef CHECK_FOR_POISON
  Check_For_Poison2(i, j, k, imin[i], imax[i], kmin[k], kmax[k], istride,
		    kstride, 1);
#endif

  for (int b=kmin[k]; b<=kmax[k]; ++b)
    {
      CCTK_REAL dval = 0.0;
      for (int a=imin[i]; a<=imax[i]; ++a)
        dval += iq[a+i*istride] * fn[CCTK_GFINDEX3D(cctkGH, a, j, b)];
      ddval += kq[b+k*kstride] * dval;
    }
  return hh*ddval;
}

static inline __attribute__((always_inline)) CCTK_REAL
calc_dbdb(const cGH* restrict const cctkGH, const CCTK_REAL* const fn, 
          const CCTK_INT i, const CCTK_INT j, const CCTK_INT k,
          const CCTK_INT stride, const CCTK_INT* restrict const jmin,
          const CCTK_INT* restrict const jmax,
          const CCTK_REAL* restrict const q, const CCTK_REAL hh)
{
  CCTK_REAL ddval=0.0;

#ifdef CHECK_FOR_POISON
    Check_For_Poison(i, j, k, jmin[j], jmax[j], stride, 1);
#endif

  for (int a=jmin[j]; a<=jmax[j]; ++a)
    ddval += q[a+j*stride] * fn[CCTK_GFINDEX3D(cctkGH, i, a, k)];
  return hh*ddval;
}

static inline __attribute__((always_inline)) CCTK_REAL
calc_dbdc(const cGH* restrict const cctkGH, const CCTK_REAL* const fn, 
          const CCTK_INT i, const CCTK_INT j, const CCTK_INT k,
          const CCTK_INT jstride, const CCTK_INT kstride,
          const CCTK_INT* restrict const jmin,
          const CCTK_INT* restrict const jmax,
          const CCTK_INT* restrict const kmin,
          const CCTK_INT* restrict const kmax,
          const CCTK_REAL* restrict const jq,
          const CCTK_REAL* restrict const kq, const CCTK_REAL hh)
{
  CCTK_REAL ddval=0.0;

#ifdef CHECK_FOR_POISON
  Check_For_Poison2(i, j, k, jmin[j], jmax[j], kmin[k], kmax[k], jstride,
		    kstride, 0);
#endif

  for (int b=kmin[k]; b<=kmax[k]; ++b)
    {
      CCTK_REAL dval = 0.0;
      for (int a=jmin[j]; a<=jmax[j]; ++a)
        dval += jq[a+j*jstride] * fn[CCTK_GFINDEX3D(cctkGH, i, a, b)];
      ddval += kq[b+k*kstride] * dval;
    }
  return hh*ddval;
}

static inline __attribute__((always_inline)) CCTK_REAL
calc_dcdc(const cGH* restrict const cctkGH, const CCTK_REAL* const fn, 
          const CCTK_INT i, const CCTK_INT j, const CCTK_INT k,
          const CCTK_INT stride, const CCTK_INT* restrict const kmin,
          const CCTK_INT* restrict const kmax,
          const CCTK_REAL* restrict const q, const CCTK_REAL hh)
{
  CCTK_REAL ddval=0.0;

#ifdef CHECK_FOR_POISON
    Check_For_Poison(i, j, k, kmin[k], kmax[k], stride, 2);
#endif

  for (int a=kmin[k]; a<=kmax[k]; ++a)
    ddval += q[a+k*stride] * fn[CCTK_GFINDEX3D(cctkGH, i, j, a)];
  return hh*ddval;
}

#endif




#ifdef HARDCODED_STENCILS_2

static inline __attribute__((always_inline)) CCTK_REAL
calc_da(const cGH* restrict const cctkGH, const CCTK_REAL* const fn,
        const CCTK_INT i, const CCTK_INT j, CCTK_INT k, const CCTK_INT stride,
        const CCTK_INT* restrict const imin,
        const CCTK_INT* restrict const imax, const CCTK_REAL* restrict const q,
        const CCTK_REAL h)
{
  if (i < 1 || i >= cctkGH->cctk_lsh[0]-1) return 0.0;
  return 0.5*h*(fn[CCTK_GFINDEX3D(cctkGH, i+1, j, k)]
                - fn[CCTK_GFINDEX3D(cctkGH, i-1, j, k)]);
}



static inline __attribute__((always_inline)) CCTK_REAL
calc_db(const cGH* restrict const cctkGH, const CCTK_REAL* const fn,
        const CCTK_INT i, const CCTK_INT j, const CCTK_INT k,
        const CCTK_INT stride, const CCTK_INT* restrict const jmin,
        const CCTK_INT* restrict const jmax, const CCTK_REAL* restrict const q,
        const CCTK_REAL h)
{
  if (j < 1 || j >= cctkGH->cctk_lsh[1]-1) return 0.0;
  return 0.5*h*(fn[CCTK_GFINDEX3D(cctkGH, i, j+1, k)]
                - fn[CCTK_GFINDEX3D(cctkGH, i, j-1, k)]);
}


static inline __attribute__((always_inline)) CCTK_REAL
calc_dc(const cGH* restrict const cctkGH, const CCTK_REAL* const fn,
        const CCTK_INT i, const CCTK_INT j, const CCTK_INT k,
        const CCTK_INT stride, const CCTK_INT* restrict const kmin,
        const CCTK_INT* restrict const kmax, const CCTK_REAL* restrict const q,
        const CCTK_REAL h)
{
  if (k < 1 || k >= cctkGH->cctk_lsh[2]-1) return 0.0;
  return 0.5*h*(fn[CCTK_GFINDEX3D(cctkGH, i, j, k+1)]
                - fn[CCTK_GFINDEX3D(cctkGH, i, j, k-1)]);
}


static inline __attribute__((always_inline)) CCTK_REAL
calc_dada(const cGH* restrict const cctkGH, const CCTK_REAL* const fn,
          const CCTK_INT i, const CCTK_INT j, const CCTK_INT k,
          const CCTK_INT stride, const CCTK_INT* restrict const imin,
          const CCTK_INT* restrict const imax,
          const CCTK_REAL* restrict const q, const CCTK_REAL hh)
{
  if (i < 1 || i >= cctkGH->cctk_lsh[0]-1) return 0.0;
  return hh*(fn[CCTK_GFINDEX3D(cctkGH, i+1, j, k)]
             + fn[CCTK_GFINDEX3D(cctkGH, i-1, j, k)]
             - 2.0*fn[CCTK_GFINDEX3D(cctkGH, i, j, k)]);
}

static inline __attribute__((always_inline)) CCTK_REAL
calc_dadb(const cGH* restrict const cctkGH, const CCTK_REAL* const fn,
          const CCTK_INT i, const CCTK_INT j, const CCTK_INT k,
          const CCTK_INT istride, const CCTK_INT jstride,
          const CCTK_INT* restrict const imin,
          const CCTK_INT* restrict const imax,
          const CCTK_INT* restrict const jmin,
          const CCTK_INT* restrict const jmax,
          const CCTK_REAL* restrict const iq,
          const CCTK_REAL* restrict const jq, const CCTK_REAL hh)
{
  if (i < 1 || i >= cctkGH->cctk_lsh[0]-1 ||
      j < 1 || j >= cctkGH->cctk_lsh[1]-1) return 0.0;
  return 0.25*hh*(fn[CCTK_GFINDEX3D(cctkGH, i+1, j+1, k)]
                  - fn[CCTK_GFINDEX3D(cctkGH, i+1, j-1, k)]
                  - fn[CCTK_GFINDEX3D(cctkGH, i-1, j+1, k)]
                  + fn[CCTK_GFINDEX3D(cctkGH, i-1, j-1, k)]);
}




static inline __attribute__((always_inline)) CCTK_REAL
calc_dadc(const cGH* restrict const cctkGH, const CCTK_REAL* const fn,
          const CCTK_INT i, const CCTK_INT j, const CCTK_INT k,
          const CCTK_INT istride, const CCTK_INT kstride,
          const CCTK_INT* restrict const imin,
          const CCTK_INT* restrict const imax,
          const CCTK_INT* restrict const kmin,
          const CCTK_INT* restrict const kmax,
          const CCTK_REAL* restrict const iq,
          const CCTK_REAL* restrict const kq, const CCTK_REAL hh)
{
  if (i < 1 || i >= cctkGH->cctk_lsh[0]-1 ||
      k < 1 || k >= cctkGH->cctk_lsh[2]-1) return 0.0;
  return 0.25*hh*(fn[CCTK_GFINDEX3D(cctkGH, i+1, j, k+1)]
                  - fn[CCTK_GFINDEX3D(cctkGH, i+1, j, k-1)]
                  - fn[CCTK_GFINDEX3D(cctkGH, i-1, j, k+1)]
                  + fn[CCTK_GFINDEX3D(cctkGH, i-1, j, k-1)]);
}




static inline __attribute__((always_inline)) CCTK_REAL
calc_dbdb(const cGH* restrict const cctkGH, const CCTK_REAL* const fn,
          const CCTK_INT i, const CCTK_INT j, const CCTK_INT k,
          const CCTK_INT stride, const CCTK_INT* restrict const jmin,
          const CCTK_INT* restrict const jmax,
          const CCTK_REAL* restrict const q, const CCTK_REAL hh)
{
  if (j < 1 || j >= cctkGH->cctk_lsh[1]-1) return 0.0;
  return hh*(fn[CCTK_GFINDEX3D(cctkGH, i, j+1, k)]
             + fn[CCTK_GFINDEX3D(cctkGH, i, j-1, k)]
             - 2.0*fn[CCTK_GFINDEX3D(cctkGH, i, j, k)]);
}



static inline __attribute__((always_inline)) CCTK_REAL
calc_dbdc(const cGH* restrict const cctkGH, const CCTK_REAL* const fn,
          const CCTK_INT i, const CCTK_INT j, const CCTK_INT k,
          const CCTK_INT jstride, const CCTK_INT kstride,
          const CCTK_INT* restrict const jmin,
          const CCTK_INT* restrict const jmax,
          const CCTK_INT* restrict const kmin,
          const CCTK_INT* restrict const kmax,
          const CCTK_REAL* restrict const jq,
          const CCTK_REAL* restrict const kq, const CCTK_REAL hh)
{
  if (j < 1 || j >= cctkGH->cctk_lsh[1]-1 ||
      k < 1 || k >= cctkGH->cctk_lsh[2]-1) return 0.0;
  return 0.25*hh*(fn[CCTK_GFINDEX3D(cctkGH, i, j+1, k+1)]
                  - fn[CCTK_GFINDEX3D(cctkGH, i, j+1, k-1)]
                  - fn[CCTK_GFINDEX3D(cctkGH, i, j-1, k+1)]
                  + fn[CCTK_GFINDEX3D(cctkGH, i, j-1, k-1)]);
}



static inline __attribute__((always_inline)) CCTK_REAL
calc_dcdc(const cGH* restrict const cctkGH, const CCTK_REAL* const fn,
          const CCTK_INT i, const CCTK_INT j, const CCTK_INT k,
          const CCTK_INT stride, const CCTK_INT* restrict const kmin,
          const CCTK_INT* restrict const kmax,
          const CCTK_REAL* restrict const q, const CCTK_REAL hh)
{
  if (k < 1 || k >= cctkGH->cctk_lsh[2]-1) return 0.0;
  return hh*(fn[CCTK_GFINDEX3D(cctkGH, i, j, k+1)]
             + fn[CCTK_GFINDEX3D(cctkGH, i, j, k-1)]
             - 2.0*fn[CCTK_GFINDEX3D(cctkGH, i, j, k)]);
}


#endif

#ifdef HARDCODED_STENCILS_4

static inline __attribute__((always_inline)) CCTK_REAL
calc_da(const cGH* restrict const cctkGH, const CCTK_REAL* const fn,
        const CCTK_INT i, const CCTK_INT j, const CCTK_INT k,
        const CCTK_INT stride, const CCTK_INT* restrict const imin,
        const CCTK_INT* restrict const imax, const CCTK_REAL* restrict const q,
        const CCTK_REAL h)
{
  if (i < 2 || i >= cctkGH->cctk_lsh[0]-2) return 0.0;
  return (h*(1.0/12.0))*(-fn[CCTK_GFINDEX3D(cctkGH, i+2, j, k)]
                         + fn[CCTK_GFINDEX3D(cctkGH, i-2, j, k)]
                         + 8.0*(fn[CCTK_GFINDEX3D(cctkGH, i+1, j, k)]
                                - fn[CCTK_GFINDEX3D(cctkGH, i-1, j, k)]));
}



static inline __attribute__((always_inline)) CCTK_REAL
calc_db(const cGH* restrict const cctkGH, const CCTK_REAL* const fn,
        const CCTK_INT i, const CCTK_INT j, const CCTK_INT k,
        const CCTK_INT stride, const CCTK_INT* restrict const jmin,
        const CCTK_INT* restrict const jmax, const CCTK_REAL* restrict const q,
        const CCTK_REAL h)
{
  if (j < 2 || j >= cctkGH->cctk_lsh[1]-2) return 0.0;
  return (h*(1.0/12.0))*(-fn[CCTK_GFINDEX3D(cctkGH, i, j+2, k)]
                         + fn[CCTK_GFINDEX3D(cctkGH, i, j-2, k)]
                         + 8.0*(fn[CCTK_GFINDEX3D(cctkGH, i, j+1, k)]
                                - fn[CCTK_GFINDEX3D(cctkGH, i, j-1, k)]));
}


static inline __attribute__((always_inline)) CCTK_REAL
calc_dc(const cGH* restrict const cctkGH, const CCTK_REAL* const fn,
        const CCTK_INT i, const CCTK_INT j, const CCTK_INT k,
        const CCTK_INT stride, const CCTK_INT* restrict const kmin,
        const CCTK_INT* restrict const kmax, const CCTK_REAL* restrict const q,
        const CCTK_REAL h)
{
  if (k < 2 || k >= cctkGH->cctk_lsh[2]-2) return 0.0;
  return (h*(1.0/12.0))*(-fn[CCTK_GFINDEX3D(cctkGH, i, j, k+2)]
                         + fn[CCTK_GFINDEX3D(cctkGH, i, j, k-2)]
                         + 8.0*(fn[CCTK_GFINDEX3D(cctkGH, i, j, k+1)]
                                - fn[CCTK_GFINDEX3D(cctkGH, i, j, k-1)]));
}


static inline __attribute__((always_inline)) CCTK_REAL
calc_dada(const cGH* restrict const cctkGH, const CCTK_REAL* const fn,
          const CCTK_INT i, const CCTK_INT j, const CCTK_INT k,
          const CCTK_INT stride, const CCTK_INT* restrict const imin,
          const CCTK_INT* restrict const imax,
          const CCTK_REAL* restrict const q, const CCTK_REAL hh)
{
  if (i < 2 || i >= cctkGH->cctk_lsh[0]-2) return 0.0;
  return (hh*(1.0/12.0))*(-(fn[CCTK_GFINDEX3D(cctkGH, i+2, j, k)]
                            + fn[CCTK_GFINDEX3D(cctkGH, i-2, j, k)])
                          + 16.0*(fn[CCTK_GFINDEX3D(cctkGH, i+1, j, k)]
                                  + fn[CCTK_GFINDEX3D(cctkGH, i-1, j, k)])
                          - 30.0*fn[CCTK_GFINDEX3D(cctkGH, i, j, k)]);
}

static inline __attribute__((always_inline)) CCTK_REAL
calc_dadb(const cGH* restrict const cctkGH, const CCTK_REAL* const fn,
          const CCTK_INT i, const CCTK_INT j, const CCTK_INT k,
          const CCTK_INT istride, const CCTK_INT jstride,
          const CCTK_INT* restrict const imin,
          const CCTK_INT* restrict const imax,
          const CCTK_INT* restrict const jmin,
          const CCTK_INT* restrict const jmax,
          const CCTK_REAL* restrict const iq,
          const CCTK_REAL* restrict const jq, const CCTK_REAL hh)
{
  if (i < 2 || i >= cctkGH->cctk_lsh[0]-2 ||
      j < 2 || j >= cctkGH->cctk_lsh[1]-2) return 0.0;
  return (hh*(1.0/144.0))*( (fn[CCTK_GFINDEX3D(cctkGH, i+2, j+2, k)]
                             - fn[CCTK_GFINDEX3D(cctkGH, i+2, j-2, k)]
                             - fn[CCTK_GFINDEX3D(cctkGH, i-2, j+2, k)]
                             + fn[CCTK_GFINDEX3D(cctkGH, i-2, j-2, k)])
                            + 8.0*(-fn[CCTK_GFINDEX3D(cctkGH, i+2, j+1, k)]
                                   +fn[CCTK_GFINDEX3D(cctkGH, i+2, j-1, k)]
                                   -fn[CCTK_GFINDEX3D(cctkGH, i+1, j+2, k)]
                                   +fn[CCTK_GFINDEX3D(cctkGH, i+1, j-2, k)]
                                   -fn[CCTK_GFINDEX3D(cctkGH, i-1, j-2, k)]
                                   +fn[CCTK_GFINDEX3D(cctkGH, i-1, j+2, k)]
                                   -fn[CCTK_GFINDEX3D(cctkGH, i-2, j-1, k)]
                                   +fn[CCTK_GFINDEX3D(cctkGH, i-2, j+1, k)])
                            + 64.0*(fn[CCTK_GFINDEX3D(cctkGH, i+1, j+1, k)]
                                    - fn[CCTK_GFINDEX3D(cctkGH, i+1, j-1, k)]
                                    - fn[CCTK_GFINDEX3D(cctkGH, i-1, j+1, k)]
                                    + fn[CCTK_GFINDEX3D(cctkGH, i-1, j-1, k)])
                            );
}




static inline __attribute__((always_inline)) CCTK_REAL
calc_dadc(const cGH* restrict const cctkGH, const CCTK_REAL* const fn,
          const CCTK_INT i, const CCTK_INT j, const CCTK_INT k,
          const CCTK_INT istride, const CCTK_INT kstride,
          const CCTK_INT* restrict const imin,
          const CCTK_INT* restrict const imax,
          const CCTK_INT* restrict const kmin,
          const CCTK_INT* restrict const kmax,
          const CCTK_REAL* restrict const iq,
          const CCTK_REAL* restrict const kq, const CCTK_REAL hh)
{
  if (i < 2 || i >= cctkGH->cctk_lsh[0]-2 ||
      k < 2 || k >= cctkGH->cctk_lsh[2]-2) return 0.0;
  return (hh*(1.0/144.0))*( (fn[CCTK_GFINDEX3D(cctkGH, i+2, j, k+2)]
                             - fn[CCTK_GFINDEX3D(cctkGH, i+2, j, k-2)]
                             - fn[CCTK_GFINDEX3D(cctkGH, i-2, j, k+2)]
                             + fn[CCTK_GFINDEX3D(cctkGH, i-2, j, k-2)])
                            + 8.0*(-fn[CCTK_GFINDEX3D(cctkGH, i+2, j, k+1)]
                                   +fn[CCTK_GFINDEX3D(cctkGH, i+2, j, k-1)]
                                   -fn[CCTK_GFINDEX3D(cctkGH, i+1, j, k+2)]
                                   +fn[CCTK_GFINDEX3D(cctkGH, i+1, j, k-2)]
                                   -fn[CCTK_GFINDEX3D(cctkGH, i-1, j, k-2)]
                                   +fn[CCTK_GFINDEX3D(cctkGH, i-1, j, k+2)]
                                   -fn[CCTK_GFINDEX3D(cctkGH, i-2, j, k-1)]
                                   +fn[CCTK_GFINDEX3D(cctkGH, i-2, j, k+1)])
                            + 64.0*(fn[CCTK_GFINDEX3D(cctkGH, i+1, j, k+1)]
                                    - fn[CCTK_GFINDEX3D(cctkGH, i+1, j, k-1)]
                                    - fn[CCTK_GFINDEX3D(cctkGH, i-1, j, k+1)]
                                    + fn[CCTK_GFINDEX3D(cctkGH, i-1, j, k-1)])
                            );
}




static inline __attribute__((always_inline)) CCTK_REAL
calc_dbdb(const cGH* restrict const cctkGH, const CCTK_REAL* const fn,
          const CCTK_INT i, const CCTK_INT j, const CCTK_INT k,
          const CCTK_INT stride, const CCTK_INT* restrict const jmin,
          const CCTK_INT* restrict const jmax,
          const CCTK_REAL* restrict const q, const CCTK_REAL hh)
{
  if (j < 2 || j >= cctkGH->cctk_lsh[1]-2) return 0.0;
  return (hh*(1.0/12.0))*(-(fn[CCTK_GFINDEX3D(cctkGH, i, j+2, k)]
                            + fn[CCTK_GFINDEX3D(cctkGH, i, j-2, k)])
                          + 16.0*(fn[CCTK_GFINDEX3D(cctkGH, i, j+1, k)]
                                  + fn[CCTK_GFINDEX3D(cctkGH, i, j-1, k)])
                          - 30.0*fn[CCTK_GFINDEX3D(cctkGH, i, j, k)]);
}



static inline __attribute__((always_inline)) CCTK_REAL
calc_dbdc(const cGH* restrict const cctkGH, const CCTK_REAL* const fn,
          const CCTK_INT i, const CCTK_INT j, const CCTK_INT k,
          const CCTK_INT jstride, const CCTK_INT kstride,
          const CCTK_INT* restrict const jmin,
          const CCTK_INT* restrict const jmax,
          const CCTK_INT* restrict const kmin,
          const CCTK_INT* restrict const kmax,
          const CCTK_REAL* restrict const jq,
          const CCTK_REAL* restrict const kq, const CCTK_REAL hh)
{
  if (j < 2 || j >= cctkGH->cctk_lsh[1]-2 ||
      k < 2 || k >= cctkGH->cctk_lsh[2]-2) return 0.0;
  return (hh*(1.0/144.0))*( (fn[CCTK_GFINDEX3D(cctkGH, i, j+2, k+2)]
                             - fn[CCTK_GFINDEX3D(cctkGH, i, j+2, k-2)]
                             - fn[CCTK_GFINDEX3D(cctkGH, i, j-2, k+2)]
                             + fn[CCTK_GFINDEX3D(cctkGH, i, j-2, k-2)])
                            + 8.0*(-fn[CCTK_GFINDEX3D(cctkGH, i, j+2, k+1)]
                                   +fn[CCTK_GFINDEX3D(cctkGH, i, j+2, k-1)]
                                   -fn[CCTK_GFINDEX3D(cctkGH, i, j+1, k+2)]
                                   +fn[CCTK_GFINDEX3D(cctkGH, i, j+1, k-2)]
                                   -fn[CCTK_GFINDEX3D(cctkGH, i, j-1, k-2)]
                                   +fn[CCTK_GFINDEX3D(cctkGH, i, j-1, k+2)]
                                   -fn[CCTK_GFINDEX3D(cctkGH, i, j-2, k-1)]
                                   +fn[CCTK_GFINDEX3D(cctkGH, i, j-2, k+1)])
                            + 64.0*(fn[CCTK_GFINDEX3D(cctkGH, i, j+1, k+1)]
                                    - fn[CCTK_GFINDEX3D(cctkGH, i, j+1, k-1)]
                                    - fn[CCTK_GFINDEX3D(cctkGH, i, j-1, k+1)]
                                    + fn[CCTK_GFINDEX3D(cctkGH, i, j-1, k-1)])
                            );
}



static inline __attribute__((always_inline)) CCTK_REAL
calc_dcdc(const cGH* restrict const cctkGH, const CCTK_REAL* const fn,
          const CCTK_INT i, const CCTK_INT j, const CCTK_INT k,
          const CCTK_INT stride, const CCTK_INT* restrict const kmin,
          const CCTK_INT* restrict const kmax,
          const CCTK_REAL* restrict const q, const CCTK_REAL hh)
{
  if (k < 2 || k >= cctkGH->cctk_lsh[2]-2) return 0.0;
  return (hh*(1.0/12.0))*(-(fn[CCTK_GFINDEX3D(cctkGH, i, j, k+2)]
                            + fn[CCTK_GFINDEX3D(cctkGH, i, j, k-2)])
                          + 16.0*(fn[CCTK_GFINDEX3D(cctkGH, i, j, k+1)]
                                  + fn[CCTK_GFINDEX3D(cctkGH, i, j, k-1)])
                          - 30.0*fn[CCTK_GFINDEX3D(cctkGH, i, j, k)]);
}


#endif


//---------------------------------------------------------------------------
//
// Global first derivatives as a linear combination of local derivatives,
// i.e.  \hat{\partial}_i = J^j_i \partial_j
//  in global coordinates (x,y,z)
//
//---------------------------------------------------------------------------



static inline __attribute__((always_inline)) CCTK_REAL
g_diff_dx(cGH const * restrict const cctkGH, const CCTK_REAL* const fn,
          const CCTK_INT general_coordinates,
          const CCTK_REAL dadx, const CCTK_REAL dbdx, const CCTK_REAL dcdx,
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
  CCTK_REAL dval;

  if (general_coordinates)
    {
      dval = dadx * calc_da(cctkGH, fn, i, j, k, ni, imin, imax, qa, ha)
        + dbdx * calc_db(cctkGH, fn, i, j, k, nj, jmin, jmax, qb, hb)
        + dcdx * calc_dc(cctkGH, fn, i, j, k, nk, kmin, kmax, qc, hc);
    }
  else
    {
      dval = calc_da(cctkGH, fn, i, j, k, ni, imin, imax, qa, ha);
    }
  
  return dval;
}



static inline __attribute__((always_inline)) CCTK_REAL
g_diff_dy(cGH const * restrict const cctkGH, const CCTK_REAL* const fn, 
          const CCTK_INT general_coordinates,
          const CCTK_REAL dady, const CCTK_REAL dbdy, const CCTK_REAL dcdy,
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
  CCTK_REAL dval;
  
  if (general_coordinates)
    {
      dval = dady * calc_da(cctkGH, fn, i, j, k, ni, imin, imax, qa, ha)
        + dbdy * calc_db(cctkGH, fn, i, j, k, nj, jmin, jmax, qb, hb)
        + dcdy * calc_dc(cctkGH, fn, i, j, k, nk, kmin, kmax, qc, hc);
    }
  else
    {
      dval = calc_db(cctkGH, fn, i, j, k, nj, jmin, jmax, qb, hb);
    }
  
  return dval;
}



static inline __attribute__((always_inline)) CCTK_REAL
g_diff_dz(cGH const * restrict const cctkGH, const CCTK_REAL* const fn, 
          const CCTK_INT general_coordinates,
          const CCTK_REAL dadz, const CCTK_REAL dbdz, const CCTK_REAL dcdz,
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
  CCTK_REAL dval;
  
  if (general_coordinates)
    {
      dval = dadz * calc_da(cctkGH, fn, i, j, k, ni, imin, imax, qa, ha)
        + dbdz * calc_db(cctkGH, fn, i, j, k, nj, jmin, jmax, qb, hb)
        + dcdz * calc_dc(cctkGH, fn, i, j, k, nk, kmin, kmax, qc, hc);
    }
  else
    {
      dval = calc_dc(cctkGH, fn, i, j, k, nk, kmin, kmax, qc, hc);
    }
  
  return dval;
}




/*
 * Returns the derivative of all directions with one function call
 */

static inline __attribute__((always_inline)) void
g_diff(cGH const * restrict const cctkGH, const CCTK_REAL* const fn, 
       CCTK_REAL* restrict const dfn_1,
       CCTK_REAL* restrict const dfn_2,
       CCTK_REAL* restrict const dfn_3,
       const CCTK_INT general_coordinates,
       const CCTK_REAL dadx, const CCTK_REAL dbdx, const CCTK_REAL dcdx,
       const CCTK_REAL dady, const CCTK_REAL dbdy, const CCTK_REAL dcdy,
       const CCTK_REAL dadz, const CCTK_REAL dbdz, const CCTK_REAL dcdz,
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
  CCTK_REAL dtemp[3];
  
  (*dfn_1) = calc_da(cctkGH, fn, i, j, k, ni, imin, imax, qa, ha);
  (*dfn_2) = calc_db(cctkGH, fn, i, j, k, nj, jmin, jmax, qb, hb);
  (*dfn_3) = calc_dc(cctkGH, fn, i, j, k, nk, kmin, kmax, qc, hc);

  if (general_coordinates)
    {
      dtemp[0] = dadx * (*dfn_1)
        + dbdx * (*dfn_2)
        + dcdx * (*dfn_3);
                  
      dtemp[1] = dady * (*dfn_1)
        + dbdy * (*dfn_2)
        + dcdy * (*dfn_3);
      
      dtemp[2] = dadz * (*dfn_1)
        + dbdz * (*dfn_2)
        + dcdz * (*dfn_3);
      
      (*dfn_1) = dtemp[0];
      (*dfn_2) = dtemp[1];
      (*dfn_3) = dtemp[2];
    }
  
  return;
}




//--------------------------------------------------------------------------
//
// Global second derivatives as a linear combination of local derivatives,
// i.e.  \hat{\partial}_i\hat{\partial}_j = 
//                      J^k_i \partial_k ( J^l_j ) \partial_l 
//                                          + J^k_i J^l_j \partial_k \partial_l
//  in global coordinates (x,y,z)
//
//--------------------------------------------------------------------------



static inline __attribute__((always_inline)) CCTK_REAL
g_diff_dxdx(cGH const * restrict const cctkGH, const CCTK_REAL* const fn,
            const CCTK_INT general_coordinates,
            const CCTK_REAL dadx, const CCTK_REAL dbdx, const CCTK_REAL dcdx,
            const CCTK_REAL ddadxdx, const CCTK_REAL ddbdxdx,
            const CCTK_REAL ddcdxdx,
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
            const CCTK_REAL ha, const CCTK_REAL hb, const CCTK_REAL hc,
            const CCTK_REAL hhaa, const CCTK_REAL hhab, const CCTK_REAL hhac,
            const CCTK_REAL hhbb,const CCTK_REAL hhbc, const CCTK_REAL hhcc)
{
  CCTK_REAL ddval;

  if (general_coordinates)
    {
      const CCTK_REAL daval = calc_da(cctkGH, fn, i, j, k, ni, imin, imax,
                                      qa, ha);
      const CCTK_REAL dbval = calc_db(cctkGH, fn, i, j, k, nj, jmin, jmax,
                                      qb, hb);
      const CCTK_REAL dcval = calc_dc(cctkGH, fn, i, j, k, nk, kmin, kmax,
                                      qc, hc);
      
      const CCTK_REAL dadaval = calc_dada(cctkGH, fn, i, j, k, ni,
                                          imin2, imax2, qa2, hhaa);
      const CCTK_REAL dadbval = calc_dadb(cctkGH, fn, i, j, k, ni, nj, imin,
                                          imax, jmin, jmax, qa, qb, hhab);
      const CCTK_REAL dadcval = calc_dadc(cctkGH, fn, i, j, k, ni, nk, imin,
                                          imax, kmin, kmax, qa, qc, hhac);
      
      const CCTK_REAL dbdbval = calc_dbdb(cctkGH, fn, i, j, k, nj, jmin2,
                                          jmax2, qb2, hhbb);
      const CCTK_REAL dbdcval = calc_dbdc(cctkGH, fn, i, j, k, nj, nk, jmin,
                                          jmax, kmin, kmax, qb, qc, hhbc);
      
      const CCTK_REAL dcdcval = calc_dcdc(cctkGH, fn, i, j, k, nk, kmin2,
                                          kmax2, qc2, hhcc);
      
      ddval =   ddadxdx * daval 
        + ddbdxdx * dbval 
        + ddcdxdx * dcval;

      ddval +=     dadx * dadx * dadaval
        + dbdx * dbdx * dbdbval
        + dcdx * dcdx * dcdcval
        + 2 * dadx * dbdx * dadbval
        + 2 * dadx * dcdx * dadcval
        + 2 * dbdx * dcdx * dbdcval;
    }
  else 
    {
      ddval = calc_dada(cctkGH, fn, i, j, k, ni, imin2, imax2, qa2, hhaa);
    }
  
  return ddval;
}






static inline __attribute__((always_inline)) CCTK_REAL
g_diff_dxdy(cGH const * restrict const cctkGH, const CCTK_REAL* const fn,
            const CCTK_INT general_coordinates,
            const CCTK_REAL dadx, const CCTK_REAL dbdx, const CCTK_REAL dcdx,
            const CCTK_REAL dady, const CCTK_REAL dbdy, const CCTK_REAL dcdy,
            const CCTK_REAL ddadxdy, const CCTK_REAL ddbdxdy,
            const CCTK_REAL ddcdxdy,
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
            const CCTK_REAL ha, const CCTK_REAL hb, const CCTK_REAL hc,
            const CCTK_REAL hhaa, const CCTK_REAL hhab, const CCTK_REAL hhac,
            const CCTK_REAL hhbb,const CCTK_REAL hhbc, const CCTK_REAL hhcc)
{
  CCTK_REAL ddval;
  
  if (general_coordinates)
    {
      const CCTK_REAL daval = calc_da(cctkGH, fn, i, j, k, ni, imin, imax,
                                      qa, ha);
      const CCTK_REAL dbval = calc_db(cctkGH, fn, i, j, k, nj, jmin, jmax,
                                      qb, hb);
      const CCTK_REAL dcval = calc_dc(cctkGH, fn, i, j, k, nk, kmin, kmax,
                                      qc, hc);
      
      const CCTK_REAL dadaval = calc_dada(cctkGH, fn, i, j, k, ni, imin2,
                                          imax2, qa2, hhaa);
      const CCTK_REAL dadbval = calc_dadb(cctkGH, fn, i, j, k, ni, nj, imin,
                                          imax, jmin, jmax, qa, qb, hhab);
      const CCTK_REAL dadcval = calc_dadc(cctkGH, fn, i, j, k, ni, nk, imin,
                                          imax, kmin, kmax, qa, qc, hhac);
      
      const CCTK_REAL dbdbval = calc_dbdb(cctkGH, fn, i, j, k, nj, jmin2,
                                          jmax2, qb2, hhbb);
      const CCTK_REAL dbdcval = calc_dbdc(cctkGH, fn, i, j, k, nj, nk, jmin,
                                          jmax, kmin, kmax, qb, qc, hhbc);
      
      const CCTK_REAL dcdcval = calc_dcdc(cctkGH, fn, i, j, k, nk, kmin2,
                                          kmax2, qc2, hhcc);
      
      ddval =   ddadxdy * daval 
        + ddbdxdy * dbval 
        + ddcdxdy * dcval;
      
      ddval +=     dadx * dady * dadaval
        + dbdx * dbdy * dbdbval
        + dcdx * dcdy * dcdcval
        + (dadx * dbdy + dady * dbdx) * dadbval
        + (dadx * dcdy + dady * dcdx) * dadcval
        + (dbdx * dcdy + dbdy * dcdx) * dbdcval;
    }
  else
    {
      ddval = calc_dadb(cctkGH, fn, i, j, k, ni, nj, imin, imax, jmin, jmax, qa, qb, hhab);
    }
  
  return ddval;
}





static inline __attribute__((always_inline)) CCTK_REAL
g_diff_dxdz(cGH const * restrict const cctkGH, const CCTK_REAL* const fn,
            const CCTK_INT general_coordinates,
            const CCTK_REAL dadx, const CCTK_REAL dbdx, const CCTK_REAL dcdx,
            const CCTK_REAL dadz, const CCTK_REAL dbdz, const CCTK_REAL dcdz,
            const CCTK_REAL ddadxdz, const CCTK_REAL ddbdxdz,
            const CCTK_REAL ddcdxdz,
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
            const CCTK_REAL ha, const CCTK_REAL hb, const CCTK_REAL hc,
            const CCTK_REAL hhaa, const CCTK_REAL hhab, const CCTK_REAL hhac,
            const CCTK_REAL hhbb,const CCTK_REAL hhbc, const CCTK_REAL hhcc)
{
  CCTK_REAL ddval;
  
  if (general_coordinates)
    {
      const CCTK_REAL daval = calc_da(cctkGH, fn, i, j, k, ni, imin, imax,
                                      qa, ha);
      const CCTK_REAL dbval = calc_db(cctkGH, fn, i, j, k, nj, jmin, jmax,
                                      qb, hb);
      const CCTK_REAL dcval = calc_dc(cctkGH, fn, i, j, k, nk, kmin, kmax,
                                      qc, hc);
      
      const CCTK_REAL dadaval = calc_dada(cctkGH, fn, i, j, k, ni, imin2,
                                          imax2, qa2, hhaa);
      const CCTK_REAL dadbval = calc_dadb(cctkGH, fn, i, j, k, ni, nj, imin,
                                          imax, jmin, jmax, qa, qb, hhab);
      const CCTK_REAL dadcval = calc_dadc(cctkGH, fn, i, j, k, ni, nk, imin,
                                          imax, kmin, kmax, qa, qc, hhac);
      
      const CCTK_REAL dbdbval = calc_dbdb(cctkGH, fn, i, j, k, nj, jmin2,
                                          jmax2, qb2, hhbb);
      const CCTK_REAL dbdcval = calc_dbdc(cctkGH, fn, i, j, k, nj, nk, jmin,
                                          jmax, kmin, kmax, qb, qc, hhbc);
      
      const CCTK_REAL dcdcval = calc_dcdc(cctkGH, fn, i, j, k, nk, kmin2,
                                          kmax2, qc2, hhcc);
      
      // first term
      ddval =   ddadxdz * daval 
        + ddbdxdz * dbval 
        + ddcdxdz * dcval;
      
      // second term
      ddval +=     dadx * dadz * dadaval
        + dbdx * dbdz * dbdbval
        + dcdx * dcdz * dcdcval
        + (dadx * dbdz + dadz * dbdx) * dadbval
        + (dadx * dcdz + dadz * dcdx) * dadcval
        + (dbdx * dcdz + dbdz * dcdx) * dbdcval;
    }
  else
    {
      ddval = calc_dadc(cctkGH, fn, i, j, k, ni, nk, imin, imax, kmin, kmax,
                        qa, qc, hhac);
    }
  
  return ddval;
}









static inline __attribute__((always_inline)) CCTK_REAL
g_diff_dydy(cGH const * restrict const cctkGH, const CCTK_REAL* const fn,
            const CCTK_INT general_coordinates,
            const CCTK_REAL dady, const CCTK_REAL dbdy, const CCTK_REAL dcdy,
            const CCTK_REAL ddadydy, const CCTK_REAL ddbdydy,
            const CCTK_REAL ddcdydy,
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
            const CCTK_REAL ha, const CCTK_REAL hb, const CCTK_REAL hc,
            const CCTK_REAL hhaa, const CCTK_REAL hhab,
            const CCTK_REAL hhac, const CCTK_REAL hhbb,const CCTK_REAL hhbc,
            const CCTK_REAL hhcc)
{
  CCTK_REAL ddval;
  
  if (general_coordinates)
    {
      const CCTK_REAL daval = calc_da(cctkGH, fn, i, j, k, ni, imin, imax,
                                      qa, ha);
      const CCTK_REAL dbval = calc_db(cctkGH, fn, i, j, k, nj, jmin, jmax,
                                      qb, hb);
      const CCTK_REAL dcval = calc_dc(cctkGH, fn, i, j, k, nk, kmin, kmax,
                                      qc, hc);
      
      const CCTK_REAL dadaval = calc_dada(cctkGH, fn, i, j, k, ni, imin2,
                                          imax2, qa2, hhaa);
      const CCTK_REAL dadbval = calc_dadb(cctkGH, fn, i, j, k, ni, nj, imin,
                                          imax, jmin, jmax, qa, qb, hhab);
      const CCTK_REAL dadcval = calc_dadc(cctkGH, fn, i, j, k, ni, nk, imin,
                                          imax, kmin, kmax, qa, qc, hhac);
      
      const CCTK_REAL dbdbval = calc_dbdb(cctkGH, fn, i, j, k, nj, jmin2,
                                          jmax2, qb2, hhbb);
      const CCTK_REAL dbdcval = calc_dbdc(cctkGH, fn, i, j, k, nj, nk, jmin,
                                          jmax, kmin, kmax, qb, qc, hhbc);
      
      const CCTK_REAL dcdcval = calc_dcdc(cctkGH, fn, i, j, k, nk, kmin2,
                                          kmax2, qc2, hhcc);
      
      ddval =   ddadydy * daval 
        + ddbdydy * dbval 
        + ddcdydy * dcval;
      
      ddval +=     dady * dady * dadaval
        + dbdy * dbdy * dbdbval
        + dcdy * dcdy * dcdcval
        + 2 * dady * dbdy * dadbval
        + 2 * dady * dcdy * dadcval
        + 2 * dbdy * dcdy * dbdcval;
    }
  else
    {
      ddval = calc_dbdb(cctkGH, fn, i, j, k, nj, jmin2, jmax2, qb2, hhbb);
    }
  
  return ddval;
}










static inline __attribute__((always_inline)) CCTK_REAL
g_diff_dydz(cGH const * restrict const cctkGH, const CCTK_REAL* const fn,
            const CCTK_INT general_coordinates,
            const CCTK_REAL dady, const CCTK_REAL dbdy, const CCTK_REAL dcdy,
            const CCTK_REAL dadz, const CCTK_REAL dbdz, const CCTK_REAL dcdz,
            const CCTK_REAL ddadydz, const CCTK_REAL ddbdydz,
            const CCTK_REAL ddcdydz,
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
            const CCTK_REAL ha, const CCTK_REAL hb, const CCTK_REAL hc,
            const CCTK_REAL hhaa, const CCTK_REAL hhab, const CCTK_REAL hhac,
            const CCTK_REAL hhbb,const CCTK_REAL hhbc, const CCTK_REAL hhcc)
{
  CCTK_REAL ddval;
  
  if (general_coordinates)
    {
      const CCTK_REAL daval = calc_da(cctkGH, fn, i, j, k, ni, imin, imax,
                                      qa, ha);
      const CCTK_REAL dbval = calc_db(cctkGH, fn, i, j, k, nj, jmin, jmax,
                                      qb, hb);
      const CCTK_REAL dcval = calc_dc(cctkGH, fn, i, j, k, nk, kmin, kmax,
                                      qc, hc);
      
      const CCTK_REAL dadaval = calc_dada(cctkGH, fn, i, j, k, ni, imin2,
                                          imax2, qa2, hhaa);
      const CCTK_REAL dadbval = calc_dadb(cctkGH, fn, i, j, k, ni, nj, imin,
                                          imax, jmin, jmax, qa, qb, hhab);
      const CCTK_REAL dadcval = calc_dadc(cctkGH, fn, i, j, k, ni, nk, imin,
                                          imax, kmin, kmax, qa, qc, hhac);
      
      const CCTK_REAL dbdbval = calc_dbdb(cctkGH, fn, i, j, k, nj, jmin2,
                                          jmax2, qb2, hhbb);
      const CCTK_REAL dbdcval = calc_dbdc(cctkGH, fn, i, j, k, nj, nk, jmin, 
                                          jmax, kmin, kmax, qb, qc, hhbc);
      
      const CCTK_REAL dcdcval = calc_dcdc(cctkGH, fn, i, j, k, nk, kmin2,
                                          kmax2, qc2, hhcc);

      ddval =   ddadydz * daval 
        + ddbdydz * dbval 
        + ddcdydz * dcval;
      
      ddval +=     dady * dadz * dadaval
        + dbdy * dbdz * dbdbval
        + dcdy * dcdz * dcdcval
        + (dady * dbdz + dadz * dbdy) * dadbval
        + (dady * dcdz + dadz * dcdy) * dadcval
        + (dbdy * dcdz + dbdz * dcdy) * dbdcval;
    }
  else
    {
      ddval = calc_dbdc(cctkGH, fn, i, j, k, ni, nk, imin, imax, kmin, kmax,
                        qa, qc, hhac);
    }
  
  return ddval;
}







static inline __attribute__((always_inline)) CCTK_REAL
g_diff_dzdz(cGH const * restrict const cctkGH, const CCTK_REAL* const fn,
            const CCTK_INT general_coordinates,
            const CCTK_REAL dadz, const CCTK_REAL dbdz, const CCTK_REAL dcdz,
            const CCTK_REAL ddadzdz, const CCTK_REAL ddbdzdz,
            const CCTK_REAL ddcdzdz,
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
            const CCTK_REAL ha, const CCTK_REAL hb, const CCTK_REAL hc,
            const CCTK_REAL hhaa, const CCTK_REAL hhab, const CCTK_REAL hhac,
            const CCTK_REAL hhbb,const CCTK_REAL hhbc, const CCTK_REAL hhcc)
{
  CCTK_REAL ddval;
  
  if (general_coordinates)
    {
      const CCTK_REAL daval = calc_da(cctkGH, fn, i, j, k, ni, imin, imax,
                                      qa, ha);
      const CCTK_REAL dbval = calc_db(cctkGH, fn, i, j, k, nj, jmin, jmax,
                                      qb, hb);
      const CCTK_REAL dcval = calc_dc(cctkGH, fn, i, j, k, nk, kmin, kmax,
                                      qc, hc);
      
      const CCTK_REAL dadaval = calc_dada(cctkGH, fn, i, j, k, ni, imin2,
                                          imax2, qa2, hhaa);
      const CCTK_REAL dadbval = calc_dadb(cctkGH, fn, i, j, k, ni, nj, imin,
                                          imax, jmin, jmax, qa, qb, hhab);
      const CCTK_REAL dadcval = calc_dadc(cctkGH, fn, i, j, k, ni, nk, imin,
                                          imax, kmin, kmax, qa, qc, hhac);
      
      const CCTK_REAL dbdbval = calc_dbdb(cctkGH, fn, i, j, k, nj, jmin2,
                                          jmax2, qb2, hhbb);
      const CCTK_REAL dbdcval = calc_dbdc(cctkGH, fn, i, j, k, nj, nk, jmin,
                                          jmax, kmin, kmax, qb, qc, hhbc);
      
      const CCTK_REAL dcdcval = calc_dcdc(cctkGH, fn, i, j, k, nk, kmin2,
                                          kmax2, qc2, hhcc);
      
      ddval =   ddadzdz * daval 
        + ddbdzdz * dbval 
        + ddcdzdz * dcval;
      
      ddval +=     dadz * dadz * dadaval
        + dbdz * dbdz * dbdbval
        + dcdz * dcdz * dcdcval
        + 2 * dadz * dbdz * dadbval
        + 2 * dadz * dcdz * dadcval
        + 2 * dbdz * dcdz * dbdcval;
    }
  else
    {
      ddval = calc_dcdc(cctkGH, fn, i, j, k, nk, kmin2, kmax2, qc2, hhcc);
    }
  
  return ddval;
}






/*
 * Returns the second derivative of all combinations of all directions
 * with one function call
 */

static inline __attribute__((always_inline)) void
g_diff2(cGH const * restrict const cctkGH, const CCTK_REAL* const fn, 
        CCTK_REAL* restrict const dfn_11, CCTK_REAL* restrict const dfn_12,
        CCTK_REAL* restrict const dfn_13, CCTK_REAL* restrict const dfn_22,
        CCTK_REAL* restrict const dfn_23, CCTK_REAL* restrict const dfn_33,
        const CCTK_INT general_coordinates,
        const CCTK_REAL dadx, const CCTK_REAL dbdx,
        const CCTK_REAL dcdx, const CCTK_REAL dady,
        const CCTK_REAL dbdy, const CCTK_REAL dcdy,
        const CCTK_REAL dadz, const CCTK_REAL dbdz,
        const CCTK_REAL dcdz, const CCTK_REAL ddadxdx,
        const CCTK_REAL ddbdxdx, const CCTK_REAL ddcdxdx,
        const CCTK_REAL ddadxdy, const CCTK_REAL ddbdxdy,
        const CCTK_REAL ddcdxdy, const CCTK_REAL ddadxdz,
        const CCTK_REAL ddbdxdz, const CCTK_REAL ddcdxdz,
        const CCTK_REAL ddadydy, const CCTK_REAL ddbdydy,
        const CCTK_REAL ddcdydy, const CCTK_REAL ddadydz,
        const CCTK_REAL ddbdydz, const CCTK_REAL ddcdydz,
        const CCTK_REAL ddadzdz, const CCTK_REAL ddbdzdz,
        const CCTK_REAL ddcdzdz,
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
        const CCTK_REAL ha, const CCTK_REAL hb, const CCTK_REAL hc,
        const CCTK_REAL hhaa, const CCTK_REAL hhab, const CCTK_REAL hhac,
        const CCTK_REAL hhbb,const CCTK_REAL hhbc, const CCTK_REAL hhcc)
{
  CCTK_REAL dtemp[6];
  
  // calculate local derivatives
  (*dfn_11) = calc_dada(cctkGH, fn, i, j, k, ni, imin2, imax2, qa2, hhaa);
  (*dfn_12) = calc_dadb(cctkGH, fn, i, j, k, ni, nj, imin, imax, jmin, jmax,
                        qa, qb, hhab);
  (*dfn_13) = calc_dadc(cctkGH, fn, i, j, k, ni, nk, imin, imax, kmin, kmax,
                        qa, qc, hhac);
  (*dfn_22) = calc_dbdb(cctkGH, fn, i, j, k, nj, jmin2, jmax2, qb2, hhbb);
  (*dfn_23) = calc_dbdc(cctkGH, fn, i, j, k, nj, nk, jmin, jmax, kmin, kmax,
                        qb, qc, hhbc);
  (*dfn_33) = calc_dcdc(cctkGH, fn, i, j, k, nk, kmin2, kmax2, qc2, hhcc);

  // if Jacobian was given, calculate global derivatives as linear
  // combinations of local derivatives along grid-coordinates
  if (general_coordinates)
    {
      const CCTK_REAL daval = calc_da(cctkGH, fn, i, j, k, ni, imin, imax, qa,
                                      ha);
      const CCTK_REAL dbval = calc_db(cctkGH, fn, i, j, k, nj, jmin, jmax, qb,
                                      hb);
      const CCTK_REAL dcval = calc_dc(cctkGH, fn, i, j, k, nk, kmin, kmax, qc,
                                      hc);
      
      // dxdx
      // first term
      dtemp[0] =   ddadxdx * daval 
        + ddbdxdx * dbval 
        + ddcdxdx * dcval;
      
      // second term
      dtemp[0] +=   dadx * dadx * (*dfn_11)
        + dbdx * dbdx * (*dfn_22)
        + dcdx * dcdx * (*dfn_33)
        + 2 * dadx * dbdx * (*dfn_12)
        + 2 * dadx * dcdx * (*dfn_13)
        + 2 * dbdx * dcdx * (*dfn_23);
      
      // dxdy
      // first term
      dtemp[1] =   ddadxdy * daval 
        + ddbdxdy * dbval 
        + ddcdxdy * dcval;
      
      // second term
      dtemp[1] +=   dadx * dady * (*dfn_11)
        + dbdx * dbdy * (*dfn_22)
        + dcdx * dcdy * (*dfn_33)
        + (dadx * dbdy + dady * dbdx) * (*dfn_12)
        + (dadx * dcdy + dady * dcdx) * (*dfn_13)
        + (dbdx * dcdy + dbdy * dcdx) * (*dfn_23);
      
      // dxdz
      // first term
      dtemp[2] =   ddadxdz * daval 
        + ddbdxdz * dbval 
        + ddcdxdz * dcval;
      
      // second term
      dtemp[2] +=   dadx * dadz * (*dfn_11)
        + dbdx * dbdz * (*dfn_22)
        + dcdx * dcdz * (*dfn_33)
        + (dadx * dbdz + dadz * dbdx) * (*dfn_12)
        + (dadx * dcdz + dadz * dcdx) * (*dfn_13)
        + (dbdx * dcdz + dbdz * dcdx) * (*dfn_23);
      
      // dydy
      // first term
      dtemp[3] =   ddadydy * daval 
        + ddbdydy * dbval 
        + ddcdydy * dcval;
      
      // second term
      dtemp[3] +=   dady * dady * (*dfn_11)
        + dbdy * dbdy * (*dfn_22)
        + dcdy * dcdy * (*dfn_33)
        + 2 * dady * dbdy * (*dfn_12)
        + 2 * dady * dcdy * (*dfn_13)
        + 2 * dbdy * dcdy * (*dfn_23);
      
      // dydz
      // first term
      dtemp[4] =   ddadydz * daval 
        + ddbdydz * dbval 
        + ddcdydz * dcval;
      
      // second term
      dtemp[4] +=   dady * dadz * (*dfn_11)
        + dbdy * dbdz * (*dfn_22)
        + dcdy * dcdz * (*dfn_33)
        + (dady * dbdz + dadz * dbdy) * (*dfn_12)
        + (dady * dcdz + dadz * dcdy) * (*dfn_13)
        + (dbdy * dcdz + dbdz * dcdy) * (*dfn_23);
      
      // dzdz
      // first term
      dtemp[5] =   ddadzdz * daval 
        + ddbdzdz * dbval 
        + ddcdzdz * dcval;
      
      // second term
      dtemp[5] +=   dadz * dadz * (*dfn_11)
        + dbdz * dbdz * (*dfn_22)
        + dcdz * dcdz * (*dfn_33)
        + 2 * dadz * dbdz * (*dfn_12)
        + 2 * dadz * dcdz * (*dfn_13)
        + 2 * dbdz * dcdz * (*dfn_23);
      
      
      (*dfn_11) = dtemp[0];
      (*dfn_12) = dtemp[1];
      (*dfn_13) = dtemp[2];
      (*dfn_22) = dtemp[3];
      (*dfn_23) = dtemp[4];
      (*dfn_33) = dtemp[5];
    }
  
  return;
}




#endif
