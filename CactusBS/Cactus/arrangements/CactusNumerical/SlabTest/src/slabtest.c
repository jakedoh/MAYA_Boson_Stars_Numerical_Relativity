/* $Header$ */

#include <assert.h>
#include <stdio.h>

#include "cctk.h"
#include "cctk_Arguments.h"

#include "Slab.h"



static const char * rcsid = "$Header$";

CCTK_FILEVERSION(TAT_SlabTest_slabtest_c);
  


static void
runtest (CCTK_ARGUMENTS,
	 struct xferinfo const * restrict const info)
{
  DECLARE_CCTK_ARGUMENTS;
  
  int i, j, k;
  
  int mismatch;

# define nvars 2
  const int vartypes[nvars] = {CCTK_VARIABLE_REAL, CCTK_VARIABLE_INT};
  const void *srcdataptrs[nvars] = {gfx, igfx};
  void *dstdataptrs[nvars] = {gfy, igfy};
  
  int ierr;
  
  ierr = Slab_MultiTransfer
    (cctkGH, cctk_dim, info, -1,
     nvars, vartypes, srcdataptrs, vartypes, dstdataptrs);
  assert (!ierr);
  
  mismatch = 0;
  for (k=0; k<cctk_lsh[2]; ++k) {
    for (j=0; j<cctk_lsh[1]; ++j) {
      for (i=0; i<cctk_lsh[0]; ++i) {
	const int ind = CCTK_GFINDEX3D(cctkGH, i,j,k);
	if (gfy[ind] != gfz[ind]) {
	  printf
	    ("Real valued mismatch at [%d,%d,%d]: should be %g, is %g\n",
	     cctk_lbnd[0] + i, cctk_lbnd[1] + j, cctk_lbnd[2] + k,
	     (double) gfz[ind], (double) gfy[ind]);
	  mismatch = 1;
	}
	if (igfy[ind] != igfz[ind]) {
	  printf
	    ("Integer valued mismatch at [%d,%d,%d]: should be %d, is %d\n",
	     cctk_lbnd[0] + i, cctk_lbnd[1] + j, cctk_lbnd[2] + k,
	     (int) igfz[ind], (int) igfy[ind]);
	  mismatch = 1;
	}
      }
    }
  }
  assert (! mismatch);
}



void SlabTest_Test (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  
  struct xferinfo info[3];
  
  int i, j, k;
  int d;
  
  CCTK_INFO ("Testing slabs...");
  
  assert (cctk_dim <= 3);
  for (d=0; d<cctk_dim; ++d) {
    info[d].src.gsh = cctk_gsh[d];
    info[d].src.lbnd = cctk_lbnd[d];
    info[d].src.lsh = cctk_lsh[d];
    info[d].src.ash = cctk_ash[d];
    info[d].src.lbbox = cctk_bbox[2*d];
    info[d].src.ubbox = cctk_bbox[2*d+1];
    info[d].src.nghostzones = cctk_nghostzones[d];
    info[d].dst.gsh = cctk_gsh[d];
    info[d].dst.lbnd = cctk_lbnd[d];
    info[d].dst.lsh = cctk_lsh[d];
    info[d].dst.ash = cctk_ash[d];
    info[d].dst.lbbox = cctk_bbox[2*d];
    info[d].dst.ubbox = cctk_bbox[2*d+1];
    info[d].dst.nghostzones = cctk_nghostzones[d];
  }
  
  
  
  CCTK_INFO ("   Identity");
  
  for (k=0; k<cctk_lsh[2]; ++k) {
    for (j=0; j<cctk_lsh[1]; ++j) {
      for (i=0; i<cctk_lsh[0]; ++i) {
	const int ind = CCTK_GFINDEX3D(cctkGH, i,j,k);
	gfx[ind] = igfx[ind] = 10000 * (cctk_lbnd[0] + i) + 100 * (cctk_lbnd[1] + j) + (cctk_lbnd[2] + k);
	gfy[ind] = igfy[ind] = - (i + 3 * j + 5 * k);
	gfz[ind] = igfz[ind] = 10000 * (cctk_lbnd[0] + i) + 100 * (cctk_lbnd[1] + j) + (cctk_lbnd[2] + k);
      }
    }
  }
  
  assert (cctk_dim <= 3);
  for (d=0; d<cctk_dim; ++d) {
    info[d].src.off = 0;
    info[d].src.len = cctk_gsh[d];
    info[d].src.str = 1;
    info[d].dst.off = 0;
    info[d].dst.len = cctk_gsh[d];
    info[d].dst.str = 1;
    info[d].xpose = d;
    info[d].flip = 0;
  }
  
  runtest (CCTK_PASS_CTOC, info);
  
  
  
  CCTK_INFO ("   Invert x direction");
  
  for (k=0; k<cctk_lsh[2]; ++k) {
    for (j=0; j<cctk_lsh[1]; ++j) {
      for (i=0; i<cctk_lsh[0]; ++i) {
	const int ind = CCTK_GFINDEX3D(cctkGH, i,j,k);
	gfx[ind] = igfx[ind] = 10000 * (cctk_lbnd[0] + i) + 100 * (cctk_lbnd[1] + j) + (cctk_lbnd[2] + k);
	gfy[ind] = igfy[ind] = - (i + 3 * j + 5 * k);
	gfz[ind] = igfz[ind] = 10000 * (cctk_gsh[0] - 1 - cctk_lbnd[0] - i) + 100 * (cctk_lbnd[1] + j) + (cctk_lbnd[2] + k);
      }
    }
  }
  
  assert (cctk_dim <= 3);
  for (d=0; d<cctk_dim; ++d) {
    info[d].src.off = 0;
    info[d].src.len = cctk_gsh[d];
    info[d].src.str = 1;
    info[d].dst.off = 0;
    info[d].dst.len = cctk_gsh[d];
    info[d].dst.str = 1;
    info[d].xpose = d;
    info[d].flip = d==0;
  }
  
  runtest (CCTK_PASS_CTOC, info);
  
  
  
  CCTK_INFO ("   Invert all directions");
  
  for (k=0; k<cctk_lsh[2]; ++k) {
    for (j=0; j<cctk_lsh[1]; ++j) {
      for (i=0; i<cctk_lsh[0]; ++i) {
	const int ind = CCTK_GFINDEX3D(cctkGH, i,j,k);
	gfx[ind] = igfx[ind] = 10000 * (cctk_lbnd[0] + i) + 100 * (cctk_lbnd[1] + j) + (cctk_lbnd[2] + k);
	gfy[ind] = igfy[ind] = - (i + 3 * j + 5 * k);
	gfz[ind] = igfz[ind] = 10000 * (cctk_gsh[0] - 1 - cctk_lbnd[0] - i) + 100 * (cctk_gsh[1] - 1 - cctk_lbnd[1] - j) + (cctk_gsh[2] - 1 - cctk_lbnd[2] - k);
      }
    }
  }
  
  assert (cctk_dim <= 3);
  for (d=0; d<cctk_dim; ++d) {
    info[d].src.off = 0;
    info[d].src.len = cctk_gsh[d];
    info[d].src.str = 1;
    info[d].dst.off = 0;
    info[d].dst.len = cctk_gsh[d];
    info[d].dst.str = 1;
    info[d].xpose = d;
    info[d].flip = 1;
  }
  
  runtest (CCTK_PASS_CTOC, info);
  
  
  
  CCTK_INFO ("   Transfer smaller slab");
  
  for (k=0; k<cctk_lsh[2]; ++k) {
    for (j=0; j<cctk_lsh[1]; ++j) {
      for (i=0; i<cctk_lsh[0]; ++i) {
	const int ind = CCTK_GFINDEX3D(cctkGH, i,j,k);
	gfx[ind] = igfx[ind] = 10000 * (cctk_lbnd[0] + i) + 100 * (cctk_lbnd[1] + j) + (cctk_lbnd[2] + k);
	gfy[ind] = igfy[ind] = - (i + 3 * j + 5 * k);
	if (cctk_lbnd[0] + i >= 2 && cctk_lbnd[0] + i < 5
	    && cctk_lbnd[1] + j >= 3 && cctk_lbnd[1] + j < 5
	    && cctk_lbnd[2] + k >= 4 && cctk_lbnd[2] + k < 6) {
	  gfz[ind] = igfz[ind] = 10000 * (cctk_lbnd[0] + i + 2) + 100 * (cctk_lbnd[1] + j + 1) + (cctk_lbnd[2] + k - 1);
	} else {
	  gfz[ind] = gfy[ind];
	  igfz[ind] = igfy[ind];
	}
      }
    }
  }
  
  assert (cctk_dim <= 3);
  for (d=0; d<cctk_dim; ++d) {
    info[d].src.str = 1;
    info[d].dst.str = 1;
    info[d].xpose = d;
    info[d].flip = 0;
  }
  info[0].src.off = 4;
  info[1].src.off = 4;
  info[2].src.off = 3;
  info[0].src.len = 3;
  info[1].src.len = 2;
  info[2].src.len = 2;
  info[0].dst.off = 2;
  info[1].dst.off = 3;
  info[2].dst.off = 4;
  info[0].dst.len = 3;
  info[1].dst.len = 2;
  info[2].dst.len = 2;
  
  runtest (CCTK_PASS_CTOC, info);
  
  
  
  CCTK_INFO ("   Transpose x and z");
  
  for (k=0; k<cctk_lsh[2]; ++k) {
    for (j=0; j<cctk_lsh[1]; ++j) {
      for (i=0; i<cctk_lsh[0]; ++i) {
	const int ind = CCTK_GFINDEX3D(cctkGH, i,j,k);
	gfx[ind] = igfx[ind] = 10000 * (cctk_lbnd[0] + i) + 100 * (cctk_lbnd[1] + j) + (cctk_lbnd[2] + k);
	gfy[ind] = igfy[ind] = - (i + 3 * j + 5 * k);
	if (cctk_lbnd[0] + i >= 2 && cctk_lbnd[0] + i < 4
	    && cctk_lbnd[1] + j >= 3 && cctk_lbnd[1] + j < 6
	    && cctk_lbnd[2] + k >= 4 && cctk_lbnd[2] + k < 6) {
	  gfz[ind] = igfz[ind] = 10000 * (cctk_lbnd[2] + k + 1) + 100 * (cctk_lbnd[1] + j + 2) + (cctk_lbnd[0] + i - 1);
	} else {
	  gfz[ind] = gfy[ind];
	  igfz[ind] = igfy[ind];
	}
      }
    }
  }
  
  assert (cctk_dim <= 3);
  for (d=0; d<cctk_dim; ++d) {
    info[d].src.str = 1;
    info[d].dst.str = 1;
    info[d].flip = 0;
  }
  info[0].src.off = 5;
  info[1].src.off = 5;
  info[2].src.off = 1;
  info[0].src.len = 2;
  info[1].src.len = 3;
  info[2].src.len = 2;
  info[0].dst.off = 2;
  info[1].dst.off = 3;
  info[2].dst.off = 4;
  info[0].dst.len = 2;
  info[1].dst.len = 3;
  info[2].dst.len = 2;
  info[0].xpose = 2;
  info[1].xpose = 1;
  info[2].xpose = 0;
  
  runtest (CCTK_PASS_CTOC, info);
  
  
  
  CCTK_INFO ("   Transpose");
  
  for (k=0; k<cctk_lsh[2]; ++k) {
    for (j=0; j<cctk_lsh[1]; ++j) {
      for (i=0; i<cctk_lsh[0]; ++i) {
	const int ind = CCTK_GFINDEX3D(cctkGH, i,j,k);
	gfx[ind] = igfx[ind] = 10000 * (cctk_lbnd[0] + i) + 100 * (cctk_lbnd[1] + j) + (cctk_lbnd[2] + k);
	gfy[ind] = igfy[ind] = - (i + 3 * j + 5 * k);
	if (cctk_lbnd[0] + i >= 2 && cctk_lbnd[0] + i < 5
	    && cctk_lbnd[1] + j >= 3 && cctk_lbnd[1] + j < 5
	    && cctk_lbnd[2] + k >= 4 && cctk_lbnd[2] + k < 6) {
	  gfz[ind] = igfz[ind] = 10000 * (cctk_lbnd[2] + k + 1) + 100 * (cctk_lbnd[0] + i + 3) + (cctk_lbnd[1] + j - 2);
	} else {
	  gfz[ind] = gfy[ind];
	  igfz[ind] = igfy[ind];
	}
      }
    }
  }
  
  assert (cctk_dim <= 3);
  for (d=0; d<cctk_dim; ++d) {
    info[d].src.str = 1;
    info[d].dst.str = 1;
    info[d].flip = 0;
  }
  info[0].src.off = 5;
  info[1].src.off = 5;
  info[2].src.off = 1;
  info[0].src.len = 2;
  info[1].src.len = 3;
  info[2].src.len = 2;
  info[0].dst.off = 2;
  info[1].dst.off = 3;
  info[2].dst.off = 4;
  info[0].dst.len = 3;
  info[1].dst.len = 2;
  info[2].dst.len = 2;
  info[0].xpose = 1;
  info[1].xpose = 2;
  info[2].xpose = 0;
  
  runtest (CCTK_PASS_CTOC, info);
  
  
  
  CCTK_INFO ("   Transpose x and z and invert x");
  
  for (k=0; k<cctk_lsh[2]; ++k) {
    for (j=0; j<cctk_lsh[1]; ++j) {
      for (i=0; i<cctk_lsh[0]; ++i) {
	const int ind = CCTK_GFINDEX3D(cctkGH, i,j,k);
	gfx[ind] = igfx[ind] = 10000 * (cctk_lbnd[0] + i) + 100 * (cctk_lbnd[1] + j) + (cctk_lbnd[2] + k);
	gfy[ind] = igfy[ind] = - (i + 3 * j + 5 * k);
	if (cctk_lbnd[0] + i >= 2 && cctk_lbnd[0] + i < 4
	    && cctk_lbnd[1] + j >= 3 && cctk_lbnd[1] + j < 6
	    && cctk_lbnd[2] + k >= 4 && cctk_lbnd[2] + k < 6) {
	  gfz[ind] = igfz[ind] = 10000 * (cctk_lbnd[2] + k + 1) + 100 * (cctk_lbnd[1] + j + 2) + (4 - cctk_lbnd[0] - i);
	} else {
	  gfz[ind] = gfy[ind];
	  igfz[ind] = igfy[ind];
	}
      }
    }
  }
  
  assert (cctk_dim <= 3);
  for (d=0; d<cctk_dim; ++d) {
    info[d].src.str = 1;
    info[d].dst.str = 1;
    info[d].flip = d==0;
  }
  info[0].src.off = 5;
  info[1].src.off = 5;
  info[2].src.off = 1;
  info[0].src.len = 2;
  info[1].src.len = 3;
  info[2].src.len = 2;
  info[0].dst.off = 2;
  info[1].dst.off = 3;
  info[2].dst.off = 4;
  info[0].dst.len = 2;
  info[1].dst.len = 3;
  info[2].dst.len = 2;
  info[0].xpose = 2;
  info[1].xpose = 1;
  info[2].xpose = 0;
  
  runtest (CCTK_PASS_CTOC, info);
  
  
  
  CCTK_INFO ("   Skip every second cell along y direction");
  
  for (k=0; k<cctk_lsh[2]; ++k) {
    for (j=0; j<cctk_lsh[1]; ++j) {
      for (i=0; i<cctk_lsh[0]; ++i) {
	const int ind = CCTK_GFINDEX3D(cctkGH, i,j,k);
	gfx[ind] = igfx[ind] = 10000 * (cctk_lbnd[0] + i) + 100 * (cctk_lbnd[1] + j) + (cctk_lbnd[2] + k);
	gfy[ind] = igfy[ind] = - (i + 3 * j + 5 * k);
	if (cctk_lbnd[1] + j < cctk_gsh[1] / 2) {
	  gfz[ind] = igfz[ind] = 10000 * (cctk_lbnd[0] + i) + 2 * 100 * (cctk_lbnd[1] + j) + (cctk_lbnd[2] + k);
	} else {
	  gfz[ind] = gfy[ind];
	  igfz[ind] = igfy[ind];
	}
      }
    }
  }
  
  assert (cctk_dim <= 3);
  for (d=0; d<cctk_dim; ++d) {
    info[d].src.off = 0;
    info[d].src.len = d == 1 ? cctk_gsh[d] / 2 : cctk_gsh[d];
    info[d].src.str = d == 1 ? 2 : 1;
    info[d].dst.off = 0;
    info[d].dst.len = d == 1 ? cctk_gsh[d] / 2 : cctk_gsh[d];
    info[d].dst.str = 1;
    info[d].xpose = d;
    info[d].flip = 0;
  }
  
  runtest (CCTK_PASS_CTOC, info);
  
  
  
  /* TODO: more strides, into array, out of array */
  
  
  
  * success = 1;
}
