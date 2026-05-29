#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <util_ErrorCodes.h>
#include <util_Table.h>

#include <TATelliptic.h>

#define DIM 3



int TATelliptic_monitor (const cGH * const cctkGH,
			 const int * const var,
			 const int * const res,
			 const int nvars,
			 const int options_table,
			 const calcfunc calcres,
			 const calcfunc applybnds,
			 void * const userdata)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  /* Domain descriptors */
  cGroup groupdata;
  cGroupDynamicData groupdyndata;
  int dim, gsh[DIM], lsh[DIM], nghostzones[DIM], bbox[2*DIM];
  CCTK_INT nboundaryzones[2*DIM];
  int imin[DIM], imax[DIM];
  
  /* Variables */
  CCTK_REAL const * restrict resptr;
  
  /* Statistics */
  CCTK_REAL sum, sum2;
  CCTK_REAL norm1, norm2, norm_inf, norm_count;
  CCTK_REAL tmp;
  
  /* Reduction handles */
  int sum_handle, norm_inf_handle;
  CCTK_REAL locals[2], globals[2];
  
  /* Indices */
  int n, nn;
  int i, j, k;
  int ind;
  int d;
  
  /* Error control */
  int nelems;
  int ierr;
  
  if (! CCTK_IsThornActive(CCTK_THORNSTRING)) {
    CCTK_ERROR ("Thorn " CCTK_THORNSTRING " has not been activated.  It is therefore not possible to call TATelliptic_monitor.");
  }
  
  /* Check arguments */
  assert (cctkGH);
  
  assert (nvars > 0);
  assert (var);
  assert (res);
  for (n=0; n<nvars; ++n) {
    assert (var[n] >= 0 && var[n] < CCTK_NumVars());
    assert (CCTK_VarTypeI(var[n]) == CCTK_VARIABLE_REAL);
    assert (CCTK_QueryGroupStorageI(cctkGH, CCTK_GroupIndexFromVarI(var[n])));
  }
  for (n=0; n<nvars; ++n) {
    assert (res[n] >= 0 && res[n] < CCTK_NumVars());
    assert (CCTK_VarTypeI(res[n]) == CCTK_VARIABLE_REAL);
    assert (CCTK_QueryGroupStorageI(cctkGH, CCTK_GroupIndexFromVarI(res[n])));
  }
  for (n=0; n<nvars; ++n) {
    assert (var[n] != res[n]);
    for (nn=0; nn<n; ++nn) {
      assert (var[nn] != var[n]);
      assert (var[nn] != res[n]);
      assert (res[nn] != var[n]);
      assert (res[nn] != res[n]);
    }
  }
  
  /* Get domain description */
  ierr = CCTK_GroupData (CCTK_GroupIndexFromVarI(var[0]), &groupdata);
  assert (!ierr);
  ierr = CCTK_GroupDynamicData (cctkGH, CCTK_GroupIndexFromVarI(var[0]),
				&groupdyndata);
  assert (!ierr);
  dim = groupdata.dim;
  assert (dim>=0 && dim<=DIM);
  assert (dim == groupdyndata.dim);
  for (d=0; d<dim; ++d) {
    gsh[d] = groupdyndata.gsh[d];
    lsh[d] = groupdyndata.lsh[d];
    nghostzones[d] = groupdyndata.nghostzones[d];
    bbox[2*d  ] = groupdyndata.bbox[2*d  ];
    bbox[2*d+1] = groupdyndata.bbox[2*d+1];
  }
  
  assert (options_table >= 0);
  
  nelems = Util_TableGetIntArray
    (options_table, 2*DIM, nboundaryzones, "nboundaryzones");
  if (nelems == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    for (d=0; d<dim; ++d) {
      nboundaryzones[2*d  ] = nghostzones[d];
      nboundaryzones[2*d+1] = nghostzones[d];
    }
  } else if (nelems != 2*dim) {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Options table key \"nboundaryzones\" is not an integer array of length %d", 2*dim);
    return -1;
  }
  
  for (d=dim; d<DIM; ++d) {
    gsh[d] = 1;
    lsh[d] = 1;
    nghostzones[d] = 0;
    bbox[2*d  ] = 0;
    bbox[2*d+1] = 0;
    nboundaryzones[2*d  ] = 0;
    nboundaryzones[2*d+1] = 0;
  }
  for (d=0; d<DIM; ++d) {
    assert (gsh[d]>=0);
    assert (lsh[d]>=0);
    assert (nghostzones[d] >= 0 && 2*nghostzones[d] <= lsh[d]);
    assert (nboundaryzones[2*d] >= 0 && nboundaryzones[2*d+1] >= 0);
    assert (nboundaryzones[2*d] + nboundaryzones[2*d+1] <= lsh[d]);
    assert (lsh[d] <= gsh[d]);
  }
  
  /* Check all variables */
  for (n=0; n<nvars; ++n) {
    ierr = CCTK_GroupData (CCTK_GroupIndexFromVarI(var[n]), &groupdata);
    assert (!ierr);
    ierr = CCTK_GroupDynamicData
      (cctkGH, CCTK_GroupIndexFromVarI(var[n]), &groupdyndata);
    assert (!ierr);
    assert (groupdata.dim == dim);
    assert (groupdyndata.dim == dim);
    for (d=0; d<dim; ++d) {
      assert (groupdyndata.lsh[d] == lsh[d]);
      assert (groupdyndata.nghostzones[d] == nghostzones[d]);
    }
  }
  
  assert (calcres);
  assert (applybnds);
  
  for (d=0; d<DIM; ++d) {
    imin[d] = bbox[2*d] ? nboundaryzones[2*d] : nghostzones[d];
    imax[d] = cctk_lsh[d] - (bbox[2*d+1] ? nboundaryzones[2*d+1] : nghostzones[d]);
  }
  
  /* Calculate residual */
  ierr = calcres (cctkGH, -1, userdata);
  if (ierr != 0) {
    CCTK_WARN (CCTK_WARN_COMPLAIN,
               "Residual calculation reported error; aborting monitor.");
    return ierr;
  }
  
  /* Calculate norms */
  sum = 0;
  sum2 = 0;
  norm_inf = 0;
  for (n=0; n<nvars; ++n) {
    resptr = CCTK_VarDataPtrI (cctkGH, 0, res[n]);
    assert (resptr);
    for (k=imin[2]; k<imax[2]; ++k) {
      for (j=imin[1]; j<imax[1]; ++j) {
        for (i=imin[0]; i<imax[0]; ++i) {
          ind = i + lsh[0] * (j + lsh[1] * k);
          tmp = resptr[ind];
          sum += tmp;
          sum2 += tmp*tmp;
          norm_inf = fabs(tmp) > norm_inf ? fabs(tmp) : norm_inf;
        }
      }
    }
  } /* for n */
  
  norm_count = nvars;
  for (d=0; d<dim; ++d) {
    norm_count *= gsh[d] - (nboundaryzones[2*d] + nboundaryzones[2*d+1]);
  }
  
  sum_handle = CCTK_ReductionArrayHandle ("sum");
  assert (sum_handle >= 0);
  locals[0] = sum;
  locals[1] = sum2;
  ierr = CCTK_ReduceLocArrayToArray1D
    (cctkGH, -1, sum_handle, &locals, &globals, 2, CCTK_VARIABLE_REAL);
  assert (!ierr);
  sum = globals[0];
  sum2 = globals[1];
  
  norm_inf_handle = CCTK_ReductionHandle ("norm_inf");
  assert (norm_inf_handle >= 0);
  locals[0] = norm_inf;
  ierr = CCTK_ReduceLocArrayToArray1D
    (cctkGH, -1, norm_inf_handle, &locals, &globals, 1, CCTK_VARIABLE_REAL);
  assert (!ierr);
  norm_inf = globals[0];
  
  norm1 = sum / norm_count;
  norm2 = sqrt(sum2 / norm_count);
  
  /* Log output */
  CCTK_VInfo (CCTK_THORNSTRING,
	      "Monitor: residual norms: count %g, L1 %g, L2 %g, Linf %g",
	      (double)norm_count,
	      (double)norm1, (double)norm2, (double)norm_inf);
  
  return 0;
}



int TATelliptic_register_monitor (void)
{
  int ierr;
  ierr = TATelliptic_RegisterSolver (TATelliptic_monitor, "TATmonitor");
  assert (!ierr);
  return 0;
}
