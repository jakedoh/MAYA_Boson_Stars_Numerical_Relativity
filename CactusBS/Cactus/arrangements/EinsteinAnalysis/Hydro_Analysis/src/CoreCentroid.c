#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Loop.h"

// math help macros
#define SQR(x) ((x)*(x))

// declare an array for the timelevels of a given variable
// we substitute the current level if a previous level does not exist
#define DECLARE_MY_VARIABLE(var) CCTK_REAL *my##var[3] = {var, var##_p ? var##_p : var, var##_p_p ? var##_p_p : var}

// Storage state of the temp. grid when we entered the loop
#define HYDRO_ANALYSIS_CORE_RHO_CENTROID_GF CCTK_THORNSTRING"::Hydro_Analysis_core_rho_centroid_gf"
static int StorageOnEntry = -1;

// get storage for the temp grid function
void Hydro_Analysis_CompCoreRhoCentroid_GetStorage(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  const int group = CCTK_GroupIndex(HYDRO_ANALYSIS_CORE_RHO_CENTROID_GF);
  assert(group >= 0);

  assert(StorageOnEntry == -1 && "We only have one level of memory, don't call us recursively");
  const int tls = timelevels;
  CCTK_GroupStorageIncrease(cctkGH, 1, &group, &tls, &StorageOnEntry);
}

// compute x^i*rho in a region r_core around the densest point
void Hydro_Analysis_CompCoreRhoCentroid_PrepareReduction(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  // this is a hack since the timelevels parameter really does not have to be
  // the number of levels that are active. It likely is however.
  DECLARE_MY_VARIABLE(rho);
  DECLARE_MY_VARIABLE(Hydro_Analysis_core_rho_centroid_gf);

  // core is centered around densest point
  // this is why we need to be scheduld after Hydro_Analysis_LocationSearch
  const CCTK_REAL corex = Hydro_Analysis_rho_max_loc[0];
  const CCTK_REAL corey = Hydro_Analysis_rho_max_loc[1];
  const CCTK_REAL corez = Hydro_Analysis_rho_max_loc[2];
  const CCTK_REAL corerhomin = Hydro_Analysis_core_rho_rel_min * (*Hydro_Analysis_rho_max);

  // need to loop over all timelevels since we do not keep the grid functions around
  // even if we did, we would have to be called in EVOL to get the step-ahead for Berger-Oliger right
  const int num_tl = CCTK_ActiveTimeLevelsGN(cctkGH, HYDRO_ANALYSIS_CORE_RHO_CENTROID_GF );
  assert(num_tl <= 3);
  for (int tl = 0 ; tl < num_tl ; tl++) {
#pragma omp parallel
    CCTK_LOOP3_ALL(CompCoreRhoCentroid_PrepareReduction, cctkGH, i,j,k) {
      const int ind = CCTK_GFINDEX3D (cctkGH, i, j, k);
      const int indr = CCTK_VECTGFINDEX3D (cctkGH, i, j, k, 0);
      const int indx = CCTK_VECTGFINDEX3D (cctkGH, i, j, k, 1);
      const int indy = CCTK_VECTGFINDEX3D (cctkGH, i, j, k, 2);
      const int indz = CCTK_VECTGFINDEX3D (cctkGH, i, j, k, 3);

      const double r2 = SQR(x[ind]-corex) + SQR(y[ind]-corey) + SQR(z[ind]-corez);
      if (r2 < SQR(Hydro_Analysis_r_core) && myrho[tl][ind] >= corerhomin) {
        myHydro_Analysis_core_rho_centroid_gf[tl][indx] = myrho[tl][ind] * x[ind];
        myHydro_Analysis_core_rho_centroid_gf[tl][indy] = myrho[tl][ind] * y[ind];
        myHydro_Analysis_core_rho_centroid_gf[tl][indz] = myrho[tl][ind] * z[ind];
        myHydro_Analysis_core_rho_centroid_gf[tl][indr] = myrho[tl][ind];
      } else {
        myHydro_Analysis_core_rho_centroid_gf[tl][indx] = 0.;
        myHydro_Analysis_core_rho_centroid_gf[tl][indy] = 0.;
        myHydro_Analysis_core_rho_centroid_gf[tl][indz] = 0.;
        myHydro_Analysis_core_rho_centroid_gf[tl][indr] = 0.;
      }
      
    } CCTK_ENDLOOP3_ALL(CompCoreRhoCentroid_PrepareReduction);
  }
}

// compute \sum rho*x^i over all processors and compute x_c^i = x^i / \sum rho
void Hydro_Analysis_CompCoreRhoCentroid_Reduction (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  const CCTK_INT handle_sum = CCTK_ReductionHandle("sum");
  if (handle_sum < 0) {
    CCTK_WARN(0, "Unable to get reduction handle 'sum'.");
    return; /* NOTREACHED */
  }

  // total "mass"
  {
    *Hydro_Analysis_core_rho_sum = -1.0;
    const int ierr = CCTK_Reduce(cctkGH, -1, handle_sum, 1, CCTK_VARIABLE_REAL,
                                 Hydro_Analysis_core_rho_sum, 1, CCTK_VarIndex(HYDRO_ANALYSIS_CORE_RHO_CENTROID_GF "[0]"));
    if (ierr) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error while reducing " HYDRO_ANALYSIS_CORE_RHO_CENTROID_GF "[0]: %d", ierr);
      return; /* NOTREACHED */
    }
  }

  // sometimes we might run out of star (though this should not happen in well
  // behaved simulations)
  static int warn_about_no_mass = 1;
  if (*Hydro_Analysis_core_rho_sum == 0.) {
    if(CCTK_MyProc(cctkGH) == 0 && warn_about_no_mass) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "No mass in core region. Not changing CoM.");
      // only trigger warning once per run (or until we find and loose mass again)
      warn_about_no_mass = 0;
    }
    return;
  } else {
    // start warning again
    warn_about_no_mass = 1;
  }

  // center of mass
  assert(*Hydro_Analysis_core_rho_sum > 0.);
  for (int i = 0 ; i < 3 ; i++) {
    char buf[] = HYDRO_ANALYSIS_CORE_RHO_CENTROID_GF"[0]";
    const char *fmt = HYDRO_ANALYSIS_CORE_RHO_CENTROID_GF"[%1d]";
    size_t len_written;

    len_written = snprintf(buf, sizeof(buf), fmt, i+1); // rho is [0], rho*x is [1], ...
    if (len_written >= sizeof(buf)) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error constructing reduction variable name. Buffer size %d too short for string of length %d.",
                 (int)sizeof(buf), (int)len_written);
      return; /* NOTREACHED */
    }

    // rho*x^i
    const int ierr = CCTK_Reduce(cctkGH, -1, handle_sum, 1, CCTK_VARIABLE_REAL,
                                 &(Hydro_Analysis_core_rho_centroid[i]), 1,
                                 CCTK_VarIndex(buf));
    if (ierr) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Error while reducing " HYDRO_ANALYSIS_CORE_RHO_CENTROID_GF "[%d]: %d",
               i+1, ierr);
      return; /* NOTREACHED */
    }

    Hydro_Analysis_core_rho_centroid[i] /= *Hydro_Analysis_core_rho_sum;
  }
}

// free storage allocated in GetStorage
void Hydro_Analysis_CompCoreRhoCentroid_FreeStorage(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  const int group = CCTK_GroupIndex(HYDRO_ANALYSIS_CORE_RHO_CENTROID_GF);
  assert(group >= 0);

  assert(StorageOnEntry != -1 && "No storage state was saved.");
  CCTK_GroupStorageDecrease(cctkGH, 1, &group, &StorageOnEntry, NULL);
  StorageOnEntry = -1;
}
