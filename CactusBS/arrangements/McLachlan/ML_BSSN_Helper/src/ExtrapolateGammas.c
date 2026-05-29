#include <cctk.h>
#include <cctk_Arguments.h>

#include <stdbool.h>

#ifndef DECLARE_CCTK_ARGUMENTS_CHECKED
#define DECLARE_CCTK_ARGUMENTS_CHECKED(func) DECLARE_CCTK_ARGUMENTS
#endif

static void extrap(const cGH *cctkGH, CCTK_REAL *var);

void ML_BSSN_ExtrapolateGammas(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_BSSN_ExtrapolateGammas);

  extrap(cctkGH, Xt1);
  extrap(cctkGH, Xt2);
  extrap(cctkGH, Xt3);

  static bool have_A_index = false;
  static int A_index;
  if (!have_A_index) {
    A_index = CCTK_VarIndex("ML_BSSN::A");
    have_A_index = true;
  }
  if (A_index >= 0) {
    if(CCTK_IsFunctionAliased("Driver_RequireValidData")) {
      CCTK_INT variables[1] = {A_index};
      CCTK_INT tls[1] = {0};
      CCTK_INT where[1] = {CCTK_VALID_INTERIOR};
      CCTK_INT ierr = Driver_RequireValidData(cctkGH, variables, tls, 1, where);
      assert(!ierr);
    }
    CCTK_REAL *A_ptr = CCTK_VarDataPtrI(cctkGH, 0, A_index);
    extrap(cctkGH, A_ptr);
    if(CCTK_IsFunctionAliased("Driver_NotifyDataModified")) {
      CCTK_INT variables[1] = {A_index};
      CCTK_INT tls[1] = {0};
      CCTK_INT where[1] = {CCTK_VALID_BOUNDARY};
      CCTK_INT ierr = Driver_NotifyDataModified(cctkGH, variables, tls, 1,
                                                where);
      assert(!ierr);
    }
  }

  static bool have_B_index = false;
  static int B_index;
  if (!have_B_index) {
    B_index = CCTK_VarIndex("ML_BSSN::B1");
    have_B_index = true;
  }
  if (B_index >= 0) {
    if(CCTK_IsFunctionAliased("Driver_RequireValidData")) {
      CCTK_INT variables[3] = {B_index, B_index + 1, B_index + 2};
      CCTK_INT tls[3] = {0, 0 ,0};
      CCTK_INT where[3] = {CCTK_VALID_INTERIOR, CCTK_VALID_INTERIOR,
                           CCTK_VALID_INTERIOR};
      CCTK_INT ierr = Driver_RequireValidData(cctkGH, variables, tls, 3, where);
      assert(!ierr);
    }
    CCTK_REAL *B1_ptr = CCTK_VarDataPtrI(cctkGH, 0, B_index);
    CCTK_REAL *B2_ptr = CCTK_VarDataPtrI(cctkGH, 0, B_index + 1);
    CCTK_REAL *B3_ptr = CCTK_VarDataPtrI(cctkGH, 0, B_index + 2);
    extrap(cctkGH, B1_ptr);
    extrap(cctkGH, B2_ptr);
    extrap(cctkGH, B3_ptr);
    if(CCTK_IsFunctionAliased("Driver_NotifyDataModified")) {
      CCTK_INT variables[3] = {B_index, B_index + 1, B_index + 2};
      CCTK_INT tls[3] = {0, 0 ,0};
      CCTK_INT where[3] = {CCTK_VALID_INTERIOR, CCTK_VALID_INTERIOR,
                           CCTK_VALID_INTERIOR};
      CCTK_INT ierr = Driver_NotifyDataModified(cctkGH, variables, tls, 3,
                                                where);
      assert(!ierr);
    }
  }
}

static void extrap(const cGH *cctkGH, CCTK_REAL *var) {
  ExtrapolateGammas(cctkGH, var);
}
