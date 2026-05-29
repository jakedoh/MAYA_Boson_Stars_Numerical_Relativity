#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

void ML_CCZ4_ParamCompat(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(my_initial_boundary_condition, "extrapolate-gammas")) {
    // Override initial_boundary_condition
    CCTK_WARN(CCTK_WARN_ALERT, "Forcing "
                               "ML_CCZ4::initial_boundary_condition="
                               "\"extrapolate-gammas\" because "
                               "ML_CCZ4::my_initial_boundary_condition="
                               "\"extrapolate-gammas\"");
    int ierr = CCTK_ParameterSet("initial_boundary_condition", "ML_CCZ4",
                                 "extrapolate-gammas");
    if (ierr) {
      CCTK_ERROR("Could not set "
                 "ML_CCZ4::initial_boundary_condition=\"extrapolate-gammas\"");
    }
  }

  if (CCTK_EQUALS(my_rhs_boundary_condition, "NewRad")) {
    // Override rhs_boundary_condition
    CCTK_WARN(CCTK_WARN_ALERT,
              "Forcing ML_CCZ4::rhs_boundary_condition=\"NewRad\" because "
              "ML_CCZ4::my_rhs_boundary_condition=\"NewRad\"");
    int ierr = CCTK_ParameterSet("rhs_boundary_condition", "ML_CCZ4", "NewRad");
    if (ierr) {
      CCTK_ERROR("Could not set ML_CCZ4::rhs_boundary_condition=\"NewRad\"");
    }
  } else if (CCTK_EQUALS(my_rhs_boundary_condition, "static")) {
    // Override rhs_boundary_condition
    CCTK_WARN(CCTK_WARN_ALERT,
              "Forcing ML_CCZ4::rhs_boundary_condition=\"scalar\" because "
              "ML_CCZ4::my_rhs_boundary_condition=\"static\" (sic!)");
    int ierr = CCTK_ParameterSet("rhs_boundary_condition", "ML_CCZ4", "scalar");
    if (ierr) {
      CCTK_ERROR("Could not set ML_CCZ4::rhs_boundary_condition=\"scalar\"");
    }
  }

  if (CCTK_EQUALS(apply_dissipation, "never")) {
    // Override epsDiss to disable dissipation
    CCTK_WARN(CCTK_WARN_ALERT, "Forcing ML_CCZ4::epsDiss=0.0 because "
                               "ML_CCZ4::apply_dissipation=\"never\"");
    int ierr = CCTK_ParameterSet("epsDiss", "ML_CCZ4", "0.0");
    if (ierr) {
      CCTK_ERROR("Could not set ML_CCZ4::epsDiss=0.0");
    }
  }

  if (LapseACoeff == 0.0) {
    // Override evolveA
    CCTK_WARN(CCTK_WARN_ALERT, "Forcing ML_CCZ4::evolveA=0 because "
                               "ML_CCZ4::LapseACoeff=0.0");
    int ierr = CCTK_ParameterSet("evolveA", "ML_CCZ4", "0");
    if (ierr) {
      CCTK_ERROR("Could not set ML_CCZ4::evolveA=0");
    }
  } else if (LapseACoeff == 1.0) {
    // Override evolveA
    CCTK_WARN(CCTK_WARN_ALERT, "Forcing ML_CCZ4::evolveA=1 because "
                               "ML_CCZ4::LapseACoeff=1.0");
    int ierr = CCTK_ParameterSet("evolveA", "ML_CCZ4", "1");
    if (ierr) {
      CCTK_ERROR("Could not set ML_CCZ4::evolveA=1");
    }
  }

  // Make a modifiable copy
  CCTK_INT my_evolveB = evolveB;
  if (ShiftBCoeff == 0.0) {
    // Override evolveB
    CCTK_WARN(CCTK_WARN_ALERT, "Forcing ML_CCZ4::evolveB=0 because "
                               "ML_CCZ4::ShiftBCoeff=0.0");
    int ierr = CCTK_ParameterSet("evolveB", "ML_CCZ4", "0");
    if (ierr) {
      CCTK_ERROR("Could not set ML_CCZ4::evolveB=0");
    }
    my_evolveB = 0;
  } else if (ShiftBCoeff == 1.0) {
    // Override evolveB
    CCTK_WARN(CCTK_WARN_ALERT, "Forcing ML_CCZ4::evolveB=1 because "
                               "ML_CCZ4::ShiftBCoeff=1.0");
    int ierr = CCTK_ParameterSet("evolveB", "ML_CCZ4", "1");
    if (ierr) {
      CCTK_ERROR("Could not set ML_CCZ4::evolveB=1");
    }
    my_evolveB = 1;
  }

  if (my_evolveB && shiftGammaCoeff == 0.0) {
    // Override evolveB
    CCTK_WARN(
        CCTK_WARN_ALERT,
        "Forcing ML_CCZ4::evolveB=0 because ML_CCZ4::shiftGammaCoeff=0.0");
    int ierr = CCTK_ParameterSet("evolveB", "ML_CCZ4", "0");
    if (ierr) {
      CCTK_ERROR("Could not set ML_CCZ4::evolveB=0");
    }
  }

  if (LapseAdvectionCoeff == 0.0) {
    // Override advectLapse
    CCTK_WARN(CCTK_WARN_ALERT, "Forcing ML_CCZ4::advectLapse=0 because "
                               "ML_CCZ4::LapseAdvectionCoeff=0.0");
    int ierr = CCTK_ParameterSet("advectLapse", "ML_CCZ4", "0");
    if (ierr) {
      CCTK_ERROR("Could not set ML_CCZ4::advectLapse=0");
    }
  } else if (LapseAdvectionCoeff == 1.0) {
    // Override advectLapse
    CCTK_WARN(CCTK_WARN_ALERT, "Forcing ML_CCZ4::advectLapse=1 because "
                               "ML_CCZ4::LapseAdvectionCoeff=1.0");
    int ierr = CCTK_ParameterSet("advectLapse", "ML_CCZ4", "1");
    if (ierr) {
      CCTK_ERROR("Could not set ML_CCZ4::advectLapse=1");
    }
  }

  if (ShiftAdvectionCoeff == 0.0) {
    // Override advectShift
    CCTK_WARN(CCTK_WARN_ALERT, "Forcing ML_CCZ4::advectShift=0 because "
                               "ML_CCZ4::ShiftAdvectionCoeff=0.0");
    int ierr = CCTK_ParameterSet("advectShift", "ML_CCZ4", "0");
    if (ierr) {
      CCTK_ERROR("Could not set ML_CCZ4::advectShift=0");
    }
  } else if (ShiftAdvectionCoeff == 1.0) {
    // Override advectShift
    CCTK_WARN(CCTK_WARN_ALERT, "Forcing ML_CCZ4::advectShift=1 because "
                               "ML_CCZ4::ShiftAdvectionCoeff=1.0");
    int ierr = CCTK_ParameterSet("advectShift", "ML_CCZ4", "1");
    if (ierr) {
      CCTK_ERROR("Could not set ML_CCZ4::advectShift=1");
    }
  }
}

void ML_CCZ4_ParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  if (!CCTK_EQUALS(my_initial_data, "default")) {
    CCTK_WARN(CCTK_WARN_ALERT, "Parameter ML_CCZ4::my_initial_data is "
                               "outdated; please update the parameter file. Do "
                               "not use this parameter, and set up initial "
                               "conditions via ADMBase as usual.");
  }

  if (!CCTK_EQUALS(my_initial_boundary_condition, "default")) {
    CCTK_WARN(CCTK_WARN_ALERT,
              "Parameter ML_CCZ4::my_initial_boundary_condition is outdated; "
              "please update the parameter file. Do not use this parameter, "
              "and set up initial boundary conditions as usual.");
  }

  if (!CCTK_EQUALS(my_rhs_boundary_condition, "default")) {
    CCTK_WARN(CCTK_WARN_ALERT, "Parameter ML_CCZ4::my_rhs_boundary_condition "
                               "is outdated; please update the parameter file. "
                               "Do not use this parameter, and set up RHS "
                               "boundary conditions as usual.");
  }

  if (!CCTK_EQUALS(my_boundary_condition, "default")) {
    CCTK_WARN(CCTK_WARN_ALERT, "Parameter ML_CCZ4::my_boundary_condition is "
                               "outdated; please update the parameter file. Do "
                               "not use this parameter, and set up RHS "
                               "boundary conditions as usual.");
  }

  if (evolveB && shiftGammaCoeff == 0.0) {
    CCTK_WARN(CCTK_WARN_ALERT, "Parameter ML_CCZ4::shiftGammaCoeff should not "
                               "be set to 0.0 any more to disable evolving "
                               "B^i. Instead, set ML_CCZ4::evolveB=0.");
  }

  if (LapseACoeff != -1.0) {
    CCTK_WARN(CCTK_WARN_ALERT, "Parameter ML_CCZ4::LapseACoeff is "
                               "outdated; please update the parameter file. "
                               "Instead of using this parameter, you should "
                               "set ML_CCZ4::evolveA.");
  }

  if (ShiftBCoeff != -1.0) {
    CCTK_WARN(CCTK_WARN_ALERT, "Parameter ML_CCZ4::ShiftBCoeff is "
                               "outdated; please update the parameter file. "
                               "Instead of using this parameter, you should "
                               "set ML_CCZ4::evolveB.");
  }

  if (LapseAdvectionCoeff != -1.0) {
    CCTK_WARN(CCTK_WARN_ALERT, "Parameter ML_CCZ4::LapseAdvectionCoeff is "
                               "outdated; please update the parameter file. "
                               "Instead of using this parameter, you should "
                               "set ML_CCZ4::advectLapse.");
  }
  if (ShiftAdvectionCoeff != -1.0) {
    CCTK_WARN(CCTK_WARN_ALERT, "Parameter ML_CCZ4::ShiftAdvectionCoeff is "
                               "outdated; please update the parameter file. "
                               "Instead of using this parameter, you should "
                               "set ML_CCZ4::advectShift.");
  }
}
