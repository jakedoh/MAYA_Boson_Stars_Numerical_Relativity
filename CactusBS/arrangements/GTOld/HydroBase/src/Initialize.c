/* Provide simple initialization, taken from public version of HydroBase */

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

void HydroBase_Bvec_Zero (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

    int const np = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];

  #pragma omp parallel for
  for (int i=0; i<np; ++i) {
    Bvec[i]      = 0.0;
    Bvec[i+  np] = 0.0;
    Bvec[i+2*np] = 0.0;
  }

  if (CCTK_EQUALS (initial_data_setup_method, "init_some_levels") ||
      CCTK_EQUALS (initial_data_setup_method, "init_single_levels"))
  {
    /* do nothing */
  }
  else if (CCTK_EQUALS (initial_data_setup_method, "init_all_levels"))
  {

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::Bvec") >= 2) {
      #pragma omp parallel for
      for (int i=0; i<np; ++i) {
        Bvec_p[i]      = 0.0;
        Bvec_p[i+  np] = 0.0;
        Bvec_p[i+2*np] = 0.0;
      }
    }

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::Bvec") >= 3) {
      #pragma omp parallel for
      for (int i=0; i<np; ++i) {
        Bvec_p_p[i]      = 0.0;
        Bvec_p_p[i+  np] = 0.0;
        Bvec_p_p[i+2*np] = 0.0;
      }
    }

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::Bvec") >= 4) {
      CCTK_WARN (CCTK_WARN_ABORT,
                 "Too many active time levels for HydroBase::Bvec");
    }

  }
  else
  {
    CCTK_WARN (CCTK_WARN_ABORT,
               "Unsupported parameter value for InitBase::initial_data_setup_method");
  }

}

void HydroBase_Evec_Zero (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

    int const np = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];

  #pragma omp parallel for
  for (int i=0; i<np; ++i) {
    Evec[i]      = 0.0;
    Evec[i+  np] = 0.0;
    Evec[i+2*np] = 0.0;
  }

  if (CCTK_EQUALS (initial_data_setup_method, "init_some_levels") ||
      CCTK_EQUALS (initial_data_setup_method, "init_single_levels"))
  {
    /* do nothing */
  }
  else if (CCTK_EQUALS (initial_data_setup_method, "init_all_levels"))
  {

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::Evec") >= 2) {
      #pragma omp parallel for
      for (int i=0; i<np; ++i) {
        Evec_p[i]      = 0.0;
        Evec_p[i+  np] = 0.0;
        Evec_p[i+2*np] = 0.0;
      }
    }

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::Evec") >= 3) {
      #pragma omp parallel for
      for (int i=0; i<np; ++i) {
        Evec_p_p[i]      = 0.0;
        Evec_p_p[i+  np] = 0.0;
        Evec_p_p[i+2*np] = 0.0;
      }
    }

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::Evec") >= 4) {
      CCTK_WARN (CCTK_WARN_ABORT,
                 "Too many active time levels for HydroBase::Evec");
    }

  }
  else
  {
    CCTK_WARN (CCTK_WARN_ABORT,
               "Unsupported parameter value for InitBase::initial_data_setup_method");
  }

}

void HydroBase_Avec_Zero (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

    int const np = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];

  #pragma omp parallel for
  for (int i=0; i<np; ++i) {
    Avec[i]      = 0.0;
    Avec[i+  np] = 0.0;
    Avec[i+2*np] = 0.0;
  }

  if (CCTK_EQUALS (initial_data_setup_method, "init_some_levels") ||
      CCTK_EQUALS (initial_data_setup_method, "init_single_levels"))
  {
    /* do nothing */
  }
  else if (CCTK_EQUALS (initial_data_setup_method, "init_all_levels"))
  {

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::Avec") >= 2) {
      #pragma omp parallel for
      for (int i=0; i<np; ++i) {
        Avec_p[i]      = 0.0;
        Avec_p[i+  np] = 0.0;
        Avec_p[i+2*np] = 0.0;
      }
    }

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::Avec") >= 3) {
      #pragma omp parallel for
      for (int i=0; i<np; ++i) {
        Avec_p_p[i]      = 0.0;
        Avec_p_p[i+  np] = 0.0;
        Avec_p_p[i+2*np] = 0.0;
      }
    }

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::Avec") >= 4) {
      CCTK_WARN (CCTK_WARN_ABORT,
                 "Too many active time levels for HydroBase::Avec");
    }

  }
  else
  {
    CCTK_WARN (CCTK_WARN_ABORT,
               "Unsupported parameter value for InitBase::initial_data_setup_method");
  }

}

void HydroBase_Aphi_Zero (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

    int const np = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];

  #pragma omp parallel for
  for (int i=0; i<np; ++i) {
    Aphi[i]      = 0.0;
  }

  if (CCTK_EQUALS (initial_data_setup_method, "init_some_levels") ||
      CCTK_EQUALS (initial_data_setup_method, "init_single_levels"))
  {
    /* do nothing */
  }
  else if (CCTK_EQUALS (initial_data_setup_method, "init_all_levels"))
  {

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::Aphi") >= 2) {
      #pragma omp parallel for
      for (int i=0; i<np; ++i) {
        Aphi_p[i]      = 0.0;
      }
    }

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::Aphi") >= 3) {
      #pragma omp parallel for
      for (int i=0; i<np; ++i) {
        Aphi_p_p[i]      = 0.0;
      }
    }

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::Aphi") >= 4) {
      CCTK_WARN (CCTK_WARN_ABORT,
                 "Too many active time levels for HydroBase::Aphi");
    }

  }
  else
  {
    CCTK_WARN (CCTK_WARN_ABORT,
               "Unsupported parameter value for InitBase::initial_data_setup_method");
  }

}

