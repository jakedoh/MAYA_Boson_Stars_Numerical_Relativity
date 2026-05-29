/* $Header$ */

#include <assert.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static void callback (int idx, const char * optstring, void * callback_arg);

/** Ensure that all ADMBase initial data that are supposed to be read
    from a file are actually scheduled for the file reader.  */
void IDFileADM_ParamCheck (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  char * variable_is_read;
  int i;
  int nvars;
  
  variable_is_read = malloc (CCTK_NumVars());
  assert (variable_is_read);
  for (i=0; i<CCTK_NumVars(); ++i) {
    variable_is_read[i] = 0;
  }
  
  nvars = CCTK_TraverseString
    (filereader_ID_vars, callback, variable_is_read, CCTK_GROUP_OR_VAR);
  assert (nvars >= 0);
  
  if (CCTK_EQUALS(initial_lapse, "read from file")) {
    int const ialp = CCTK_VarIndex ("ADMBase::alp");
    assert (ialp >= 0);
    if (! variable_is_read[ialp]) {
      CCTK_PARAMWARN ("The lapse is initialised using the file reader, but the group ADMBase::lapse has not been scheduled to be read.  Please set the parameter \"IO::filereader_ID_vars\" accordingly.");
    }
  }
  
  if (CCTK_EQUALS(initial_shift, "read from file")) {
    int const ibetax = CCTK_VarIndex ("ADMBase::betax");
    int const ibetay = CCTK_VarIndex ("ADMBase::betay");
    int const ibetaz = CCTK_VarIndex ("ADMBase::betaz");
    assert (ibetax >= 0);
    assert (ibetay >= 0);
    assert (ibetaz >= 0);
    if (! variable_is_read[ibetax]
        || ! variable_is_read[ibetay]
        || ! variable_is_read[ibetaz]) {
      CCTK_PARAMWARN ("The shift is initialised using the file reader, but the group ADMBase::shift has not been scheduled to be read.  Please set the parameter \"IO::filereader_ID_vars\" accordingly.");
    }
  }
  
  if (CCTK_EQUALS(initial_data, "read from file")) {
    int const igxx = CCTK_VarIndex ("ADMBase::gxx");
    int const igxy = CCTK_VarIndex ("ADMBase::gxy");
    int const igxz = CCTK_VarIndex ("ADMBase::gxz");
    int const igyy = CCTK_VarIndex ("ADMBase::gyy");
    int const igyz = CCTK_VarIndex ("ADMBase::gyz");
    int const igzz = CCTK_VarIndex ("ADMBase::gzz");
    int const ikxx = CCTK_VarIndex ("ADMBase::kxx");
    int const ikxy = CCTK_VarIndex ("ADMBase::kxy");
    int const ikxz = CCTK_VarIndex ("ADMBase::kxz");
    int const ikyy = CCTK_VarIndex ("ADMBase::kyy");
    int const ikyz = CCTK_VarIndex ("ADMBase::kyz");
    int const ikzz = CCTK_VarIndex ("ADMBase::kzz");
    int const ipsi = CCTK_VarIndex ("StaticConformal::psi");
    int const ipsix = CCTK_VarIndex ("StaticConformal::psix");
    int const ipsiy = CCTK_VarIndex ("StaticConformal::psiy");
    int const ipsiz = CCTK_VarIndex ("StaticConformal::psiz");
    int const ipsixx = CCTK_VarIndex ("StaticConformal::psixx");
    int const ipsixy = CCTK_VarIndex ("StaticConformal::psixy");
    int const ipsixz = CCTK_VarIndex ("StaticConformal::psixz");
    int const ipsiyy = CCTK_VarIndex ("StaticConformal::psiyy");
    int const ipsiyz = CCTK_VarIndex ("StaticConformal::psiyz");
    int const ipsizz = CCTK_VarIndex ("StaticConformal::psizz");
    assert (igxx >= 0);
    assert (igxy >= 0);
    assert (igxz >= 0);
    assert (igyy >= 0);
    assert (igyz >= 0);
    assert (igzz >= 0);
    if (! variable_is_read[igxx]
        || ! variable_is_read[igxy]
        || ! variable_is_read[igxz]
        || ! variable_is_read[igyy]
        || ! variable_is_read[igyz]
        || ! variable_is_read[igzz]) {
      CCTK_PARAMWARN ("The metric is initialised using the file reader, but the group ADMBase::metric has not been scheduled to be read.  Please set the parameter \"IO::filereader_ID_vars\" accordingly.");
    }
    assert (ikxx >= 0);
    assert (ikxy >= 0);
    assert (ikxz >= 0);
    assert (ikyy >= 0);
    assert (ikyz >= 0);
    assert (ikzz >= 0);
    if (! variable_is_read[ikxx]
        || ! variable_is_read[ikxy]
        || ! variable_is_read[ikxz]
        || ! variable_is_read[ikyy]
        || ! variable_is_read[ikyz]
        || ! variable_is_read[ikzz]) {
      CCTK_PARAMWARN ("The metric is initialised using the file reader, but the group ADMBase::curv has not been scheduled to be read.  Please set the parameter \"IO::filereader_ID_vars\" accordingly.");
    }
    if (CCTK_EQUALS (metric_type, "physical")) {
      /* do nothing */
    } else if (CCTK_EQUALS (metric_type, "static conformal")) {
      if (CCTK_EQUALS (conformal_storage, "factor")) {
        assert (ipsi >= 0);
        if (! variable_is_read[ipsi]) {
          CCTK_PARAMWARN ("The metric is initialised using the file reader, and the metric_type is \"static conformal\", and the conformal_storage is \"factor\", but the group StaticConformal::confac has not been scheduled to be read.  Please set the parameter \"IO::filereader_ID_vars\" accordingly.");
        }
      } else if (CCTK_EQUALS (conformal_storage, "factor+derivs")) {
        assert (ipsi >= 0);
        if (! variable_is_read[ipsi]) {
          CCTK_PARAMWARN ("The metric is initialised using the file reader, and the metric_type is \"static conformal\", and the conformal_storage is \"factor+derivs\", but the group StaticConformal::confac has not been scheduled to be read.  Please set the parameter \"IO::filereader_ID_vars\" accordingly.");
        }
        if (! variable_is_read[ipsix]
            || ! variable_is_read[ipsiy]
            || ! variable_is_read[ipsiz]) {
          CCTK_PARAMWARN ("The metric is initialised using the file reader, and the metric_type is \"static conformal\", and the conformal_storage is \"factor+derivs\", but the group StaticConformal::confac_1derivs has not been scheduled to be read.  Please set the parameter \"IO::filereader_ID_vars\" accordingly.");
        }
      } else if (CCTK_EQUALS (conformal_storage, "factor+derivs+2nd derivs")) {
        assert (ipsi >= 0);
        if (! variable_is_read[ipsi]) {
          CCTK_PARAMWARN ("The metric is initialised using the file reader, and the metric_type is \"static conformal\", and the conformal_storage is \"factor+derivs+2nd derivs\", but the group StaticConformal::confac has not been scheduled to be read.  Please set the parameter \"IO::filereader_ID_vars\" accordingly.");
        }
        if (! variable_is_read[ipsix]
            || ! variable_is_read[ipsiy]
            || ! variable_is_read[ipsiz]) {
          CCTK_PARAMWARN ("The metric is initialised using the file reader, and the metric_type is \"static conformal\", and the conformal_storage is \"factor+derivs+2nd derivs\", but the group StaticConformal::confac_1derivs has not been scheduled to be read.  Please set the parameter \"IO::filereader_ID_vars\" accordingly.");
        }
        if (! variable_is_read[ipsixx]
            || ! variable_is_read[ipsixy]
            || ! variable_is_read[ipsixz]
            || ! variable_is_read[ipsiyy]
            || ! variable_is_read[ipsiyz]
            || ! variable_is_read[ipsizz]) {
          CCTK_PARAMWARN ("The metric is initialised using the file reader, and the metric_type is \"static conformal\", and the conformal_storage is \"factor+derivs+2nd derivs\", but the group StaticConformal::confac_2derivs has not been scheduled to be read.  Please set the parameter \"IO::filereader_ID_vars\" accordingly.");
        }
      } else {
        CCTK_PARAMWARN ("Unknown conformal_storage type");
      }
    } else {
      CCTK_PARAMWARN ("Unknown metric type");
    }
  }
  
  free (variable_is_read);
}

/** Mark a variable as to be read from the file reader.  */
static void callback (int idx, const char * optstring, void * callback_arg)
{
  assert (idx>=0 && idx<CCTK_NumVars());
  ((char *)callback_arg)[idx] = 1;
}
