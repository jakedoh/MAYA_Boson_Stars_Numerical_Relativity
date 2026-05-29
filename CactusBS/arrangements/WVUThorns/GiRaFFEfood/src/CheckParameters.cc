// Please consult GiRaFFEfood/README for the list of authors and current maintainers //
/*
 * Check the parameter usage
 * */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void GiRaFFEfood_CheckParameters(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INFO ("Setting Initial Data for GiRaFFE");

  if (CCTK_Equals(test_case,"FastWave")) {
     CCTK_INFO(" Feeding the GiRaFFE... with the Fast Wave Solution");
  }
  else if (CCTK_Equals(test_case,"AlfvenWave")) {
     CCTK_INFO(" Feeding the GiRaFFE... with the Alfven Wave Solution");
  }
  else if (CCTK_Equals(test_case,"DegenAlfvenWave")) {
     CCTK_INFO(" Feeding the GiRaFFE... with the Degenerate Alfven Wave Solution");
  }
  else if (CCTK_Equals(test_case,"ThreeAlfvenWave")) {
     CCTK_INFO(" Feeding the GiRaFFE... with the Three Alfven Waves Solution");
  }
  else if(CCTK_Equals(test_case,"FFEBreakdown")) {
     CCTK_INFO(" Feeding the GiRaFFE... with the Force Free Breakdown Solution");
  }
  else if(CCTK_Equals(test_case,"SplitMonopole")) {
     CCTK_INFO(" Feeding the GiRaFFE... with the Split Monopole Solution");
  }
  else if(CCTK_Equals(test_case,"ExactWald")) {
     CCTK_INFO(" Feeding the GiRaFFE... with the Exact Wald Solution");
  }
  else if(CCTK_Equals(test_case,"MagnetoWald")) {
     CCTK_INFO(" Feeding the GiRaFFE... with the Magnetic Wald Solution");
  }
  else if(CCTK_Equals(test_case,"AlignedRotator")) {
     CCTK_INFO(" Feeding the GiRaFFE... with the Aligned Rotator Solution");
  }
  else if(CCTK_Equals(test_case,"VerticalBfield")) {
     CCTK_INFO(" Feeding the GiRaFFE... with the Vertical B field");
  }
  else {
    CCTK_VParamWarn(CCTK_THORNSTRING, "Test '%s' is not implemented", test_case);
  }

}
