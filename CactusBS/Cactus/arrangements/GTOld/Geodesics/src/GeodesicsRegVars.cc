#include "Geodesics.hh"

extern "C" void  Geodesics_RegVars(CCTK_ARGUMENTS);

static void register_evolved(string v1, string v2)
{
  MoLRegisterEvolvedGroup(CCTK_GroupIndex(v1.c_str()),CCTK_GroupIndex(v2.c_str()));
}

static void register_saveandrestore(string v1)
{
  MoLRegisterSaveAndRestoreGroup(CCTK_GroupIndex(v1.c_str()));
}

static void register_constrained(string v1)
{
  MoLRegisterConstrainedGroup(CCTK_GroupIndex(v1.c_str()));
}


extern "C" void Geodesics_RegVars(CCTK_ARGUMENTS) 
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Register the evolved groups for MoL
  register_evolved("Geodesics::geo_pos","Geodesics::geo_pos_rhs");
  register_evolved("Geodesics::geo_u","Geodesics::geo_u_rhs");
  if(CCTK_Equals(geo_moltest,"yes")) register_evolved("Geodesics::geo_tg","Geodesics::geo_t_rhsg");

  // This should be registered...
  // for some reason MoL wont accumulate to max_num_constrained_array_vars (MoL bug?)
//  register_constrained("Geodesics::geo_iu");
//  if(CCTK_Equals(em_output,"yes")) register_constrained("Geodesics::geo_em");
//  if(CCTK_Equals(geo_test,"yes")) register_constrained("Geodesics::geo_test");

  // Register SnR vars
  if(CCTK_IsThornActive("Kranc2BSSNChi") != 0) {
    register_saveandrestore("Kranc2BSSNChi::alpha_group");
    register_saveandrestore("Kranc2BSSNChi::beta_group");
    register_saveandrestore("Kranc2BSSNChi::chi_group");
    register_saveandrestore("Kranc2BSSNChi::h_group");
  }

  if(CCTK_IsThornActive("Kranc2BSSNChiMatter") != 0) {  
    register_saveandrestore("Kranc2BSSNChiMatter::alpha_group");
    register_saveandrestore("Kranc2BSSNChiMatter::beta_group");
    register_saveandrestore("Kranc2BSSNChiMatter::chi_group");
    register_saveandrestore("Kranc2BSSNChiMatter::h_group");
  }

  if(CCTK_IsThornActive("EMevoPlasmaChi") != 0) {
    register_saveandrestore("EMevoPlasmaChi::Bff_group");
    register_saveandrestore("EMevoPlasmaChi::Eff_group");
  }

  if(CCTK_IsThornActive("Pic") != 0) {
    register_saveandrestore("Pic::pic_rhs");
  }
}
