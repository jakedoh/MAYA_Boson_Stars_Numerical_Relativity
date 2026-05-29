
void TOVSolverC_export_local_variables(CCTK_REAL **exported_TOV_press_1d, CCTK_REAL **exported_TOV_m_1d, CCTK_REAL **exported_TOV_phi_1d, CCTK_REAL **exported_TOV_rbar_1d, CCTK_REAL **exported_TOV_r_1d)
{
  int i;

  *exported_TOV_press_1d = TOV_press_1d;
  *exported_TOV_m_1d     = TOV_m_1d;
  *exported_TOV_phi_1d   = TOV_phi_1d;
  *exported_TOV_rbar_1d  = TOV_rbar_1d;
  *exported_TOV_r_1d     = TOV_r_1d;
}

