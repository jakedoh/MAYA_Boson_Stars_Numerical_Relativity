
void TOV_C_AllocateMemory(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  assert(TOV_Surface==0);
  assert(TOV_R_Surface==0);
  assert(TOV_Atmosphere==0);

  assert(TOV_r_1d==0);
  assert(TOV_rbar_1d==0);
  assert(TOV_press_1d==0);
  assert(TOV_phi_1d==0);
  assert(TOV_m_1d==0);

  TOV_Surface = malloc(TOV_Num_TOVs * sizeof(CCTK_REAL));
  TOV_R_Surface = malloc(TOV_Num_TOVs * sizeof(CCTK_REAL));
  TOV_Atmosphere = malloc(TOV_Num_TOVs * sizeof(CCTK_REAL));

  TOV_r_1d = malloc(TOV_Num_Radial * TOV_Num_TOVs * sizeof(CCTK_REAL));
  TOV_rbar_1d = malloc(TOV_Num_Radial * TOV_Num_TOVs * sizeof(CCTK_REAL));
  TOV_press_1d = malloc(TOV_Num_Radial * TOV_Num_TOVs * sizeof(CCTK_REAL));
  TOV_phi_1d = malloc(TOV_Num_Radial * TOV_Num_TOVs * sizeof(CCTK_REAL));
  TOV_m_1d = malloc(TOV_Num_Radial * TOV_Num_TOVs * sizeof(CCTK_REAL));
}

void TOV_C_FreeMemory (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  assert(TOV_Surface!=0);
  assert(TOV_R_Surface!=0);
  assert(TOV_Atmosphere!=0);

  assert(TOV_r_1d!=0);
  assert(TOV_rbar_1d!=0);
  assert(TOV_press_1d!=0);
  assert(TOV_phi_1d!=0);
  assert(TOV_m_1d!=0);

  free(TOV_Surface);
  free(TOV_R_Surface);
  free(TOV_Atmosphere);

  free(TOV_r_1d);
  free(TOV_rbar_1d);
  free(TOV_press_1d);
  free(TOV_phi_1d);
  free(TOV_m_1d);

  TOV_Surface=0;
  TOV_R_Surface=0;
  TOV_Atmosphere=0;

  TOV_r_1d=0;
  TOV_rbar_1d=0;
  TOV_press_1d=0;
  TOV_phi_1d=0;
  TOV_m_1d=0;
}

/* - utility routine
   - fills an real-array 'var' of size 'i' with value 'r' */
void TOV_C_fill(CCTK_REAL *var, CCTK_INT i, CCTK_REAL r)
{
  for (i-- ;i >= 0; i--)
    var[i]=r;
}

void TOV_Copy(CCTK_INT size, CCTK_REAL *var_p, CCTK_REAL *var)
{
    for(; size; )
    {
        --size;
        var_p[size] = var[size];
    }
}
