//-----------------------------------------------------------------------------
// Compute the flux for advecting S_i .
//-----------------------------------------------------------------------------
static inline void GRFFE__S_i__flux(const int i,const int j,const int k,const int flux_dirn, CCTK_REAL *Ul, CCTK_REAL *Ur, CCTK_REAL *FACEVAL,const CCTK_REAL *FACEVAL_LAPSE_PSI4,
                                    CCTK_REAL &cmax,CCTK_REAL &cmin,
                                    CCTK_REAL &st_x_flux,CCTK_REAL &st_y_flux,CCTK_REAL &st_z_flux) {
  const CCTK_REAL psi4 = FACEVAL_LAPSE_PSI4[PSI4];
  const CCTK_REAL psi6 = FACEVAL_LAPSE_PSI4[PSI4]*FACEVAL_LAPSE_PSI4[PSI2];
  const CCTK_REAL psim4 = 1.0/(psi4);

  const CCTK_REAL alpha_sqrt_gamma = FACEVAL_LAPSE_PSI4[LAPSE]*psi6;
  const CCTK_REAL ONE_OVER_LAPSE = 1.0/FACEVAL_LAPSE_PSI4[LAPSE];
  const CCTK_REAL ONE_OVER_LAPSE_SQUARED=SQR(ONE_OVER_LAPSE);

  //Compute face velocities
  // Begin by computing u0
  output_stats stats; stats.failure_checker=0;
  CCTK_REAL u0_r,u0_l;
  impose_speed_limit_output_u0(FACEVAL,Ur,psi4,ONE_OVER_LAPSE,stats,u0_r);
  impose_speed_limit_output_u0(FACEVAL,Ul,psi4,ONE_OVER_LAPSE,stats,u0_l);

  //Next compute b^{\mu}, the magnetic field measured in the comoving fluid frame:
  const CCTK_REAL ONE_OVER_LAPSE_SQRT_4PI = ONE_OVER_LAPSE*ONE_OVER_SQRT_4PI;
  /***********************************************************/
  /********** RIGHT FACE ************/
  // Note that smallbr[4] = b^a defined in Gammie's paper, on the right face.
  CCTK_REAL u_x_over_u0_psi4r,u_y_over_u0_psi4r,u_z_over_u0_psi4r;
  CCTK_REAL smallbr[NUMVARS_SMALLB];
  // Compute b^{a}, b^2, and u_i over u^0
  compute_smallba_b2_and_u_i_over_u0_psi4(FACEVAL,FACEVAL_LAPSE_PSI4,Ur,u0_r,ONE_OVER_LAPSE_SQRT_4PI,
                                          u_x_over_u0_psi4r,u_y_over_u0_psi4r,u_z_over_u0_psi4r,smallbr);
  // Then compute u_xr,u_yr, and u_zr. We need to set the zeroth component so we can specify U_LOWER{r,l}[{UX,UY,UZ}] (UX=1,UY=2,UZ=3).
  const CCTK_REAL U_LOWERr[4] = { 0.0, u_x_over_u0_psi4r*u0_r*FACEVAL_LAPSE_PSI4[PSI4], u_y_over_u0_psi4r*u0_r*FACEVAL_LAPSE_PSI4[PSI4],
                                  u_z_over_u0_psi4r*u0_r*FACEVAL_LAPSE_PSI4[PSI4] };
  /********** LEFT FACE ************/
  // Note that smallbl[4] = b^a defined in Gammie's paper, on the left face.
  CCTK_REAL u_x_over_u0_psi4l,u_y_over_u0_psi4l,u_z_over_u0_psi4l;
  CCTK_REAL smallbl[NUMVARS_SMALLB];
  // Compute b^{a}, b^2, and u_i over u^0
  compute_smallba_b2_and_u_i_over_u0_psi4(FACEVAL,FACEVAL_LAPSE_PSI4,Ul,u0_l,ONE_OVER_LAPSE_SQRT_4PI,
                                          u_x_over_u0_psi4l,u_y_over_u0_psi4l,u_z_over_u0_psi4l,smallbl);
  // Then compute u_xr,u_yr, and u_zr. We need to set the zeroth component so we can specify U_LOWER{r,l}[{UX,UY,UZ}]
  const CCTK_REAL U_LOWERl[4] = { 0.0, u_x_over_u0_psi4l*u0_l*FACEVAL_LAPSE_PSI4[PSI4], u_y_over_u0_psi4l*u0_l*FACEVAL_LAPSE_PSI4[PSI4],
                                  u_z_over_u0_psi4l*u0_l*FACEVAL_LAPSE_PSI4[PSI4] };
  /***********************************************************/

  // v02 is the speed of light in FFE, which is 1.0 in G=c=1.0.
  const CCTK_REAL v02r=1.0,v02l=1.0;

  int offset=flux_dirn-1;

  CCTK_REAL cplusr,cminusr,cplusl,cminusl;
  find_cp_cm(cplusr,cminusr,v02r,u0_r,
             Ur[VX+offset],ONE_OVER_LAPSE_SQUARED,FACEVAL[SHIFTX+offset],psim4,FACEVAL[GUPXX+offset]);
  find_cp_cm(cplusl,cminusl,v02l,u0_l,
             Ul[VX+offset],ONE_OVER_LAPSE_SQUARED,FACEVAL[SHIFTX+offset],psim4,FACEVAL[GUPXX+offset]);

  // Then compute cmax, cmin. This is required for the HLL flux.
  const CCTK_REAL cmaxL =  MAX(0.0,MAX(cplusl,cplusr));
  const CCTK_REAL cminL = -MIN(0.0,MIN(cminusl,cminusr));


  //*********************************************************************
  // momentum flux = \alpha \sqrt{\gamma} T^m_j, where m is the current flux direction (the m index)
  //*********************************************************************
  // b_j = g_{ij} (b^i + b^t shift^i), g_{ij} = physical metric
  //CCTK_REAL sbtr=0,sbtl=0;
  CCTK_REAL smallb_lowerr[NUMVARS_SMALLB],smallb_lowerl[NUMVARS_SMALLB];
  lower_4vector_output_spatial_part(psi4,FACEVAL,smallbr,smallb_lowerr);
  lower_4vector_output_spatial_part(psi4,FACEVAL,smallbl,smallb_lowerl);

  /********** Flux for S_x **********/
  // [S_x flux] = \alpha \sqrt{\gamma} T^m_x, where m is the current flux direction (the m index)
  //    Again, offset = 0 for reconstruction in x direction, 1 for y, and 2 for z
  //    Note that kronecker_delta[flux_dirn][0] = { 1 if flux_dirn==1, 0 otherwise }.
  /********** RIGHT FACE ************/
  // Compute a couple useful hydro quantities:
  const CCTK_REAL rho0_h_plus_b2_r = smallbr[SMALLB2];     // Since rho=0 in GRFFE
  const CCTK_REAL P_plus_half_b2_r = 0.5*smallbr[SMALLB2]; // Since P=0 in GRFFE
  /********** LEFT FACE *************/
  // Compute a couple useful hydro quantities:
  const CCTK_REAL rho0_h_plus_b2_l = smallbl[SMALLB2];     // Since rho=0 in GRFFE
  const CCTK_REAL P_plus_half_b2_l = 0.5*smallbl[SMALLB2]; // Since P=0 in GRFFE
  CCTK_REAL Fr = alpha_sqrt_gamma*( rho0_h_plus_b2_r*(u0_r*Ur[VX+offset])*U_LOWERr[UX]
                                    + P_plus_half_b2_r*kronecker_delta[flux_dirn][0] - smallbr[SMALLBX+offset]*smallb_lowerr[SMALLBX] );
  CCTK_REAL Fl = alpha_sqrt_gamma*( rho0_h_plus_b2_l*(u0_l*Ul[VX+offset])*U_LOWERl[UX]
                                    + P_plus_half_b2_l*kronecker_delta[flux_dirn][0] - smallbl[SMALLBX+offset]*smallb_lowerl[SMALLBX] );

  //        S_x =\alpha\sqrt{\gamma}( T^0_x )
  const CCTK_REAL st_x_r = alpha_sqrt_gamma*( rho0_h_plus_b2_r*u0_r*U_LOWERr[UX] - smallbr[SMALLBT]*smallb_lowerr[SMALLBX] );
  const CCTK_REAL st_x_l = alpha_sqrt_gamma*( rho0_h_plus_b2_l*u0_l*U_LOWERl[UX] - smallbl[SMALLBT]*smallb_lowerl[SMALLBX] );

  // HLL step for Sx:
  st_x_flux = (cminL*Fr + cmaxL*Fl - cminL*cmaxL*(st_x_r-st_x_l) )/(cmaxL + cminL);

  /********** Flux for S_y **********/
  // [S_y flux] = \alpha \sqrt{\gamma} T^m_y, where m is the current flux direction (the m index)
  //    Again, offset = 1 for reconstruction in x direction, 2 for y, and 3 for z
  //    Note that kronecker_delta[flux_dirn][1] = { 1 if flux_dirn==2, 0 otherwise }.
  Fr = alpha_sqrt_gamma*( rho0_h_plus_b2_r*(u0_r*Ur[VX+offset])*U_LOWERr[UY] + P_plus_half_b2_r*kronecker_delta[flux_dirn][1]
                          - smallbr[SMALLBX+offset]*smallb_lowerr[SMALLBY] );
  Fl = alpha_sqrt_gamma*( rho0_h_plus_b2_l*(u0_l*Ul[VX+offset])*U_LOWERl[UY] + P_plus_half_b2_l*kronecker_delta[flux_dirn][1]
                          - smallbl[SMALLBX+offset]*smallb_lowerl[SMALLBY] );

  //        S_y =\alpha\sqrt{\gamma}( T^0_y )
  const CCTK_REAL st_y_r = alpha_sqrt_gamma*( rho0_h_plus_b2_r*u0_r*U_LOWERr[UY] - smallbr[SMALLBT]*smallb_lowerr[SMALLBY] );
  const CCTK_REAL st_y_l = alpha_sqrt_gamma*( rho0_h_plus_b2_l*u0_l*U_LOWERl[UY] - smallbl[SMALLBT]*smallb_lowerl[SMALLBY] );

  // HLL step for Sy:
  st_y_flux = (cminL*Fr + cmaxL*Fl - cminL*cmaxL*(st_y_r-st_y_l) )/(cmaxL + cminL);
  /********** Flux for S_z **********/
  // [S_z flux] = \alpha \sqrt{\gamma} T^m_z, where m is the current flux direction (the m index)
  //    Again, offset = 1 for reconstruction in x direction, 2 for y, and 3 for z
  //    Note that kronecker_delta[flux_dirn][2] = { 1 if flux_dirn==3, 0 otherwise }.
  Fr = alpha_sqrt_gamma*( rho0_h_plus_b2_r*(u0_r*Ur[VX+offset])*U_LOWERr[UZ] + P_plus_half_b2_r*kronecker_delta[flux_dirn][2]
                          - smallbr[SMALLBX+offset]*smallb_lowerr[SMALLBZ] );
  Fl = alpha_sqrt_gamma*( rho0_h_plus_b2_l*(u0_l*Ul[VX+offset])*U_LOWERl[UZ] + P_plus_half_b2_l*kronecker_delta[flux_dirn][2]
                          - smallbl[SMALLBX+offset]*smallb_lowerl[SMALLBZ] );

  //        S_z =\alpha\sqrt{\gamma}( T^0_z )
  const CCTK_REAL st_z_r = alpha_sqrt_gamma*( rho0_h_plus_b2_r*u0_r*U_LOWERr[UZ] - smallbr[SMALLBT]*smallb_lowerr[SMALLBZ] );
  const CCTK_REAL st_z_l = alpha_sqrt_gamma*( rho0_h_plus_b2_l*u0_l*U_LOWERl[UZ] - smallbl[SMALLBT]*smallb_lowerl[SMALLBZ] );

  // HLL step for Sz:
  st_z_flux = (cminL*Fr + cmaxL*Fl - cminL*cmaxL*(st_z_r-st_z_l) )/(cmaxL + cminL);

  cmax = cmaxL;
  cmin = cminL;
}
