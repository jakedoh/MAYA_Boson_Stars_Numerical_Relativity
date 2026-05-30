static void compute_TUPmunu
(const cGH *cctkGH,const int *cctk_lsh,const int *cctk_nghostzones,const CCTK_REAL *dX,CCTK_REAL **metric,const gf_and_gz_struct *prims,
 const CCTK_REAL *gupxy,const CCTK_REAL *gupxz,const CCTK_REAL *gupyz,
 CCTK_REAL **TUPmunu) {

  // These loop extents must be consistent with add_fluxes_and_source_terms_to_hydro_rhss(), since we use TUPmunu there as well.
#pragma omp parallel for
  for(int k=cctk_nghostzones[2];k<cctk_lsh[2]-(cctk_nghostzones[2]-1);k++) for(int j=cctk_nghostzones[1];j<cctk_lsh[1]-(cctk_nghostzones[1]-1);j++) for(int i=cctk_nghostzones[0];i<cctk_lsh[0]-(cctk_nghostzones[0]-1);i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        // First we pull in needed hydrodynamic and metric variables from memory: PART 1.
        // Reading from main memory is a SLOW operation, usually resulting in
        //   cache misses, which
        //   will waste precious runtime. Note that cache misses will often not
        //   show up when using, e.g., gprof. The slowdown due to cache misses
        //   can more than double the amount of time in this routine, so instead
        //   of reading in variables from main memory multiple times, the below
        //   forces us to read in only once, storing in local variables
        //   PRIMS{p,m}{1,2,3}, METRIC{p,m}{1,2}, PRIMS, and METRIC.

        //-----------------------------------------------------------------------------
        // Compute T^{\mu \nu}

        CCTK_REAL METRIC[NUMVARS_FOR_METRIC_FACEVALS]; for(int ii=0;ii<NUMVARS_FOR_METRIC_FACEVALS;ii++) METRIC[ii] = metric[ii][index];
        CCTK_REAL METRIC_LAP_PSI4[NUMVARS_METRIC_AUX]; SET_LAPSE_PSI4(METRIC_LAP_PSI4,METRIC);

        // The "vector" PRIMS represents the primitive variables: vx, vy, vz, Bx, By, and Bz.
        CCTK_REAL PRIMS[6]; // 6 primitives in the set: {vx,vy,vz,Bx,By,Bz}
        for(int ii=0;ii<6;ii++) PRIMS[ii] = prims[ii].gf[index];

        struct output_stats stats; stats.failure_checker=0;
        CCTK_REAL u0L;
        impose_speed_limit_output_u0(METRIC,PRIMS,METRIC_LAP_PSI4[PSI4],METRIC_LAP_PSI4[LAPSEINV],stats, u0L);

        /***********************************************************/
        // Compute b^{\mu} and b^2
        const CCTK_REAL ONE_OVER_LAPSE = 1.0/METRIC_LAP_PSI4[LAPSE];
        const CCTK_REAL ONE_OVER_LAPSE_SQRT_4PI = ONE_OVER_LAPSE*ONE_OVER_SQRT_4PI;
        CCTK_REAL u_x_over_u0_psi4,u_y_over_u0_psi4,u_z_over_u0_psi4;
        CCTK_REAL smallb[NUMVARS_SMALLB];
        compute_smallba_b2_and_u_i_over_u0_psi4(METRIC,METRIC_LAP_PSI4,PRIMS,u0L,ONE_OVER_LAPSE_SQRT_4PI,
                                                u_x_over_u0_psi4,u_y_over_u0_psi4,u_z_over_u0_psi4,smallb);

        /***********************************************************/
        // Next compute T^{\mu \nu} in the FFE limit:
        const CCTK_REAL Psim4 = 1.0/METRIC_LAP_PSI4[PSI4];

        const CCTK_REAL rho0_h_plus_b2 = smallb[SMALLB2];     // Since rho=0 in FFE.
        const CCTK_REAL P_plus_half_b2 = 0.5*smallb[SMALLB2]; // Since P=0 in FFE.

        const CCTK_REAL uUP[4] = { u0L, u0L*PRIMS[VX],u0L*PRIMS[VY],u0L*PRIMS[VZ] };
        // If you like, see Eq 2.119 in Numerical Relativity, by Baumgarte & Shapiro:
        const CCTK_REAL ONE_OVER_LAPSE_SQUARED = SQR(ONE_OVER_LAPSE);
        CCTK_REAL g4up[4][4];

        g4up[0][0] = -ONE_OVER_LAPSE_SQUARED;
        g4up[0][1] = ONE_OVER_LAPSE_SQUARED*METRIC[SHIFTX];
        g4up[0][2] = ONE_OVER_LAPSE_SQUARED*METRIC[SHIFTY];
        g4up[0][3] = ONE_OVER_LAPSE_SQUARED*METRIC[SHIFTZ];
        g4up[1][1] = METRIC[GUPXX]*Psim4 - ONE_OVER_LAPSE_SQUARED*METRIC[SHIFTX]*METRIC[SHIFTX];
        // Note that for i!=j, gupij is not stored in METRIC, since we don't need it in face value calculations.
        g4up[1][2] = gupxy[index] *Psim4 - ONE_OVER_LAPSE_SQUARED*METRIC[SHIFTX]*METRIC[SHIFTY];
        g4up[1][3] = gupxz[index] *Psim4 - ONE_OVER_LAPSE_SQUARED*METRIC[SHIFTX]*METRIC[SHIFTZ];
        g4up[2][2] = METRIC[GUPYY]*Psim4 - ONE_OVER_LAPSE_SQUARED*METRIC[SHIFTY]*METRIC[SHIFTY];
        g4up[2][3] = gupyz[index] *Psim4 - ONE_OVER_LAPSE_SQUARED*METRIC[SHIFTY]*METRIC[SHIFTZ];
        g4up[3][3] = METRIC[GUPZZ]*Psim4 - ONE_OVER_LAPSE_SQUARED*METRIC[SHIFTZ]*METRIC[SHIFTZ];

        // Next compute T^{\mu \nu}:
        // (Eq. 33 in http://arxiv.org/pdf/astro-ph/0503420.pdf):
        // T^{mn} = (rho_0 h + b^2) u^m u^n + (P + 0.5 b^2) g^{mn} - b^m b^n, where m and n both run from 0 to 3.
        //-----------------------------------------------------------------------------
        // Set the T^{\mu \nu} gridfunction, since computing T^{\mu \nu} is expensive
        int counter=0;
        for(int ii=0;ii<4;ii++) for(int jj=ii;jj<4;jj++) { TUPmunu[counter][index] = rho0_h_plus_b2*uUP[ii]*uUP[jj] + P_plus_half_b2*g4up[ii][jj] - smallb[SMALLBT+ii]*smallb[SMALLBT+jj]; counter++; }
      }
}
