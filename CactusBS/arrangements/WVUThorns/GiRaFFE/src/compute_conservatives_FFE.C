void GiRaFFE_compute_conservatives(const CCTK_REAL *PRIMS,  const CCTK_REAL *METRIC, CCTK_REAL *CONSERVS) {
  const CCTK_REAL psi_bssnL = exp(METRIC[PHI]);
  const CCTK_REAL psi2 = psi_bssnL*psi_bssnL;
  const CCTK_REAL psi4 = psi2*psi2;
  const CCTK_REAL sqrtg = psi4*psi2;

  // \gamma_{ij}, computed from \tilde{\gamma}_{ij}
  const CCTK_REAL gxxL = psi4*METRIC[GXX];
  const CCTK_REAL gxyL = psi4*METRIC[GXY];
  const CCTK_REAL gxzL = psi4*METRIC[GXZ];
  const CCTK_REAL gyyL = psi4*METRIC[GYY];
  const CCTK_REAL gyzL = psi4*METRIC[GYZ];
  const CCTK_REAL gzzL = psi4*METRIC[GZZ];

  // Read in magnetic field and momentum variables once from memory, since memory access is expensive:
  const CCTK_REAL BxL = PRIMS[BX_CENTER];
  const CCTK_REAL ByL = PRIMS[BY_CENTER];
  const CCTK_REAL BzL = PRIMS[BZ_CENTER];

  const CCTK_REAL vxL = PRIMS[VX];
  const CCTK_REAL vyL = PRIMS[VY];
  const CCTK_REAL vzL = PRIMS[VZ];

  const CCTK_REAL fourpialpha_inv = 1.0/( 4.0*M_PI*(METRIC[LAPM1] + 1.0) );

  const CCTK_REAL betaxL = METRIC[SHIFTX];
  const CCTK_REAL betayL = METRIC[SHIFTY];
  const CCTK_REAL betazL = METRIC[SHIFTZ];

  const CCTK_REAL B2 = gxxL*BxL*BxL + gyyL*ByL*ByL + gzzL*BzL*BzL
    + 2.0*(gxyL*BxL*ByL + gxzL*BxL*BzL + gyzL*ByL*BzL);


  // NOTE: SIGNIFICANTLY MODIFIED FROM ILLINOISGRMHD VERSION:
  //       velocities in GiRaFFE are defined to be "drift" velocity.
  //       cf. Eqs 47 and 85 in http://arxiv.org/pdf/1310.3274.pdf

  const CCTK_REAL vplusbetaxL = vxL + betaxL;
  const CCTK_REAL vplusbetayL = vyL + betayL;
  const CCTK_REAL vplusbetazL = vzL + betazL;

  const CCTK_REAL vplusbeta_xL = gxxL*vplusbetaxL + gxyL*vplusbetayL + gxzL*vplusbetazL;
  const CCTK_REAL vplusbeta_yL = gxyL*vplusbetaxL + gyyL*vplusbetayL + gyzL*vplusbetazL;
  const CCTK_REAL vplusbeta_zL = gxzL*vplusbetaxL + gyzL*vplusbetayL + gzzL*vplusbetazL;

  /*
   * Comments:
   * Eq. 85 in https://arxiv.org/pdf/1310.3274.pdf:
   * v^i = 4 pi alpha * (gamma^{ij} tilde{S}_j) / (sqrtgamma * B^2) - beta^i
   * which implies that
   * (v^i + beta^i)*(sqrtgamma * B^2)/(4 pi alpha) = gamma^{ij} tilde{S}_j
   * Multiply both sides by gamma_{ik}:
   * gamma_{ik} (v^i + beta^i)*(sqrtgamma * B^2)/(4 pi alpha) = gamma_{ik} gamma^{ij} tilde{S}_j
   *
   * -> tilde{S}_k = gamma_{ik} (v^i + beta^i)*(sqrtgamma * B^2)/(4 pi alpha)
   */

  CONSERVS[STILDEX] = vplusbeta_xL * sqrtg * B2 * fourpialpha_inv;
  CONSERVS[STILDEY] = vplusbeta_yL * sqrtg * B2 * fourpialpha_inv;
  CONSERVS[STILDEZ] = vplusbeta_zL * sqrtg * B2 * fourpialpha_inv;
}
