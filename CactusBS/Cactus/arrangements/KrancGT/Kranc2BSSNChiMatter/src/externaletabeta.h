const CCTK_REAL *externalEtaBeta = (CCTK_REAL *)CCTK_VarDataPtr(cctkGH, 0, "ExternalEtaBeta::etaBeta");
if(externalEtaBeta == NULL)
{
  CCTK_WARN(CCTK_WARN_COMPLAIN,"Could not get data pointer for ExternalEtaBeta::etaBeta. Is ExternalEtaBeta active?");
  /* NOTREACHED */
}
