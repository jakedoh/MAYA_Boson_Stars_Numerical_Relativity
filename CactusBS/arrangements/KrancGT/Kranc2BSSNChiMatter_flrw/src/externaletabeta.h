const CCTK_REAL *externalEtaBeta = (CCTK_REAL *)CCTK_VarDataPtr(cctkGH, 0, "ExternalEtaBeta::exetaBeta");
if(externalEtaBeta == NULL)
{
  CCTK_WARN(CCTK_WARN_COMPLAIN,"Could not get data pointer for ExternalEtaBeta::exetaBeta. Is ExternalEtaBeta active?");
  /* NOTREACHED */
}
