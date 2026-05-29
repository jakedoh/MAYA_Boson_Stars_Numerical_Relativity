  
  // define pointer array for past timelevels!
  vector<CCTK_REAL*> uBondiP(max_timelevels, (CCTK_REAL*)NULL);
  vector<CCTK_COMPLEX*> Psi4P(max_timelevels, (CCTK_COMPLEX*)NULL);
  vector<CCTK_COMPLEX*> NewsP(max_timelevels, (CCTK_COMPLEX*)NULL);
  vector<CCTK_COMPLEX*> NewsBP(max_timelevels, (CCTK_COMPLEX*)NULL);
  vector<CCTK_COMPLEX*> linStrainP(max_timelevels, (CCTK_COMPLEX*)NULL);

  // Get current timelevel...
  // We assume here that the two vars uBondi[0] amd uBondi[1] are right next to each
  // other in memory, i.e. the pointer to uBondi[0] can be used to access uBondi[1] as well!
  const int varindex_uBondi = CCTK_VarIndex("NullNews::uBondi[0]");
  const int varindex_Psi4 = CCTK_VarIndex("NullNews::Psi4[0]");
  const int varindex_News = CCTK_VarIndex("NullNews::News[0]");
  const int varindex_NewsB = CCTK_VarIndex("NullNews::NewsB[0]");
  const int varindex_linStrain = CCTK_VarIndex("NullNews::linStrain[0]");

  if (varindex_uBondi < 0 || varindex_Psi4 < 0 || varindex_News < 0 || varindex_NewsB < 0)
     CCTK_WARN(0, "Error getting variable index!");

  if (compute_lin_strain && varindex_linStrain < 0)
     CCTK_WARN(0, "Error getting variable index!");

  uBondiP[0] = (CCTK_REAL*) CCTK_VarDataPtrI(cctkGH, 0, varindex_uBondi);
  Psi4P[0] = (CCTK_COMPLEX*) CCTK_VarDataPtrI(cctkGH, 0, varindex_Psi4);
  NewsP[0] = (CCTK_COMPLEX*) CCTK_VarDataPtrI(cctkGH, 0, varindex_News);
  NewsBP[0] = (CCTK_COMPLEX*) CCTK_VarDataPtrI(cctkGH, 0, varindex_NewsB);

  if (compute_lin_strain)
     linStrainP[0] = (CCTK_COMPLEX*) CCTK_VarDataPtrI(cctkGH, 0, varindex_linStrain);

  // Get past timelevels...
  // We assume here that the two vars uBondi[0] amd uBondi[1] are right next to each
  // other in memory, i.e. the pointer to uBondi[0] can be used to access uBondi[1] as well!
  const int varindex_uBondi_past = CCTK_VarIndex("NullNews::uBondi_past[0]");
  const int varindex_Psi4_past = CCTK_VarIndex("NullNews::Psi4_past[0]");
  const int varindex_News_past = CCTK_VarIndex("NullNews::News_past[0]");
  const int varindex_NewsB_past = CCTK_VarIndex("NullNews::NewsB_past[0]");
  const int varindex_linStrain_past = CCTK_VarIndex("NullNews::linStrain_past[0]");

  if (varindex_uBondi_past < 0 || varindex_Psi4_past < 0 || varindex_News_past < 0 || varindex_NewsB_past < 0)
     CCTK_WARN(0, "Error getting variable index!");

  if (compute_lin_strain && varindex_linStrain_past < 0)
     CCTK_WARN(0, "Error getting variable index!");

  const int offset = 2*null_lsh[0]*null_lsh[1];

  for (int tl=1; tl < max_timelevels; ++tl)
  {
     uBondiP[tl] = &((CCTK_REAL*) CCTK_VarDataPtrI(cctkGH, 0, varindex_uBondi_past))[offset*(tl-1)];
     Psi4P[tl] = &((CCTK_COMPLEX*) CCTK_VarDataPtrI(cctkGH, 0, varindex_Psi4_past))[offset*(tl-1)];
     NewsP[tl] = &((CCTK_COMPLEX*) CCTK_VarDataPtrI(cctkGH, 0, varindex_News_past))[offset*(tl-1)];
     NewsBP[tl] = &((CCTK_COMPLEX*) CCTK_VarDataPtrI(cctkGH, 0, varindex_NewsB_past))[offset*(tl-1)];
     
     if (compute_lin_strain)
        linStrainP[tl] = &((CCTK_COMPLEX*) CCTK_VarDataPtrI(cctkGH, 0, varindex_linStrain_past))[offset*(tl-1)];
  }
