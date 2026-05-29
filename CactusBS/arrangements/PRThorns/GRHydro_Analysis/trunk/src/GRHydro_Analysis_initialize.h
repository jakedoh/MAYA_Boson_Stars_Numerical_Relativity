
struct GlobalGRHydroAnalysis {
  CCTK_REAL UnitVolume;
  CCTK_INT isPUGH;
  CCTK_INT symXref;
  CCTK_INT symYref;
  CCTK_INT symZref;
  CCTK_INT symRot90;
  CCTK_INT symRot180;
  CCTK_INT  masktype;
  CCTK_REAL maskvalue;
};

extern struct GlobalGRHydroAnalysis GlobalGRHydroAnalysisDetect; 

