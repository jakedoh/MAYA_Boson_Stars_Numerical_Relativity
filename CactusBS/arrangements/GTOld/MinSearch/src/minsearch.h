// MinSearch: minsearch.h

#ifdef __cplusplus
extern "C" {
#endif

// we require parameter_type == 0
enum target_types {invalid_type = -1, parameter_type=0, grid_scalar_type};

//*******************************//
//****** External Routines ******//
//*******************************//
extern void MinSearch_ParamCheck(CCTK_ARGUMENTS);
extern void MinSearch_SetupOutputFiles(CCTK_ARGUMENTS);
extern void MinSearch_Init(CCTK_ARGUMENTS);
extern void MinSearch_SetSphericalSurface(CCTK_ARGUMENTS);
extern void MinSearch_ReduceGlobalMinimum(CCTK_ARGUMENTS);
extern void MinSearch_Output(CCTK_ARGUMENTS);
extern void MinSearch_FindGlobalMinimum(CCTK_ARGUMENTS);
extern void MinSearch_InterpolateMinimum(CCTK_ARGUMENTS);

extern void MinSearch_FindLocalizedMinimum(CCTK_ARGUMENTS);
extern void MinSearch_ReduceLocalizedMinimum(CCTK_ARGUMENTS);

/* Internal use only */
extern CCTK_INT MinSearch_TraverseTargetString(cGH *cctkGH, CCTK_STRING targetstring, 
        CCTK_INT (*callback)(CCTK_POINTER_TO_CONST target, const CCTK_INT target_type, CCTK_STRING targetname, CCTK_POINTER callback_arg),
        CCTK_POINTER callback_arg);

#ifdef __cplusplus
}
#endif


