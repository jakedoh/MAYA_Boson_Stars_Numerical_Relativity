// MassDistribution: matterlibs.h

#ifdef __cplusplus
extern "C" {
#endif



/* known special variable names */
#define NUM_SPECIAL_VARIABLE_TYPES 5
#define SPECIAL_VARIABLE_TYPE_NAMES {"gridfunction", "specific energy", "u0", "E_Newton", "E_SR"}
enum bin_variable_type {binvar_gridfunction=0, binvar_specific_energy, binvar_u0, binvar_ENewton, binvar_ESR};
extern const char *special_variables[NUM_SPECIAL_VARIABLE_TYPES];

/* known mass variable names */
#define NUM_MASS_VARIABLE_TYPES 3
#define MASS_VARIABLE_TYPE_NAMES {"gridfunction", "Mass", "debug1"}
enum mass_variable_type {massvar_gridfunction=0, massvar_relativistic_rest_mass_density, massvar_debug1};
extern const char *mass_variables[NUM_MASS_VARIABLE_TYPES];

/* known scalings */
#define NUM_BIN_SCALINGS 2
#define BIN_SCALING_NAMES {"linear", "logarithmic"}
enum bin_scaling_type {bin_scaling_linear=0, bin_scaling_logarithmic};
extern const char *bin_scalings[NUM_BIN_SCALINGS];

/* translation of Cactus KEYWORD arguments to ints for speedup */
enum operation {find_minmax, integrate};
enum bin_mode {automatic, automatic_volume, manual, by_mass_fraction};

/* extra bins for total mass etc. */
#define TOTAL_NUMBER_OF_BINS (number_of_bins+3)
#define TOTAL_MASS_BIN (number_of_bins)
#define BELOW_MASS_BIN (number_of_bins+1)
#define ABOVE_MASS_BIN (number_of_bins+2)
#define NUMBER_OF_LABELS (number_of_bins+1)

/*******************************
 ****** External Routines ******
 *******************************/
extern void MassDistribution_InitMinMax(CCTK_ARGUMENTS);
extern void MassDistribution_ZeroMasses(CCTK_ARGUMENTS);
extern void MassDistribution_SumLocalRestMass(CCTK_ARGUMENTS);
extern void MassDistribution_SumGlobalRestMass(CCTK_ARGUMENTS);
extern void MassDistribution_FindMinMax(CCTK_ARGUMENTS);
extern void MassDistribution_ReduceMinMax(CCTK_ARGUMENTS);
extern void MassDistribution_Output(CCTK_ARGUMENTS);
extern void MassDistribution_RescaleBins(CCTK_ARGUMENTS);

/* Internal use only */
extern void MassDistribution_get_centre(CCTK_ARGUMENTS, int c, CCTK_REAL *xc, CCTK_REAL *yc, CCTK_REAL *zc);
extern void MassDistribution_sumLocalMass(CCTK_ARGUMENTS, enum operation opcode);
extern void MassDistribution_get_do_level_every(const cGH * cctkGH, CCTK_INT * do_every);
extern void MassDistribution_resolve_special_variable(cGH *cctkGH, const char *varstring, int *vartype, CCTK_REAL **varptr, int num_vartypes, const char *vartypes[]);
extern int MassDistribution_validate_special_variable(cGH *cctkGH, const char *varstring, int num_vartypes, const char *vartypes[]);
extern int MassDistribution_resolve_keyword_parameter(const char *name, const char *kw, int num_kws, const char *kws[]);

extern void query_symmetries(const char *applied_symmetries, int *ncopies, double symsigns[8][3]);

#ifdef __cplusplus
}
#endif

