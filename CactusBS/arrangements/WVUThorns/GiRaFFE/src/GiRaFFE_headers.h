// To safeguard against double-including this header file:
#ifndef GIRAFFE_HEADERS_H_
#define GIRAFFE_HEADERS_H_

#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )
#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )
#define SQR(x) ((x) * (x))
#define ONE_OVER_SQRT_4PI 0.282094791773878143474039725780

#define VERR_DEF_PARAMS __LINE__, __FILE__, CCTK_THORNSTRING

// The order here MATTERS, as we assume that GUPXX+1=GUPYY, etc.
static const int PHI=0,PSI=1,GXX=2,GXY=3,GXZ=4,GYY=5,GYZ=6,GZZ=7,
  LAPM1=8,SHIFTX=9,SHIFTY=10,SHIFTZ=11,GUPXX=12,GUPYY=13,GUPZZ=14,
  NUMVARS_FOR_METRIC_FACEVALS=15; //<-- Be _sure_ to set this correctly, or you'll have memory access bugs!

// These are not used for facevals in the reconstruction step, but boy are they useful anyway.
static const int GUPXY=15,GUPXZ=16,GUPYZ=17,
  NUMVARS_FOR_METRIC=18; //<-- Be _sure_ to set this correctly, or you'll have memory access bugs!

// The order here MATTERS, and must be consistent with the order in the in_prims[] array in driver_evaluate_FFE_rhs.C.
static const int VX=0,VY=1,VZ=2,
  BX_CENTER=3,BY_CENTER=4,BZ_CENTER=5,BX_STAGGER=6,BY_STAGGER=7,BZ_STAGGER=8,
  VXR=9,VYR=10,VZR=11,VXL=12,VYL=13,VZL=14,MAXNUMVARS=15;  //<-- Be _sure_ to define MAXNUMVARS appropriately!

static const int UT=0,UX=1,UY=2,UZ=3;

// The "I" suffix denotes interpolation. In other words, these
//    definitions are used for interpolation ONLY. The order here
//    matters as well!
static const int SHIFTXI=0,SHIFTYI=1,SHIFTZI=2,GUPXXI=3,GUPXYI=4,GUPXZI=5,GUPYYI=6,GUPYZI=7,GUPZZI=8,
  PSII=9,LAPM1I=10,A_XI=11,A_YI=12,A_ZI=13,LAPSE_PSI2I=14,LAPSE_OVER_PSI6I=15,MAXNUMINTERP=16;

// Again, the order here MATTERS, since we assume in the code that, e.g., smallb[0]=b^t, smallb[3]=b^z, etc.
static const int SMALLBT=0,SMALLBX=1,SMALLBY=2,SMALLBZ=3,SMALLB2=4,NUMVARS_SMALLB=5;

// Again, the order here MATTERS, since we assume in the code that, CONSERV[STILDEX+1] = \tilde{S}_y
static const int STILDEX=0,STILDEY=1,STILDEZ=2,NUM_CONSERVS=3;

static const int LAPSE=0,PSI2=1,PSI4=2,PSI6=3,PSIM4=4,LAPSEINV=5,NUMVARS_METRIC_AUX=6;
#define SET_LAPSE_PSI4(array_name,METRIC)   {                   \
      array_name[LAPSE] = METRIC[LAPM1]+1.0;                    \
      array_name[PSI2]  = exp(2.0*METRIC[PHI]);                 \
      array_name[PSI4]  = SQR(array_name[PSI2]);                \
      array_name[PSI6]  = array_name[PSI4]*array_name[PSI2];    \
      array_name[PSIM4]  = 1.0/array_name[PSI4];                \
      array_name[LAPSEINV]  = 1.0/array_name[LAPSE];            \
  }

// Keeping track of ghostzones between routines is a nightmare, so
//   we instead attach ghostzone info to each gridfunction and set
//   the ghostzone information correctly within each routine.
struct gf_and_gz_struct {
  CCTK_REAL *gf;
  int gz_lo[4],gz_hi[4];
};

struct output_stats {
  int font_fixed,vel_limited,failure_checker;
  long n_iter;
};


// FIXME: For cosmetic purposes, we might want to make everything either zero-offset or one-offset, instead of a mixture.
const int kronecker_delta[4][3] = { { 0,0,0 },
                                    { 1,0,0 },
                                    { 0,1,0 },
                                    { 0,0,1 } };

/* PUBLIC FUNCTIONS, USED OUTSIDE GiRaFFE AS WELL */
void GiRaFFE_compute_conservatives(const CCTK_REAL *PRIMS,  const CCTK_REAL *METRIC, CCTK_REAL *CONSERVS);

void GiRaFFE_convert_ADM_to_BSSN__enforce_detgtij_eq_1__and_compute_gtupij
(const cGH *cctkGH,const int *cctk_lsh,
 CCTK_REAL *gxx,CCTK_REAL *gxy,CCTK_REAL *gxz,CCTK_REAL *gyy,CCTK_REAL *gyz,CCTK_REAL *gzz,const CCTK_REAL *alp,
 CCTK_REAL *gtxx,CCTK_REAL *gtxy,CCTK_REAL *gtxz,CCTK_REAL *gtyy,CCTK_REAL *gtyz,CCTK_REAL *gtzz,
 CCTK_REAL *gtupxx,CCTK_REAL *gtupxy,CCTK_REAL *gtupxz,CCTK_REAL *gtupyy,CCTK_REAL *gtupyz,CCTK_REAL *gtupzz,
 CCTK_REAL *phi,CCTK_REAL *psi,CCTK_REAL *lapm1);

void GiRaFFE_set_symmetry_gzs_staggered(const cGH *cctkGH, const int *cctk_lsh,const CCTK_REAL *X,const CCTK_REAL *Y,const CCTK_REAL *Z, CCTK_REAL *gridfunc,
                                        const CCTK_REAL *gridfunc_syms,const int stagger_x,const int stagger_y,const int stagger_z);

#endif // GIRAFFE_HEADERS_H
