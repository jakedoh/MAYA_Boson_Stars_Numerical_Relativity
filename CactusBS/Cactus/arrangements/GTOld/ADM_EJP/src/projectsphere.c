#include "math.h"
#include "stdio.h"

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "util_Table.h"

#include "sphereintegrate.h"

#define NUM_INPUT_ARRAYS 16
#define NUM_OUTPUT_ARRAYS 43
#define DIM 3

void ADM_EJP_ProjectSphere(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  CCTK_REAL th=0,ph=0,ct=0,st=0,cp=0,sp=0;
  CCTK_INT ierr=0;

  if (verbose>2) {
    CCTK_INFO("projecting metric and extrinsic curvature");
  }

  if (verbose>2) {
    CCTK_VInfo(CCTK_THORNSTRING,"current_detector=%d",*current_detector);
  }

  if (*do_nothing==1) {
    return;
  }

  // current detector radius
  const CCTK_REAL rad=detector_radius[*current_detector];

  // XXX we actually want the areal radius here
  *rdet=rad;
  if (rad<=0) {
    if (verbose>1) {
      CCTK_VInfo(CCTK_THORNSTRING,"detector out of range: rad=%g\n",rad);
    }

    return;
  }

  // setup interpolation coords
  const CCTK_REAL pi=acos(-1.0);
  const CCTK_REAL dth=pi/(ntheta-1);
  const CCTK_REAL dph=2.0*pi/(nphi-1);
  const CCTK_INT npoints=ntheta*nphi;

  CCTK_INT interp_npoints=npoints;

  // uni-processor code - only work on CPU 0
  const CCTK_INT myproc= CCTK_MyProc(cctkGH);
  if ( myproc != 0 ) {
    interp_npoints=0;
  }

  // coordinates setup
  static CCTK_REAL *p4_x, *p4_y, *p4_z;
  static CCTK_INT have_memory=0;
  if (have_memory==0) {
    p4_x=(CCTK_REAL *) malloc(npoints*sizeof(CCTK_REAL));
    if (p4_x==NULL) {
      CCTK_INFO("failed to allocate p4_x");
      return;
    }
    p4_y=(CCTK_REAL *) malloc(npoints*sizeof(CCTK_REAL));
    if (p4_y==NULL) {
      CCTK_INFO("failed to allocate p4_y");
      return;
    }
    p4_z=(CCTK_REAL *) malloc(npoints*sizeof(CCTK_REAL));
    if (p4_z==NULL) {
      CCTK_INFO("failed to allocate p4_z");
      return;
    }

    have_memory=1;
  }
  ierr=SphereIntegrate_setup_coords_open_theta(ntheta, nphi, rad, 0., 0., 0.,
					       p4_x, p4_y, p4_z);
  if (ierr<0) {
    CCTK_WARN(1,"failed to setup coords");
  }

  const void* interp_coords[3] 
    = { (const void *) p4_x,
        (const void *) p4_y,
        (const void *) p4_z };


  // 3d input arrays
  const CCTK_INT input_array_indices[NUM_INPUT_ARRAYS]
    = { CCTK_VarIndex("ADMBase::gxx"),
        CCTK_VarIndex("ADMBase::gxy"),
        CCTK_VarIndex("ADMBase::gxz"),
        CCTK_VarIndex("ADMBase::gyy"),
        CCTK_VarIndex("ADMBase::gyz"),
        CCTK_VarIndex("ADMBase::gzz"),
        CCTK_VarIndex("ADMBase::kxx"),
        CCTK_VarIndex("ADMBase::kxy"),
        CCTK_VarIndex("ADMBase::kxz"),
        CCTK_VarIndex("ADMBase::kyy"),
        CCTK_VarIndex("ADMBase::kyz"),
        CCTK_VarIndex("ADMBase::kzz"),
        CCTK_VarIndex("ADMBase::alp"),
        CCTK_VarIndex("ADMBase::betax"),
        CCTK_VarIndex("ADMBase::betay"),
        CCTK_VarIndex("ADMBase::betaz")  };

  const CCTK_INT output_array_types[NUM_OUTPUT_ARRAYS]
    = { CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL
 };

  // 2d output arrays
  void * output_arrays[NUM_OUTPUT_ARRAYS]
    = { (void *) gxx2d, 
        (void *) dx_gxx,
        (void *) dy_gxx,
        (void *) dz_gxx,
        (void *) gxy2d, 
        (void *) dx_gxy,
        (void *) dy_gxy,
        (void *) dz_gxy,
        (void *) gxz2d, 
        (void *) dx_gxz,
        (void *) dy_gxz,
        (void *) dz_gxz,
        (void *) gyy2d, 
        (void *) dx_gyy,
        (void *) dy_gyy,
        (void *) dz_gyy,
        (void *) gyz2d, 
        (void *) dx_gyz,
        (void *) dy_gyz,
        (void *) dz_gyz,
        (void *) gzz2d, 
        (void *) dx_gzz,
        (void *) dy_gzz,
        (void *) dz_gzz,
        (void *) kxx2d,
        (void *) kxy2d,
        (void *) kxz2d,
        (void *) kyy2d,
        (void *) kyz2d,
        (void *) kzz2d,
        (void *) alpha2d,
        (void *) betax2d,
        (void *) dx_betax2d,
        (void *) dy_betax2d,
        (void *) dz_betax2d,
        (void *) betay2d,
        (void *) dx_betay2d,
        (void *) dy_betay2d,
        (void *) dz_betay2d,
        (void *) betaz2d,
        (void *) dx_betaz2d,
        (void *) dy_betaz2d,
        (void *) dz_betaz2d };

  const CCTK_INT operand_indices[NUM_OUTPUT_ARRAYS]
    = { 0, 0, 0, 0, //gxx
        1, 1, 1, 1, //gxy
        2, 2, 2, 2, //gxz
        3, 3, 3, 3, //gyy
        4, 4, 4, 4, //gyz
        5, 5, 5, 5, //gzz
        6,          //kxx
        7,          //kxy
        8,          //kxz
        9,          //kyy
        10,         //kyz
        11,         //kzz
        12,         //alpha
        13, 13, 13, 13, //betax
        14, 14, 14, 14, //betay
        15, 15, 15, 15  //betaz
  };

  // derivative definitions
#define DERIV(x) x
  const CCTK_INT opcodes[NUM_OUTPUT_ARRAYS]
    = { 0, DERIV(1), DERIV(2), DERIV(3), //gxx
        0, DERIV(1), DERIV(2), DERIV(3), //gxy
        0, DERIV(1), DERIV(2), DERIV(3), //gxz
        0, DERIV(1), DERIV(2), DERIV(3), //gyy
        0, DERIV(1), DERIV(2), DERIV(3), //gyz
        0, DERIV(1), DERIV(2), DERIV(3), //gzz
        0, //kxx
        0, //kxy
        0, //kxz
        0, //kyy
        0, //kyz
        0, //kzz
        0, //alpha
        0, DERIV(1), DERIV(2), DERIV(3), //betax
        0, DERIV(1), DERIV(2), DERIV(3), //betay
        0, DERIV(1), DERIV(2), DERIV(3)  //betaz
  };

  // handles setup
  const int operator_handle = CCTK_InterpHandle(interpolator_name);
  if (operator_handle < 0)
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "couldn't find interpolator \"%s\"!",
               interpolator_name);

  int param_table_handle = Util_TableCreateFromString(interpolator_pars);

  Util_TableSetIntArray(param_table_handle, NUM_OUTPUT_ARRAYS,
                        operand_indices, "operand_indices");

  Util_TableSetIntArray(param_table_handle, NUM_OUTPUT_ARRAYS, 
                        opcodes, "operation_codes");

  if (param_table_handle < 0)
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "bad interpolator parameter(s) \"%s\"!",
               interpolator_pars);

  const int coord_system_handle = CCTK_CoordSystemHandle(coord_system);
  if (coord_system_handle < 0)
    CCTK_VWarn(-1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "can't get coordinate system handle for coordinate system \"%s\"!",
               coord_system);

  // actual interpolation call
  ierr = CCTK_InterpGridArrays(cctkGH,
                               DIM, // number of dimensions 
                               operator_handle,
                               param_table_handle,
                               coord_system_handle,
                               interp_npoints,
                               CCTK_VARIABLE_REAL,
                               interp_coords,
                               NUM_INPUT_ARRAYS, // Number of input arrays
                               input_array_indices,
                               NUM_OUTPUT_ARRAYS, // Number of output arrays
                               output_array_types,
                               output_arrays);

  if (ierr<0) {
    CCTK_WARN(1,"interpolation screwed up");
    *rdet=-1; // flag detector as out of range
  }

}
