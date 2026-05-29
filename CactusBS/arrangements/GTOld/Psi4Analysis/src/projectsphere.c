#include "math.h"
#include "stdio.h"

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "util_Table.h"


#define NUM_INPUT_ARRAYS 2
#define NUM_OUTPUT_ARRAYS 2
#define DIM 3

void Psi4Analysis_ProjectSphere(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  CCTK_REAL th=0,ph=0,ct=0,st=0,cp=0,sp=0;
  CCTK_INT ind2d=0;

  if (verbose>2) {
    CCTK_INFO("projecting psi4");
  }

  if (verbose>2) {
    CCTK_VInfo(CCTK_THORNSTRING,"current_detector=%d",*current_detector);
  }

  if (*do_nothing==1) {
    return;
  }

  // current detector radius
  CCTK_REAL rad=detector_radius[*current_detector];

  // XXX we actually want the areal radius here
  *rdet=rad;
  if (rad<=0) {
    if (verbose>1) {
      CCTK_VInfo(CCTK_THORNSTRING,"detector out of range: rad=%g detector number: %d\n",
		 rad,*current_detector);
    }

    return;
  }

  CCTK_INT ierr=0;

  // setup interpolation coords
  const CCTK_REAL pi=acos(-1.0);
  const CCTK_REAL dth=pi/ (*ntheta-1);
  const CCTK_REAL dph=2.0*pi/(*nphi-1);
  const CCTK_INT npoints=maxntheta*maxnphi; // XXX would be cool if we could use ntheta and nphi itself!

  CCTK_INT interp_npoints=npoints;

  // uni-processor code - only work on CPU 0
  const CCTK_INT myproc= CCTK_MyProc(cctkGH);
  if ( myproc != 0 ) {
    interp_npoints=0;
  }

  // origin setup
  CCTK_REAL origin_x=0,origin_y=0,origin_z=0;
  const CCTK_INT sno=track_origin_surface_number[*current_detector];
  const CCTK_REAL cox=constant_origin_x[*current_detector];
  const CCTK_REAL coy=constant_origin_y[*current_detector];
  const CCTK_REAL coz=constant_origin_z[*current_detector];
  const CCTK_REAL eps=1e-15;
  if (track_origin[*current_detector]) {
    if (fabs(cox)>eps && fabs(coy)>eps && fabs(coz)>eps) {
      origin_x=cox;
      origin_y=coy;
      origin_z=coz;
      if (verbose>1) {
	CCTK_VInfo(CCTK_THORNSTRING,"Setting origin for detector %d to fixed value (%g,%g,%g)",
		   *current_detector, cox,coy,coz);
      }
    }
    // spherical surface gives center of mass location
    else {
      if (sno<0) {
      CCTK_VInfo(CCTK_THORNSTRING,"invalid track_origin_surface_number[%d]=%d",
		 *current_detector, sno);
      }
      else if (!sf_valid[sno]) {
	CCTK_VInfo(CCTK_THORNSTRING,"invalid surface: detector=%d sno=%d",
		   *current_detector, sno);
      }
      else {
	origin_x=sf_centroid_x[sno];
	origin_y=sf_centroid_y[sno];
	origin_z=sf_centroid_z[sno];
	if (verbose>1) {
	  CCTK_VInfo(CCTK_THORNSTRING,
		     "change origin for detector %d from surface %d to (%g,%g,%g)",
		     *current_detector, sno, origin_x,origin_y,origin_z);
	}
      }
    }
  }

  // coordinates setup
  CCTK_REAL p4_x[npoints], p4_y[npoints], p4_z[npoints];
  for (int i=0;i<npoints;i++) {
    p4_x[i]=0;p4_y[i]=0;p4_z[i]=0;
  }

  if (x_axis_up[*current_detector] && verbose>1) {
    CCTK_VInfo(CCTK_THORNSTRING,"x axis points up for detector %d",*current_detector);
  }

  if (y_axis_up[*current_detector] && verbose>1) {
    CCTK_VInfo(CCTK_THORNSTRING,"y axis points up for detector %d",*current_detector);
  }

  for (int i=0;i<*ntheta;i++) {
    th=i*dth;
    ct=cos(th);
    st=sin(th);
    for (int j=0;j<*nphi;j++) {
      ph=j*dph;
      cp=cos(ph);
      sp=sin(ph);
      ind2d=i + maxntheta*j;
      p4_x[ind2d]=origin_x+ rad*cp*st;
      p4_y[ind2d]=origin_y+ rad*sp*st;
      p4_z[ind2d]=origin_z+ rad*   ct;
      // allow different axis directions
      if (x_axis_up[*current_detector]) {
	p4_x[ind2d]=origin_x+  rad*   ct;
	p4_y[ind2d]=origin_y+  rad*sp*st;
	p4_z[ind2d]=origin_z+ -rad*cp*st;
      }
      if (y_axis_up[*current_detector]) {
	p4_x[ind2d]=origin_x+  rad*cp*st;
	p4_y[ind2d]=origin_y+  rad*   ct;
	p4_z[ind2d]=origin_z+ -rad*sp*st;
      }
    }
  }

  if (*ntheta<maxntheta) {
    for (int i=*ntheta;i<maxntheta;i++) {
      for (int j=0;j<maxnphi;j++) {
	ind2d=i + maxntheta*j;
	p4_x[ind2d]=rad;
	p4_y[ind2d]=0;
	p4_z[ind2d]=0;
      }
    }
  }

  if (*nphi<maxnphi) {
    for (int i=0;i<maxntheta;i++) {
      for (int j=*nphi;j<maxnphi;j++) {
	ind2d=i + maxntheta*j;
	p4_x[ind2d]=rad;
	p4_y[ind2d]=0;
	p4_z[ind2d]=0;
      }
    }
  }


  const void* interp_coords[3] 
    = { (const void *) p4_x,
        (const void *) p4_y,
        (const void *) p4_z };


  // 3d input arrays
  const CCTK_INT input_array_indices[NUM_INPUT_ARRAYS]
    = { CCTK_VarIndex("WeylScal4::Psi4r"),
        CCTK_VarIndex("WeylScal4::Psi4i") };

  const CCTK_INT output_array_types[NUM_OUTPUT_ARRAYS]
    = { CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL };

  // 2d output arrays
  void * output_arrays[NUM_OUTPUT_ARRAYS]
    = { (void *) psi4_2d_re,
        (void *) psi4_2d_im };

  const CCTK_INT operand_indices[NUM_OUTPUT_ARRAYS]
    = { 0, 1 };

  const CCTK_INT opcodes[NUM_OUTPUT_ARRAYS]
    = { 0, 0 };

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
