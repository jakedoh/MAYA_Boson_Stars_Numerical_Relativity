#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "util_Table.h"

// uses theory in PHYSICAL REVIEW D 67, 024018 2003

#define DIM 3
#define NUM_INPUT_ARRAYS 12
#define NUM_OUTPUT_ARRAYS 12
#define NGHOSTS 2

#define MAX_NUMBER_DETECTORS 100
static CCTK_INT file_created[MAX_NUMBER_DETECTORS];

/* External routines */
void ihspin_init(CCTK_ARGUMENTS);
void ihspin(CCTK_ARGUMENTS);


/* IO */
static int IHSpin_write_output(CCTK_ARGUMENTS, CCTK_INT hn, 
			CCTK_REAL sum_x, CCTK_REAL sum_y, CCTK_REAL sum_z,
			CCTK_REAL sum_px, CCTK_REAL sum_py, CCTK_REAL sum_pz)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  const char *fmode;
  char *filename;
  char varname[1024]; // XXX fixed size
  char file_extension[5]=".asc";
  char format_str_real[2048]; // XXX fixed size
  FILE *file;


  if (verbose>3) {
    CCTK_VInfo(CCTK_THORNSTRING, "writing output");
  }

  if (hn>=MAX_NUMBER_DETECTORS) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "warn: hn=%d, but MAX_NUMBER_DETECTORS=%d, increase",
               hn,MAX_NUMBER_DETECTORS);
  }

  // file mode: append if already written
  fmode = (file_created[hn]>0) ? "a" : "w";

  // filename
  sprintf(varname, "ihspin_hn_%d",hn);

  filename = (char *) malloc (strlen (out_dir) + strlen (varname) +
                              strlen (file_extension) +2);
  assert(filename);
  sprintf (filename, "%s/%s%s", out_dir, varname, file_extension);

  // open file
  file = fopen (filename, fmode);
  if (!file) {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "write_ihspin: Could not open scalar output file '%s'",
                filename);
    return -1;
  }

  // write header on startup
  if (file_created[hn]<=0) {
    fprintf(file,"# IHSpin\n");
    fprintf(file,"# horizon no.=%d\n",hn);
    fprintf(file,"# gnuplot column index:\n");
    fprintf(file,"# 1:t 2:Sx 3:Sy 4:Sz 5:Px 6:Py 7:Pz\n");
  }

  // write data
  sprintf (format_str_real,
           "%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\n", 
           out_format,out_format,out_format,out_format,
           out_format,out_format,out_format);
  fprintf(file, format_str_real, cctk_time, 
	  sum_x,sum_y,sum_z,sum_px,sum_py,sum_pz);

  fclose(file); 
  free(filename);

  if (file_created[hn]==0) {
    file_created[hn]=1;
  }
  return 1;
}





/* fills g11...g33, k11...k33 with the interpolated numbers */
static int get_gab_kab_onto_horizon(CCTK_ARGUMENTS,
                             CCTK_INT sn,
                             CCTK_REAL *g11, CCTK_REAL *g12, CCTK_REAL *g13,
                             CCTK_REAL *g22, CCTK_REAL *g23, CCTK_REAL *g33,
                             CCTK_REAL *k11, CCTK_REAL *k12, CCTK_REAL *k13,
                             CCTK_REAL *k22, CCTK_REAL *k23, CCTK_REAL *k33)
{
  DECLARE_CCTK_ARGUMENTS; 
  DECLARE_CCTK_PARAMETERS; 
  int ierr=-1;
  CCTK_INT ind,ind2d;
  CCTK_REAL th,ph,ct,st,cp,sp;
  CCTK_INT ntheta,nphi;

  ntheta=sf_ntheta[sn]-2*nghoststheta[sn];
  nphi=sf_nphi[sn]-2*nghostsphi[sn];

  if (verbose>1) {
    CCTK_VInfo(CCTK_THORNSTRING,"surface %d (%g,%g,%g) nth,nph (%d,%d)",
	       sn,sf_centroid_x[sn],sf_centroid_y[sn],sf_centroid_z[sn],
	       ntheta,nphi);
  }

  const CCTK_INT npoints=ntheta*nphi;

  CCTK_INT interp_npoints=npoints;
  // uni-processor code - only work on CPU 0
  const CCTK_INT myproc= CCTK_MyProc(cctkGH);
  if ( myproc != 0 ) {
    interp_npoints=0;
  }

  // coordinates setup
  const CCTK_INT imin=nghoststheta[sn], imax=sf_ntheta[sn]-nghoststheta[sn]-1;
  const CCTK_INT jmin=nghostsphi[sn], jmax=sf_nphi[sn]-nghostsphi[sn]-1;
  const CCTK_REAL oth=sf_origin_theta[sn];
  const CCTK_REAL oph=sf_origin_phi[sn];
  const CCTK_REAL dth=sf_delta_theta[sn];
  const CCTK_REAL dph=sf_delta_phi[sn];
  CCTK_REAL ah_x[npoints], ah_y[npoints], ah_z[npoints];
  for (int i=imin,n=0;i<=imax;i++,n++) { // theta in [0.5 delta_th, pi-0.5 delta_th]
    th=oth + i * dth;
    ct=cos(th);
    st=sin(th);
    for (int j=jmin,m=0;j<=jmax;j++,m++) { // phi in [0,2pi-delta_phi]
      ind=i + maxntheta *(j+maxnphi*sn); // XXX not sf_ntheta!
      ph=oph + j * dph;
      cp=cos(ph);
      sp=sin(ph);
      ind2d=n + ntheta*m;
      ah_x[ind2d]=sf_centroid_x[sn]+sf_radius[ind]*cp*st;
      ah_y[ind2d]=sf_centroid_y[sn]+sf_radius[ind]*sp*st;
      ah_z[ind2d]=sf_centroid_z[sn]+sf_radius[ind]*ct;
    }
  }

  const void* interp_coords[3] 
    = { (const void *) ah_x,
        (const void *) ah_y,
        (const void *) ah_z };

  // 3d input arrays
  CCTK_STRING input_array_names[NUM_INPUT_ARRAYS]
    = { "ADMBase::gxx", "ADMBase::gxy", "ADMBase::gxz",
        "ADMBase::gyy", "ADMBase::gyz", "ADMBase::gzz",

        "ADMBase::kxx", "ADMBase::kxy", "ADMBase::kxz",
        "ADMBase::kyy", "ADMBase::kyz", "ADMBase::kzz" };
  CCTK_INT input_array_indices[NUM_INPUT_ARRAYS];
  for(unsigned int i = 0 ; i < sizeof(input_array_indices)/sizeof(input_array_indices[0]) ; i++) {
    input_array_indices[i] = CCTK_VarIndex(input_array_names[i]);
    if(input_array_indices[i] < 0) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
        "couldn't find variable '%s'", input_array_names[i]);
        return -1; /*NOTREACHED*/
    }
  }

  const CCTK_INT output_array_types[NUM_OUTPUT_ARRAYS]
    = { CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL  };

  // 2d output arrays
  void * output_arrays[NUM_OUTPUT_ARRAYS]
    = { (void *) g11, 
        (void *) g12,
        (void *) g13,
        (void *) g22,
        (void *) g23,
        (void *) g33,

        (void *) k11, 
        (void *) k12,
        (void *) k13,
        (void *) k22,
        (void *) k23,
        (void *) k33 };

  const CCTK_INT operand_indices[NUM_OUTPUT_ARRAYS]
    = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };

  const CCTK_INT opcodes[NUM_OUTPUT_ARRAYS]
    = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

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
                        opcodes, "opcodes");
  
  if (param_table_handle < 0)
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "bad interpolator parameter(s) \"%s\"!",
               interpolator_pars);
  
  const int coord_system_handle = CCTK_CoordSystemHandle(coord_system);
  if (coord_system_handle < 0)
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
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
    return -1;
  }
  else {
    return 1;
  }

  ierr = Util_TableDestroy(param_table_handle);

  if (ierr != 0) {
    CCTK_WARN(1,"Could not destroy table");
    return -1;
  }
  else {
    return 1;
  }

}

/* invert 3-metric to get g^ij - see TAT/TGRtensor/src/tensor.F90*/
static int upstairs_metric(CCTK_REAL gg[3][3],CCTK_REAL gu[3][3])
{
  const CCTK_REAL dtg=(gg[0][0]*gg[1][1]-gg[0][1]*gg[0][1])*gg[2][2]
    -gg[0][0]*gg[1][2]*gg[1][2]
    +2.0*gg[0][1]*gg[0][2]*gg[1][2]
    -gg[0][2]*gg[0][2]*gg[1][1];

  if (dtg==0) {
      return -1;
  }
  if (dtg<0) {
      return -2;
  }

  gu[0][0] = (gg[1][1] * gg[2][2] - gg[1][2] *gg[1][2] ) / dtg;
  gu[0][1] = (gg[0][2] * gg[1][2] - gg[0][1] *gg[2][2] ) / dtg;
  gu[0][2] = (gg[0][1] * gg[1][2] - gg[0][2] *gg[1][1] ) / dtg;
  gu[1][1] = (gg[0][0] * gg[2][2] - gg[0][2] *gg[0][2] ) / dtg;
  gu[1][2] = (gg[0][2] * gg[0][1] - gg[1][2] *gg[0][0] ) / dtg;
  gu[2][2] = (gg[0][0] * gg[1][1] - gg[0][1] *gg[0][1] ) / dtg;

  gu[1][0] = gu[0][1];
  gu[2][0] = gu[0][2];
  gu[2][1] = gu[1][2];

#if 0
  CCTK_REAL tmp;
  for (int i=0;i<3;i++) {
    for (int j=0;j<3;j++) {
      tmp=0;
      for (int a=0;a<3;a++) {
        tmp+=gu[i][a]*gg[a][j];
      }
      fprintf(stderr,"test: %d %d: %g\n",i,j,tmp);
    }
  }
#endif

  return 1;
}



static int get_g_k_local(int i, int j, int ntheta,
                  CCTK_REAL *g11_ah,CCTK_REAL *g12_ah,CCTK_REAL *g13_ah,
                  CCTK_REAL *g22_ah,CCTK_REAL *g23_ah,CCTK_REAL *g33_ah,
                  CCTK_REAL *k11_ah,CCTK_REAL *k12_ah,CCTK_REAL *k13_ah,
                  CCTK_REAL *k22_ah,CCTK_REAL *k23_ah,CCTK_REAL *k33_ah,
                  CCTK_REAL gloc[3][3],CCTK_REAL kloc[3][3])
{
  /* gloc_ij = g_ij, kloc_ij=K_ij - downstairs index */
  CCTK_INT ind2d=i + ntheta*j;
  gloc[0][0]=g11_ah[ind2d];
  gloc[0][1]=g12_ah[ind2d];
  gloc[0][2]=g13_ah[ind2d];
  gloc[1][0]=g12_ah[ind2d];
  gloc[1][1]=g22_ah[ind2d];
  gloc[1][2]=g23_ah[ind2d];
  gloc[2][0]=g13_ah[ind2d];
  gloc[2][1]=g23_ah[ind2d];
  gloc[2][2]=g33_ah[ind2d];
  
  kloc[0][0]=k11_ah[ind2d];
  kloc[0][1]=k12_ah[ind2d];
  kloc[0][2]=k13_ah[ind2d];
  kloc[1][0]=k12_ah[ind2d];
  kloc[1][1]=k22_ah[ind2d];
  kloc[1][2]=k23_ah[ind2d];
  kloc[2][0]=k13_ah[ind2d];
  kloc[2][1]=k23_ah[ind2d];
  kloc[2][2]=k33_ah[ind2d];
  return 1;
}

static int drdth_drdph(int i, int j,
                int sn,
                CCTK_REAL dth, CCTK_REAL dph,
                CCTK_INT verbose,
                CCTK_INT maxntheta, CCTK_INT maxnphi,
                CCTK_REAL *sf_radius,
                CCTK_REAL *ht, CCTK_REAL *hp)
{
  CCTK_INT ind;
  CCTK_REAL htp1,htm1,hpp1,hpm1;
  CCTK_REAL htp2,htm2,hpp2,hpm2;

  /* dr/dth dr/dph */
  ind=(i+1) + maxntheta *(j+maxnphi*sn); // XXX not sf_ntheta!
  htp1=sf_radius[ind];
  ind=(i-1) + maxntheta *(j+maxnphi*sn); // XXX not sf_ntheta!
  htm1=sf_radius[ind];
  ind=(i+2) + maxntheta *(j+maxnphi*sn); // XXX not sf_ntheta!
  htp2=sf_radius[ind];
  ind=(i-2) + maxntheta *(j+maxnphi*sn); // XXX not sf_ntheta!
  htm2=sf_radius[ind];
  *ht = (1./12.*htm2-2./3.*htm1+2./3.*htp1-1./12.*htp2)/dth;
  if (verbose>5) {
    fprintf(stderr,"  normal : i=%d j=%d ht=%g\n",i,j,*ht);
  }

  ind=i + maxntheta *((j+1)+maxnphi*sn); // XXX not sf_ntheta!
  hpp1=sf_radius[ind];
  ind=i + maxntheta *((j-1)+maxnphi*sn); // XXX not sf_ntheta!
  hpm1=sf_radius[ind];
  ind=i + maxntheta *((j+2)+maxnphi*sn); // XXX not sf_ntheta!
  hpp2=sf_radius[ind];
  ind=i + maxntheta *((j-2)+maxnphi*sn); // XXX not sf_ntheta!
  hpm2=sf_radius[ind];
  *hp = (1./12.*hpm2-2./3.*hpm1+2./3.*hpp1-1./12.*hpp2)/dph;
  if (verbose>5) { 
    fprintf(stderr,"  normal : i=%d j=%d hp=%g\n",i,j,*hp);            
  }
  return 1;
}


void ihspin_init(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
 
  if (verbose>0) {
    CCTK_INFO("initialize ihspin stuff");
  }

  if (out_offset !=0) {
    CCTK_WARN(1,"out_offset is deprecated - assuming 0");
  }
}

static CCTK_REAL *ihspin_allocate_array(CCTK_INT npoints, char const *name)
{
  DECLARE_CCTK_PARAMETERS;

  if (npoints<=0) {
    CCTK_WARN(0,"can't allocate array with npoints <=0");
  }
  if (name==NULL) {
    CCTK_WARN(0,"give a name");
  }

  CCTK_REAL *res;
  res=(CCTK_REAL *) malloc(sizeof(CCTK_REAL)*npoints);
  if (res==NULL) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "%s allocation for npoints=%d failed",name,npoints);
  }
  for (int i=0;i<npoints;i++) res[i]=0;
  const CCTK_REAL mbyt=npoints*8.0/(1024.0*1024.0);
  if (verbose>0) {
    CCTK_VInfo(CCTK_THORNSTRING,"allocated array %s with %d elements -> %g MB",name,npoints,mbyt);
  }
  return res;
}


static CCTK_INT have_integrand_memory=0;
static CCTK_REAL *g11_ah, *g12_ah, *g13_ah, *g22_ah, *g23_ah, *g33_ah;
static CCTK_REAL *k11_ah, *k12_ah, *k13_ah, *k22_ah, *k23_ah, *k33_ah;
static CCTK_INT ihspin_get_local_memory(CCTK_INT npoints)
{
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose>1) CCTK_INFO("in allocate_memory");

  if (have_integrand_memory==0) {
    if (verbose>0) CCTK_INFO("allocating new memory");
    // metric on horizon
    g11_ah=ihspin_allocate_array(npoints,"g11_ah");
    g12_ah=ihspin_allocate_array(npoints,"g12_ah");
    g13_ah=ihspin_allocate_array(npoints,"g13_ah");
    g22_ah=ihspin_allocate_array(npoints,"g22_ah");
    g23_ah=ihspin_allocate_array(npoints,"g23_ah");
    g33_ah=ihspin_allocate_array(npoints,"g33_ah");
    // extrinsic curvature on horizon
    k11_ah=ihspin_allocate_array(npoints,"k11_ah");
    k12_ah=ihspin_allocate_array(npoints,"k12_ah");
    k13_ah=ihspin_allocate_array(npoints,"k13_ah");
    k22_ah=ihspin_allocate_array(npoints,"k22_ah");
    k23_ah=ihspin_allocate_array(npoints,"k23_ah");
    k33_ah=ihspin_allocate_array(npoints,"k33_ah");
    // update memory allocation flag
    have_integrand_memory=1;
  }
  else {
    if (verbose>1) CCTK_INFO("already allocated memory");
    return 2;
  }

  return 1;
}



void ihspin (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  CCTK_INT sn, ind, ierr;
  CCTK_REAL ht,hp;
  CCTK_REAL sint,sinp,cost,cosp,rp;

  /* local memory allocation */
  CCTK_INT npoints=maxntheta*maxnphi;
  ierr=ihspin_get_local_memory(npoints);
  if (ierr<0) {
    CCTK_WARN(1,"failed to allocate memory");
    return;
  }

  /* loop over horizons */
  for (int hn=0;hn<num_horizons;hn++)
  {

    CCTK_INT my_compute_every;

    /* check parameters and decide if we have to do anythin */
    if ( compute_every_individual[hn] >= 0 ) {
        my_compute_every = compute_every_individual[hn];
    } else {
        my_compute_every = compute_every;
    }
    if ( my_compute_every == 0 || cctk_iteration % my_compute_every != 0 
		|| cctk_time < find_after_individual_time[hn] || cctk_iteration < find_after_individual[hn] ) {
      continue;
    }


    sn=surface_index[hn];
    if (sn>=sphericalsurfaces_nsurfaces) {
      CCTK_VInfo(CCTK_THORNSTRING,
                 "surface number sn=%d too large, increase SphericalSurface::nsurfaces from its current value %d",
                 sn,sphericalsurfaces_nsurfaces);
      continue;
    }
    if (sf_valid[sn]<=0) {
      CCTK_VInfo(CCTK_THORNSTRING,
                 "didn't find valid horizon surface for sn=%d, hn=%d",sn,hn);
      continue;
    }
    if(nghoststheta[sn] < NGHOSTS || nghostsphi[sn] < NGHOSTS) { // we need at least NGHOSTS ghost zones 
      CCTK_VInfo(CCTK_THORNSTRING,
                 "number of ghost zones for spherical surface %d must be at least %d.",sn,NGHOSTS);
      continue;
    }
    const CCTK_REAL pi=acos(-1.0);
    const CCTK_INT imin=nghoststheta[sn], imax=sf_ntheta[sn]-nghoststheta[sn]-1;
    const CCTK_INT jmin=nghostsphi[sn], jmax=sf_nphi[sn]-nghostsphi[sn]-1;
    const CCTK_INT ntheta = imax-imin+1, nphi = jmax-jmin+1;
    const CCTK_REAL oth=sf_origin_theta[sn];
    const CCTK_REAL oph=sf_origin_phi[sn];
    const CCTK_REAL dth=sf_delta_theta[sn];
    const CCTK_REAL dph=sf_delta_phi[sn];
    const CCTK_REAL dtp=dth*dph;

    if (verbose>2) {
      CCTK_VInfo(CCTK_THORNSTRING,"ntheta=%d nphi=%d dth=%g dph=%g",
                 ntheta,nphi,dth,dph);
    }

    ierr=get_gab_kab_onto_horizon(CCTK_PASS_CTOC, sn,
                              g11_ah, g12_ah, g13_ah,
                              g22_ah, g23_ah, g33_ah,
                              k11_ah, k12_ah, k13_ah,
                              k22_ah, k23_ah, k33_ah);
    if (ierr<0) {
      CCTK_WARN(1,"unable to get g_ab and K_ab on the horizon. not doing anything.");
      continue;
    }

    CCTK_REAL phi[3],rvec[3],rdn[3],rup[3];
    CCTK_REAL phi_x[3], phi_y[3], phi_z[3];
    CCTK_REAL gloc[3][3],kloc[3][3],gup[3][3];
    CCTK_REAL th,ph,dA;
    CCTK_REAL sum=0; // the value of the spin integral
    CCTK_REAL sum_x=0,sum_y=0,sum_z=0;
    CCTK_REAL sum_px=0,sum_py=0,sum_pz=0;

    const CCTK_INT myproc= CCTK_MyProc(cctkGH);

    CCTK_REAL iwtheta,iwphi,intweight;
    /* loop over horizon surface */
    for (int i=imin,n=0;i<=imax;i++,n++) // theta in [0.5 delta_th, pi-0.5 delta_th]
    {
      th=oth + i * dth;
      cost=cos(th);
      sint=sin(th);
      // weigths from NR(C++,2) 4.1.14 plus extrapolated 1/2 integral (see Maple worksheet)
      if      (i==imin+0 || i==imax-0) iwtheta=13.0/12.0;
      else if (i==imin+1 || i==imax-1) iwtheta= 7.0/ 8.0;
      else if (i==imin+2 || i==imax-2) iwtheta=25.0/24.0;
      else iwtheta=1;

      for (int j=jmin,m=0;j<=jmax;j++,m++) // phi in [0,2pi-delta_phi]
      {
        ph=oph + j * dph;
        cosp=cos(ph);
        sinp=sin(ph);
	iwphi=1; // trapezoid rule
	intweight=iwphi*iwtheta;

        ind=i + maxntheta *(j+maxnphi*sn); // XXX not sf_ntheta!
        rp=sf_radius[ind];

        get_g_k_local(n,m,ntheta,
                      g11_ah,g12_ah,g13_ah,
                      g22_ah,g23_ah,g33_ah,
                      k11_ah,k12_ah,k13_ah,
                      k22_ah,k23_ah,k33_ah,
                      gloc,kloc);

        /* invert 3-metric to get g^ij */
        ierr=upstairs_metric(gloc,gup);
        if (ierr<0) {
          if (verbose>2) {
            CCTK_VInfo(CCTK_THORNSTRING,"unable to compute inverse metric at th=%g ph=%g",
                       th,ph);
            if (ierr==-1) {
                CCTK_INFO("detg=0!");
            }
            if (ierr==-2) {
                CCTK_INFO("metric negative");
            }
          }
          continue; // skip this point
        }

        /* get derivatives of r in theta and phi direction */
        ierr=drdth_drdph(i, j, 
                         sn,
                         dth,dph,
                         verbose,
                         maxntheta, maxnphi,
                         sf_radius,
                         &ht, &hp);
        if (ierr<0) {
          CCTK_WARN(1,"derivative computation failed");
          continue;
        }

        /* STEP 1 : COMPUTATION OF APPROXIMATE KILLING VECTOR */
        /* dx^a/dr^b */
        CCTK_REAL dxdr=sint*cosp;
        CCTK_REAL dydr=sint*sinp;
        CCTK_REAL dzdr=cost;
        CCTK_REAL dxdp=rp*sint*(-sinp);
        CCTK_REAL dydp=rp*sint*cosp;
        CCTK_REAL dzdp=0;
        CCTK_REAL dxdt=rp*cosp*cost;
        CCTK_REAL dydt=rp*sinp*cost;
        CCTK_REAL dzdt=rp*(-sint);

        /* phi^x - not killing, simply coordinate, but this vector has to be in the surface */
        phi[0]=dxdr*hp+dxdp;
        phi[1]=dydr*hp+dydp;
        phi[2]=dzdr*hp+dzdp;

        /* STEP 2 : COMPUTATION OF SURFACE NORMAL VECTOR */
        /* d_i theta, d_i phi */
        CCTK_REAL ah_x=rp*cosp*sint;
        CCTK_REAL ah_y=rp*sinp*sint;
        CCTK_REAL ah_z=rp*cost;
        /* stay away from ah_z-axis */
        if (ah_x*ah_x+ah_y*ah_y == 0) {
          ah_x=0;
          ah_y=1e-8;
        }

	/* flat space rotational Killing vectors */
	phi_x[0]=0;
	phi_x[1]=-ah_z;
	phi_x[2]=ah_y;

	phi_y[0]=ah_z;
	phi_y[1]=0;
	phi_y[2]=-ah_x;

	phi_z[0]=-ah_y;
	phi_z[1]=ah_x;
	phi_z[2]=0;

        CCTK_REAL dtdx=ah_x*1.0*ah_z*1.0/(sqrt(ah_y*ah_y+ah_x*ah_x)*(ah_z*ah_z+ah_y*ah_y+ah_x*ah_x));
        CCTK_REAL dtdy=ah_y*1.0*ah_z*1.0/(sqrt(ah_y*ah_y+ah_x*ah_x)*(ah_z*ah_z+ah_y*ah_y+ah_x*ah_x));
        CCTK_REAL dtdz=-sqrt(ah_y*ah_y+ah_x*ah_x)*1.0/(ah_z*ah_z+ah_y*ah_y+ah_x*ah_x);
        CCTK_REAL dpdx=-ah_y*1.0/(ah_y*ah_y+ah_x*ah_x);
        CCTK_REAL dpdy=ah_x*1.0/(ah_y*ah_y+ah_x*ah_x);
        CCTK_REAL dpdz=0;
        CCTK_REAL dhdi[3];
        dhdi[0]=ht*dtdx+hp*dpdx;
        dhdi[1]=ht*dtdy+hp*dpdy;
        dhdi[2]=ht*dtdz+hp*dpdz;
        CCTK_REAL ni[3],fi[3];
        /* ni: unit normal vector */
        ni[0]=cosp*sint;
        ni[1]=sinp*sint;
        ni[2]=     cost;
        fi[0]=ni[0];
        fi[1]=ni[1];
        fi[2]=ni[2];
        for (int idir=0;idir<3;idir++) {
          fi[idir]=fi[idir]-dhdi[idir];
        }

        CCTK_REAL metnorm=0.;
        for (int idir=0;idir<3;idir++) {
          for (int jdir=0;jdir<3;jdir++) {
            metnorm = metnorm+fi[idir]*fi[jdir]*gup[idir][jdir];
          }
        }
        metnorm=sqrt(metnorm);
        for (int idir=0;idir<3;idir++) {
          rdn[idir]=fi[idir]/metnorm;
        }
        for (int idir=0;idir<3;idir++) {
          rup[idir]=0;
          for (int jdir=0;jdir<3;jdir++) {
            rup[idir]=rup[idir]+gup[idir][jdir]*rdn[jdir];
          }
        }
        rvec[0]=rup[0];
        rvec[1]=rup[1];
        rvec[2]=rup[2];

        /* check rvec[] and phi[] are orthogonal */
        if (verbose>5) {
          CCTK_REAL tmp=0;
          for (int a=0;a<3;a++) {
            for (int b=0;b<3;b++) {
              tmp+=gloc[a][b]*rvec[a]*phi[b];
            }
          }
          fprintf(stderr,"i=%d j=%d normalization test=%g ==0?\n",i,j,tmp);
        }

        /* STEP 3 : COMPUTATION OF AREA ELEMENT dA */
        /* two metric tdn_ij: g_ab e_i^a e_j^b
           e^a_b=dx^a/dth^b */
        CCTK_REAL ee[3][2];
        ee[0][0] = dxdt;
        ee[0][1] = dxdp;
        ee[1][0] = dydt;
        ee[1][1] = dydp;
        ee[2][0] = dzdt;
        ee[2][1] = dzdp;

        CCTK_REAL tdn[2][2];
        for (int a=0;a<2;a++) {
          for (int b=0;b<2;b++) {
            tdn[a][b]=0;
            for (int c=0;c<3;c++) {
              for (int d=0;d<3;d++) {
                tdn[a][b]+=gloc[c][d]*ee[c][a]*ee[d][b];
              }
            }
          }
        }

        /* determinant */
        CCTK_REAL tdetg;
        tdetg=tdn[0][0]*tdn[1][1]-tdn[1][0]*tdn[0][1];

        if (tdetg>0) {
          dA=sqrt(tdetg)*dtp*intweight; // area element  (weighted for 4th order)
        }
        else {
          dA=0;
        }

        /* trK for linear momentum (not needed for spin) */
        CCTK_REAL trK=0;
        for (int a=0;a<3;a++) {
          for (int b=0;b<3;b++) {
            trK += gup[a][b]*kloc[a][b];
          }
        }

        CCTK_REAL xvec[3],yvec[3],zvec[3];
        xvec[0]=1;xvec[1]=0;xvec[2]=0;
        yvec[0]=0;yvec[1]=1;yvec[2]=0;
        zvec[0]=0;zvec[1]=0;zvec[2]=1;

        /* STEP 4 : SPIN FORMULA */
        for (int a=0;a<3;a++) {
          for (int b=0;b<3;b++) {
	    CCTK_REAL rtmp=rvec[b]*kloc[a][b]*dA;
            sum   += phi[a]  *rtmp;
	    sum_x += phi_x[a]*rtmp;
	    sum_y += phi_y[a]*rtmp;
	    sum_z += phi_z[a]*rtmp;

            CCTK_REAL ptmp=rtmp-trK*gloc[a][b]*rvec[b]*dA;
            sum_px += xvec[a]*ptmp;
            sum_py += yvec[a]*ptmp;
            sum_pz += zvec[a]*ptmp;
          }
        }

        if (verbose>4) {
          fprintf(stderr,"sum=%g [%g,%g,%g]\n",sum,sum_x,sum_y,sum_z);
          fprintf(stderr,"sum=[%g,%g,%g]\n",sum_px,sum_py,sum_pz);
        }
      } // j : phi
    } // i : theta

    /* spin */
    const CCTK_REAL prefac=1.0/(8.0*pi);
    sum  =prefac* sum; // normalization of spin integral
    sum_x=prefac* sum_x; // normalization of spin integral
    sum_y=prefac* sum_y; // normalization of spin integral
    sum_z=prefac* sum_z; // normalization of spin integral

    if (verbose>0) {
      CCTK_VInfo(CCTK_THORNSTRING,"spin value=%g S^i=[%g,%g,%g] on horizon %d",
	         sum,sum_x,sum_y,sum_z,hn);
    }

    ihspin_spin[hn]=sum;
    ihspin_spin_x[hn]=sum_x;
    ihspin_spin_y[hn]=sum_y;
    ihspin_spin_z[hn]=sum_z;

    /* linear momentum */
    sum_px=prefac* sum_px; // normalization of momentum integral
    sum_py=prefac* sum_py; // normalization of momentum integral
    sum_pz=prefac* sum_pz; // normalization of momentum integral

    if (verbose>0) {
      CCTK_VInfo(CCTK_THORNSTRING,"P^i=[%g,%g,%g] on horizon %d",
	         sum_px,sum_py,sum_pz,hn);
    }

    ihspin_mom_x[hn]=sum_px;
    ihspin_mom_y[hn]=sum_py;
    ihspin_mom_z[hn]=sum_pz;

    /* IO on CPU 0 */
    if (myproc == 0) {
      ierr=IHSpin_write_output(CCTK_PASS_CTOC,hn, sum_x,sum_y,sum_z,
                                                  sum_px,sum_py,sum_pz);
      if (ierr<0) {
	CCTK_WARN(1,"writing of information to files failed");
      }
    }
  } // hn loop over horizon number
}
