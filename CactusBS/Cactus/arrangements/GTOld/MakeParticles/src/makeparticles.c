
#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <string.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#define KRANC_C
#include "GenericFD.h"

#define SQR(x) ((x)*(x))
#define DIM(x) (sizeof(x)/sizeof((x)[0]))

/* does the particle sit on an edge? If so, which one? */
#define IS_LEFT_EDGE 0x01   /* x direction */
#define IS_RIGHT_EDGE 0x02
#define IS_FRONT_EDGE 0x04  /* y direction */
#define IS_BACK_EDGE 0x08
#define IS_TOP_EDGE 0x10    /* z direction */
#define IS_BOTTOM_EDGE 0x20

void MakeParticles(CCTK_ARGUMENTS);
void MakeParticles_ParamCheck(CCTK_ARGUMENTS);

/* information type encodings */
#include "makeparticles_types.h"

/* from GenericFD: obtain range of indices to cover the interior domain (plus
 * the outer physical boundary) */
static void get_interior_domain_specs(CCTK_ARGUMENTS, CCTK_INT imin[3], CCTK_INT imax[3])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /* Summation limits */
  CCTK_INT is_symbnd[6], is_physbnd[6], is_ipbnd[6];
  CCTK_INT nboundaryzones[6], is_internal[6], is_staggered[6], shiftout[6];
  CCTK_INT npoints;

  GenericFD_GetBoundaryInfo(cctkGH, cctk_ash, cctk_lsh, cctk_bbox,
                            cctk_nghostzones,
                            imin, imax, is_symbnd, is_physbnd, is_ipbnd);

  if (CCTK_IsFunctionAliased ("GetBoundarySpecification")) {
    int ierr = GetBoundarySpecification(6, nboundaryzones, is_internal, is_staggered, shiftout);
    if (ierr != 0)
      CCTK_WARN(0, "Could not obtain boundary specification");
  } else {
    /* assume that nboundaryzones = nghostzones
     * leave is_internal, is_staggered, shiftout undefined
     */
    for (int dir = 0; dir<3; dir++)
    {
      for (int face = 0; face < 2; face++)
      {
         CCTK_INT indx = dir*2 + face;
         nboundaryzones[indx] = cctk_nghostzones[dir];
      }
    }
  }

  for (int dir = 0; dir<3; dir++)
  {
    for (int face = 0; face < 2; face++)
    {

       CCTK_INT indx = dir*2 + face;
       if (is_physbnd[indx]) {
          npoints=0;
       } else if (is_symbnd[indx]) {
          npoints=nboundaryzones[indx];
       } else {
          npoints=cctk_nghostzones[dir];
       }

       switch(face)
       {
       case 0: /*lower*/
          imin[dir]=npoints;
          break;
       case 1: /*upper*/
          imax[dir]=cctk_lsh[dir] - npoints;
          break;
       default:
          CCTK_WARN(0,"internal error");
          return; /* NOTREACHED */
       }

    }
  }

  if (verbose >= 3) { 
        CCTK_VInfo(CCTK_THORNSTRING,"Extent: (%d,%d,%d) to (%d,%d,%d)",imin[0],imin[1],imin[2],imax[0],imax[1],imax[2]);
  }
}

/* output a set of nvals CCTK_REALs to fh. If binary!=0 then output them as
 * binaries, otherwise as ASCII numbers. */
static void output(FILE *fh, int binary, int nvals, const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);

  if ( binary ) {
    CCTK_REAL vals[128];
    size_t written;

    assert(nvals <= DIM(vals));
    for(int i = 0 ; i < nvals ; i++)
      vals[i] = va_arg(ap, CCTK_REAL);

    written = fwrite(vals, sizeof(vals[0]), nvals, fh);
    if ( written != nvals ) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING, "Error writing binary data: %s", strerror(errno));
      return; /* NOTREACHED */
    }
  } else {
    int written;

    written = vfprintf(fh, fmt, ap);
    if ( written < 0 ) {
      CCTK_VWarn(0,  __LINE__, __FILE__, CCTK_THORNSTRING, "Error writing ASCII data: %s", strerror(errno));
      return; /* NOTREACHED */
    }
  }

  va_end(ap);
}

static FILE *open(const cGH *cctkGH, CCTK_INT cctk_iteration, const char *hdr)
{
  DECLARE_CCTK_PARAMETERS;

  char fn[256];
  FILE *fh;
  CCTK_INT myproc, nprocs;
  int sz;
  static CCTK_INT last_iteration_seen = -1; // write header into new file?

  // build filename for output
  myproc = CCTK_MyProc(cctkGH);
  nprocs = CCTK_nProcs(cctkGH);
  sz = snprintf(fn, DIM(fn), "%s/%s/particles.it_%d", out_dir, 
                particle_subdir, cctk_iteration);
  assert(sz < DIM(fn));
  if ( nprocs > 1 )  {
    sz += snprintf(fn+sz, DIM(fn)-sz, ".file_%d", myproc);
    assert(sz < DIM(fn));
  }
  sz += snprintf(fn+sz, DIM(fn)-sz, ".%s", binary_output ? "dat" : "asc");
  assert(sz < DIM(fn));

  // open and write header
  if ( cctk_iteration > last_iteration_seen )  {
    fh = fopen(fn, binary_output ? "wb" : "w");
    assert(fh);
    if ( !binary_output )
      fputs(hdr, fh);
    if ( verbose > 4 )
      CCTK_VInfo(CCTK_THORNSTRING, "Opened file '%s' in write mode", fn);
  } else {
    fh = fopen(fn, binary_output ? "ab" : "a");
    assert(fh);
    if ( verbose > 4 )
      CCTK_VInfo(CCTK_THORNSTRING, "Opened file '%s' in append mode", fn);
  }
  if ( verbose > 5 )
    CCTK_VInfo(CCTK_THORNSTRING, "it: %d", cctk_iteration);

  last_iteration_seen = cctk_iteration;

  return fh;
}

void MakeParticles(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  FILE *fh;                                 // output file and output filename
  CCTK_INT imin[3], imax[3];                // domain to loop over
  CCTK_REAL *weight;                        // weight to use to avoid AH and refined points
  int num_negative_mass_points;             // just that
  enum info_type_t info_type;               // parameter encoding
  CCTK_INT reflevel;                        // refinement level a point lives on

  if ( calc_every == 0 || cctk_iteration % calc_every != 0 )
    return;

  if ( verbose > 0 )
     CCTK_INFO("Entering makeparticles");

  /* parse parameters into numbers */
  info_type = info_invalid;
  for (int t = 0 ; t < DIM(info) ; t++) {
    if ( CCTK_EQUALS(information, info[t].name) ) {
      info_type = info[t].type;
      break;
    }
  }
  assert(info_type != info_invalid);

  /* open output file, write header */
  fh = open(cctkGH, cctk_iteration, info[info_type].header);

  reflevel = GetRefinementLevel(cctkGH);

  /* loop over grid and output particles */
  get_interior_domain_specs(CCTK_PASS_CTOC, imin, imax);
  weight = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH,0,weight_variable));
  assert(weight);
  num_negative_mass_points = 0;
#pragma omp parallel for schedule(dynamic)
  for (int k = imin[2] ; k < imax[2] ; k++) {
    for (int j = imin[1] ; j < imax[1] ; j++) {
      for (int i = imin[0] ; i < imax[0] ; i++) {
        CCTK_REAL rhoL,epsL,pressL,vxL,vyL,vzL,w2L, wL;
        CCTK_REAL gxxL,gyyL,gzzL,gxyL,gxzL,gyzL,detg;
        CCTK_REAL vlowxL,vlowyL,vlowzL;
        CCTK_REAL vcosobsY, vcosobsZ;
        CCTK_REAL Df4Y, Df4Z;
        CCTK_REAL alphaL, betaxL, betayL, betazL;
        CCTK_REAL weightL,MassL,EnergyL;
        CCTK_REAL dx,dy,dz,dV,wdV;
        CCTK_REAL xL,yL,zL;
        CCTK_INT is_on_edge; /* bitmask describing which edge we site on, see IS_*_EDGE above */

        /* Information from the grid */
        int idx = CCTK_GFINDEX3D(cctkGH, i,j,k);
        
        /* anything that is too low is treated as vacuum and does not 
         * contribute to the mass */
        if ( rho[idx] <= rho_min || weight[idx] == 0.0 )
          continue;
        
        rhoL = rho[idx];
        epsL = eps[idx];
        pressL = press[idx];
        
        xL = x[idx];
        yL = y[idx];
        zL = z[idx];
        
        gxxL = gxx[idx];
        gxyL = gxy[idx];
        gxzL = gxz[idx];
        gyyL = gyy[idx];
        gyzL = gyz[idx];
        gzzL = gzz[idx];
        
        vxL = velx[idx];
        vyL = vely[idx];
        vzL = velz[idx];

        alphaL = alp[idx];
        betaxL = betax[idx];
        betayL = betay[idx];
        betazL = betaz[idx];

        weightL = weight[idx];

        vlowxL = gxxL*vxL + gxyL*vyL + gxzL*vzL;
        vlowyL = gxyL*vxL + gyyL*vyL + gyzL*vzL;
        vlowzL = gxzL*vxL + gyzL*vyL + gzzL*vzL;

        w2L = 1. / ( 1. - (vlowxL*vxL + vlowyL*vyL + vlowzL*vzL) );
        wL = sqrt(w2L);

        /* Get volume element */
        dx = CCTK_DELTA_SPACE(0);
        dy = CCTK_DELTA_SPACE(1);
        dz = CCTK_DELTA_SPACE(2);
        dV = dx*dy*dz;
        wdV = weightL*dV;
        
        detg = 2*gxyL*gxzL*gyzL + gzzL*(gxxL*gyyL - SQR(gxyL)) - gyyL*SQR(gxzL)
               - gxxL*SQR(gyzL);
                      
        // compute the actual payload
        MassL = sqrt(detg) * wL * rhoL * dV * weightL;
        EnergyL = wL*(alphaL - (betaxL*vlowxL+betayL*vlowyL+betazL*vlowzL));

        /* Doppler factor, observer at -\hat{y} and \hat{z} */
        vcosobsY = -vyL;
        vcosobsZ = vzL;
        Df4Y = SQR(1./(w2L * SQR(1. - vcosobsY)));
        Df4Z = SQR(1./(w2L * SQR(1. - vcosobsZ)));

        // try to find out if we sit on an edge and on which one
        // this assumes that at refinement and symmetry boundaries there is
        // only one layer on the coarse/fine levels with non-unity weight
        is_on_edge = 0x0;
        if ( weightL < 1.0 ) {
          // the logic goes as follows:
          // For each cell with non-unit weight I check if I am either at the
          // left/right edge of a patch [as seen in imin/imax] or my neighbour
          // to the left/right has weight zero. If either of those is true then
          // I am a left/right edge.
          // This should hopefully catch both symmetry boundaries (via edges)
          // and fine-level edges (via edges) and coarse-level edges (via the
          // weight).

          if ( i == imin[0] || weight[CCTK_GFINDEX3D(cctkGH, i-1,j,k)] == 0.0 ) 
            is_on_edge |= IS_LEFT_EDGE;
          else if ( i == imax[0] && weight[CCTK_GFINDEX3D(cctkGH, i+1,j,k)] == 0.0 ) 
            is_on_edge |= IS_RIGHT_EDGE;
          
          if ( j == imin[1] || weight[CCTK_GFINDEX3D(cctkGH, i,j-1,k)] == 0.0 ) 
            is_on_edge |= IS_FRONT_EDGE;
          else if ( j == imax[1] && weight[CCTK_GFINDEX3D(cctkGH, i,j+1,k)] == 0.0 ) 
            is_on_edge |= IS_BACK_EDGE;

          if ( k == imin[2] || weight[CCTK_GFINDEX3D(cctkGH, i,j,k-1)] == 0.0 ) 
            is_on_edge |= IS_TOP_EDGE;
          else if ( k == imax[2] && weight[CCTK_GFINDEX3D(cctkGH, i,j,k+1)] == 0.0 ) 
            is_on_edge |= IS_BOTTOM_EDGE;
        }

        if ( MassL < 0 ) {
#pragma omp atomic
          num_negative_mass_points += 1;
          if (verbose >= 4) {
#pragma omp critical
            CCTK_VInfo(CCTK_THORNSTRING,
                       "Negative MassL found in MakeParticles at (%g,%g,%g).  MassL=%g.", 
                       xL, yL, zL, MassL);
          }
        }

        if ( MassL > 0 ) {
          switch ( info_type ) {
            case info_fallback:
              output(fh, binary_output, info[info_fallback].nvals, 
                     info[info_fallback].fmt, xL, yL, zL, 
                     vxL, vyL, vzL, gxxL, gxyL, gxzL, gyyL, gyzL, gzzL,
                     alphaL, betaxL, betayL, betazL, MassL, EnergyL, rhoL, epsL);
              break;
            case info_luminosity: //FIX: Currently assumes dx=dy=dz to minimize output
              // TODO: include is_on_edge in output?
              output(fh, binary_output, info[info_luminosity].nvals, 
                     info[info_luminosity].fmt, xL, yL, zL, 
                     rhoL, epsL, MassL, detg, wL, Df4Y, Df4Z, dx, wdV,
                     vxL, vyL, vzL);
              break;
            case info_sph:
              if(weightL == 1.0) { /* the SPH code does not like two particles at the same location */
                output(fh, binary_output, info[info_fallback].nvals, 
                       info[info_sph].fmt, xL, yL, zL, 
                       vxL, vyL, vzL, rhoL, epsL, MassL);
              }
            case info_everything:
              output(fh, binary_output, info[info_everything].nvals, 
                     info[info_everything].fmt, xL, yL, zL, 
                     vxL, vyL, vzL, gxxL, gxyL, gxzL, gyyL, gyzL, gzzL,
                     alphaL, betaxL, betayL, betazL, rhoL, epsL, pressL,
                     (double)is_on_edge, weightL, dx, (double)reflevel);
              break;
              break;
            case info_invalid:
            default:
              assert(0 == "Internal error");
              break;
          } /* switch */
        } /* MassL > 0 */
      } /* i */
    } /* j */
  } /* k */

  if ( num_negative_mass_points > 0 ) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Local mass integrand in MassDistribution is negative in %d points.  This implies there's something wrong.", 
               num_negative_mass_points);
  }

  fclose(fh);
  if(verbose > 4)
    CCTK_VInfo(CCTK_THORNSTRING, "Closed file.");

} 

void MakeParticles_ParamCheck(CCTK_ARGUMENTS) 
{
  DECLARE_CCTK_PARAMETERS;

  // Create particle_subdir if it doesn't exist.
  char fulldir[256];
  int sz;
  sz = snprintf(fulldir, DIM(fulldir), "%s/%s",out_dir,particle_subdir);
  assert(sz < DIM(fulldir));
  int ierr=CCTK_CreateDirectory(0777,fulldir);
  if ( ierr < 0 ) {
     CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
                "Error in creating directory %s. Error code %d.",
                fulldir,ierr);
  }

}
