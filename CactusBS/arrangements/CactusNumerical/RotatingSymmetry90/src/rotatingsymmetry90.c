#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <util_ErrorCodes.h>
#include <util_Table.h>

#include <Slab.h>

#include "rotatingsymmetry90.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>



/* This is pretty hard coded into all the tensor types and cannot be
   changed easily.  */
#define DIM 3 /* spatial dimension */



static int convert_index (int const step,
                          int const index,
                          int const * restrict const alldirs,
                          int * restrict const parity)
{
  int srcindex;
  assert (index>=0 && index<DIM);
  assert (alldirs);
  assert (alldirs[0]>=0 && alldirs[0]<DIM);
  assert (alldirs[1]>=0 && alldirs[1]<DIM);
  assert (alldirs[1] != alldirs[0]);
  assert (parity);
  assert (abs(*parity) == 1);
  switch (step) {
  case 0:
    /* x face: counterclockwise, 90 degrees */
    if      (index == alldirs[0]) srcindex = alldirs[1], *parity *= -1;
    else if (index == alldirs[1]) srcindex = alldirs[0], *parity *= +1;
    else                          srcindex = index     , *parity *= +1;
    break;
  case 1:
    /* y face: clockwise, 90 degrees */
    if      (index == alldirs[0]) srcindex = alldirs[1], *parity *= +1;
    else if (index == alldirs[1]) srcindex = alldirs[0], *parity *= -1;
    else                          srcindex = index     , *parity *= +1;
    break;
  case 2:
    /* xy edge: 180 degrees */
    if      (index == alldirs[0]) srcindex = alldirs[0], *parity *= -1;
    else if (index == alldirs[1]) srcindex = alldirs[1], *parity *= -1;
    else                          srcindex = index     , *parity *= +1;
    break;
  default:
    assert (0);
  }
  return srcindex;
}

/* bbox, lbnd and ubnd of combined level box */
static int global_bbox[2*DIM];
static int global_lbnd[DIM], global_ubnd[DIM];
static int extent_valid_for_iteration = -1;

int BndRot90VI (cGH const * restrict const cctkGH,
                int const nvars,
                CCTK_INT const * restrict const vis)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  int * restrict gis;
  cGroup group;
  cGroupDynamicData data;
  char * restrict fullname;
  void * restrict * restrict varptrs;
  int * restrict vartypes;
  int * restrict vectorlengths;
  int firstvar;
  
  char tensortypealias[100];
  int numvars;
  int index;
  int srcvi;
  void const * restrict * restrict srcptrs;
  int * restrict parities;
  
  int fake_bbox[2*DIM];
  
  CCTK_REAL x0[DIM], dx[DIM];
  CCTK_REAL origin[DIM], dorigin[DIM];
  int avoid_origin[DIM];
  int offset[DIM];              /* offset 0..1 due to avoid_origin */
  
  struct xferinfo * restrict xferinfo;
  int options;
  
  int have_global_bbox, have_local_bbox;
  
  int var;
  
  int alldirs[2];
  int step;
  int ndirs;                    /* 1 for face, 2 for edge */
  int dir[DIM];                 /* current face or edge */
  int otherdir[DIM];            /* the other direction of the rotation */
  int q;                        /* 0..ndirs-1 */
  int d;                        /* 0..group.dim-1 */
  int ierr;
  
  /* Check arguments */
  assert (cctkGH);
  assert (nvars==0 || vis);
  for (var=0; var<nvars; ++var) {
    assert (vis[var]>=0 && vis[var]<CCTK_NumVars());
  }
  
  if (verbose) {
    for (var=0; var<nvars; ++var) {
      fullname = CCTK_FullName(vis[var]);
      assert (fullname);
      CCTK_VInfo (CCTK_THORNSTRING,
                  "Applying 90 degree rotating symmetry boundary conditions to \"%s\"",
                  fullname);
      free (fullname);
    }
  }
  
  /* Return early if there is nothing to do */
  if (nvars == 0) return 0;
  
  /* Get and check group info */
  assert (nvars>0);
  gis = malloc (nvars * sizeof *gis);
  assert (nvars==0 || gis);
  varptrs = malloc (nvars * sizeof *varptrs);
  assert (nvars==0 || varptrs);
  vartypes = malloc (nvars * sizeof  *vartypes);
  assert (nvars==0 || vartypes);
  vectorlengths = malloc (nvars * sizeof  *vectorlengths);
  assert (nvars==0 || vectorlengths);
  for (var=0; var<nvars; ++var) {
    gis[var] = CCTK_GroupIndexFromVarI (vis[var]);
    assert (gis[var]>=0 && gis[var]<CCTK_NumGroups());
    
    ierr = CCTK_GroupData (gis[var], &group);
    assert (!ierr);
    assert (group.grouptype == CCTK_GF);
    assert (group.disttype == CCTK_DISTRIB_DEFAULT);
    
    ierr = CCTK_GroupDynamicData (cctkGH, gis[var], &data);
    assert (!ierr);
    
    varptrs[var] = CCTK_VarDataPtrI (cctkGH, 0, vis[var]);
    assert (varptrs[var]);
    
    vartypes[var] = group.vartype;
    assert (vartypes[var] >= 0);

    vectorlengths[var] = group.vectorlength;
    assert (vectorlengths[var]>=0);
    assert (vectorlengths[var]==1 || group.vectorgroup);
  }
  
  for (d=0; d<group.dim; ++d) {
    x0[d] = CCTK_ORIGIN_SPACE(d);
    dx[d] = CCTK_DELTA_SPACE(d);
  }
  
  {
    assert(extent_valid_for_iteration == cctk_iteration);
  
    for (d=0; d<DIM; ++d) {
      fake_bbox[2*d  ] = data.lbnd[d] == global_lbnd[d];
      fake_bbox[2*d+1] = data.ubnd[d] == global_ubnd[d];
    }
  }
  
  /* directions */
  alldirs[0] = 0;
  alldirs[1] = 1;
  
  /* Find grid point that corresponds to the origin */
  for (q=0; q<2; ++q) {
    d = alldirs[q];
    /* x0 + dx * origin == 0 */
    origin[d] = - x0[d] / dx[d];
    dorigin[d] = origin[d] - floor(origin[d]);
    if (fabs(dorigin[d]) < 1.0e-6 || fabs(dorigin[d] - 1.0) < 1.0e-6) {
      avoid_origin[d] = 0;
    } else if (fabs(dorigin[d] - 0.5) < 1.0e-6) {
      avoid_origin[d] = 1;
    } else {
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "The coordinate origin in the %c-direction falls neither onto a grid point nor into the middle between two grid points.",
                  "xyz"[d]);
    }
    offset[d] = avoid_origin[d] ? 0 : 1;
  }
  
  parities = malloc (nvars * sizeof *parities);
  assert (nvars==0 || parities);
  srcptrs = malloc (nvars * sizeof *srcptrs);
  assert (nvars==0 || srcptrs);
  
  for (step=0; step<3; ++step) {
    switch (step) {
    case 0:
      /* x face */
      ndirs = 1;
      dir[0] = 0;
      otherdir[0] = 1;
      break;
    case 1:
      /* y face */
      ndirs = 1;
      dir[0] = 1;
      otherdir[0] = 0;
      break;
    case 2:
      /* xy edge */
      ndirs = 2;
      dir[0] = 0;
      dir[1] = 1;
      otherdir[0] = 1;
      otherdir[1] = 0;
      break;
    default:
      assert (0);
    }
    
    for (var=0; var<nvars; ++var) {
      
      /* Cactus lays out variables of a group
       * CCTK_REAL vel[3] {vx,vy,vz}
       * like this:
       * vx[0], vx[1], vx[2], vy[0], vy[1], ...
       * so the component x,y,z is given by index / vectorlength and the
       * vector index vx[i] by index % vectorlength.
       */
      const int vectorlength = vectorlengths[var];

      /* Find tensor type, source variable, and parity */
      {
        int table;

        numvars = CCTK_NumVarsInGroupI(gis[var]);
        assert (numvars>0);
        firstvar = CCTK_FirstVarIndexI(gis[var]);
        assert (firstvar>=0);
        index = vis[var] - firstvar;
        assert (index>=0 && index<numvars);
        table = CCTK_GroupTagsTableI(gis[var]);
        assert (table>=0);

        ierr = Util_TableGetString
          (table, sizeof tensortypealias, tensortypealias, "tensortypealias");
        if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
          /* assume a scalar */
          if (numvars != 1) {
            static int * restrict didwarn = 0;
            if (! didwarn) {
              didwarn = calloc (CCTK_NumGroups(), sizeof *didwarn);
            }
            if (! didwarn[gis[var]]) {
              char * groupname = CCTK_GroupName(gis[var]);
              assert (groupname);
              CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                          "Group \"%s\" has no tensor type and contains more than one element -- treating these as \"scalar\"",
                          groupname);
              free (groupname);
              didwarn[gis[var]] = 1;
            }
          }
          strcpy (tensortypealias, "scalar");
        } else if (ierr<0) {
          char * groupname = CCTK_GroupName(gis[var]);
          assert (groupname);
          CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "Error in tensor type alias declaration for group \"%s\"",
                      groupname);
          free (groupname);
        }
        
        if (CCTK_EQUALS (tensortypealias, "scalar")) {
          /* scalar */
        } else if (CCTK_EQUALS (tensortypealias, "4scalar")) {
          /* 4-scalar */
        } else if (CCTK_EQUALS (tensortypealias, "u")
                   || CCTK_EQUALS (tensortypealias, "d")) {
          /* vector */
          assert ((numvars == DIM*vectorlength) ||
                  (numvars == DIM && vectorlength == DIM));
        } else if (CCTK_EQUALS (tensortypealias, "4u")
                   || CCTK_EQUALS (tensortypealias, "4d")) {
          /* 4-vector */
          assert ((numvars == (DIM+1)*vectorlength) ||
                  (numvars == DIM+1 && vectorlength == DIM+1));
        } else if (CCTK_EQUALS (tensortypealias, "uu")
                   || CCTK_EQUALS (tensortypealias, "ud")
                   || CCTK_EQUALS (tensortypealias, "du")
                   || CCTK_EQUALS (tensortypealias, "dd")) {
          /* tensor */
          assert (numvars == DIM*DIM*vectorlength);
        } else if (CCTK_EQUALS (tensortypealias, "uu_sym")
                   || CCTK_EQUALS (tensortypealias, "dd_sym")) {
          /* symmetric tensor */
          assert (numvars == DIM*(DIM+1)/2*vectorlength);
        } else if (CCTK_EQUALS (tensortypealias, "ddd_sym")) {
          /* 3rd rank tensor, symmetric in last 2 indices */
          assert (numvars == DIM*DIM*(DIM+1)/2*vectorlength);
        } else if (CCTK_EQUALS (tensortypealias, "4uu_sym")
                   || CCTK_EQUALS (tensortypealias, "4dd_sym")) {
          /* symmetric 4-tensor */
          assert (numvars == (DIM+1)*(DIM+2)/2*vectorlength);
        } else if (CCTK_EQUALS (tensortypealias, "weylscalars_real")) {
          /* Weyl scalars, stored as 10 real numbers */
          assert (numvars == 10*vectorlength);
        } else {
          char * groupname = CCTK_GroupName(gis[var]);
          assert (groupname);
          CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "Illegal tensor type alias for group \"%s\"",
                      groupname);
          free (groupname);
        }
      }
      
      parities[var] = +1;
      if (CCTK_EQUALS (tensortypealias, "scalar")) {
        /* scalar */
        srcvi = vis[var];
        /* do nothing */
      } else if (CCTK_EQUALS (tensortypealias, "4scalar")) {
        /* 4-scalar */
        srcvi = vis[var];
        /* do nothing */
      } else if (CCTK_EQUALS (tensortypealias, "u")
                 || CCTK_EQUALS (tensortypealias, "d")) {
        /* vector */
        int srcindex;
        /* special case for vel[i] */
        assert (index>=0 && index<DIM*vectorlength);
        if(numvars == DIM && vectorlength == DIM) {
          srcindex = convert_index (step, index, alldirs, &parities[var]);
        } else {
          srcindex = index % vectorlength + vectorlength *
                     convert_index (step, index/vectorlength, alldirs, &parities[var]);
        }
        srcvi = firstvar + srcindex;
      } else if (CCTK_EQUALS (tensortypealias, "4u")
                 || CCTK_EQUALS (tensortypealias, "4d")) {
        /* 4-vector */
        assert (index>=0 && index<(DIM+1)*vectorlength);
        if(numvars == DIM+1 && vectorlength == DIM+1) {
          if (index==0) {
            /* temporal component */
            srcvi = firstvar;
          } else {
            /* spatial components */
            int srcindex;
            int const off = 1;
            srcindex = off + convert_index (step, index-off, alldirs, &parities[var]);
            srcvi = firstvar + srcindex;
          }
        } else {
          if (index/vectorlength==0) {
            /* temporal component */
            srcvi = firstvar;
          } else {
            /* spatial components */
            int srcindex;
            int const off = 1;
            srcindex = index % vectorlength + vectorlength * (off +
                       convert_index (step, index/vectorlength-off, alldirs, &parities[var]));
            srcvi = firstvar + srcindex;
          }
        }
      } else if (CCTK_EQUALS (tensortypealias, "uu")
                 || CCTK_EQUALS (tensortypealias, "ud")
                 || CCTK_EQUALS (tensortypealias, "du")
                 || CCTK_EQUALS (tensortypealias, "dd")) {
        /* tensor */
        assert (numvars == DIM*DIM);
        int index1, index2;
        int srcindex1, srcindex2;
        int srcindex;
        assert (index>=0 && index<DIM*DIM*vectorlength);
        index1 = (index / vectorlength) / DIM;
        index2 = (index / vectorlength) % DIM;
        srcindex1 = convert_index (step, index1, alldirs, &parities[var]);
        srcindex2 = convert_index (step, index2, alldirs, &parities[var]);
        srcindex = index % vectorlength +
                   (srcindex1 * DIM + srcindex2) * vectorlength;
        srcvi = firstvar + srcindex;
      } else if (CCTK_EQUALS (tensortypealias, "uu_sym")
                 || CCTK_EQUALS (tensortypealias, "dd_sym")) {
        /* symmetric tensor */
        int index1, index2;
        int srcindex1, srcindex2;
        int srcindex;
        assert (index>=0 && index<DIM*(DIM+1)/2*vectorlength);
        {
          int const expand1[DIM*(DIM+1)/2] = { 0,0,0,1,1,2 };
          int const expand2[DIM*(DIM+1)/2] = { 0,1,2,1,2,2 };
          index1 = expand1[index/vectorlength];
          index2 = expand2[index/vectorlength];
        }
        srcindex1 = convert_index (step, index1, alldirs, &parities[var]);
        srcindex2 = convert_index (step, index2, alldirs, &parities[var]);
        {
          int const compact[DIM][DIM] = { { 0,1,2 }, { 1,3,4 }, { 2,4,5 } };
          srcindex = index % vectorlength +
                     compact[srcindex1][srcindex2]*vectorlength;
        }
        srcvi = firstvar + srcindex;
      } else if (CCTK_EQUALS (tensortypealias, "ddd_sym")) {
        /* 3rd rank tensor, symmetric in last 2 indices */
        int index1, index2, index3;
        int srcindex1, srcindex2, srcindex3;
        int srcindex;
        assert (index>=0 && index<DIM*DIM*(DIM+1)/2*vectorlength);
        {
          int const expand1[DIM*DIM*(DIM+1)/2] =
            { 0,0,0,0,0,0, 1,1,1,1,1,1, 2,2,2,2,2,2 };
          int const expand2[DIM*DIM*(DIM+1)/2] =
            { 0,0,0,1,1,2, 0,0,0,1,1,2, 0,0,0,1,1,2 };
          int const expand3[DIM*DIM*(DIM+1)/2] =
            { 0,1,2,1,2,2, 0,1,2,1,2,2, 0,1,2,1,2,2 };
          index1 = expand1[index/vectorlength];
          index2 = expand2[index/vectorlength];
          index3 = expand3[index/vectorlength];
        }
        srcindex1 = convert_index (step, index1, alldirs, &parities[var]);
        srcindex2 = convert_index (step, index2, alldirs, &parities[var]);
        srcindex3 = convert_index (step, index3, alldirs, &parities[var]);
        {
          int const compact[DIM][DIM][DIM] = 
            { { {  0, 1, 2 }, {  1, 3, 4 }, {  2, 4, 5 } },
              { {  6, 7, 8 }, {  7, 9,10 }, {  8,10,11 } },
              { { 12,13,14 }, { 13,15,16 }, { 14,16,17 } } };
          srcindex = index % vectorlength +
                     compact[srcindex1][srcindex2][srcindex3]*vectorlength;
        }
        srcvi = firstvar + srcindex;
      } else if (CCTK_EQUALS (tensortypealias, "4uu_sym")
                 || CCTK_EQUALS (tensortypealias, "4dd_sym")) {
        /* symmetric 4-tensor */
        assert (index>=0 && index<(DIM+1)*(DIM+2)/2*vectorlength);
        if (index/vectorlength==0) {
          /* temporal-temporal component */
          srcvi = firstvar;
        } else if (index<(DIM+1)*vectorlength) {
          /* temporal-spatial components */
          int srcindex;
          int const off = 1;
          srcindex = index % vectorlength + vectorlength * (off +
                     convert_index (step, index/vectorlength-off, alldirs, &parities[var]));
          srcvi = firstvar + srcindex;
        } else {
          /* spatial-spatial components */
          int index1, index2;
          int srcindex1, srcindex2;
          int srcindex;
          int const off = DIM+1;
          {
            int const expand1[DIM*(DIM+1)/2] = { 0,0,0,1,1,2 };
            int const expand2[DIM*(DIM+1)/2] = { 0,1,2,1,2,2 };
            index1 = expand1[index/vectorlength-off];
            index2 = expand2[index/vectorlength-off];
          }
          srcindex1 = convert_index (step, index1, alldirs, &parities[var]);
          srcindex2 = convert_index (step, index2, alldirs, &parities[var]);
          {
            int const compact[DIM][DIM] = { { 0,1,2 }, { 1,3,4 }, { 2,4,5 } };
            srcindex = index % vectorlength + vectorlength *
                       (off + compact[srcindex1][srcindex2]);
          }
          srcvi = firstvar + srcindex;
        }
      } else if (CCTK_EQUALS (tensortypealias, "weylscalars_real")) {
        /* Weyl scalars, stored as 10 reals */
        assert (index>=0 && index<10*vectorlength);
        srcvi = vis[var];
      } else {
        assert (0);
      }
      
      srcptrs[var] = CCTK_VarDataPtrI (cctkGH, 0, srcvi);
      assert (srcptrs[var]);
      
    } /* for var */
    
    have_global_bbox = 1;
    for (q=0; q<ndirs; ++q) {
      have_global_bbox = have_global_bbox && global_bbox[2*dir[q]];
    }
    if (have_global_bbox) {
      
      for (q=0; q<ndirs; ++q) {
        assert (data.nghostzones[dir[q]] >= 0);
        
        if (data.gsh[otherdir[q]] < data.nghostzones[dir[q]] + data.nghostzones[otherdir[q]] + offset[otherdir[q]]) {
          assert (nvars > 0);
          var = 0;
          CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "The group \"%s\" has in the %c-direction only %d grid points.  "
                      "This is not large enough for a 90 degree rotating symmetry boundary that is %d grid points wide in the %c-direction.  "
                        "The group needs to have at least %d grid points in the %c-direction.",
                      CCTK_GroupNameFromVarI(vis[var]), "xyz"[dir[q]], data.gsh[dir[q]],
                      data.nghostzones[dir[q]],  "xyz"[otherdir[q]],
                      data.nghostzones[dir[q]] + data.nghostzones[otherdir[q]] + offset[otherdir[q]], "xyz"[dir[q]]);
        }
      }
      
      if (poison_boundaries) {
        /* poison destination grid points */
        int have_bnd = 1;
        for (q=0; q<ndirs; ++q) {
          have_bnd = have_bnd && cctkGH->cctk_bbox[2*dir[q]];
        }
        if (have_bnd) {
          
          int imin[3], imax[3];
          int i, j, k;
          for (d=0; d<3; ++d) {
            imin[d] = 0;
            imax[d] = cctk_lsh[d];
          }
          for (q=0; q<ndirs; ++q) {
            imax[dir[q]] = cctk_nghostzones[dir[q]];
          }
          
          var = 0;
          
          assert (group.dim == 3);
          switch (group.vartype) {
          case CCTK_VARIABLE_INT:
            /* do nothing */
            break;
          case CCTK_VARIABLE_REAL: {
            CCTK_REAL * restrict const varptr = varptrs[var];
            for (k=imin[0]; k<imax[2]; ++k) {
              for (j=imin[1]; j<imax[1]; ++j) {
                for (i=imin[2]; i<imax[0]; ++i) {
                  const int ind = CCTK_GFINDEX3D(cctkGH, i, j, k);
                  memset (&varptr[ind], poison_value, sizeof varptr[ind]);
                }
              }
            }
            break;
          }
          case CCTK_VARIABLE_COMPLEX:
            /* do nothing */
            break;
          default:
            assert (0);
          } /* switch grouptype */
        } /* if bbox */
      } /* if poison_boundaries */
      
      /* Allocate slab transfer description */
      xferinfo = malloc(group.dim * sizeof *xferinfo);
      assert (xferinfo);
      
      /* Fill in slab transfer description */
      for (d=0; d<group.dim; ++d) {
        xferinfo[d].src.gsh         = data.gsh        [d];
        xferinfo[d].src.lbnd        = data.lbnd       [d];
        xferinfo[d].src.lsh         = data.lsh        [d];
        xferinfo[d].src.ash         = data.ash        [d];
        xferinfo[d].src.lbbox       = fake_bbox       [2*d  ];
        xferinfo[d].src.ubbox       = fake_bbox       [2*d+1];
        xferinfo[d].src.nghostzones = data.nghostzones[d];
        xferinfo[d].src.off         = 0;
        xferinfo[d].src.str         = 1;
        xferinfo[d].src.len         = data.gsh        [d];
        
        xferinfo[d].dst.gsh         = data.gsh        [d];
        xferinfo[d].dst.lbnd        = data.lbnd       [d];
        xferinfo[d].dst.lsh         = data.lsh        [d];
        xferinfo[d].dst.ash         = data.ash        [d];
        xferinfo[d].dst.lbbox       = fake_bbox       [2*d  ];
        xferinfo[d].dst.ubbox       = fake_bbox       [2*d+1];
        xferinfo[d].dst.nghostzones = data.nghostzones[d];
        xferinfo[d].dst.off         = 0;
        xferinfo[d].dst.str         = 1;
        xferinfo[d].dst.len         = data.gsh        [d];
        
        xferinfo[d].xpose = d;
        xferinfo[d].flip  = 0;
      }
      
      switch (ndirs) {
      case 1:
        /* face */
        for (q=0; q<ndirs; ++q) {
          xferinfo[otherdir[q]].src.off = data.nghostzones[otherdir[q]] + offset[otherdir[q]];
          xferinfo[otherdir[q]].src.len = data.nghostzones[     dir[q]];
          xferinfo[     dir[q]].dst.off = 0;
          xferinfo[     dir[q]].dst.len = data.nghostzones[     dir[q]];
          
          xferinfo[     dir[q]].xpose = otherdir[q];
          xferinfo[otherdir[q]].xpose =      dir[q];
          /* Note: For rotations other than in the xy plane, the
             flipping might be different */
          xferinfo[     dir[q]].flip = 1;
        }
        break;
      case 2:
        /* edge */
        for (q=0; q<ndirs; ++q) {
          xferinfo[otherdir[q]].src.off = data.nghostzones[otherdir[q]] + offset[otherdir[q]];
          xferinfo[otherdir[q]].src.len = data.nghostzones[     dir[q]];
          xferinfo[     dir[q]].dst.off = 0;
          xferinfo[     dir[q]].dst.len = data.nghostzones[     dir[q]];
          
          xferinfo[dir[q]].flip = 1;
        }
        break;
      default:
        assert (0);
      }
      
      if (CCTK_IsFunctionAliased("GetLocalComponents")) {
        int const num_local_components = GetLocalComponents(cctkGH);
        if (num_local_components > 1) {
          CCTK_WARN (CCTK_WARN_ABORT,
                     "TAT/Slab can only be used if there is a single local component per MPI process");
        }
      }
      
      options = Util_TableCreateFromString ("useghosts=1");
      assert (options>=0);
      
      /* Can we use the more efficient interface? */
      if (CCTK_IsFunctionAliased ("GetRegriddingEpoch")) {
        static int epoch = -1;
        static struct slabsetup *restrict *restrict slab_setups[3] =
          { NULL, NULL, NULL};
        static int num_reflevels = 0;
        int const new_epoch = GetRegriddingEpoch (cctkGH);
        if (new_epoch > epoch) {
          epoch = new_epoch;
          /* Delete old slabbing setups */
          for (int s=0; s<3; ++s) {
            for (int rl=0; rl<num_reflevels; ++rl) {
              if (slab_setups[s][rl]) {
                Slab_MultiTransfer_Finalize (cctkGH, slab_setups[s][rl]);
              }
            }
            free (slab_setups[s]);
          }
          /* Allocate space for new slabbing setups */
          num_reflevels = GetRefinementLevels (cctkGH);
          for (int s=0; s<3; ++s) {
            slab_setups[s] = malloc (num_reflevels * sizeof *slab_setups[s]);
            assert (slab_setups[s]);
            for (int rl=0; rl<num_reflevels; ++rl) {
              slab_setups[s][rl] = NULL;
            }
          }
        }
        /* Use existing slabbing setup */
        assert (slab_setups[step]);
        int const reflevel = GetRefinementLevel (cctkGH);
        assert (reflevel>=0 && reflevel<num_reflevels);
        if (!slab_setups[step][reflevel]) {
          slab_setups[step][reflevel] = 
            Slab_MultiTransfer_Init (cctkGH, group.dim, xferinfo, options);
        }
        assert (slab_setups[step][reflevel]);
        ierr = Slab_MultiTransfer_Apply
          (cctkGH, slab_setups[step][reflevel],
           nvars, vartypes, srcptrs, vartypes, varptrs);
        assert (!ierr);
      } else {
        /* We can't use the more efficient interface, so fall back to
           the old one */
        ierr = Slab_MultiTransfer
          (cctkGH, group.dim, xferinfo, options,
           nvars, vartypes, srcptrs, vartypes, varptrs);
        assert (!ierr);
      }
      
      ierr = Util_TableDestroy (options);
      assert (!ierr);
      
      if (check_boundaries) {
        /* check destination grid points for poison */
        int have_bnd = 1;
        for (q=0; q<ndirs; ++q) {
          have_bnd = have_bnd && cctkGH->cctk_bbox[2*dir[q]];
        }
        if (have_bnd) {
          
          int imin[3], imax[3];
          int i, j, k;
          for (d=0; d<3; ++d) {
            imin[d] = 0;
            imax[d] = cctk_lsh[d];
          }
          for (q=0; q<ndirs; ++q) {
            imax[dir[q]] = cctk_nghostzones[dir[q]];
          }
          
          var = 0;
          
          assert (group.dim == 3);
          switch (group.vartype) {
          case CCTK_VARIABLE_INT:
            /* do nothing */
            break;
          case CCTK_VARIABLE_REAL: {
            int poison_found = 0;
            int nonfinite_found = 0;
            CCTK_REAL const * restrict const varptr = varptrs[var];
            CCTK_REAL poison;
            memset (&poison, poison_value, sizeof poison);
            for (k=imin[0]; k<imax[2]; ++k) {
              for (j=imin[1]; j<imax[1]; ++j) {
                for (i=imin[2]; i<imax[0]; ++i) {
                  const int ind = CCTK_GFINDEX3D(cctkGH, i, j, k);
                  if (memcmp (&varptr[ind], &poison, sizeof poison) == 0) {
                    printf ("   poison ijk=[%d,%d,%d] val=%g\n",
                            i, j, k, (double)varptr[ind]);
                    poison_found = 1;
                  }
                  if (! isfinite(varptr[ind])) {
                    printf ("   nonfinite ijk=[%d,%d,%d] val=%g\n",
                            i, j, k, (double)varptr[ind]);
                    nonfinite_found = 1;
                  }
                }
              }
            }
            if (poison_found || nonfinite_found) {
              if (poison_found) {
                printf ("Poison found:\n");
              }
              if (nonfinite_found) {
                printf ("Non-finite number found:\n");
              }
              printf ("   levfac=[%d,%d,%d]\n", cctk_levfac[0], cctk_levfac[1], cctk_levfac[2]);
              printf ("   origin_space=[%g,%g,%g]\n", cctk_origin_space[0], cctk_origin_space[1], cctk_origin_space[2]);
              printf ("   delta_space=[%g,%g,%g]\n", cctk_delta_space[0], cctk_delta_space[1], cctk_delta_space[2]);
              printf ("   lbnd=[%d,%d,%d]\n", cctk_lbnd[0], cctk_lbnd[1], cctk_lbnd[2]);
              printf ("   lsh=[%d,%d,%d]\n", cctk_lsh[0], cctk_lsh[1], cctk_lsh[2]);
              printf ("   gsh=[%d,%d,%d]\n", cctk_gsh[0], cctk_gsh[1], cctk_gsh[2]);
              printf ("   bbox=[%d,%d,%d,%d,%d,%d]\n", cctk_bbox[0], cctk_bbox[1], cctk_bbox[2], cctk_bbox[3], cctk_bbox[4], cctk_bbox[5]);
              if (poison_found) {
                CCTK_WARN (CCTK_WARN_ABORT, "Poison found in symmetry regions -- there is an error in this thorn");
              }
              if (nonfinite_found) {
                CCTK_WARN (CCTK_WARN_ALERT, "Non-finite number found in symmetry regions");
              }
            }
            break;
          }
          case CCTK_VARIABLE_COMPLEX:
            /* do nothing */
            break;
          default:
            assert (0);
          } /* switch grouptype */
        } /* if bbox */
      } /* if check_boundaries */
      
      /* take parity into account */
      have_local_bbox = 1;
      for (q=0; q<ndirs; ++q) {
        have_local_bbox = have_local_bbox && cctk_bbox[2*dir[q]];
      }
      if (have_local_bbox) {
        for (var=0; var<nvars; ++var) {
          assert (abs(parities[var]) == 1);
          if (parities[var] == -1) {
            int imin[DIM], imax[DIM];
            int i, j, k;
            for (d=0; d<DIM; ++d) {
              imin[d] = xferinfo[d].dst.off;
              imax[d] = xferinfo[d].dst.off + xferinfo[d].dst.len;
              imin[d] -= cctk_lbnd[d];
              imax[d] -= cctk_lbnd[d];
              if (imin[d] < 0) imin[d] = 0;
              if (imax[d] >= cctk_lsh[d]) imax[d] = cctk_lsh[d];
            }
            assert (group.dim == DIM);
            switch (group.vartype) {
            case CCTK_VARIABLE_INT: {
              CCTK_INT * restrict const varptr = varptrs[var];
              for (k=imin[2]; k<imax[2]; ++k) {
                for (j=imin[1]; j<imax[1]; ++j) {
                  for (i=imin[0]; i<imax[0]; ++i) {
                    const int ind = CCTK_GFINDEX3D(cctkGH, i, j, k);
                    varptr[ind] *= -1;
                  }
                }
              }
              break;
            }
            case CCTK_VARIABLE_REAL: {
              CCTK_REAL * restrict const varptr = varptrs[var];
              for (k=imin[2]; k<imax[2]; ++k) {
                for (j=imin[1]; j<imax[1]; ++j) {
                  for (i=imin[0]; i<imax[0]; ++i) {
                    const int ind = CCTK_GFINDEX3D(cctkGH, i, j, k);
                    varptr[ind] *= -1;
                  }
                }
              }
              break;
            }
            case CCTK_VARIABLE_COMPLEX: {
              CCTK_COMPLEX const czero = CCTK_Cmplx (0.0, 0.0);
              CCTK_COMPLEX * restrict const varptr = varptrs[var];
              for (k=imin[2]; k<imax[2]; ++k) {
                for (j=imin[1]; j<imax[1]; ++j) {
                  for (i=imin[0]; i<imax[0]; ++i) {
                    const int ind = CCTK_GFINDEX3D(cctkGH, i, j, k);
                    CCTK_COMPLEX * restrict const ptr = & varptr[ind];
                    * ptr = CCTK_CmplxSub (czero, * ptr);
                  }
                }
              }
              break;
            }
            default:
              assert (0);
            } /* switch grouptype */
          }
        } /* for var */
      
      } /* if bbox */
      
      free (xferinfo);
      
    } /* if have_global_bbox */
    
  } /* for step */
  
  free (gis);
  free (varptrs);
  free (vectorlengths);
  free (vartypes);
  free (parities);
  free (srcptrs);
  
  return 0;
}



void Rot90_ComputeLevelExtent (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  
  int nvars;

  assert (cctkGH);
  
  nvars = Boundary_SelectedGVs (cctkGH, 0, 0, 0, 0, 0, 0);
  assert (nvars>=0);
  
  if (nvars==0) return;

  {
    int max_handle;
    int d;                        /* 0..group.dim-1 */
    int ierr;
    CCTK_INT local[4*DIM], global[4*DIM];
    max_handle = CCTK_ReductionArrayHandle ("maximum");
    if (max_handle<0) CCTK_WARN (0, "Could not obtain reduction handle");
    
    for (d=0; d<2*DIM; ++d) local[      d] =  cctkGH->cctk_bbox[d];
    for (d=0; d<  DIM; ++d) local[2*DIM+d] = -cctkGH->cctk_lbnd[d];
    for (d=0; d<  DIM; ++d) local[3*DIM+d] =  cctkGH->cctk_ubnd[d];
    ierr = CCTK_ReduceLocArrayToArray1D
      (cctkGH, -1, max_handle, local, global, 4*DIM, CCTK_VARIABLE_INT);
    assert(!ierr);
    for (d=0; d<2*DIM; ++d) global_bbox[d] =  global[      d];
    for (d=0; d<  DIM; ++d) global_lbnd[d] = -global[2*DIM+d];
    for (d=0; d<  DIM; ++d) global_ubnd[d] =  global[3*DIM+d];

    /* record when we ran to have some sanity check against using old data in
     * global variables */
    extent_valid_for_iteration = cctk_iteration;
  }
}



void Rot90_ApplyBC (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  
  int nvars;
  CCTK_INT * restrict indices;
  CCTK_INT * restrict faces;
  CCTK_INT * restrict widths;
  CCTK_INT * restrict tables;
  int var;
  int ierr;
  
  assert (cctkGH);
  
  nvars = Boundary_SelectedGVs (cctkGH, 0, 0, 0, 0, 0, 0);
  assert (nvars>=0);
  
  if (nvars==0) return;
  
  indices = malloc (nvars * sizeof *indices);
  assert (indices);
  faces = malloc (nvars * sizeof *faces);
  assert (faces);
  widths = malloc (nvars * sizeof *widths);
  assert (widths);
  tables = malloc (nvars * sizeof *tables);
  assert (tables);
  
  ierr =  Boundary_SelectedGVs
    (cctkGH, nvars, indices, faces, widths, tables, 0);
  assert (ierr == nvars);
  
  for (var=0; var<nvars; ++var) {
    assert (indices[var]>=0 && indices[var]<CCTK_NumVars());
    assert (widths[var]>=0);
  }
  
  ierr = BndRot90VI (cctkGH, nvars, indices);
  assert (!ierr);
  
  free (indices);
  free (faces);
  free (widths);
  free (tables);
}
