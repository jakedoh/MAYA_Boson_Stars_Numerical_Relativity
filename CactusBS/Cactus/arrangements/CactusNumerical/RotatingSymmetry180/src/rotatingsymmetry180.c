#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <util_ErrorCodes.h>
#include <util_Table.h>

#include <Slab.h>

#include "rotatingsymmetry180.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* bbox, lbnd and ubnd of combined level box */
static int global_bbox[6];
static int global_lbnd[3], global_ubnd[3];
static int extent_valid_for_iteration = -1;


int BndRot180VI (cGH const * restrict const cctkGH,
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
  
  int (* restrict paritiess)[3];
  
  int fake_bbox[6];
  
  CCTK_REAL x0[3], dx[3];
  CCTK_REAL symbnd[3];          /* location of symmetry boundary */
  CCTK_REAL origin[3], dorigin[3];
  int avoid_origin[3], iorigin[3];
  int offset[3];                /* offset 0..1 due to avoid_origin */
  
  struct xferinfo * restrict xferinfo;
  
  int var;
  
  int alldirs[2];
  int dir;                      /* direction of the symmetry face */
  int otherdir;                 /* the other direction of the rotation */
  int q;
  int d;
  int icnt;
  int ierr;
  
  /* Check arguments */
  assert (cctkGH);
  assert (nvars >= 0);
  assert (nvars==0 || vis);
  for (var=0; var<nvars; ++var) {
    assert (vis[var]>=0 && vis[var]<CCTK_NumVars());
  }

  if (verbose) {
    for (var=0; var<nvars; ++var) {
      fullname = CCTK_FullName(vis[var]);
      assert (fullname);
      CCTK_VInfo (CCTK_THORNSTRING,
                  "Applying 180 degree rotating symmetry boundary conditions to \"%s\"",
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
  
  /* find parity */
  paritiess = malloc (nvars * sizeof *paritiess);
  assert (nvars==0 || paritiess);
  
  for (var=0; var<nvars; ++var) {
    
    char tensortypealias[100];
    int firstvar, numvars;
    int table;
    int index;

    /* Cactus lays out variables of a group
     * CCTK_REAL vel[3] {vx,vy,vz}
     * like this:
     * vx[0], vx[1], vx[2], vy[0], vy[1], ...
     * so the component x,y,z is given by index / vectorlength and the
     * vector index vx[i] by index % vectorlength.
     */
    const int vectorlength = vectorlengths[var];
    
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
      if (numvars != vectorlength) {
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
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "Error in tensor type alias declaration for group \"%s\"",
                  groupname);
    }
    
    if (CCTK_EQUALS (tensortypealias, "scalar")) {
      paritiess[var][0] = paritiess[var][1] = paritiess[var][2] = +1;
    } else if (CCTK_EQUALS (tensortypealias, "4scalar")) {
      paritiess[var][0] = paritiess[var][1] = paritiess[var][2] = +1;
    } else if (CCTK_EQUALS (tensortypealias, "u")
               || CCTK_EQUALS (tensortypealias, "d")) {
      assert ((numvars == 3*vectorlength) ||
              (numvars == 3 && vectorlength ==3));
      paritiess[var][0] = paritiess[var][1] = paritiess[var][2] = +1;
      /* special case for vel[i] */
      if(numvars == vectorlength) {
        paritiess[var][index] = -1;
      } else {
        paritiess[var][index/vectorlength] = -1;
      }
    } else if (CCTK_EQUALS (tensortypealias, "4u")
               || CCTK_EQUALS (tensortypealias, "4d")) {
      assert ((numvars == 4*vectorlength) ||
              (numvars == 4 && vectorlength ==4));
      if(numvars == vectorlength) {
        if (index == 0) {
          paritiess[var][0] = paritiess[var][1] = paritiess[var][2] = +1;
        } else {
          paritiess[var][0] = paritiess[var][1] = paritiess[var][2] = +1;
          paritiess[var][index-1] = -1;
        }
      } else {
        if (index/vectorlength == 0) {
          paritiess[var][0] = paritiess[var][1] = paritiess[var][2] = +1;
        } else {
          paritiess[var][0] = paritiess[var][1] = paritiess[var][2] = +1;
          paritiess[var][index/vectorlength-1] = -1;
        }
      }
    } else if (CCTK_EQUALS (tensortypealias, "uu_sym")
               || CCTK_EQUALS (tensortypealias, "dd_sym")) {
      assert (numvars == 6*vectorlength);
      paritiess[var][0] = paritiess[var][1] = paritiess[var][2] = +1;
      switch (index/vectorlength) {
      case 0: break;
      case 1: paritiess[var][0] = paritiess[var][1] = -1; break;
      case 2: paritiess[var][0] = paritiess[var][2] = -1; break;
      case 3: break;
      case 4: paritiess[var][1] = paritiess[var][2] = -1; break;
      case 5: break;
      default: assert(0);
      }
    } else if (CCTK_EQUALS (tensortypealias, "uu")
               || CCTK_EQUALS (tensortypealias, "ud")
               || CCTK_EQUALS (tensortypealias, "du")
               || CCTK_EQUALS (tensortypealias, "dd")) {
      assert (numvars == 9*vectorlength);
      paritiess[var][0] = paritiess[var][1] = paritiess[var][2] = +1;
      int const d1 = (index / vectorlength) % 3;
      int const d2 = (index / vectorlength) / 3;
      paritiess[var][d1] *= -1;
      paritiess[var][d2] *= -1;
    } else if (CCTK_EQUALS (tensortypealias, "dd_sym_d")) {
      assert (numvars == 18*vectorlength);
      paritiess[var][0] = paritiess[var][1] = paritiess[var][2] = +1;
      switch ((index / vectorlength) / 3) {
      case 0: break;
      case 1: paritiess[var][0] = paritiess[var][1] = -1; break;
      case 2: paritiess[var][0] = paritiess[var][2] = -1; break;
      case 3: break;
      case 4: paritiess[var][1] = paritiess[var][2] = -1; break;
      case 5: break;
      default: assert(0);
      }
      paritiess[var][(index / vectorlength) % 3] *= -1;
    } else if (CCTK_EQUALS (tensortypealias, "ddd_sym")) {
      assert (numvars == 18*vectorlength);
      paritiess[var][0] = paritiess[var][1] = paritiess[var][2] = +1;
      switch ((index / vectorlength) % 6) {
      case 0: break;
      case 1: paritiess[var][0] = paritiess[var][1] = -1; break;
      case 2: paritiess[var][0] = paritiess[var][2] = -1; break;
      case 3: break;
      case 4: paritiess[var][1] = paritiess[var][2] = -1; break;
      case 5: break;
      default: assert(0);
      }
      paritiess[var][(index / vectorlength) / 6] *= -1;
    } else if (CCTK_EQUALS (tensortypealias, "4uu_sym")
               || CCTK_EQUALS (tensortypealias, "4dd_sym")) {
      assert (numvars == 10*vectorlength);
      if (index / vectorlength == 0) {
        paritiess[var][0] = paritiess[var][1] = paritiess[var][2] = +1;
      } else if (index / vectorlength < 4) {
        paritiess[var][0] = paritiess[var][1] = paritiess[var][2] = +1;
        paritiess[var][index/vectorlength-1] = -1;
      } else {
        paritiess[var][0] = paritiess[var][1] = paritiess[var][2] = +1;
        switch (index/vectorlength-4) {
        case 0: break;
        case 1: paritiess[var][0] = paritiess[var][1] = -1; break;
        case 2: paritiess[var][0] = paritiess[var][2] = -1; break;
        case 3: break;
        case 4: paritiess[var][1] = paritiess[var][2] = -1; break;
        case 5: break;
        default: assert(0);
        }
      }
    } else if (CCTK_EQUALS (tensortypealias, "weylscalars_real")) {
      assert (numvars == 10*vectorlength);
      {
        static int const weylparities[10][3] =
          {{+1,+1,+1},
           {-1,-1,-1},
           {+1,+1,+1},
           {-1,-1,-1},
           {+1,+1,+1},
           {-1,-1,-1},
           {+1,+1,+1},
           {-1,-1,-1},
           {+1,+1,+1},
           {-1,-1,-1}};
        for (d=0; d<3; ++d) {
          paritiess[var][d] = weylparities[index/vectorlength][d];
        }
      }
    } else if (CCTK_EQUALS (tensortypealias, "ManualCartesian")) {
	RotatingSymmetry180_GetManualParities(table, gis[var], paritiess[var]);
    } else {
      char * groupname = CCTK_GroupName(gis[var]);
      assert (groupname);
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "Illegal tensor type alias for group \"%s\"",
                  groupname);
    }
    
    assert (abs(paritiess[var][0])==1 && abs(paritiess[var][1])==1
            && abs(paritiess[var][2])==1);
    
  } /* for var */
  
  {
    assert(extent_valid_for_iteration == cctk_iteration);
  
    for (d=0; d<3; ++d) {
      fake_bbox[2*d  ] = data.lbnd[d] == global_lbnd[d];
      fake_bbox[2*d+1] = data.ubnd[d] == global_ubnd[d];
    }
  }
  
  /* directions */
  alldirs[0] = 0;
  alldirs[1] = 1;
  
  /* Find location of symmetry boundary */
  if (use_coordbase) {
    CCTK_INT const size = 3;
    CCTK_REAL physical_min[3];
    CCTK_REAL physical_max[3];
    CCTK_REAL interior_min[3];
    CCTK_REAL interior_max[3];
    CCTK_REAL exterior_min[3];
    CCTK_REAL exterior_max[3];
    CCTK_REAL spacing;
    GetDomainSpecification
      (size,
       physical_min, physical_max,
       interior_min, interior_max,
       exterior_min, exterior_max,
       & spacing);
    symbnd[0] = physical_min[0];
    symbnd[1] = physical_min[1];
    symbnd[2] = physical_min[2];
  } else {
    symbnd[0] = symmetry_boundary_x;
    symbnd[1] = symmetry_boundary_y;
    symbnd[2] = 0;              /* unused */
  }
  
  /* Find grid point that corresponds to the origin */
  for (q=0; q<2; ++q) {
    d = alldirs[q];
    /* x0 + dx * origin == symbnd */
    origin[d] = (symbnd[d] - x0[d]) / dx[d];
    dorigin[d] = origin[d] - floor(origin[d]);
    if (fabs(dorigin[d]) < 1.0e-6 || fabs(dorigin[d] - 1.0) < 1.0e-6) {
      avoid_origin[d] = 0;
    } else if (fabs(dorigin[d] - 0.5) < 1.0e-6) {
      avoid_origin[d] = 1;
    } else {
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "The coordinate origin in the %c-direction falls neither onto a grid point nor into the middle between two grid points.",
                  "xyz"[d]);
    }
    iorigin[d] = lrint(origin[d] + (avoid_origin[d] ? 0.5 : 0.0));
    offset[d] = avoid_origin[d] ? 0 : 1;
  }
  
  /* x direction */
  dir = 0;
  otherdir = 1;
  
  assert (data.nghostzones[dir] >= 0);
  
  if (global_bbox[2*dir]) {
    
    if (2*iorigin[otherdir] + offset[otherdir] != cctk_gsh[otherdir]) {
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "The coordinate origin in the %c-direction is not in the centre of the domain.  The boundary condition cannot be applied.",
                  "xyz"[otherdir]);
    }
    
    if (iorigin[dir] != data.nghostzones[dir]) {
      assert (nvars > 0);
      var = 0;
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "The group \"%s\" has in the %c-direction %d symmetry points (grid points outside of the symmetry boundary).  "
                  "This is not equal to the number of ghost zones, which is %d.  "
                  "The number of symmetry points must be equal to the number of ghost zones.",
                  CCTK_GroupNameFromVarI(vis[var]), "xyz"[dir],
                  iorigin[dir], data.nghostzones[dir]);
    }
    
    if (data.gsh[dir] < 2*data.nghostzones[dir] + offset[dir]) {
      assert (nvars > 0);
      var = 0;
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "The group \"%s\" has in the %c-direction only %d grid points.  "
                  "This is not large enough for a 180 degree rotating symmetry boundary that is %d grid points wide.  "
                  "The group needs to have at least %d grid points in that direction.",
                  CCTK_GroupNameFromVarI(vis[var]), "xyz"[dir], data.gsh[dir],
                  data.nghostzones[dir],
                  2*data.nghostzones[dir] + offset[dir]);
    }
    
    if (poison_boundaries) {
      /* poison destination grid points */
      if (cctkGH->cctk_bbox[2*dir]) {
        
        int imin[3], imax[3];
        int i, j, k;
        for (d=0; d<3; ++d) {
          imin[d] = 0;
          imax[d] = cctk_lsh[d];
        }
        imax[dir] = cctk_nghostzones[dir];
        
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
          assert(0);
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
    
    xferinfo[dir].src.off = data.nghostzones[dir] + offset[dir];
    xferinfo[dir].src.len = data.nghostzones[dir];
    xferinfo[dir].dst.off = 0;
    xferinfo[dir].dst.len = data.nghostzones[dir];
    
    xferinfo[     dir].flip = 1;
    xferinfo[otherdir].flip = 1;
    
    if (CCTK_EQUALS (hyperslabber, "TAT/Slab")) {
      
      if (CCTK_IsFunctionAliased("GetLocalComponents")) {
        int const num_local_components = GetLocalComponents(cctkGH);
        if (num_local_components > 1) {
          CCTK_ERROR("TAT/Slab can only be used if there is a single local component per MPI process");
        }
      }
      
      int options;
      options = Util_TableCreateFromString ("useghosts=1");
      assert (options>=0);
      
      /* Can we use the more efficient interface? */
      if (CCTK_IsFunctionAliased ("GetRegriddingEpoch")) {
        static int old_epoch = -1;
        static struct slabsetup *restrict *restrict slab_setups = NULL;
        static int num_slab_setups = 0;
        int const epoch = GetRegriddingEpoch (cctkGH);
        if (epoch > old_epoch) {
          old_epoch = epoch;
          /* Delete old slabbing setups */
          for (int rl=0; rl<num_slab_setups; ++rl) {
            if (slab_setups[rl]) {
              Slab_MultiTransfer_Finalize (cctkGH, slab_setups[rl]);
            }
          }
          free ((void*)slab_setups);
          num_slab_setups = GetRefinementLevels (cctkGH);
          slab_setups = malloc (num_slab_setups * sizeof *slab_setups);
          assert (slab_setups);
          for (int rl=0; rl<num_slab_setups; ++rl) {
            slab_setups[rl] = NULL;
          }
        }
        /* Use existing slabbing setup */
        assert (slab_setups);
        int const reflevel = GetRefinementLevel (cctkGH);
        assert (reflevel>=0 && reflevel<num_slab_setups);
        if (!slab_setups[reflevel]) {
          slab_setups[reflevel] = 
            Slab_MultiTransfer_Init (cctkGH, group.dim, xferinfo, options);
        }
        assert (slab_setups[reflevel]);
        ierr = Slab_MultiTransfer_Apply
          (cctkGH, slab_setups[reflevel],
           nvars, vartypes, varptrs, vartypes, varptrs);
        assert (!ierr);
      } else {
        /* We can't use the more efficient interface, so fall back to
           the old one */
        ierr = Slab_MultiTransfer
          (cctkGH, group.dim, xferinfo, options,
           nvars, vartypes, varptrs, vartypes, varptrs);
        assert (!ierr);
      }
      
      ierr = Util_TableDestroy (options);
      assert (!ierr);
      
    } else if (CCTK_EQUALS (hyperslabber, "GetHyperslab")) {
      
      CCTK_INT const direction[3*3] = {1, 0, 0,   0, 1, 0,   0, 0, 1};
      CCTK_INT const origin[3] =
        {xferinfo[0].src.off, xferinfo[1].src.off, xferinfo[2].src.off};
      CCTK_INT const extent[3] =
        {xferinfo[0].src.len, xferinfo[1].src.len, xferinfo[2].src.len};
      CCTK_INT const downsample[3] =
        {xferinfo[0].src.str, xferinfo[1].src.str, xferinfo[2].src.str};
      int mapping;
      CCTK_INT hsize[3], total_hsize;
      CCTK_POINTER hdata[nvars];
      
      assert (nvars > 0);
      mapping = Hyperslab_GlobalMappingByIndex
        (cctkGH, vis[0], 3, direction, origin, extent, downsample,
         -1, NULL, hsize);
      assert (mapping>=0);
      
      total_hsize = hsize[0] * hsize[1] * hsize[2];
      assert (total_hsize >= 0);
      for (var=0; var<nvars; ++var) {
        hdata[var] = malloc(total_hsize * CCTK_VarTypeSize(CCTK_VarTypeI(var)));
        assert (total_hsize==0 || hdata[var]);
      }
      
      icnt = Hyperslab_GetList
        (cctkGH, mapping, nvars, NULL, vis, NULL, NULL, hdata, NULL);
      assert (icnt == nvars);
      
      // copy hyperslabs into grid functions
      assert (0);
      
      for (var=0; var<nvars; ++var) {
        free (hdata[var]);
      }
      
      ierr = Hyperslab_FreeMapping (mapping);
      assert (!ierr);
      
    } else {
      assert(0);
    }
    
    if (check_boundaries) {
      /* check destination grid points for poison */
      if (cctkGH->cctk_bbox[2*dir]) {
        
        int imin[3], imax[3];
        int i, j, k;
        for (d=0; d<3; ++d) {
          imin[d] = 0;
          imax[d] = cctk_lsh[d];
        }
        imax[dir] = cctk_nghostzones[dir];
        
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
              CCTK_ERROR("Poison found in symmetry regions -- there is an error in this thorn");
            }
            if (nonfinite_found) {
              CCTK_WARN(CCTK_WARN_ALERT, "Non-finite number found in symmetry regions");
            }
          }
          break;
        }
        case CCTK_VARIABLE_COMPLEX:
          /* do nothing */
          break;
        default:
          assert(0);
        } /* switch grouptype */
      } /* if bbox */
    } /* if check_boundaries */
    
    /* take parity into account */
    if (cctkGH->cctk_bbox[2*dir]) {
      for (var=0; var<nvars; ++var) {
        int parity;
        
        parity = paritiess[var][0] * paritiess[var][1];
        assert (abs(parity) == 1);
        if (parity == -1) {
          int imin[3], imax[3];
          int i, j, k;
          for (d=0; d<3; ++d) {
            imin[d] = xferinfo[d].dst.off;
            imax[d] = xferinfo[d].dst.off + xferinfo[d].dst.len;
            imin[d] -= cctk_lbnd[d];
            imax[d] -= cctk_lbnd[d];
            if (imin[d] < 0) imin[d] = 0;
            if (imax[d] >= cctk_lsh[d]) imax[d] = cctk_lsh[d];
          }
          assert (group.dim == 3);
          switch (group.vartype) {
          case CCTK_VARIABLE_INT: {
            CCTK_INT * restrict const varptr = varptrs[var];
            for (k=imin[0]; k<imax[2]; ++k) {
              for (j=imin[1]; j<imax[1]; ++j) {
                for (i=imin[2]; i<imax[0]; ++i) {
                  const int ind = CCTK_GFINDEX3D(cctkGH, i, j, k);
                  varptr[ind] *= -1;
                }
              }
            }
            break;
          }
          case CCTK_VARIABLE_REAL: {
            CCTK_REAL * restrict const varptr = varptrs[var];
            for (k=imin[0]; k<imax[2]; ++k) {
              for (j=imin[1]; j<imax[1]; ++j) {
                for (i=imin[2]; i<imax[0]; ++i) {
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
            for (k=imin[0]; k<imax[2]; ++k) {
              for (j=imin[1]; j<imax[1]; ++j) {
                for (i=imin[2]; i<imax[0]; ++i) {
                  const int ind = CCTK_GFINDEX3D(cctkGH, i, j, k);
                  CCTK_COMPLEX * restrict const ptr = & varptr[ind];
                  * ptr = CCTK_CmplxSub (czero, * ptr);
                }
              }
            }
            break;
          }
          default:
            assert(0);
          } /* switch grouptype */
        }
      } /* for var */
    } /* if bbox */
    
    free (xferinfo);
    
  } /* if global_bbox */
  
  free (gis);
  free ((void*)varptrs);
  free (vectorlengths);
  free (vartypes);
  free (paritiess);
  
  return 0;
}



void Rot180_ComputeLevelExtent (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  
  int nvars;

  assert (cctkGH);
  
  nvars = Boundary_SelectedGVs (cctkGH, 0, 0, 0, 0, 0, 0);
  assert (nvars>=0);
  
  if (nvars==0) return;
  
  {
#if 0
    int min_handle, max_handle;
    int d;                        /* 0..group.dim-1 */
    int ierr;
    CCTK_REAL local[6], global[6];
    min_handle = CCTK_ReductionArrayHandle ("minimum");
    if (min_handle<0) CCTK_ERROR("Could not obtain reduction handle");
    max_handle = CCTK_ReductionArrayHandle ("maximum");
    if (max_handle<0) CCTK_ERROR("Could not obtain reduction handle");
    
    for (d=0; d<6; ++d) local[d] = cctkGH->cctk_bbox[d];
    ierr = CCTK_ReduceLocArrayToArray1D
      (cctkGH, -1, max_handle, local, global, 6, CCTK_VARIABLE_REAL);
    for (d=0; d<6; ++d) global_bbox[d] = (int)global[d];
    
    for (d=0; d<3; ++d) local[d] = cctkGH->cctk_lbnd[d];
    ierr = CCTK_ReduceLocArrayToArray1D
      (cctkGH, -1, min_handle, local, global, 3, CCTK_VARIABLE_REAL);
    for (d=0; d<3; ++d) global_lbnd[d] = (int)global[d];
    
    for (d=0; d<3; ++d) local[d] = cctkGH->cctk_ubnd[d];
    ierr = CCTK_ReduceLocArrayToArray1D
      (cctkGH, -1, max_handle, local, global, 3, CCTK_VARIABLE_REAL);
    for (d=0; d<3; ++d) global_ubnd[d] = (int)global[d];
#else
    int max_handle;
    int d;                        /* 0..group.dim-1 */
    int ierr;
    CCTK_INT local[12], global[12];
    max_handle = CCTK_ReductionArrayHandle ("maximum");
    if (max_handle<0) CCTK_ERROR("Could not obtain reduction handle");
    
    for (d=0; d<6; ++d) local[  d] =  cctkGH->cctk_bbox[d];
    for (d=0; d<3; ++d) local[6+d] = -cctkGH->cctk_lbnd[d];
    for (d=0; d<3; ++d) local[9+d] =  cctkGH->cctk_ubnd[d];
    ierr = CCTK_ReduceLocArrayToArray1D
      (cctkGH, -1, max_handle, local, global, 12, CCTK_VARIABLE_INT);
    assert(!ierr);
    for (d=0; d<6; ++d) global_bbox[d] =  global[  d];
    for (d=0; d<3; ++d) global_lbnd[d] = -global[6+d];
    for (d=0; d<3; ++d) global_ubnd[d] =  global[9+d];

    /* record when we ran to have some sanity check against using old data in
     * global variables */
    extent_valid_for_iteration = cctk_iteration;
#endif
  }
}



void Rot180_ApplyBC (CCTK_ARGUMENTS)
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
    assert (widths[var] >= 0);
  }
  
  ierr = BndRot180VI (cctkGH, nvars, indices);
  assert (!ierr);
  
  free (indices);
  free (faces);
  free (widths);
  free (tables);
}
