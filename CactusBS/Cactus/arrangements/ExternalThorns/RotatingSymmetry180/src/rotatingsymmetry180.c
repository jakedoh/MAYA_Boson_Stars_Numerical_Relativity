/* $Header$ */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

#include "Slab.h"
#include <stdio.h>
#include "rotatingsymmetry180.h"

void *pmalloc(size_t size)
{
  void *p = malloc(size);
  memset(p, 255, size);
  return p;
}

int BndRot180VI (cGH const * restrict const cctkGH,
                 int const nvars,
                 CCTK_INT const * restrict const vis)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  int * restrict gis = NULL;
  cGroup group;
  cGroupDynamicData data;
  char * restrict fullname = NULL;
  void * restrict * restrict varptrs = NULL;
  int * restrict vartypes = NULL;
  
  //int (*paritiess)[3] = {NULL}; // Not restricted in Ian's version
  int (* restrict paritiess)[3];
  
  int global_bbox[6] = {0,0,0,0,0,0};
  int global_lbnd[3] = {0,0,0}, global_ubnd[3] = {0,0,0};
  int fake_bbox[6] = {0,0,0,0,0,0};
  
  CCTK_REAL x0[3] = {0,0,0}, dx[3] = {0,0,0};
  CCTK_REAL symbnd[3] = {0,0,0};          /* location of symmetry boundary */
  CCTK_REAL origin[3] = {0,0,0}, dorigin[3] = {0,0,0};
  int avoid_origin[3] = {0,0,0}, iorigin[3] = {0,0,0};
  int offset[3] = {0,0,0};                /* offset 0..1 due to avoid_origin */
  
  struct xferinfo * restrict xferinfo = NULL;
  int options = 0;
  
  int var = 0;
  
  int alldirs[2] = {0,0};
  int dir = 0;                      /* direction of the symmetry face */
  int otherdir = 0;                 /* the other direction of the rotation */
  int q = 0;
  int d = 0;
  int ierr = 0;
  
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
  gis = pmalloc (nvars * sizeof *gis);
  assert (nvars==0 || gis);
  varptrs = pmalloc (nvars * sizeof *varptrs);
  assert (nvars==0 || varptrs);
  vartypes = pmalloc (nvars * sizeof  *vartypes);
  assert (nvars==0 || vartypes);
  for (var=0; var<nvars; ++var) {
    gis[var] = CCTK_GroupIndexFromVarI (vis[var]);
    assert (gis[var]>=0 && gis[var]<CCTK_NumGroups());
    
    ierr = CCTK_GroupData (gis[var], &group);
    assert (!ierr);
    assert (group.grouptype == CCTK_GF);
    assert (group.disttype == CCTK_DISTRIB_DEFAULT);
    assert (group.stagtype == 0);
    
    ierr = CCTK_GroupDynamicData (cctkGH, gis[var], &data);
    assert (!ierr);
    
    varptrs[var] = CCTK_VarDataPtrI (cctkGH, 0, vis[var]);
    assert (varptrs[var]);
    
    vartypes[var] = group.vartype;
    assert (vartypes[var] >= 0);
  }
  
  for (d=0; d<group.dim; ++d) {
    x0[d] = CCTK_ORIGIN_SPACE(d);
    dx[d] = CCTK_DELTA_SPACE(d);
  }
  
  /* find parity */
  paritiess = pmalloc (nvars * sizeof *paritiess);
  assert (nvars==0 || paritiess);
  
  for (var=0; var<nvars; ++var) {
    
    char tensortypealias[100];
    int firstvar=0, numvars=0;
    int table=0;
    int index=0;
    
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
      paritiess[var][0] = paritiess[var][1] = paritiess[var][2] = +1;
    } else if (CCTK_EQUALS (tensortypealias, "4scalar")) {
      paritiess[var][0] = paritiess[var][1] = paritiess[var][2] = +1;
    } else if (CCTK_EQUALS (tensortypealias, "u")
               || CCTK_EQUALS (tensortypealias, "d")) {
      assert (numvars == 3);
      paritiess[var][0] = paritiess[var][1] = paritiess[var][2] = +1;
      paritiess[var][index] = -1;
    } else if (CCTK_EQUALS (tensortypealias, "4u")
               || CCTK_EQUALS (tensortypealias, "4d")) {
      assert (numvars == 4);
      if (index == 0) {
        paritiess[var][0] = paritiess[var][1] = paritiess[var][2] = +1;
      } else {
        paritiess[var][0] = paritiess[var][1] = paritiess[var][2] = +1;
        paritiess[var][index-1] = -1;
      }
    } else if (CCTK_EQUALS (tensortypealias, "uu_sym")
               || CCTK_EQUALS (tensortypealias, "dd_sym")) {
      assert (numvars == 6);
      paritiess[var][0] = paritiess[var][1] = paritiess[var][2] = +1;
      switch (index) {
      case 0: break;
      case 1: paritiess[var][0] = paritiess[var][1] = -1; break;
      case 2: paritiess[var][0] = paritiess[var][2] = -1; break;
      case 3: break;
      case 4: paritiess[var][1] = paritiess[var][2] = -1; break;
      case 5: break;
      default: assert(0);
      }
    } else if (CCTK_EQUALS (tensortypealias, "ManualCartesian")) {
	RotatingSymmetry180_GetManualParities(table, gis[var], paritiess[var]);
	/*    printf("RotatingSymmetry180: Parities for group %s are %d, %d, %d\n", 
	 CCTK_GroupName(gis[var]), paritiess[var][0], 
	 paritiess[var][1], paritiess[var][2]);*/
    } else if (CCTK_EQUALS (tensortypealias, "4uu_sym")
               || CCTK_EQUALS (tensortypealias, "4dd_sym")) {
      assert (numvars == 10);
      if (index == 0) {
        paritiess[var][0] = paritiess[var][1] = paritiess[var][2] = +1;
      } else if (index < 4) {
        paritiess[var][0] = paritiess[var][1] = paritiess[var][2] = +1;
        paritiess[var][index-1] = -1;
      } else {
        paritiess[var][0] = paritiess[var][1] = paritiess[var][2] = +1;
        switch (index-4) {
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
      assert (numvars == 10);
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
          paritiess[var][d] = weylparities[index][d];
        }
      }
    } else {
      char * groupname = CCTK_GroupName(gis[var]);
      assert (groupname);
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Illegal tensor type alias for group \"%s\"",
                  groupname);
      free (groupname);
    }
    
    assert (abs(paritiess[var][0])==1 && abs(paritiess[var][1])==1
            && abs(paritiess[var][2])==1);
    
  } /* for var */
  
  {
#if 0
    int min_handle, max_handle;
    CCTK_REAL local[6], global[6];
    min_handle = CCTK_ReductionArrayHandle ("minimum");
    if (min_handle<0) CCTK_WARN (0, "Could not obtain reduction handle");
    max_handle = CCTK_ReductionArrayHandle ("maximum");
    if (max_handle<0) CCTK_WARN (0, "Could not obtain reduction handle");
    
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
    int max_handle = 0;
    CCTK_INT local[12]={0,0,0,0,0,0,0,0,0,0,0,0}, global[12]={0,0,0,0,0,0,0,0,0,0,0,0};
    max_handle = CCTK_ReductionArrayHandle ("maximum");
    if (max_handle<0) CCTK_WARN (0, "Could not obtain reduction handle");
    
    for (d=0; d<6; ++d) local[  d] =  cctkGH->cctk_bbox[d];
    for (d=0; d<3; ++d) local[6+d] = -cctkGH->cctk_lbnd[d];
    for (d=0; d<3; ++d) local[9+d] =  cctkGH->cctk_ubnd[d];
    ierr = CCTK_ReduceLocArrayToArray1D
      (cctkGH, -1, max_handle, local, global, 12, CCTK_VARIABLE_INT);
    for (d=0; d<6; ++d) global_bbox[d] =  global[  d];
    for (d=0; d<3; ++d) global_lbnd[d] = -global[6+d];
    for (d=0; d<3; ++d) global_ubnd[d] =  global[9+d];
#endif
    
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
    CCTK_REAL physical_min[3] = {0,0,0};
    CCTK_REAL physical_max[3] = {0,0,0};
    CCTK_REAL interior_min[3] = {0,0,0};
    CCTK_REAL interior_max[3] = {0,0,0};
    CCTK_REAL exterior_min[3] = {0,0,0};
    CCTK_REAL exterior_max[3] = {0,0,0};
    CCTK_REAL spacing = 0;
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
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "The coordinate origin in the %c-direction falls neither onto a grid point nor into the middle between two grid points.",
                  "xyz"[d]);
    }
    iorigin[d] = floor(origin[d] + (avoid_origin[d] ? 0.5 : 0.0) + 0.5);
    offset[d] = avoid_origin[d] ? 0 : 1;
  }
  
  /* x direction */
  dir = 0;
  otherdir = 1;
  
  assert (data.nghostzones[dir] >= 0);
  
  if (global_bbox[2*dir]) {

/*         printf("iorigin[otherdir] == %d\n", iorigin[otherdir]); */
/*     printf("offset[otherdir] == %d\n", offset[otherdir]); */
/*     printf("cctk_gsh[otherdir] == %d\n", cctk_gsh[otherdir]);  */

    if (2*iorigin[otherdir] + offset[otherdir] != cctk_gsh[otherdir]) {
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "The coordinate origin in the %c-direction is not in the centre of the domain.  The boundary condition cannot be applied.",
                  "xyz"[otherdir]);
    }
    
    if (iorigin[dir] != data.nghostzones[dir]) {
      assert (nvars > 0);
      var = 0;
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "The group \"%s\" has in the %c-direction %d symmetry points (grid points outside of the symmetry boundary).  "
                  "This is not equal to the number of ghost zones, which is %d.  "
                  "The number of symmetry points must be equal to the number of ghost zones.",
                  CCTK_GroupNameFromVarI(vis[var]), "xyz"[dir],
                  iorigin[dir], data.nghostzones[dir]);
    }
    
    if (data.gsh[dir] < 2*data.nghostzones[dir] + offset[dir]) {
      assert (nvars > 0);
      var = 0;
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
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
        
        int imin[3] = {0,0,0}, imax[3] = {0,0,0};
        int i=0, j=0, k=0;
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
              for (i=imin[0]; i<imax[0]; ++i) {
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
    xferinfo = pmalloc(group.dim * sizeof *xferinfo);
    assert (xferinfo);
    
    /* Fill in slab transfer description */
    for (d=0; d<group.dim; ++d) {
      xferinfo[d].src.gsh         = data.gsh        [d];
      xferinfo[d].src.lbnd        = data.lbnd       [d];
      xferinfo[d].src.lsh         = data.lsh        [d];
      xferinfo[d].src.lbbox       = fake_bbox       [2*d  ];
      xferinfo[d].src.ubbox       = fake_bbox       [2*d+1];
      xferinfo[d].src.nghostzones = data.nghostzones[d];
      xferinfo[d].src.off         = 0;
      xferinfo[d].src.str         = 1;
      xferinfo[d].src.len         = data.gsh        [d];
      
      xferinfo[d].dst.gsh         = data.gsh        [d];
      xferinfo[d].dst.lbnd        = data.lbnd       [d];
      xferinfo[d].dst.lsh         = data.lsh        [d];
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
    
    options = Util_TableCreateFromString ("useghosts=1");
    assert (options>=0);
    
    ierr = Slab_MultiTransfer
      (cctkGH, group.dim, xferinfo, options,
       nvars, vartypes, varptrs, vartypes, varptrs);
    assert (!ierr);
    
    ierr = Util_TableDestroy (options);
    assert (!ierr);
    
    if ( poison_boundaries ) {
      /* check destination grid points for poison */
      if (cctkGH->cctk_bbox[2*dir]) {
        
        int imin[3]={0,0,0}, imax[3] = {0,0,0};
        int i = 0, j=0, k=0;
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
          CCTK_REAL const * restrict const varptr = varptrs[var];
          CCTK_REAL poison;
          memset (&poison, poison_value, sizeof poison);
          for (k=imin[0]; k<imax[2]; ++k) {
            for (j=imin[1]; j<imax[1]; ++j) {
              for (i=imin[0]; i<imax[0]; ++i) {
                const int ind = CCTK_GFINDEX3D(cctkGH, i, j, k);
                if (memcmp (&varptr[ind], &poison, sizeof poison) == 0) {
                  printf ("   ijk=[%d,%d,%d]\n", i, j, k);
                  poison_found = 1;
                }
              }
            }
          }
          if (poison_found) {
            printf ("Poison found:\n");
            printf ("   levfac=[%d,%d,%d]\n", cctk_levfac[0], cctk_levfac[1], cctk_levfac[2]);
            printf ("   origin_space=[%g,%g,%g]\n", cctk_origin_space[0], cctk_origin_space[1], cctk_origin_space[2]);
            printf ("   delta_space=[%g,%g,%g]\n", cctk_delta_space[0], cctk_delta_space[1], cctk_delta_space[2]);
            printf ("   lbnd=[%d,%d,%d]\n", cctk_lbnd[0], cctk_lbnd[1], cctk_lbnd[2]);
            printf ("   lsh=[%d,%d,%d]\n", cctk_lsh[0], cctk_lsh[1], cctk_lsh[2]);
            printf ("   gsh=[%d,%d,%d]\n", cctk_gsh[0], cctk_gsh[1], cctk_gsh[2]);
            printf ("   bbox=[%d,%d,%d,%d,%d,%d]\n", cctk_bbox[0], cctk_bbox[1], cctk_bbox[2], cctk_bbox[3], cctk_bbox[4], cctk_bbox[5]);
            CCTK_WARN (CCTK_WARN_ABORT, "Poison found in symmetry regions -- there is an error in this thorn");
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
    
    /* take parity into account */
    if (cctkGH->cctk_bbox[2*dir]) {
      for (var=0; var<nvars; ++var) {
        int parity = 0;
        
        parity = paritiess[var][0] * paritiess[var][1];
        assert (abs(parity) == 1);
        if (parity == -1) {
          int imin[3]={0,0,0}, imax[3]={0,0,0};
          int i=0, j=0, k=0;
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
            assert (0);
          } /* switch grouptype */
        }
      } /* for var */
    } /* if bbox */
    
    free (xferinfo);
    
  } /* if global_bbox */
  
  free (gis);
  free (varptrs);
  free (vartypes);
  free (paritiess);
  
  return 0;
}



void Rot180_ApplyBC (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  int nvars = 0;
  CCTK_INT * restrict indices = NULL;
  CCTK_INT * restrict faces = NULL;
  CCTK_INT * restrict widths = NULL;
  CCTK_INT * restrict tables = NULL;
  CCTK_INT * restrict vartypes = NULL;
  int var = 0;
  int ierr = 0;
  int gap, i, j, ind_type_a, ind_type_b;
  
  assert (cctkGH);
  
  nvars = Boundary_SelectedGVs (cctkGH, 0, 0, 0, 0, 0, 0);
  assert (nvars>=0);
  
  if (nvars==0) return;
  
  indices = pmalloc (nvars * sizeof *indices);
  assert (indices);
  faces = pmalloc (nvars * sizeof *faces);
  assert (faces);
  widths = pmalloc (nvars * sizeof *widths);
  assert (widths);
  tables = pmalloc (nvars * sizeof *tables);
  assert (tables);
  vartypes = pmalloc (nvars * sizeof *vartypes);
  assert(vartypes);
  
  ierr =  Boundary_SelectedGVs
    (cctkGH, nvars, indices, faces, widths, tables, 0);
  assert (ierr == nvars);
  
  for (var=0; var<nvars; ++var) {
    assert (indices[var]>=0 && indices[var]<CCTK_NumVars());
    assert (widths[var] >= 0);
  }

  if(sort_vartypes) {
  /* sort variables into bins according to their type */
#   define SWAP(x,y) temp=x,x=y,y=temp
    for (var=0; var<nvars; ++var) {
      vartypes[var] = CCTK_VarTypeI(indices[var]);
    }
    for (gap = nvars/2; gap > 0; gap /= 2) { /* a shell-sort */
      for (i = gap; i < nvars; i++) {
        for (j=i-gap; j>=0 && vartypes[j]>vartypes[j+gap]; j-=gap) {
          CCTK_INT temp;

          SWAP(vartypes[j], vartypes[j+gap]);
          SWAP(indices[j], indices[j+gap]);
        }
      }
    }
#   undef SWAP     
  
    /* process each type of variables separately */
    /* ind_type_a is the first index of one type, ind_type_b is the first index
     * of a different type */
    for (ind_type_a = 0, ind_type_b = 1; ind_type_a < nvars; ind_type_a = ind_type_b) {
      /* search for end of block of one type */ 
      while(ind_type_b < nvars && vartypes[ind_type_b] == vartypes[ind_type_a]) {
        ind_type_b += 1;
      }

      /* ind_type_b points behind the block of common type => no +1 for length */
      ierr = BndRot180VI (cctkGH, ind_type_b - ind_type_a, &indices[ind_type_a]);
      assert (!ierr);
    }
  } else {
      ierr = BndRot180VI (cctkGH, nvars, indices);
      assert (!ierr);
  }
  
  free (indices);
  free (faces);
  free (widths);
  free (tables);
  free (vartypes);
}
