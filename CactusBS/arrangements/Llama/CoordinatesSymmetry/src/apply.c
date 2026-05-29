#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

#include "reflection.h"



#define COPY_PRE(VARTYPE)                                               \
  static void                                                           \
  copy_##VARTYPE (VARTYPE const * restrict const srcvar,                \
                  VARTYPE * restrict const dstvar,                      \
                  int const ni, int const nj, int const nk,             \
                  int const imin, int const jmin, int const kmin,       \
                  int const imax, int const jmax, int const kmax,       \
                  int const ioff, int const joff, int const koff,       \
                  int const idir, int const jdir, int const kdir,       \
                  int const parity)                                     \
  {                                                                     \
    assert (abs(idir)==1);                                              \
    assert (abs(jdir)==1);                                              \
    assert (abs(kdir)==1);                                              \
    int const iioff = ioff + (1 - idir) * imin;                         \
    int const jjoff = joff + (1 - jdir) * jmin;                         \
    int const kkoff = koff + (1 - kdir) * kmin;                         \
    int const iimin = iioff + idir * imin;                              \
    int const jjmin = jjoff + jdir * jmin;                              \
    int const kkmin = kkoff + kdir * kmin;                              \
    int const iimax = iioff + idir * imax;                              \
    int const jjmax = jjoff + jdir * jmax;                              \
    int const kkmax = kkoff + kdir * kmax;                              \
    assert (imin>=0 && imax<=ni);                                       \
    assert (jmin>=0 && jmax<=nj);                                       \
    assert (kmin>=0 && kmax<=nk);                                       \
    assert (iimin>=0 && iimax<=ni);                                     \
    assert (jjmin>=0 && jjmax<=nj);                                     \
    assert (kkmin>=0 && kkmax<=nk);                                     \
    assert (iimax>=-1 && iimin<ni);                                     \
    assert (jjmax>=-1 && jjmin<nj);                                     \
    assert (kkmax>=-1 && kkmin<nk);

#define COPY_LOOP(VARTYPE)                                              \
    for (int k=kmin; k<kmax; ++k) {                                     \
      for (int j=jmin; j<jmax; ++j) {                                   \
        for (int i=imin; i<imax; ++i) {                                 \
          int const dstind = i + ni * (j + nj * k);                     \
          int const ii = iioff + idir * i;                              \
          int const jj = jjoff + jdir * j;                              \
          int const kk = kkoff + kdir * k;                              \
          int const srcind = ii + ni * (jj + nj * kk);                  \
          dstvar[dstind] = parity * srcvar[srcind];                     \
        }                                                               \
      }                                                                 \
    }                                                                   \

#define COPY_POST(VARTYPE)                                              \
  }

#ifdef HAVE_CCTK_INT1
COPY_PRE(CCTK_INT1)
#pragma omp parallel for
COPY_LOOP(CCTK_INT1)
COPY_POST(CCTK_INT1)
#endif

#ifdef HAVE_CCTK_INT2
COPY_PRE(CCTK_INT2)
#pragma omp parallel for
COPY_LOOP(CCTK_INT2)
COPY_POST(CCTK_INT2)
#endif

#ifdef HAVE_CCTK_INT4
COPY_PRE(CCTK_INT4)
#pragma omp parallel for
COPY_LOOP(CCTK_INT4)
COPY_POST(CCTK_INT4)
#endif

#ifdef HAVE_CCTK_INT8
COPY_PRE(CCTK_INT8)
#pragma omp parallel for
COPY_LOOP(CCTK_INT8)
COPY_POST(CCTK_INT8)
#endif

#ifdef HAVE_CCTK_INT16
COPY_PRE(CCTK_INT16)
#pragma omp parallel for
COPY_LOOP(CCTK_INT16)
COPY_POST(CCTK_INT16)
#endif

#ifdef HAVE_CCTK_REAL4
COPY_PRE(CCTK_REAL4)
#pragma omp parallel for
COPY_LOOP(CCTK_REAL4)
COPY_POST(CCTK_REAL4)
#endif

#ifdef HAVE_CCTK_REAL8
COPY_PRE(CCTK_REAL8)
#pragma omp parallel for
COPY_LOOP(CCTK_REAL8)
COPY_POST(CCTK_REAL8)
#endif

#ifdef HAVE_CCTK_REAL16
COPY_PRE(CCTK_REAL16)
#pragma omp parallel for
COPY_LOOP(CCTK_REAL16)
COPY_POST(CCTK_REAL16)
#endif

#ifdef HAVE_CCTK_COMPLEX8
COPY_PRE(CCTK_COMPLEX8)
#pragma omp parallel for
COPY_LOOP(CCTK_COMPLEX8)
COPY_POST(CCTK_COMPLEX8)
#endif

#ifdef HAVE_CCTK_COMPLEX16
COPY_PRE(CCTK_COMPLEX16)
#pragma omp parallel for
COPY_LOOP(CCTK_COMPLEX16)
COPY_POST(CCTK_COMPLEX16)
#endif

#ifdef HAVE_CCTK_COMPLEX32
COPY_PRE(CCTK_COMPLEX32)
#pragma omp parallel for
COPY_LOOP(CCTK_COMPLEX32)
COPY_POST(CCTK_COMPLEX32)
#endif

#undef COPY_PRE
#undef COPY_LOOP
#undef COPY_POST



static int
BndReflectVI (cGH const * restrict const cctkGH,
              int const vi)
{
  DECLARE_CCTK_PARAMETERS;
  
  int gi;
  cGroup group;
  cGroupDynamicData data;
  int firstvar, numvars, vectorlength;
  char * restrict fullname;
  
  void * restrict varptr;
  
  int table;
  char tensortypealias[1000];
  enum tensortype { TT_UNKNOWN,
                    SCALAR, VECTOR, SYMTENSOR, SYMTENSOR3, TENSOR,
                    WEYLSCALARS_REAL, MANUALCARTESIAN };
  enum tensortype ttype;
  CCTK_INT tensorparity;
  int tcomponent;
  
  /* The stagger type is defined for grid variable groups, and
     depending on this type, each variable in the group may be
     staggered differently. This implies that groups with FACE or EDGE
     staggering have well-defined tensor types. */
  /* The stagger type is defined assuming a cell-centered grid, since
     this seems to be the common case where different stagger types
     are used. If the grid is vertex centered, then we currently abort
     when staggered grid variables are encountered. (We determine the
     centering of the grid via the avoid_origin_* parameters.) */
  char staggertype[1000];
  enum staggertype { ST_UNKNOWN, CELL, FACE, EDGE, VERTEX };
  enum staggertype stype;
  
  int do_reflection[6];
  int do_stagger_grid[6];
  
  int dir, face;
  
  int ash[3], imin[3], imax[3], ioff[3], idir[3];
  
  int parity;
  int do_stagger;
  int manual_parities[3];

  int d;
  
  int ierr;
  
  
  
  /* Check arguments */
  if (! cctkGH) {
    CCTK_WARN (0, "Argument cctkGH is NULL");
  }
  if (vi < 0 || vi >= CCTK_NumVars()) {
    CCTK_WARN (0, "Illegal variable index");
  }
  
  if (verbose) {
    fullname = CCTK_FullName (vi);
    if (! fullname)   {
      CCTK_WARN (0, "Internal error in CCTK_FullName");
    }
    CCTK_VInfo (CCTK_THORNSTRING,
                "Applying reflection boundary conditions to \"%s\"",
                fullname);
    free (fullname);
  }
  
  
  
  /* Get and check group information */
  gi = CCTK_GroupIndexFromVarI (vi);
  if (gi < 0 || gi > CCTK_NumGroups()) {
    CCTK_WARN (0, "Internal error in CCTK_GroupIndexFromVarI");
  }
  
  ierr = CCTK_GroupData (gi, &group);
  assert (!ierr);
  assert (group.grouptype == CCTK_GF);
  assert (group.disttype == CCTK_DISTRIB_DEFAULT);
  
  firstvar = CCTK_FirstVarIndexI (gi);
  assert (firstvar>=0 && firstvar<CCTK_NumVars());
  numvars = CCTK_NumVarsInGroupI (gi);
  assert (numvars>=0);
  vectorlength = group.vectorlength;
  assert (vectorlength>=0);
  assert (vectorlength==1 || group.vectorgroup);
  
  ierr = CCTK_GroupDynamicData (cctkGH, gi, &data);
  assert (!ierr);
  
  table = CCTK_GroupTagsTableI(gi);
  assert (table>=0);
  
  varptr = CCTK_VarDataPtrI (cctkGH, 0, vi);
  if (!varptr) {
    fullname = CCTK_FullName (vi);
    CCTK_VInfo (CCTK_THORNSTRING,
                "assertion for \"%s\"",
                fullname);
    free (fullname);
  }
  assert (varptr);
  
  
  
  /* Get and check tensor type information */
  ierr = Util_TableGetString
    (table, sizeof tensortypealias, tensortypealias, "tensortypealias");
  if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    /* assume a scalar */
    if (numvars != 1) {
      static int * restrict didwarn = 0;
      if (! didwarn) {
        didwarn = calloc (CCTK_NumGroups(), sizeof *didwarn);
      }
      if (! didwarn[gi]) {
        didwarn[gi] = 1;
        {
          char * groupname = CCTK_GroupName(gi);
          assert (groupname);
          CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "Group \"%s\" has no tensor type and contains more than one element -- treating these as \"scalar\"",
                      groupname);
          free (groupname);
        }
      }
    }
    strcpy (tensortypealias, "scalar");
  } else if (ierr<0) {
    char * groupname = CCTK_GroupName(gi);
    assert (groupname);
    CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Error in tensor type alias declaration for group \"%s\": %d",
                groupname, ierr);
    free (groupname);
  }
  
  ttype = TT_UNKNOWN;
  tcomponent = 0;
  if (CCTK_EQUALS (tensortypealias, "scalar")) {
    /* scalar */
    ttype = SCALAR;
    tcomponent = 0;
  } else if (CCTK_EQUALS (tensortypealias, "4scalar")) {
    /* 4-scalar */
    ttype = SCALAR;
    tcomponent = 0;
  } else if (CCTK_EQUALS (tensortypealias, "u")
             || CCTK_EQUALS (tensortypealias, "d"))
  {
    /* vector */
    const int numcomps = 3;
    /* special case to handle things like vel[3] */
    assert (numvars % numcomps == 0 && 
            (numvars == numcomps * vectorlength || numvars == vectorlength));
    ttype = VECTOR;
    if(numvars == vectorlength) {
      tcomponent = (vi - firstvar);
    } else {
      tcomponent = (vi - firstvar) / vectorlength;
    }
  } else if (CCTK_EQUALS (tensortypealias, "4u")
             || CCTK_EQUALS (tensortypealias, "4d"))
  {
    /* 4-vector */
    const int numcomps = 4;
    assert (numvars % numcomps == 0 && numvars == numcomps * vectorlength);
    if ((vi - firstvar) / vectorlength == 0) {
      ttype = SCALAR;
      tcomponent = 0;
    } else {
      ttype = VECTOR;
      tcomponent = (vi - firstvar) / vectorlength - 1;
    }
  } else if (CCTK_EQUALS (tensortypealias, "uu_sym")
             || CCTK_EQUALS (tensortypealias, "dd_sym"))
  {
    /* symmetric tensor */
    const int numcomps = 6;
    assert (numvars % numcomps == 0 && numvars == numcomps * vectorlength);
    ttype = SYMTENSOR;
    tcomponent = (vi - firstvar) / vectorlength;
  } else if (CCTK_EQUALS (tensortypealias, "uu")
             || CCTK_EQUALS (tensortypealias, "ud")
             || CCTK_EQUALS (tensortypealias, "du")
             || CCTK_EQUALS (tensortypealias, "dd"))
  {
    /* non-symmetric tensor */
    const int numcomps = 9;
    assert (numvars % numcomps == 0 && numvars == numcomps * vectorlength);
    ttype = TENSOR;
    tcomponent = (vi - firstvar) / vectorlength;
  } else if (CCTK_EQUALS (tensortypealias, "4uu_sym")
             || CCTK_EQUALS (tensortypealias, "4dd_sym"))
  {
    /* symmetric 4-tensor */
    const int numcomps = 10;
    assert (numvars % numcomps == 0 && numvars == numcomps * vectorlength);
    if ((vi - firstvar) / vectorlength == 0) {
      ttype = SCALAR;
      tcomponent = 0;
    } else if ((vi - firstvar) / vectorlength <= 3) {
      ttype = VECTOR;
      tcomponent = (vi - firstvar) / vectorlength - 1;
    } else {
      ttype = SYMTENSOR;
      tcomponent = (vi - firstvar) / vectorlength - 4;
    }
  } else if (CCTK_EQUALS (tensortypealias, "ddd_sym")) {
    /* 3rd rank tensor, symmetric in last 2 indices */
    const int numcomps = 18;
    assert (numvars % numcomps == 0 && numvars == numcomps * vectorlength);
    ttype = SYMTENSOR3;
    tcomponent = (vi - firstvar) / vectorlength;
  } else if (CCTK_EQUALS (tensortypealias, "weylscalars_real")) {
    /* Weyl scalars, stored as 10 real values.  NOTE: This assumes
       that Psi_0 comes first, which is NOT the default with
       PsiKadelia.  */
    const int numcomps = 10;
    assert (numvars % numcomps == 0 && numvars == numcomps * vectorlength);
    ttype = WEYLSCALARS_REAL;
    tcomponent = (vi - firstvar) / vectorlength;
  } else if (CCTK_EQUALS (tensortypealias, "ManualCartesian")) {
      /* Reflection symmetries specified by hand */
      ttype = MANUALCARTESIAN;
      tcomponent = vi - firstvar;
  } else {
    char * groupname = CCTK_GroupName(gi);
    assert (groupname);
    CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Illegal tensor type alias \"%s\" for group \"%s\"",
                tensortypealias, groupname);
    free (groupname);
  }
  
  switch (ttype) {
  case SCALAR:
    assert (tcomponent>=0 && tcomponent<1);
    break;
  case VECTOR:
    assert (tcomponent>=0 && tcomponent<3);
    break;
  case SYMTENSOR:
    assert (tcomponent>=0 && tcomponent<6);
    break;
  case SYMTENSOR3:
    assert (tcomponent>=0 && tcomponent<18);
    break;
  case TENSOR:
    assert (tcomponent>=0 && tcomponent<9);
    break;
  case WEYLSCALARS_REAL:
    assert (tcomponent>=0 && tcomponent<10);
    break;
  case MANUALCARTESIAN:
    /* No restriction on number of components */
    CoordinatesSymmetry_GetManualParities(table, gi, manual_parities);
    break;

  default:
    assert (0);
  }
  
  ierr = Util_TableGetInt (table, & tensorparity, "tensorparity");
  if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    tensorparity = +1;
  } else if (ierr<0) {
    char * groupname = CCTK_GroupName(gi);
    assert (groupname);
    CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Error in tensor parity declaration for group \"%s\": %d",
                groupname, ierr);
    free (groupname);
  }
  
  
  
  /* Get and check stagger type information */
  ierr = Util_TableGetString
    (table, sizeof staggertype, staggertype, "staggertype");
  if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    /* assume cell centering */
    strcpy (staggertype, "cell");
  } else if (ierr<0) {
    char * groupname = CCTK_GroupName(gi);
    assert (groupname);
    CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Error in stagger type declaration for group \"%s\": %d",
                groupname, ierr);
    free (groupname);
  }
  
  stype = ST_UNKNOWN;
  if (CCTK_EQUALS (staggertype, "cell")) {
    /* cell */
    stype = CELL;
  } else if (CCTK_EQUALS (staggertype, "face")) {
    /* face */
    stype = FACE;
  } else if (CCTK_EQUALS (staggertype, "edge")) {
    /* edge */
    stype = EDGE;
  } else if (CCTK_EQUALS (staggertype, "vertex")) {
    /* vertex */
    stype = VERTEX;
  } else {
    char * groupname = CCTK_GroupName(gi);
    assert (groupname);
    CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Illegal stagger type \"%s\" for group \"%s\"",
                staggertype, groupname);
    free (groupname);
  }
  
  
  
  /* Reflection symmetry information */
  
  const int map = MultiPatch_GetMap(cctkGH);
  int type[6];
  int sym_dir[6];
  int overlap;
  ierr += MultiPatch_GetSymmetrySpecification(map, 6, type, &overlap, sym_dir);
  
  /* extract symmetry direction */
  
  for (dir=0; dir<3; ++dir) {
     for (face=0; face<2; ++face) {
        const int real_dir = sym_dir[2*dir+face];
        do_reflection[2*dir+face] = (reflection_x && real_dir == 0) || (reflection_y && real_dir == 1) || (reflection_z && real_dir == 2);
        if (do_reflection[2*dir+face] && type[2*dir+face] != 1) {
           CCTK_WARN(0, "Reflection symmetry cannot be applied! Patch system not compatible with that!");
        }
     }
  }
  
  /*
  do_reflection[0] = type[0] == 1 ? 1 : 0;
  do_reflection[1] = type[1] == 1 ? 1 : 0;
  do_reflection[2] = type[2] == 1 ? 1 : 0;
  do_reflection[3] = type[3] == 1 ? 1 : 0;
  do_reflection[4] = type[4] == 1 ? 1 : 0;
  do_reflection[5] = type[5] == 1 ? 1 : 0;
  */
  
  do_stagger_grid[0] = stagger;
  do_stagger_grid[1] = stagger;
  do_stagger_grid[2] = stagger;
  do_stagger_grid[3] = stagger;
  do_stagger_grid[4] = stagger;
  do_stagger_grid[5] = stagger;
    
  int bbox[6];
  ierr += MultiPatch_GetBbox(cctkGH, 6, bbox);
  
  /* Loop over all directions and faces */
  for (dir=0; dir<3; ++dir) {
    for (face=0; face<2; ++face) {
      /* If there is a reflection symmetry on that face */
      if (do_reflection[2*dir+face]) {
        /* If we have the outer boundary of that face */
        if (cctkGH->cctk_bbox[2*dir+face] && bbox[2*dir+face]) {
          
          const int real_dir = sym_dir[2*dir+face];
          
          /* Find parity */
          parity = tensorparity;
          switch (ttype) {
          case SCALAR:
            parity *= +1;
            break;
          case VECTOR:
            parity *= real_dir == tcomponent ? -1 : +1;
            break;
          case SYMTENSOR:
            switch (tcomponent) {
            case 0: parity *= +1; break;
            case 1: parity *= (real_dir == 0 || real_dir == 1) ? -1 : +1; break;
            case 2: parity *= (real_dir == 0 || real_dir == 2) ? -1 : +1; break;
            case 3: parity *= +1; break;
            case 4: parity *= (real_dir == 1 || real_dir == 2) ? -1 : +1; break;
            case 5: parity *= +1; break;
            default: assert (0);
            }
            break;
          case SYMTENSOR3:
            switch (tcomponent % 6) {
            case 0: parity *= +1; break;
            case 1: parity *= (real_dir == 0 || real_dir == 1) ? -1 : +1; break;
            case 2: parity *= (real_dir == 0 || real_dir == 2) ? -1 : +1; break;
            case 3: parity *= +1; break;
            case 4: parity *= (real_dir == 1 || real_dir == 2) ? -1 : +1; break;
            case 5: parity *= +1; break;
            default: assert (0);
            }
            switch (tcomponent / 6) {
            case 0: parity *= real_dir == 0 ? -1 : +1; break;
            case 1: parity *= real_dir == 1 ? -1 : +1; break;
            case 2: parity *= real_dir == 2 ? -1 : +1; break;
            default: assert (0);
            }
            break;
          case TENSOR:
            switch (tcomponent) {
            case 0: parity *= +1; break;
            case 1: parity *= (real_dir == 0 || real_dir == 1) ? -1 : +1; break;
            case 2: parity *= (real_dir == 0 || real_dir == 2) ? -1 : +1; break;
            case 3: parity *= (real_dir == 1 || real_dir == 0) ? -1 : +1; break;
            case 4: parity *= +1; break;
            case 5: parity *= (real_dir == 1 || real_dir == 2) ? -1 : +1; break;
            case 6: parity *= (real_dir == 2 || real_dir == 0) ? -1 : +1; break;
            case 7: parity *= (real_dir == 2 || real_dir == 1) ? -1 : +1; break;
            case 8: parity *= +1; break;
            default: assert (0);
            }
            break;
          case WEYLSCALARS_REAL: {
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
            parity *= weylparities[tcomponent][real_dir];
            break;
          }
          case MANUALCARTESIAN:
            parity = manual_parities[real_dir];
            break;
          default:
            assert (0);
          }
          
          /* Find staggering */
          do_stagger = do_stagger_grid[2*dir+face];
          switch (stype) {
          case CELL:
            /* do nothing */
            break;
          case FACE:
            assert (face == 0);  /* assume lower face */
            assert (do_stagger); /* assume cell-centered grid */
            assert (ttype == VECTOR); /* TODO: support other tensor types */
            if (dir == tcomponent) do_stagger = !do_stagger;
            break;
          case EDGE:
            assert (face == 0);  /* assume lower face */
            assert (do_stagger); /* assume cell-centered grid */
            assert (ttype == VECTOR); /* TODO: support other tensor types */
            if (dir != tcomponent) do_stagger = !do_stagger;
            break;
          case VERTEX:
            assert (face == 0);  /* assume lower face */
            assert (do_stagger); /* assume cell-centered grid */
            do_stagger = !do_stagger;
            break;
          default:
            assert (0);
          }
          
          /* Find region extent */
          for (d=0; d<3; ++d) {
            ash[d] = cctkGH->cctk_ash[d];
            imin[d] = 0;
            imax[d] = cctkGH->cctk_lsh[d];
            ioff[d] = 0;
            idir[d] = 1;
          }
          /* To determine the correct expressions for ioff below,
           * consider the definitions of iimin and iimax in the macro
           * COPY_PRE above:
           *
           *    iioff = ioff + (1 - idir) * imin
           *    iimin = iioff + idir * imin
           *    iimax = iioff + idir * imax
           *
           * We need idir = -1 since the loop copying into the
           * symmetry zones traverses the grid points in opposite
           * directions.
           *
           * Using the expressions above, one chooses the desired
           * values of iimin and iimax, and can then solve the
           * equations above for ioff.
           */
          if (face == 0) {
            imax[dir] = cctkGH->cctk_nghostzones[dir];
            ioff[dir] = (+ 2*(cctkGH->cctk_nghostzones[dir]) - 1 
                         + (do_stagger ? 0 : 1));
            idir[dir] = -1;
          } else {
            imin[dir] = cctkGH->cctk_lsh[dir] - (cctkGH->cctk_nghostzones[dir]);
            ioff[dir] = - 1 - (do_stagger ? 0 : 1); 
            idir[dir] = -1;
          }
          
          /* Ensure that there are sufficient interior zones, since
             this thorn does not support filling symmetry zones from
             other symmetry zones */
          {
            int const have_points = cctkGH->cctk_gsh[dir];
            int const need_points =
              3 * (cctkGH->cctk_nghostzones[dir])
              + !do_stagger_grid[2*dir] + !do_stagger_grid[2*dir+1];
            if (need_points > have_points) {
              CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                          "Cannot apply symmetry boundary zones in the %s %c direction, since there seem to be more symmetry zones than interior zones",
                          (face==0 ? "lower" : "upper"),
                          "xyz"[dir]);
            }
          }
          
          /* Copy region */
          switch (group.vartype) {
            
#define ARGS  varptr, varptr,                   \
              ash[0], ash[1], ash[2],           \
              imin[0], imin[1], imin[2],        \
              imax[0], imax[1], imax[2],        \
              ioff[0], ioff[1], ioff[2],        \
              idir[0], idir[1], idir[2],        \
              parity

#ifdef HAVE_CCTK_INT1
          case CCTK_VARIABLE_INT1:
#ifdef CCTK_INTEGER_PRECISION_1
          case CCTK_VARIABLE_INT:
#endif
            copy_CCTK_INT1 (ARGS);
            break;
#endif
            
#ifdef HAVE_CCTK_INT2
          case CCTK_VARIABLE_INT2:
#ifdef CCTK_INTEGER_PRECISION_2
          case CCTK_VARIABLE_INT:
#endif
            copy_CCTK_INT2 (ARGS);
            break;
#endif
            
#ifdef HAVE_CCTK_INT4
          case CCTK_VARIABLE_INT4:
#ifdef CCTK_INTEGER_PRECISION_4
          case CCTK_VARIABLE_INT:
#endif
            copy_CCTK_INT4 (ARGS);
            break;
#endif
            
#ifdef HAVE_CCTK_INT8
          case CCTK_VARIABLE_INT8:
#ifdef CCTK_INTEGER_PRECISION_8
          case CCTK_VARIABLE_INT:
#endif
            copy_CCTK_INT8 (ARGS);
            break;
#endif
            
#ifdef HAVE_CCTK_INT16
          case CCTK_VARIABLE_INT16:
#ifdef CCTK_INTEGER_PRECISION_16
          case CCTK_VARIABLE_INT:
#endif
            copy_CCTK_INT16 (ARGS);
            break;
#endif
            
#ifdef HAVE_CCTK_REAL4
          case CCTK_VARIABLE_REAL4:
#ifdef CCTK_REAL_PRECISION_4
          case CCTK_VARIABLE_REAL:
#endif
            copy_CCTK_REAL4 (ARGS);
            break;
#endif
            
#ifdef HAVE_CCTK_REAL8
          case CCTK_VARIABLE_REAL8:
#ifdef CCTK_REAL_PRECISION_8
          case CCTK_VARIABLE_REAL:
#endif
            copy_CCTK_REAL8 (ARGS);
            break;
#endif
            
#ifdef HAVE_CCTK_REAL16
          case CCTK_VARIABLE_REAL16:
#ifdef CCTK_REAL_PRECISION_16
          case CCTK_VARIABLE_REAL:
#endif
            copy_CCTK_REAL16 (ARGS);
            break;
#endif
            
#ifdef HAVE_CCTK_COMPLEX8
          case CCTK_VARIABLE_COMPLEX8:
#ifdef CCTK_COMPLEX_PRECISION_8
          case CCTK_VARIABLE_COMPLEX:
#endif
            copy_CCTK_COMPLEX8 (ARGS);
            break;
#endif
            
#ifdef HAVE_CCTK_COMPLEX16
          case CCTK_VARIABLE_COMPLEX16:
#ifdef CCTK_COMPLEX_PRECISION_16
          case CCTK_VARIABLE_COMPLEX:
#endif
            copy_CCTK_COMPLEX16 (ARGS);
            break;
#endif
            
#ifdef HAVE_CCTK_COMPLEX32
          case CCTK_VARIABLE_COMPLEX32:
#ifdef CCTK_COMPLEX_PRECISION_32
          case CCTK_VARIABLE_COMPLEX:
#endif
            copy_CCTK_COMPLEX32 (ARGS);
            break;
#endif
            
#undef ARGS
            
          default:
            CCTK_WARN (0, "Unsupported variable type");
          }
          
        } /* if cctk_bbox */
      } /* if do_reflection */
    } /* for face */
  } /* for dir */
  
  /* Success */
  return 0;
}



/* When CoordBase is used to specify the location of the boundary
   points, then ensure that the CoordBase parameters and this thorn's
   parameters are consistent.  */
static
void
CheckBoundaryParameters (cGH const * restrict const cctkGH,
                         int const vi,
                         int const * restrict const stencil)
{
  DECLARE_CCTK_PARAMETERS;
  
  static int did_check = 0;
  
  int type;                     /* Parameter type */
  void const * ptr;             /* Pointer to parameter value */
  char const * coordtype;       /* CartGrid3D::type */
  
  int dim;                      /* Number of dimensions of vi */
  
  CCTK_INT * restrict nboundaryzones; /* CoordBase boundary location */
  CCTK_INT * restrict is_internal;
  CCTK_INT * restrict is_staggered;
  CCTK_INT * restrict shiftout;
  
  int do_reflection[6];         /* This thorn's parameters */
  int do_stagger_grid[6];
  
  int d;
  int ierr;
  
  
  
  /* Check only once to save time */
  // TODO: This doesn't work with multiple maps!
  //if (did_check) return;
  
  /* Check only for grid functions */
  if (CCTK_GroupTypeFromVarI (vi) != CCTK_GF) return;
  
  did_check = 1;
  
  /* Check whether Coordinates is active */
  if (! CCTK_IsThornActive ("CartGrid3D")) return;
  
  /* Check whether CoordBase is used */
  ptr = CCTK_ParameterGet ("type", "CartGrid3D", & type);
  assert (ptr != 0);
  assert (type == PARAMETER_KEYWORD);
  coordtype = * (char const * const *) ptr;
  if (! CCTK_EQUALS (coordtype, "multipatch")) return;
  
  /* Get the boundary specification */
  dim = CCTK_GroupDimFromVarI (vi);
  assert (dim >= 0);
  nboundaryzones = malloc (2*dim * sizeof *nboundaryzones);
  is_internal = malloc (2*dim * sizeof *is_internal);
  is_staggered = malloc (2*dim * sizeof *is_staggered);
  shiftout = malloc (2*dim * sizeof *shiftout);
  const int map = MultiPatch_GetMap(cctkGH);
  ierr = MultiPatch_GetBoundarySpecification
    (map, 2*dim, nboundaryzones, is_internal, is_staggered, shiftout);
  assert (! ierr);
  
  /* Reflection symmetry information */
  assert(dim == 3);
  int btype[6];
  int sym_dir[6];
  int overlap;
  ierr += MultiPatch_GetSymmetrySpecification(map, 6, btype, &overlap, sym_dir);
  
  do_reflection[0] = btype[0] == 1 ? 1 : 0;
  do_reflection[1] = btype[1] == 1 ? 1 : 0;
  do_reflection[2] = btype[2] == 1 ? 1 : 0;
  do_reflection[3] = btype[3] == 1 ? 1 : 0;
  do_reflection[4] = btype[4] == 1 ? 1 : 0;
  do_reflection[5] = btype[5] == 1 ? 1 : 0;
  
  do_stagger_grid[0] = stagger;
  do_stagger_grid[1] = stagger;
  do_stagger_grid[2] = stagger;
  do_stagger_grid[3] = stagger;
  do_stagger_grid[4] = stagger;
  do_stagger_grid[5] = stagger;
  
  /* Check the boundary sizes */
  for (d=0; d<6; ++d) {
    if (do_reflection[d]) {
      char const dir = "xyz"[d/2];
      char const * const face = (d%2==0) ? "lower" : "upper";
      if (stencil[d/2] != nboundaryzones[d]) {
        CCTK_VWarn (CCTK_WARN_ABORT,
                    __LINE__, __FILE__, CCTK_THORNSTRING,
                    "The %s %c face is a symmetry boundary.  Since there are %d ghost zones in the %c direction, the corresponding Coordinates boundary width must also be %d.  The boundary width is currently %d.",
                    face, dir,
                    stencil[d/2], dir, stencil[d/2],
                    (int) nboundaryzones[d]);
      }
      if (is_internal[d]) {
        CCTK_VWarn (CCTK_WARN_ABORT,
                    __LINE__, __FILE__, CCTK_THORNSTRING,
                    "The %s %c face is a symmetry boundary.  The corresponding Coordinates boundary must not be internal.",
                    face, dir);
      }
      if (do_stagger_grid[d] != is_staggered[d]) {
        CCTK_VWarn (CCTK_WARN_ABORT,
                    __LINE__, __FILE__, CCTK_THORNSTRING,
                    "The %s %c face of map %d is a symmetry boundary.  The symmetry condition and the corresponding Coordinates boundary must either be both staggered or both not staggered.",
                    face, dir, map);
      }
      /*if ((do_stagger_grid[d] ? 0 : 1) != shiftout[d]) {
        CCTK_VWarn (CCTK_WARN_ABORT,
                    __LINE__, __FILE__, CCTK_THORNSTRING,
                    "The %s %c face is a symmetry boundary.  If the symmetry condition is staggered, then the corresponding CoordBase shiftout must be 0; otherwise it must be 1.",
                    face, dir);
      }*/
    }
  }
  
  /* Free memory */
  free (nboundaryzones);
  free (is_internal);
  free (is_staggered);
  free (shiftout);
}



void
CoordinatesSymmetry_Apply (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  int nvars;
  CCTK_INT * restrict indices;
  CCTK_INT * restrict faces;
  CCTK_INT * restrict widths;
  CCTK_INT * restrict tables;
  int vi;
  int dim;
  int * restrict stencil;
  int i;
  int istat;
  int ierr;
  
  if (verbose) {
     int reflevel = GetRefinementLevel(cctkGH);
     int map = MultiPatch_GetMap(cctkGH);
     
     printf("Applying reflecting BCs on map %d reflevel %d\n", map, reflevel);
  }
  
  nvars = Boundary_SelectedGVs (cctkGH, 0, NULL, NULL, NULL, NULL, NULL);
  if (nvars < 0) {
    CCTK_WARN (0, "Internal error in Boundary_SelectedGVs");
  }
  
  if (nvars == 0) {
    /* Nothing to do */
    return;
  }
  
  indices = malloc (nvars * sizeof *indices);
  if (! indices) {
    CCTK_WARN (0, "Out of memory");
  }
  faces = malloc (nvars * sizeof *faces);
  if (! faces) {
    CCTK_WARN (0, "Out of memory");
  }
  widths = malloc (nvars * sizeof *widths);
  if (! widths) {
    CCTK_WARN (0, "Out of memory");
  }
  tables = malloc (nvars * sizeof *tables);
  if (! tables) {
    CCTK_WARN (0, "Out of memory");
  }
  
  istat =  Boundary_SelectedGVs
    (cctkGH, nvars, indices, faces, widths, tables, 0);
  if (istat != nvars) {
    CCTK_WARN (0, "Internal error in Boundary_SelectedGVs");
  }
  
  for (i=0; i<nvars; ++i) {
    vi = indices[i];
    if (vi < 0 || vi >= CCTK_NumVars()) {
      CCTK_WARN (0, "Illegal variable index");
    }
    
    if (widths[i] < 0) {
      CCTK_WARN (0, "Illegal boundary width");
    }
    
    dim = CCTK_GroupDimFromVarI (vi);
    if (dim < 0) {
      CCTK_WARN (0, "Illegal dimension");
    }
    
    stencil = malloc (dim * sizeof *stencil);
    if (! stencil) {
      CCTK_WARN (0, "Out of memory");
    }
    ierr = CCTK_GroupnghostzonesVI (cctkGH, dim, stencil, vi);
    if (ierr) {
      CCTK_WARN (0, "Internal error in CCTK_GroupnghostzonesVI");
    }
    
    CheckBoundaryParameters (cctkGH, vi, stencil);
    
    ierr = BndReflectVI (cctkGH, vi);
    if (ierr) {
      CCTK_WARN (0, "Internal error in BndReflectVI");
    }
    
    free (stencil);
  }
  
  free (indices);
  free (faces);
  free (widths);
  free (tables);
}
