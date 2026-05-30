 /*@@
   @file      LinearExtrapBnd.c
   @date      24 Jan 2003
   @author    David Rideout
   @desc
              Function which is registered as handling Carsten Gundlach's 
              "linear_extrap_one_bndry" boundary condition
   @enddesc
   @history
   @hdate     
   @hauthor   
   @hdesc     
   @endhistory
   @version   $Id$
 @@*/

#include "cctk.h"

#include "util_Table.h"

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusExamples_SampleBoundary_LinearExtrapBnd_c);

/* #define DEBUG */

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* External Routine Prototypes ******************
 ********************************************************************/

void CCTK_FCALL CCTK_FNAME(Linear_extrap_one_bndry)(int *doBC, const int *lsh, 
                                                    CCTK_REAL *var_ptr);

int BndLinExtrap (const cGH *GH, int num_vars, int *vars, int *faces, 
                  int *widths, int *tables);

/********************************************************************
 ********************    Scheduled Routines   ***********************
 ********************************************************************/


/********************************************************************
 ********************    External Routines   ************************
 ********************************************************************/

/*@@
   @routine    BndLinExtrap
   @date       24 Jan 2003
   @author     David Rideout
   @desc
               Apply linear extrapolation boundary condition to a
               group of grid functions given by their indices.
               This routine is registered to handle the linear
               extrapolation boundary condition.
               Can only handle 3D grid functions.

               All symmetries are ignored for now -- the symmetry bcs must 
               overwrite the output of this bc where necessary
   @enddesc

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        num_vars
   @vdesc      number of variables passed in through vars[]
   @vtype      int
   @vio        in
   @endvar
   @var        var_indices
   @vdesc      array of variable indicies to which to apply this boundary 
               condition
   @vtype      int *
   @vio        in
   @endvar
   @var        faces
   @vdesc      array of set of faces to which to apply the bc
   @vtype      int
   @vio        in
   @endvar
   @var        widths
   @vdesc      array of boundary widths for each variable
   @vtype      int
   @vio        in
   @endvar
   @var        table_handles
   @vdesc      array of table handles which hold extra arguments
   @vtype      int
   @vio        in
   @endvar

   @calls      CCTK_GroupIndexFromVarI
               CCTK_GroupDimI
               CCTK_VarTypeI

               Linear_extrap_one_bndry
   @history
   @hdate      
   @hauthor    
   @hdesc                     
   @endhistory

   @returntype int
   @returndesc
      each bit of return value indicates a possible error code, as follows:
       0  for success
      -1  invalid faces specification
      -2  unsupported staggering
      -4  boundary width != 1
      -8  potentially valid table handle
      -16 dimension is not supported
      -32 possibly called with a grid scalar
   @endreturndesc
@@*/

int BndLinExtrap (const cGH *GH, int num_vars, int *vars, int *faces, 
                  int *widths, int *tables)
{
  int i, j, gi, gtype, dim, retval, err;
  int doBC[6];
  const int *lsh, *bbox;
  CCTK_INT symtable;
  CCTK_INT symbnd[6];
  CCTK_INT is_physical[6];
  CCTK_REAL *var_ptr;
  cGroupDynamicData group_data;

  retval = 0;

#ifdef DEBUG
  printf("calling BndLinExtrap at iter %d\n", GH->cctk_iteration);
#endif

  /* loop through variables, one at a time (since this is all that 
     Linear_extrap_one_bndry can handle) */
  for (i=0; i<num_vars; ++i) {

    /* Check to see if faces specification is valid */
    if (faces[i] != CCTK_ALL_FACES)
    {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Faces specification %d for LinExtrap boundary conditions on "
                 "%s is not implemented yet.  "
                 "Not applying bc.", faces[i],
                 CCTK_VarName(vars[i]));
      retval |= 1;
    }

    /* Check to see if the boundary width might be something other than one */
    if (widths[i] != 1)
    {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "LinExtrapBnd does not handle boundary widths other than one."
                 "  Assuming boundary width of one on all faces.");
      retval |= 4;
    }

    /* Ignore table handles */
    if (tables[i] >= 0) 
    {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Possibly valid table handle.  LinExtrapBnd ignores "
                  "information stored in tables.");
      retval |= 8;
    }

    /* Gather some important information about this grid variable */
    gi = CCTK_GroupIndexFromVarI(vars[i]);
    gtype = CCTK_GroupTypeI(gi);
    if (gtype==CCTK_GF)
    {
      dim = GH->cctk_dim;
      bbox = GH->cctk_bbox;
      lsh = GH->cctk_lsh;
    }
    else
    {
      err = CCTK_GroupDynamicData(GH, gi, &group_data);
      if (err)
      {
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Error getting group information for group %s.  "
                   "Perhaps it is a grid scalar?", CCTK_GroupName(gi));
        retval |= 32;
      }
      dim = group_data.dim;
      bbox = group_data.bbox;
      lsh = group_data.lsh;
    }
    var_ptr = GH->data[vars[i]][0];
#ifdef DEBUG
    printf("dim=%d bbox[0]=%d lsh[1]=%d\n", dim, bbox[0], lsh[1]);
    printf("var_ptr=%p\n", var_ptr);
#endif

    /* Check dimension */
    if (dim != 3)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Variable dimension of %d not supported.  "
		  "Not applying bc.", dim);
      retval |= 16;
    }

    /* Decide on which faces the bc should be applied */
    for (j=0; j<2*dim; ++j)
    {
      doBC[j] = lsh[j/2] > widths[i]+2 && bbox[j];
    }

    /* see if we have a physical boundary */
    symtable = SymmetryTableHandleForGrid (GH);
    if (symtable < 0) CCTK_WARN (0, "internal error");
    err = Util_TableGetIntArray (symtable, 2 * dim, symbnd, "symmetry_handle");
    if (err != 2 * dim) CCTK_WARN (0, "internal error");
    for (j = 0; j < 2 * dim; j++)
    {
      is_physical[j] = symbnd[j] < 0;
    }

    /* Only do bc on faces without a symmetry bc */
    for (j=0; j<2*dim; ++j)
    {
      doBC[j] &= is_physical[i];
    }

    /* Apply the boundary condition */
    if( !( retval & ( 1 | 16 ) ) ) /* unless particularly bad errors */
      CCTK_FNAME(Linear_extrap_one_bndry)(doBC, lsh, var_ptr);
  }

  return -retval;
}

/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
