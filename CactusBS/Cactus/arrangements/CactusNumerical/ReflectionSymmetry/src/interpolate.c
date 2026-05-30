#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

#include "reflection.h"



CCTK_INT
ReflectionSymmetry_Interpolate (CCTK_POINTER_TO_CONST restrict const cctkGH,
                                CCTK_INT const N_dims,
                                CCTK_INT const local_interp_handle,
                                CCTK_INT const param_table_handle,
                                CCTK_INT const coord_system_handle,
                                CCTK_INT const N_interp_points,
                                CCTK_INT const interp_coords_type,
                                CCTK_POINTER_TO_CONST restrict const interp_coords[],
                                CCTK_INT const N_input_arrays,
                                CCTK_INT const input_array_indices[],
                                CCTK_INT const N_output_arrays,
                                CCTK_INT const output_array_types[],
                                CCTK_POINTER restrict const output_arrays[],
                                CCTK_INT const faces)
{
  DECLARE_CCTK_PARAMETERS;
  
  int do_reflection[6];
  
  CCTK_POINTER_TO_CONST restrict new_interp_coords[3];
  CCTK_INT newfaces;
  
  CCTK_INT * restrict operand_indices;
  CCTK_INT * restrict operation_codes;
  CCTK_INT * restrict output_array_indices;
  
  int dir;
  
  int iret;
  
  int m;
  int n;
  int d;
  
  int ierr;
  
  
  
  /* Get symmetry information */
  do_reflection[0] = reflection_x;
  do_reflection[1] = reflection_upper_x;
  do_reflection[2] = reflection_y;
  do_reflection[3] = reflection_upper_y;
  do_reflection[4] = reflection_z;
  do_reflection[5] = reflection_upper_z;
  
  newfaces = faces;
  for (d=0; d<6; ++d) {
    if (do_reflection[d]) {
      assert (newfaces & (1 << d));
      newfaces &= ~ (1 << d);
    }
  }
  
  
  
  /* Fold coordinates */
  assert (interp_coords_type == CCTK_VARIABLE_REAL);
  for (dir=0; dir<3; ++dir) {
    assert (! do_reflection[2*dir+1]);
    
    if (do_reflection[2*dir]) {
      CCTK_REAL * restrict coords_in_dir
        = malloc (N_interp_points * sizeof (CCTK_REAL));
      assert (N_interp_points == 0 || coords_in_dir);
      
      for (n=0; n<N_interp_points; ++n) {
        CCTK_REAL const pos = ((CCTK_REAL const *)interp_coords[dir])[n];
        CCTK_REAL const newpos = fabs(pos);
        coords_in_dir[n] = newpos;
      }
      new_interp_coords[dir] = coords_in_dir;
    } else {
      new_interp_coords[dir] = interp_coords[dir];
    }
  }
  
  
  
  /* Recursive call */
  iret = SymmetryInterpolateFaces
    (cctkGH, 
     N_dims, local_interp_handle, param_table_handle, coord_system_handle,
     N_interp_points, interp_coords_type, (CCTK_POINTER_TO_CONST)new_interp_coords,
     N_input_arrays, input_array_indices,
     N_output_arrays, output_array_types, (CCTK_POINTER_TO_CONST)output_arrays,
     newfaces);
  
  
  
  /* Free coordinates */
  for (dir=0; dir<3; ++dir) {
    if (do_reflection[2*dir]) {
      free ((void*)new_interp_coords[dir]);
    }
  }
  
  
  
  /* Find output variable indices */
  operand_indices = malloc (N_output_arrays * sizeof *operand_indices);
  assert (operand_indices);
  ierr = Util_TableGetIntArray
    (param_table_handle, N_output_arrays, operand_indices, "operand_indices");
  if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    assert (N_output_arrays == N_input_arrays);
    for (m=0; m<N_output_arrays; ++m) {
      operand_indices[m] = m;   /* set output index to input index */
    }
  } else {
    assert (ierr == N_output_arrays);
  }
  
  operation_codes = malloc (N_output_arrays * sizeof *operation_codes);
  assert (operation_codes);
  ierr = Util_TableGetIntArray
    (param_table_handle, N_output_arrays, operation_codes, "operation_codes");
  if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    assert (N_output_arrays == N_input_arrays);
    for (m=0; m<N_output_arrays; ++m) {
      operation_codes[m] = 0;     /* do not take derivatives */
    }
  } else {
    assert (ierr == N_output_arrays);
  }
  
  output_array_indices
    = malloc (N_output_arrays * sizeof *output_array_indices);
  assert (output_array_indices);
  for (m=0; m<N_output_arrays; ++m) {
    assert (operand_indices[m]>=0 && operand_indices[m]<N_input_arrays);
    output_array_indices[m] = input_array_indices[operand_indices[m]];
     assert (output_array_indices[m]==-1
             || (output_array_indices[m]>=0
                 && output_array_indices[m]<CCTK_NumVars()));
  }
  
  
  
  /* Unfold tensor types */
  for (m=0; m<N_output_arrays; ++m) {
    if (output_array_indices[m]!=-1) {
      
      int vi, gi;
      int numvars, firstvar, vectorlength;
      cGroup group;
      char * groupname;
      
      int table;
      char tensortypealias[1000];
      enum tensortype { UNKNOWN,
                        SCALAR, VECTOR, SYMTENSOR, SYMTENSOR3, TENSOR,
                        WEYLSCALARS_REAL, MANUALCARTESIAN };
      enum tensortype ttype;
      CCTK_INT tensorparity;
      int tcomponent;
      
      int parities[3];
      int check_dir[3];
      int needs_checking;

      
      vi = output_array_indices[m];
      assert (vi>=0 && vi<CCTK_NumVars());
      gi = CCTK_GroupIndexFromVarI (vi);
      assert (gi>=0 && gi<CCTK_NumGroups());

      ierr = CCTK_GroupData (gi, &group);
      assert (!ierr);
      
      numvars = CCTK_NumVarsInGroupI(gi);
      assert (numvars>0);
      firstvar = CCTK_FirstVarIndexI(gi);
      assert (firstvar>=0);
      vectorlength = group.vectorlength;
      assert (vectorlength>=0);
      assert (vectorlength==1 || group.vectorgroup);
      table = CCTK_GroupTagsTableI(gi);
      assert (table>=0);
  
      
      
      /* Get and check tensor type information */
      ierr = Util_TableGetString
        (table, sizeof tensortypealias, tensortypealias, "tensortypealias");
      if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
        groupname = CCTK_GroupName(gi);
        assert (groupname);
        CCTK_VWarn (4, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Tensor type alias not declared for group \"%s\" -- assuming a scalar",
                    groupname);
        free (groupname);
        strcpy (tensortypealias, "scalar");
      } else if (ierr<0) {
        groupname = CCTK_GroupName(gi);
        assert (groupname);
        CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Error in tensor type alias declaration for group \"%s\"",
                    groupname);
        free (groupname);
      }
      
      ttype = UNKNOWN;
      tcomponent = 0;
      if (CCTK_EQUALS (tensortypealias, "scalar")) {
        /* scalar */
        if (numvars != vectorlength) {
          groupname = CCTK_GroupName(gi);
          assert (groupname);
          CCTK_VWarn (4, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "Group \"%s\" has the tensor type alias \"scalar\", but contains more than 1 element",
                      groupname);
          free (groupname);
        }
        ttype = SCALAR;
        tcomponent = 0;
      } else if (CCTK_EQUALS (tensortypealias, "4scalar")) {
        /* 4-scalar */
        assert (numvars == vectorlength);
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
        /* Weyl scalars, stored as 10 real numbers */
        const int numcomps = 10;
        assert (numvars % numcomps == 0 && numvars == numcomps * vectorlength);
        ttype = WEYLSCALARS_REAL;
        tcomponent = (vi - firstvar) / vectorlength;
      } else if (CCTK_EQUALS (tensortypealias, "ManualCartesian")) {
        /* Reflection symmetries specified by hand */
        ttype = MANUALCARTESIAN;
      tcomponent = vi - firstvar;
      } else {
        groupname = CCTK_GroupName(gi);
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
        break;
      default:
        assert (0);
      }
      
      ierr = Util_TableGetInt (table, & tensorparity, "tensorparity");
      if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
        tensorparity = +1;
      } else if (ierr<0) {
        groupname = CCTK_GroupName(gi);
        assert (groupname);
        CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Error in tensor parity declaration for group \"%s\"",
                    groupname);
        free (groupname);
      }
      
      
      
      /* Calculate parities */
      parities[0] = parities[1] = parities[2] = +1;
      switch (ttype) {
      case SCALAR:
        /* do nothing */
        break;
      case VECTOR:
        parities[tcomponent] = -1;
        break;
      case SYMTENSOR:
        switch (tcomponent) {
        case 0: break;
        case 1: parities[0] = parities[1] = -1; break;
        case 2: parities[0] = parities[2] = -1; break;
        case 3: break;
        case 4: parities[1] = parities[2] = -1; break;
        case 5: break;
        default: assert (0);
        }
        break;
      case SYMTENSOR3:
        switch (tcomponent % 6) {
        case 0: break;
        case 1: parities[0] = parities[1] = -1; break;
        case 2: parities[0] = parities[2] = -1; break;
        case 3: break;
        case 4: parities[1] = parities[2] = -1; break;
        case 5: break;
        default: assert (0);
        }
        switch (tcomponent / 6) {
        case 0: parities[0] *= -1; break;
        case 1: parities[1] *= -1; break;
        case 2: parities[2] *= -1; break;
        default: assert (0);
        }
        break;
      case TENSOR:
        switch (tcomponent) {
        case 0: break;
        case 1: parities[0] = parities[1] = -1; break;
        case 2: parities[0] = parities[2] = -1; break;
        case 3: parities[1] = parities[0] = -1; break;
        case 4: break;
        case 5: parities[1] = parities[2] = -1; break;
        case 6: parities[2] = parities[0] = -1; break;
        case 7: parities[2] = parities[1] = -1; break;
        case 8: break;
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
        for (dir=0; dir<3; ++dir) {
          parities[dir] = weylparities[tcomponent][dir];
        }
        break;
      }
      case MANUALCARTESIAN:
        ReflectionSymmetry_GetManualParities(table, gi, parities);
        break;
      default:
        assert (0);
      }
      
      
      
      /* Take derivatives into account */
      {
        int code = operation_codes[m];
        while (code) {
          d = code % 10 - 1;
          code /= 10;
          assert (d>=0 && d<3);
          parities[d] *= -1;
        }
      }
      
      
      
      /* Are there negative parities? */
      needs_checking = 0;
      for (dir=0; dir<3; ++dir) {
        check_dir[dir] = do_reflection[2*dir] && parities[dir]*tensorparity < 0;
        needs_checking |= check_dir[dir];
      }

      
      
      
      /* Loop over all points and unfold */
      if (needs_checking) {
        for (n=0; n<N_interp_points; ++n) {
          int parity = 1; /* all types have their "natural" sign in the +++ octant */
          /* Is the point outside the domain? */
          for (dir=0; dir<3; ++dir) {
            if (check_dir[dir]) {
              CCTK_REAL const pos = ((CCTK_REAL const *)interp_coords[dir])[n];
              if (pos < 0) {
                /* Reflect the tensor component */
                parity *= -1; /* we only get here if parities[dir]*tensorparity < 0 */
              }
            }
          }
          ((CCTK_REAL *)output_arrays[m])[n] *= parity;
        }
      }
      
    }
    } /* for m */
  
  
  
  /* Free output variable indices */
  free (operand_indices);
  free (operation_codes);
  free (output_array_indices);
  
  
  
  /* Return */
  return iret;
}

void  ReflectionSymmetry_GetManualParities(int table, int gi, int *parities)
{
  char cartsyms[100];
  char *groupname = NULL;
  int i = 0;
  int ierr = -1;

  /* Get and check tensor type information */
  ierr = Util_TableGetString
    (table, sizeof cartsyms, cartsyms, "cartesianreflectionparities");
  if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    groupname = CCTK_GroupName(gi);
    assert (groupname);
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
		"Cartesian refection parity not declared for group \"%s\" -- aborting",
		groupname);
    assert(0);
  } else if (ierr<0) {
    groupname = CCTK_GroupName(gi);
    assert (groupname);
    CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
		"Error in tensor type alias declaration for group \"%s\"",
		groupname);
    free (groupname);
  }

  if (strlen(cartsyms) != 3)
  {
    CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
		"Invalid format for cartesianreflectionparities: must be xxx where x is + or - for group %s",
		groupname);
  }

  for (i = 0; i < 3; i++)
  {
    switch (cartsyms[i])
    {
    case '+':
      parities[i] = 1;
      break;
    case '-':
      parities[i] = -1;
      break;
    default:
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "Invalid format for cartesianreflectionparities: must be xxx where x is + or - for group %s",
		  groupname);
    }
  }
}
