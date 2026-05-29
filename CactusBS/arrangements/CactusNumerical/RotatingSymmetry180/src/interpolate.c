#include "rotatingsymmetry180.h"

#include <cctk.h>
#include <util_ErrorCodes.h>
#include <util_Table.h>

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>



void
Rot180_CheckTensorTypes (CCTK_ARGUMENTS)
{
  /* Check tensor types of all groups */
  for (int gi=0; gi<CCTK_NumGroups(); ++gi) {
    
    char tensortypealias[100];
    int numvars, firstvar, vectorlength;
    int table;
    int ierr;
    cGroup group;
    
    numvars = CCTK_NumVarsInGroupI(gi);
    if (numvars == 0) continue;
    assert (numvars>0);
    firstvar = CCTK_FirstVarIndexI(gi);
    assert (firstvar>=0);
    ierr = CCTK_GroupData (gi, &group);
    assert (!ierr);
    vectorlength = group.vectorlength;
    table = CCTK_GroupTagsTableI(gi);
    assert (table>=0);
    
    ierr = Util_TableGetString
      (table, sizeof tensortypealias, tensortypealias, "tensortypealias");
    if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
      char * groupname = CCTK_GroupName(gi);
      assert (groupname);
      CCTK_VWarn (4, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Tensor type alias not declared for group \"%s\" -- assuming a scalar",
                  groupname);
      free (groupname);
      strcpy (tensortypealias, "");
    } else if (ierr<0) {
      char * groupname = CCTK_GroupName(gi);
      assert (groupname);
      CCTK_VWarn (CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Error in tensor type alias declaration for group \"%s\"",
                  groupname);
      free (groupname);
    }
    
    if (CCTK_EQUALS (tensortypealias, "")) {
      /* do nothing */
    } else if (CCTK_EQUALS (tensortypealias, "scalar")) {
      /* scalar */
      if (numvars != 1) {
        char * groupname = CCTK_GroupName(gi);
        assert (groupname);
        CCTK_VWarn (4, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Group \"%s\" has the tensor type alias \"scalar\", but contains more than 1 element",
                    groupname);
        free (groupname);
      }
    } else if (CCTK_EQUALS (tensortypealias, "4scalar")) {
      /* 4-scalar */
      if (numvars != 1) {
        char * groupname = CCTK_GroupName(gi);
        assert (groupname);
        CCTK_VWarn (4, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Group \"%s\" has the tensor type alias \"4scalar\", but contains more than 1 element",
                    groupname);
        free (groupname);
      }
    } else if (CCTK_EQUALS (tensortypealias, "u")
               || CCTK_EQUALS (tensortypealias, "d")) {
      /* vector, special case */
      if (!((numvars == 3*vectorlength) ||
            (numvars == 3 && vectorlength == 3))) {
        char * groupname = CCTK_GroupName(gi);
        assert (groupname);
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Group \"%s\" has the tensor type alias \"u\", but does not contain a multiple of 3 elements (but %d)",
                    groupname, numvars);
        free (groupname);
      }
      assert ((numvars == 3*vectorlength) ||
              (numvars == 3 && vectorlength == 3));
    } else if (CCTK_EQUALS (tensortypealias, "4u")
               || CCTK_EQUALS (tensortypealias, "4d")) {
      /* 4-vector */
      assert ((numvars == 4*vectorlength) ||
              (numvars == 4 && vectorlength == 4));
    } else if (CCTK_EQUALS (tensortypealias, "uu_sym")
               || CCTK_EQUALS (tensortypealias, "dd_sym")) {
      /* symmetric tensor */
      assert (numvars == 6*vectorlength);
    } else if (CCTK_EQUALS (tensortypealias, "uu")
               || CCTK_EQUALS (tensortypealias, "ud")
               || CCTK_EQUALS (tensortypealias, "du")
               || CCTK_EQUALS (tensortypealias, "dd")) {
      /* tensor */
      assert (numvars == 9*vectorlength);
    } else if (CCTK_EQUALS (tensortypealias, "ddd_sym")) {
      /* 3rd rank tensor, symmetric in last 2 indices */
      assert (numvars == 18*vectorlength);
    } else if (CCTK_EQUALS (tensortypealias, "4uu_sym")
               || CCTK_EQUALS (tensortypealias, "4dd_sym")) {
      /* symmetric 4-tensor */
      assert (numvars == 10*vectorlength);
    } else if (CCTK_EQUALS (tensortypealias, "weylscalars_real")) {
      /* Weyl scalars, stored as 10 real values */
      assert (numvars == 10*vectorlength);
    } else if (CCTK_EQUALS (tensortypealias, "ManualCartesian")) {
      /* No restriction */
    } else {
      char * groupname = CCTK_GroupName(gi);
      assert (groupname);
      CCTK_VWarn (CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Illegal tensor type alias \"%s\" for group \"%s\"",
                  tensortypealias, groupname);
      free (groupname);
    }
    
  }
}



/* Symmetry interpolation */
CCTK_INT
Rot180_SymmetryInterpolate (CCTK_POINTER_TO_CONST const cctkGH,
                            CCTK_INT const N_dims,
                            CCTK_INT const local_interp_handle,
                            CCTK_INT const param_table_handle,
                            CCTK_INT const coord_system_handle,
                            CCTK_INT const N_interp_points,
                            CCTK_INT const interp_coords_type,
                            CCTK_POINTER_TO_CONST const interp_coords[],
                            CCTK_INT const N_input_arrays,
                            CCTK_INT const input_array_indices[],
                            CCTK_INT const N_output_arrays,
                            CCTK_INT const output_array_types[],
                            CCTK_POINTER const output_arrays[],
                            CCTK_INT const faces)
{
  CCTK_POINTER new_interp_coords[3];
  CCTK_INT new_faces;
  CCTK_INT * restrict operand_indices;
  CCTK_INT * restrict operation_codes;
  CCTK_INT * restrict output_array_indices;
  int parities[3];              /* variable parities */
  int m;                        /* output array */
  int n;                        /* point */
  int d;                        /* dimension */
  int iret;                     /* interpolator return value */
  int ierr;
  
  /* Check arguments */
  assert (N_dims == 3);
  assert (N_interp_points >= 0);
  assert (interp_coords_type >= 0);
  for (d=0; d<3; ++d) {
    assert (N_interp_points == 0 || interp_coords[d]);
  }
  assert (N_output_arrays >= 0);
  
  /* Coordinates must be CCTK_REAL */
  assert (interp_coords_type == CCTK_VARIABLE_REAL);
  for (m=0; m<N_output_arrays; ++m) {
    assert (output_array_types[m] == CCTK_VARIABLE_REAL);
  }
  
  
  
  /* Claim faces */
  assert (faces & (1 << 0));
  new_faces = faces;
  new_faces &= ~ (1 << 0);
  
  /* Copy coordinates */
  for (d=0; d<3; ++d) {
    new_interp_coords[d] = malloc (N_interp_points * sizeof(CCTK_REAL));
    assert (N_interp_points == 0 || new_interp_coords[d]);
    memcpy (new_interp_coords[d], interp_coords[d],
            N_interp_points * sizeof(CCTK_REAL));
  }
  
  /* Fold coordinates */
  for (n=0; n<N_interp_points; ++n) {
    /* Is the point outside the domain? */
    if (((CCTK_REAL *)new_interp_coords[0])[n] < 0) {
      /* Rotate the point by 180 degrees */
      ((CCTK_REAL *)new_interp_coords[0])[n] *= -1;
      ((CCTK_REAL *)new_interp_coords[1])[n] *= -1;
    }
  }
  
  
  
  /* Recursive call */
  iret = SymmetryInterpolateFaces
    (cctkGH, N_dims,
     local_interp_handle, param_table_handle, coord_system_handle,
     N_interp_points, interp_coords_type, new_interp_coords,
     N_input_arrays, input_array_indices,
     N_output_arrays, output_array_types, output_arrays,
     new_faces);
  
  
  
  /* Free coordinates */
  for (d=0; d<3; ++d) {
    free (new_interp_coords[d]);
  }
  
  
  
  /* Find output variable indices */
  operand_indices = malloc (N_output_arrays * sizeof *operand_indices);
  assert (N_output_arrays == 0 || operand_indices);
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
  assert (N_output_arrays == 0 || operation_codes);
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
  assert (N_output_arrays == 0 || output_array_indices);
  for (m=0; m<N_output_arrays; ++m) {
    assert (operand_indices[m]>=0 && operand_indices[m]<N_input_arrays);
    output_array_indices[m] = input_array_indices[operand_indices[m]];
     assert (output_array_indices[m]==-1
             || (output_array_indices[m]>=0
                 && output_array_indices[m]<CCTK_NumVars()));
  }
  
  
  
  /* Loop over all output arrays */
  for (m=0; m<N_output_arrays; ++m) {
  if (output_array_indices[m]!=-1) {
    
    
    
    /* Find variable parity */
    {
      int vi, gi;
      char tensortypealias[100];
      int firstvar, numvars, vectorlength;
      cGroup group;
      int table;
      int index;
      
      vi = output_array_indices[m];
      assert (vi>=0 && vi<CCTK_NumVars());
      gi = CCTK_GroupIndexFromVarI (vi);
      assert (gi>=0 && gi<CCTK_NumGroups());
      numvars = CCTK_NumVarsInGroupI(gi);
      assert (numvars>0);
      firstvar = CCTK_FirstVarIndexI(gi);
      assert (firstvar>=0);
      index = vi - firstvar;
      assert (index>=0 && index<numvars);
      table = CCTK_GroupTagsTableI(gi);
      assert (table>=0);
      ierr = CCTK_GroupData (gi, &group);
      assert (!ierr);
      vectorlength = group.vectorlength;

      ierr = Util_TableGetString
        (table, sizeof tensortypealias, tensortypealias, "tensortypealias");
      if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
        char * groupname = CCTK_GroupName(gi);
        assert (groupname);
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Tensor type alias not declared for group \"%s\" -- assuming a scalar",
                    groupname);
        free (groupname);
        strcpy (tensortypealias, "scalar");
      } else if (ierr<0) {
        char * groupname = CCTK_GroupName(gi);
        assert (groupname);
        CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Error in tensor type alias declaration for group \"%s\"",
                    groupname);
        free (groupname);
      }
      
      if (CCTK_EQUALS (tensortypealias, "scalar")) {
        if (numvars != vectorlength) {
          char * groupname = CCTK_GroupName(gi);
          assert (groupname);
          CCTK_VWarn (4, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "Group \"%s\" has the tensor type alias \"scalar\", but contains more than 1 element",
                      groupname);
          free (groupname);
        }
        parities[0] = parities[1] = parities[2] = +1;
      } else if (CCTK_EQUALS (tensortypealias, "4scalar")) {
        if (numvars != vectorlength) {
          char * groupname = CCTK_GroupName(gi);
          assert (groupname);
          CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "Group \"%s\" has the tensor type alias \"scalar\", but contains more than 1 element",
                      groupname);
          free (groupname);
        }
        parities[0] = parities[1] = parities[2] = +1;
      } else if (CCTK_EQUALS (tensortypealias, "u")
                 || CCTK_EQUALS (tensortypealias, "d")) {
        assert ((numvars == 3*vectorlength) ||
                (numvars == 3 && vectorlength == 3));
        parities[0] = parities[1] = parities[2] = +1;
        if(numvars == vectorlength) {
          parities[index] = -1;
        } else {
          parities[index/vectorlength] = -1;
        }
      } else if (CCTK_EQUALS (tensortypealias, "4u")
                 || CCTK_EQUALS (tensortypealias, "4d")) {
        assert ((numvars == 4*vectorlength) ||
                (numvars == 4 && vectorlength == 4));
        if (vi == firstvar) {
          parities[0] = parities[1] = parities[2] = +1;
        } else {
          parities[0] = parities[1] = parities[2] = +1;
          parities[index/vectorlength-1] = -1;
        }
      } else if (CCTK_EQUALS (tensortypealias, "uu_sym")
                 || CCTK_EQUALS (tensortypealias, "dd_sym")) {
        assert (numvars == 6*vectorlength);
        parities[0] = parities[1] = parities[2] = +1;
        switch (index/vectorlength) {
        case 0: break;
        case 1: parities[0] = parities[1] = -1; break;
        case 2: parities[0] = parities[2] = -1; break;
        case 3: break;
        case 4: parities[1] = parities[2] = -1; break;
        case 5: break;
        default: assert(0);
        }
      } else if (CCTK_EQUALS (tensortypealias, "uu")
                 || CCTK_EQUALS (tensortypealias, "ud")
                 || CCTK_EQUALS (tensortypealias, "du")
                 || CCTK_EQUALS (tensortypealias, "dd")) {
        assert (numvars == 9*vectorlength);
        parities[0] = parities[1] = parities[2] = +1;
        int const d1 = (index / vectorlength) % 3;
        int const d2 = (index / vectorlength) / 3;
        parities[d1] *= -1;
        parities[d2] *= -1;
      } else if (CCTK_EQUALS (tensortypealias, "ddd_sym")) {
        assert (numvars == 18*vectorlength);
        parities[0] = parities[1] = parities[2] = +1;
        switch ((index/vectorlength)%6) {
        case 0: break;
        case 1: parities[0] = parities[1] = -1; break;
        case 2: parities[0] = parities[2] = -1; break;
        case 3: break;
        case 4: parities[1] = parities[2] = -1; break;
        case 5: break;
        default: assert(0);
        }
        parities[index/vectorlength/6] *= -1;
      } else if (CCTK_EQUALS (tensortypealias, "4uu_sym")
                 || CCTK_EQUALS (tensortypealias, "4dd_sym")) {
        assert (numvars == 10*vectorlength);
        if (vi == firstvar) {
          parities[0] = parities[1] = parities[2] = +1;
        } else if (vi < firstvar + 4) {
          parities[0] = parities[1] = parities[2] = +1;
          parities[index/vectorlength-1] = -1;
        } else {
          parities[0] = parities[1] = parities[2] = +1;
          switch (index/vectorlength-4) {
          case 0: break;
          case 1: parities[0] = parities[1] = -1; break;
          case 2: parities[0] = parities[2] = -1; break;
          case 3: break;
          case 4: parities[1] = parities[2] = -1; break;
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
            parities[d] = weylparities[index/vectorlength][d];
          }
        }
      } else if (CCTK_EQUALS (tensortypealias, "ManualCartesian")) {
	RotatingSymmetry180_GetManualParities(table, gi, parities);
      } else {
        char * groupname = CCTK_GroupName(gi);
        assert (groupname);
        CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Illegal tensor type alias for group \"%s\"",
                    groupname);
        free (groupname);
      }
      
      for (d=0; d<3; ++d) {
        assert (abs(parities[d])==1);
      }
    } /* Find parity */
    
    
    
    /* Modify parity */
    {
      int code = operation_codes[m];
      while (code) {
        const int dir = code % 10;
        code /= 10;
        assert (dir>=1 && dir<=3);
        parities[dir-1] *= -1;
      }
    }
    
    
    
    /* Loop over all points and unfold */
    if (parities[0] * parities[1] == -1) {
      for (n=0; n<N_interp_points; ++n) {
        /* Is the point outside the domain? */
        if (((const CCTK_REAL *)interp_coords[0])[n] < 0) {
          /* Rotate the tensor components back by 180 degrees */
          ((CCTK_REAL *)output_arrays[m])[n] *= -1;
        }
      }
    }
    
    
    
  }
  } /* for m */
  
  
  
  /* Free output variable descriptions */
  free (operand_indices);
  free (operation_codes);
  free (output_array_indices);
  
  

  return iret;
}

void  RotatingSymmetry180_GetManualParities(int table, int gi, int *parities)
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
		  "Invalid format for cartesianreflectionparities: must be of the form xxx where x is + or - for group %s",
		  groupname);
    }
  }
}
