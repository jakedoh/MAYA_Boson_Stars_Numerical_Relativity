#include <assert.h>
#include <stdlib.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <util_Table.h>



void
InterpToArray (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_InterpToArray;
  DECLARE_CCTK_PARAMETERS;
  
  
  
  int interpolator = -1;
  if (!use_carpetinterp2) {
    interpolator = CCTK_InterpHandle (interpolator_name);
    assert (interpolator >= 0);
  }
  
  int options_table = -1;
  if (!use_carpetinterp2) {
    options_table = Util_TableCreateFromString (interpolator_options);
    assert (options_table >= 0);
  }
  
  int coord_handle = -1;
  if (!use_carpetinterp2) {
    coord_handle = CCTK_CoordSystemHandle (interpolator_coordinates);
    assert (coord_handle >= 0);
  }
  
  

  /* Scalars */
  {
    int const nvars = nscalars;
    if (nvars > 0) {
      int const npoints = 1;
      
      CCTK_REAL * restrict const
        coordsx = malloc (npoints * sizeof * coordsx);
      assert (npoints==0 || coordsx);
      CCTK_REAL * restrict const
        coordsy = malloc (npoints * sizeof * coordsy);
      assert (npoints==0 || coordsy);
      CCTK_REAL * restrict const
        coordsz = malloc (npoints * sizeof * coordsz);
      assert (npoints==0 || coordsz);
      CCTK_POINTER_TO_CONST coords[3];
      coords[0] = coordsx;
      coords[1] = coordsy;
      coords[2] = coordsz;
      
      {
        int const n = 0;
        coordsx[n] = scalar_x0;
        coordsy[n] = scalar_y0;
        coordsz[n] = scalar_z0;
      }
      
      CCTK_INT * restrict const
        inputs = malloc (nvars * sizeof * inputs);
      assert (inputs);
      
      for (int n=0; n<nvars; ++n) {
        inputs[n] = CCTK_VarIndex (scalar_vars[n]);
        if (inputs[n] < 0) {
          inputs[n] = -1;
        }
      }
      
      CCTK_INT * restrict const
        output_types = malloc (nvars * sizeof * output_types);
      assert (output_types);
      
      for (int n=0; n<nvars; ++n) {
        output_types[n] = CCTK_VARIABLE_REAL;
      }
      
      CCTK_POINTER * restrict const outputs = malloc (nvars * sizeof * outputs);
      assert (outputs);
      
      for (int n=0; n<nvars; ++n) {
        outputs[n] = &scalars[npoints * n];
      }
      
      if (use_carpetinterp2)
      {
        int const ierr = InterpGridArrays
          (cctkGH, 3, carpetinterp2_interpolator_order,
           npoints, coords,
           nvars, inputs,
           nvars, outputs);
        assert (! ierr);
      }
      else
      {
        int const ierr = CCTK_InterpGridArrays
          (cctkGH, 3, interpolator, options_table, coord_handle,
           npoints, CCTK_VARIABLE_REAL, coords,
           nvars, inputs,
           nvars, output_types, outputs);
        assert (! ierr);
      }
      
      free (coordsx);
      free (coordsy);
      free (coordsz);
      free (inputs);
      free (output_types);
      free (outputs);
    }
  }
  
  
  
  /* 1D Arrays */
  {
    int const nvars = narrays1d;
    if (nvars > 0) {
      int const npoints = array1d_npoints_i;
      
      CCTK_REAL * restrict const
        coordsx = malloc (npoints * sizeof * coordsx);
      assert (npoints==0 || coordsx);
      CCTK_REAL * restrict const
        coordsy = malloc (npoints * sizeof * coordsy);
      assert (npoints==0 || coordsy);
      CCTK_REAL * restrict const
        coordsz = malloc (npoints * sizeof * coordsz);
      assert (npoints==0 || coordsz);
      CCTK_POINTER_TO_CONST coords[3];
      coords[0] = coordsx;
      coords[1] = coordsy;
      coords[2] = coordsz;
      
#pragma omp parallel for
      for (int i=0; i<array1d_npoints_i; ++i) {
        int const n = i;
        coordsx[n] = array1d_x0 + i * array1d_dx_i;
        coordsy[n] = array1d_y0 + i * array1d_dy_i;
        coordsz[n] = array1d_z0 + i * array1d_dz_i;
      }
      
      CCTK_INT * restrict const
        inputs = malloc (nvars * sizeof * inputs);
      assert (inputs);
      
      for (int n=0; n<nvars; ++n) {
        inputs[n] = CCTK_VarIndex (array1d_vars[n]);
        assert (inputs[n] >= 0);
        if (inputs[n] < 0) {
          inputs[n] = -1;
        }
      }
      
      CCTK_INT * restrict const
        output_types = malloc (nvars * sizeof * output_types);
      assert (output_types);
      
      for (int n=0; n<nvars; ++n) {
        output_types[n] = CCTK_VARIABLE_REAL;
      }
      
      CCTK_POINTER * restrict const
        outputs = malloc (nvars * sizeof * outputs);
      assert (outputs);
      
      for (int n=0; n<nvars; ++n) {
        outputs[n] = &arrays1d[npoints * n];
      }
      
      CCTK_INT * restrict const
        operation_codes = malloc (nvars * sizeof * operation_codes);
      assert (operation_codes);
      
      for (int n=0; n<nvars; ++n) {
        operation_codes[n] = array1d_spacederivs[n];
      }
      
      {
        int const ierr = Util_TableSetIntArray
          (options_table, nvars, operation_codes, "operation_codes");
        assert (! ierr);
      }
      
      CCTK_INT * restrict const
        time_deriv_order = malloc (nvars * sizeof * time_deriv_order);
      assert (time_deriv_order);
      
      for (int n=0; n<nvars; ++n) {
        time_deriv_order[n] = array1d_timederivs[n];
      }
      
      {
        int const ierr = Util_TableSetIntArray
          (options_table, nvars, time_deriv_order, "time_deriv_order");
        assert (! ierr);
      }
      
      if (use_carpetinterp2)
      {
        int const ierr = InterpGridArrays
          (cctkGH, 3, carpetinterp2_interpolator_order,
           npoints, coords,
           nvars, inputs,
           nvars, outputs);
        assert (! ierr);
      }
      else
      {
        int const ierr = CCTK_InterpGridArrays
          (cctkGH, 3, interpolator, options_table, coord_handle,
           npoints, CCTK_VARIABLE_REAL, coords,
           nvars, inputs,
           nvars, output_types, outputs);
        assert (! ierr);
      }
      
      {
        int const ierr = Util_TableDeleteKey (options_table, "operation_codes");
        assert (! ierr);
      }
      
      {
        int const ierr =
          Util_TableDeleteKey (options_table, "time_deriv_order");
        assert (! ierr);
      }
      
      free (coordsx);
      free (coordsy);
      free (coordsz);
      free (inputs);
      free (output_types);
      free (outputs);
      free (operation_codes);
      free (time_deriv_order);
    }
  }
  
  
  
  /* 2D Arrays */
  {
    int const nvars = narrays2d;
    if (nvars > 0) {
      int const npoints = array2d_npoints_i * array2d_npoints_j;
      
      CCTK_REAL * restrict const
        coordsx = malloc (npoints * sizeof * coordsx);
      assert (npoints==0 || coordsx);
      CCTK_REAL * restrict const
        coordsy = malloc (npoints * sizeof * coordsy);
      assert (npoints==0 || coordsy);
      CCTK_REAL * restrict const
        coordsz = malloc (npoints * sizeof * coordsz);
      assert (npoints==0 || coordsz);
      CCTK_POINTER_TO_CONST coords[3];
      coords[0] = coordsx;
      coords[1] = coordsy;
      coords[2] = coordsz;
      
#pragma omp parallel for
      for (int j=0; j<array2d_npoints_j; ++j) {
        for (int i=0; i<array2d_npoints_i; ++i) {
          const int n = i + array2d_npoints_i * j;
          coordsx[n] = array2d_x0 + i * array2d_dx_i + j * array2d_dx_j;
          coordsy[n] = array2d_y0 + i * array2d_dy_i + j * array2d_dy_j;
          coordsz[n] = array2d_z0 + i * array2d_dz_i + j * array2d_dz_j;
        }
      }
      
      CCTK_INT * restrict const
        inputs = malloc (nvars * sizeof * inputs);
      assert (inputs);
      
      for (int n=0; n<nvars; ++n) {
        inputs[n] = CCTK_VarIndex (array2d_vars[n]);
        if (inputs[n] < 0) {
          inputs[n] = -1;
        }
      }
      
      CCTK_INT * restrict const
        output_types = malloc (nvars * sizeof * output_types);
      assert (output_types);
      
      for (int n=0; n<nvars; ++n) {
        output_types[n] = CCTK_VARIABLE_REAL;
      }
      
      CCTK_POINTER * restrict const
        outputs = malloc (nvars * sizeof * outputs);
      assert (outputs);
      
      for (int n=0; n<nvars; ++n) {
        outputs[n] = &arrays2d[npoints * n];
      }
      
      if (use_carpetinterp2)
      {
        int const ierr = InterpGridArrays
          (cctkGH, 3, carpetinterp2_interpolator_order,
           npoints, coords,
           nvars, inputs,
           nvars, outputs);
        assert (! ierr);
      }
      else
      {
        int const ierr = CCTK_InterpGridArrays
          (cctkGH, 3, interpolator, options_table, coord_handle,
           npoints, CCTK_VARIABLE_REAL, coords,
           nvars, inputs,
           nvars, output_types, outputs);
        assert (! ierr);
      }
      
      free (coordsx);
      free (coordsy);
      free (coordsz);
      free (inputs);
      free (output_types);
      free (outputs);
    }
  }
  
  
  
  /* 3D Arrays */
  {
    int const nvars = narrays3d;
    if (nvars > 0) {
      int const npoints =
        array3d_npoints_i * array3d_npoints_j * array3d_npoints_k;
      
      CCTK_REAL * restrict const
        coordsx = malloc (npoints * sizeof * coordsx);
      assert (npoints==0 || coordsx);
      CCTK_REAL * restrict const
        coordsy = malloc (npoints * sizeof * coordsy);
      assert (npoints==0 || coordsy);
      CCTK_REAL * restrict const
        coordsz = malloc (npoints * sizeof * coordsz);
      assert (npoints==0 || coordsz);
      CCTK_POINTER_TO_CONST coords[3];
      coords[0] = coordsx;
      coords[1] = coordsy;
      coords[2] = coordsz;
      
#pragma omp parallel for
      for (int k=0; k<array3d_npoints_k; ++k) {
        for (int j=0; j<array3d_npoints_j; ++j) {
          for (int i=0; i<array3d_npoints_i; ++i) {
            int const n = i + array3d_npoints_i * (j + array3d_npoints_j * k);
            coordsx[n] = array3d_x0 + i * array3d_dx_i + j * array3d_dx_j + k * array3d_dx_k;
            coordsy[n] = array3d_y0 + i * array3d_dy_i + j * array3d_dy_j + k * array3d_dy_k;
            coordsz[n] = array3d_z0 + i * array3d_dz_i + j * array3d_dz_j + k * array3d_dz_k;
          }
        }
      }
      
      CCTK_INT * restrict const
        inputs = malloc (nvars * sizeof * inputs);
      assert (inputs);
      
      for (int n=0; n<nvars; ++n) {
        inputs[n] = CCTK_VarIndex (array3d_vars[n]);
        if (inputs[n] < 0) {
          inputs[n] = -1;
        }
      }
      
      CCTK_INT * restrict const
        output_types = malloc (nvars * sizeof * output_types);
      assert (output_types);
      
      for (int n=0; n<nvars; ++n) {
        output_types[n] = CCTK_VARIABLE_REAL;
      }
      
      CCTK_POINTER * restrict const
        outputs = malloc (nvars * sizeof * outputs);
      assert (outputs);
      
      for (int n=0; n<nvars; ++n) {
        outputs[n] = &arrays3d[npoints * n];
      }
      
      if (use_carpetinterp2)
      {
        int const ierr = InterpGridArrays
          (cctkGH, 3, carpetinterp2_interpolator_order,
           npoints, coords,
           nvars, inputs,
           nvars, outputs);
        assert (! ierr);
      }
      else
      {
        int const ierr = CCTK_InterpGridArrays
          (cctkGH, 3, interpolator, options_table, coord_handle,
           npoints, CCTK_VARIABLE_REAL, coords,
           nvars, inputs,
           nvars, output_types, outputs);
        assert (! ierr);
      }
      
      free (coordsx);
      free (coordsy);
      free (coordsz);
      free (inputs);
      free (output_types);
      free (outputs);
    }
  }
  
  
  
  /* Parallel 1D Arrays */
  {
    int const nvars = nparrays1d;
    if (nvars > 0) {
      int const group = CCTK_GroupIndex ("InterpToArray::parrays1d");
      assert (group >= 0);
      cGroupDynamicData dyndata;
      {
        int const ierr = CCTK_GroupDynamicData (cctkGH, group, &dyndata);
        assert (! ierr);
      }
      assert (dyndata.dim == 1);
      int const npoints = dyndata.lsh[0];
      
      CCTK_REAL * restrict const
        coordsx = malloc (npoints * sizeof * coordsx);
      assert (npoints==0 || coordsx);
      CCTK_REAL * restrict const
        coordsy = malloc (npoints * sizeof * coordsy);
      assert (npoints==0 || coordsy);
      CCTK_REAL * restrict const
        coordsz = malloc (npoints * sizeof * coordsz);
      assert (npoints==0 || coordsz);
      CCTK_POINTER_TO_CONST coords[3];
      coords[0] = coordsx;
      coords[1] = coordsy;
      coords[2] = coordsz;
      
#pragma omp parallel for
      for (int i=0; i<dyndata.lsh[0]; ++i) {
        int const n = i;
        coordsx[n] = (parray1d_x0 +
                      (dyndata.lbnd[0] + i) * parray1d_dx_i);
        coordsy[n] = (parray1d_y0 +
                      (dyndata.lbnd[0] + i) * parray1d_dy_i);
        coordsz[n] = (parray1d_z0 +
                      (dyndata.lbnd[0] + i) * parray1d_dz_i);
      }
      
      CCTK_INT * restrict const
        inputs = malloc (nvars * sizeof * inputs);
      assert (inputs);
      
      for (int n=0; n<nvars; ++n) {
        inputs[n] = CCTK_VarIndex (parray1d_vars[n]);
        assert (inputs[n] >= 0);
        if (inputs[n] < 0) {
          inputs[n] = -1;
        }
      }
      
      CCTK_INT * restrict const
        output_types = malloc (nvars * sizeof * output_types);
      assert (output_types);
      
      for (int n=0; n<nvars; ++n) {
        output_types[n] = CCTK_VARIABLE_REAL;
      }
      
      CCTK_POINTER * restrict const
        outputs = malloc (nvars * sizeof * outputs);
      assert (outputs);
      
      for (int n=0; n<nvars; ++n) {
        outputs[n] = &parrays1d[npoints * n];
      }
      
      CCTK_INT * restrict const
        operation_codes = malloc (nvars * sizeof * operation_codes);
      assert (operation_codes);
      
      for (int n=0; n<nvars; ++n) {
        operation_codes[n] = parray1d_spacederivs[n];
      }
      
      {
        int const ierr = Util_TableSetIntArray
          (options_table, nvars, operation_codes, "operation_codes");
        assert (! ierr);
      }
      
      CCTK_INT * restrict const
        time_deriv_order = malloc (nvars * sizeof * time_deriv_order);
      assert (time_deriv_order);
      
      for (int n=0; n<nvars; ++n) {
        time_deriv_order[n] = parray1d_timederivs[n];
      }
      
      {
        int const ierr = Util_TableSetIntArray
          (options_table, nvars, time_deriv_order, "time_deriv_order");
        assert (! ierr);
      }
      
      if (use_carpetinterp2)
      {
        int const ierr = InterpGridArrays
          (cctkGH, 3, carpetinterp2_interpolator_order,
           npoints, coords,
           nvars, inputs,
           nvars, outputs);
        assert (! ierr);
      }
      else
      {
        int const ierr = CCTK_InterpGridArrays
          (cctkGH, 3, interpolator, options_table, coord_handle,
           npoints, CCTK_VARIABLE_REAL, coords,
           nvars, inputs,
           nvars, output_types, outputs);
        assert (! ierr);
      }
      
      {
        int const ierr = Util_TableDeleteKey (options_table, "operation_codes");
        assert (! ierr);
      }
      
      {
        int const ierr =
          Util_TableDeleteKey (options_table, "time_deriv_order");
        assert (! ierr);
      }
      
      free (coordsx);
      free (coordsy);
      free (coordsz);
      free (inputs);
      free (output_types);
      free (outputs);
      free (operation_codes);
      free (time_deriv_order);
    }
  }
  
  
  
  /* Parallel 2D Arrays */
  {
    int const nvars = nparrays2d;
    if (nvars > 0) {
      int const group = CCTK_GroupIndex ("InterpToArray::parrays2d");
      assert (group >= 0);
      cGroupDynamicData dyndata;
      {
        int const ierr = CCTK_GroupDynamicData (cctkGH, group, &dyndata);
        assert (! ierr);
      }
      assert (dyndata.dim == 2);
      int const npoints = dyndata.lsh[0] * dyndata.lsh[1];
      
      CCTK_REAL * restrict const
        coordsx = malloc (npoints * sizeof * coordsx);
      assert (npoints==0 || coordsx);
      CCTK_REAL * restrict const
        coordsy = malloc (npoints * sizeof * coordsy);
      assert (npoints==0 || coordsy);
      CCTK_REAL * restrict const
        coordsz = malloc (npoints * sizeof * coordsz);
      assert (npoints==0 || coordsz);
      CCTK_POINTER_TO_CONST coords[3];
      coords[0] = coordsx;
      coords[1] = coordsy;
      coords[2] = coordsz;
      
#pragma omp parallel for
      for (int j=0; j<dyndata.lsh[1]; ++j) {
        for (int i=0; i<dyndata.lsh[0]; ++i) {
          int const n = i + dyndata.lsh[0] * j;
          coordsx[n] = (parray2d_x0 +
                        (dyndata.lbnd[0] + i) * parray2d_dx_i +
                        (dyndata.lbnd[1] + j) * parray2d_dx_j);
          coordsy[n] = (parray2d_y0 +
                        (dyndata.lbnd[0] + i) * parray2d_dy_i +
                        (dyndata.lbnd[1] + j) * parray2d_dy_j);
          coordsz[n] = (parray2d_z0 +
                        (dyndata.lbnd[0] + i) * parray2d_dz_i +
                        (dyndata.lbnd[1] + j) * parray2d_dz_j);
        }
      }
      
      CCTK_INT * restrict const
        inputs = malloc (nvars * sizeof * inputs);
      assert (inputs);
      
      for (int n=0; n<nvars; ++n) {
        inputs[n] = CCTK_VarIndex (parray2d_vars[n]);
        if (inputs[n] < 0) {
          inputs[n] = -1;
        }
      }
      
      CCTK_INT * restrict const
        output_types = malloc (nvars * sizeof * output_types);
      assert (output_types);
      
      for (int n=0; n<nvars; ++n) {
        output_types[n] = CCTK_VARIABLE_REAL;
      }
      
      CCTK_POINTER * restrict const
        outputs = malloc (nvars * sizeof * outputs);
      assert (outputs);
      
      for (int n=0; n<nvars; ++n) {
        outputs[n] = &parrays2d[npoints * n];
      }
      
      if (use_carpetinterp2)
      {
        int const ierr = InterpGridArrays
          (cctkGH, 3, carpetinterp2_interpolator_order,
           npoints, coords,
           nvars, inputs,
           nvars, outputs);
        assert (! ierr);
      }
      else
      {
        int const ierr = CCTK_InterpGridArrays
          (cctkGH, 3, interpolator, options_table, coord_handle,
           npoints, CCTK_VARIABLE_REAL, coords,
           nvars, inputs,
           nvars, output_types, outputs);
        assert (! ierr);
      }
      
      free (coordsx);
      free (coordsy);
      free (coordsz);
      free (inputs);
      free (output_types);
      free (outputs);
    }
  }
  
  
  
  /* Parallel 3D Arrays */
  {
    int const nvars = nparrays3d;
    if (nvars > 0) {
      int const group = CCTK_GroupIndex ("InterpToArray::parrays3d");
      assert (group >= 0);
      cGroupDynamicData dyndata;
      {
        int const ierr = CCTK_GroupDynamicData (cctkGH, group, &dyndata);
        assert (! ierr);
      }
      assert (dyndata.dim == 3);
      int const npoints = dyndata.lsh[0] * dyndata.lsh[1] * dyndata.lsh[2];
      
      CCTK_REAL * restrict const
        coordsx = malloc (npoints * sizeof * coordsx);
      assert (npoints==0 || coordsx);
      CCTK_REAL * restrict const
        coordsy = malloc (npoints * sizeof * coordsy);
      assert (npoints==0 || coordsy);
      CCTK_REAL * restrict const
        coordsz = malloc (npoints * sizeof * coordsz);
      assert (npoints==0 || coordsz);
      CCTK_POINTER_TO_CONST coords[3];
      coords[0] = coordsx;
      coords[1] = coordsy;
      coords[2] = coordsz;
      
#pragma omp parallel for
      for (int k=0; k<dyndata.lsh[2]; ++k) {
        for (int j=0; j<dyndata.lsh[1]; ++j) {
          for (int i=0; i<dyndata.lsh[0]; ++i) {
            int const n = i + dyndata.lsh[0] * (j + dyndata.lsh[1] * k);
            coordsx[n] = (parray3d_x0 +
                          (dyndata.lbnd[0] + i) * parray3d_dx_i +
                          (dyndata.lbnd[1] + j) * parray3d_dx_j +
                          (dyndata.lbnd[2] + k) * parray3d_dx_k);
            coordsy[n] = (parray3d_y0 +
                          (dyndata.lbnd[0] + i) * parray3d_dy_i +
                          (dyndata.lbnd[1] + j) * parray3d_dy_j +
                          (dyndata.lbnd[2] + k) * parray3d_dy_k);
            coordsz[n] = (parray3d_z0 +
                          (dyndata.lbnd[0] + i) * parray3d_dz_i +
                          (dyndata.lbnd[1] + j) * parray3d_dz_j +
                          (dyndata.lbnd[2] + k) * parray3d_dz_k);
          }
        }
      }
      
      CCTK_INT * restrict const
        inputs = malloc (nvars * sizeof * inputs);
      assert (inputs);
      
      for (int n=0; n<nvars; ++n) {
        inputs[n] = CCTK_VarIndex (parray3d_vars[n]);
        if (inputs[n] < 0) {
          inputs[n] = -1;
        }
      }
      
      CCTK_INT * restrict const
        output_types = malloc (nvars * sizeof * output_types);
      assert (output_types);
      
      for (int n=0; n<nvars; ++n) {
        output_types[n] = CCTK_VARIABLE_REAL;
      }
      
      CCTK_POINTER * restrict const
        outputs = malloc (nvars * sizeof * outputs);
      assert (outputs);
      
      for (int n=0; n<nvars; ++n) {
        outputs[n] = &parrays3d[npoints * n];
      }
      
      if (use_carpetinterp2)
      {
        int const ierr = InterpGridArrays
          (cctkGH, 3, carpetinterp2_interpolator_order,
           npoints, coords,
           nvars, inputs,
           nvars, outputs);
        assert (! ierr);
      }
      else
      {
        int const ierr = CCTK_InterpGridArrays
          (cctkGH, 3, interpolator, options_table, coord_handle,
           npoints, CCTK_VARIABLE_REAL, coords,
           nvars, inputs,
           nvars, output_types, outputs);
        assert (! ierr);
      }
      
      free (coordsx);
      free (coordsy);
      free (coordsz);
      free (inputs);
      free (output_types);
      free (outputs);
    }
  }
  
  
  
  if (options_table >= 0)
  {
    int const ierr = Util_TableDestroy (options_table);
    assert (! ierr);
  }
}
