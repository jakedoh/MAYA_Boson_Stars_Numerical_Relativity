/*@@
   @file      Write.c
   @date      Thu 18 April 2002
   @author    Thomas Radke
   @desc
              Output two-dimensional slices in Jpeg image format.
   @enddesc
   @version   $Id$
 @@*/


#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "util_Table.h"
#include "CactusBase/IOUtil/src/ioutil_AdvertisedFiles.h"
#include "ioJpegGH.h"

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Header$";
CCTK_FILEVERSION(CactusIO_IOJpeg_Write_c)


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
static void WriteData (const cGH *GH, int vindex, const char *alias, int dim,
                       int dir, CCTK_REAL min, CCTK_REAL max,
                       const CCTK_INT hsize[2], const void *hdata);
static void *RefineData (CCTK_INT input_dims[2], const void *input_data);


/*@@
   @routine   IOJpeg_Write
   @date      Thu 18 April 2002
   @author    Thomas Radke
   @desc
              Writes the slices of a variable into separate Jpeg files.
   @enddesc
   @calls     Hyperslab_GlobalMappingByIndex
              Hyperslab_FreeMapping
              Hyperslab_GetList
              OpenFile
              WriteData

   @var       GH
   @vdesc     pointer to CCTK GH
   @vtype     const cGH *
   @vio       in
   @endvar
   @var       vindex
   @vdesc     index of variable to output
   @vtype     CCTK_INT
   @vio       in
   @endvar
   @var       alias
   @vdesc     alias name of variable to output
   @vtype     const char *
   @vio       in
   @endvar

   @returntype int
   @returndesc
                0 for success, or<BR>
               -1 if variable has no storage assigned<BR>
               -2 if output file couldn't be opened<BR>
               -3 if hyperslab for coordinates and/or variable couldn't be
                  extracted
   @endreturndesc
@@*/
int IOJpeg_Write (const cGH *GH, CCTK_INT vindex, const char *alias)
{
  const ioJpegGH *myGH;
  int mapping, total_hsize;
  int dir_i, dir_j, maxdir, myproc, groupindex;
  cGroup gdata;
  char *fullname;
  CCTK_INT origin[3], direction[6], hsize[2];
  const CCTK_INT extent[2] = {-1, -1};
  void *hdata;
  CCTK_REAL min, max;
  const CCTK_INT htype = CCTK_VARIABLE_REAL;
  DECLARE_CCTK_PARAMETERS
    
    /* Reduction variables */
  CCTK_INT input_array[1];
  CCTK_INT input_array_type_codes[1];
  void*  value_min[1], *value_max[1];

  input_array_type_codes[0]  = CCTK_VARIABLE_REAL;
  input_array[0] = vindex;
  value_min[0] = &min;
  value_max[0] = &max;

  /* get the variable name and group information */
  fullname = CCTK_FullName (vindex);
  groupindex = CCTK_GroupIndexFromVarI (vindex);
  CCTK_GroupData (groupindex, &gdata);

  /* check if variable has storage assigned */
  if (! CCTK_QueryGroupStorageI (GH, groupindex))
  {
    CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                "No IOJpeg output for '%s' (no storage)", fullname);
    free (fullname);
    return (-1);
  }

  /* get the handle for IOJpeg extensions */
  myproc = CCTK_MyProc (GH);
  myGH = (const ioJpegGH *) CCTK_GHExtension (GH, "IOJpeg");

  /* set the minimum/maximum for the colormap */
  if (CCTK_Equals (colormap, "custom"))
  {
    min = colormap_min;
    max = colormap_max;
  }
  else if (CCTK_Equals (colormap, "auto"))
  {
    int i = CCTK_LocalArrayReductionHandle ("minimum");
    CCTK_ReduceGridArrays (GH, 0, i, -1, 1, input_array , 1, input_array_type_codes, value_min);
    i = CCTK_LocalArrayReductionHandle ("maximum");
    CCTK_ReduceGridArrays (GH, 0, i, -1, 1, input_array , 1, input_array_type_codes, value_max);
  }
  else if (CCTK_Equals (colormap, "auto-old"))
  {
    int i = CCTK_ReductionHandle ("minimum");
    CCTK_Reduce (GH, 0, i, 1, CCTK_VARIABLE_REAL, &min, 1, vindex);
    i = CCTK_ReductionHandle ("maximum");
    CCTK_Reduce (GH, 0, i, 1, CCTK_VARIABLE_REAL, &max, 1, vindex);
  }
  else
  {
    CCTK_WARN (CCTK_WARN_ABORT, "Setting for parameter IOpjeg::colormap not recognised");
  }

  /* get the number of slices to output */
  /* in general: maxdir = gdata.dim * (gdata.dim - 1) / 2; */
  maxdir = gdata.dim == 2 ? 1 : 3;

  if (CCTK_EQUALS(gridpoints, "hyperslab"))
  {
  
    /* now do the actual I/O looping over all directions */
    for (int dir = 0; dir < maxdir; dir++)
    {
      /* get the directions to span the hyperslab */
      if (dir == 0)
      {
        dir_i = 0; dir_j = 1;   /* xy */
        if(GH->cctk_gsh[0]==1 || GH->cctk_gsh[1]==1)
        	continue;
      }
      else if (dir == 1)
      {
        dir_i = 0; dir_j = 2;   /* xz */
        if(GH->cctk_gsh[0]==1 || GH->cctk_gsh[2]==1)
        	continue;
      }
      else
      {
        dir_i = 1; dir_j = 2;   /* yz */
        if(GH->cctk_gsh[1]==1 || GH->cctk_gsh[2]==1)
        	continue;
      }
  
      /* set the origin using the slice center from IOUtil */
      memset (origin, 0, sizeof (origin));
      if (gdata.grouptype == CCTK_GF)
      {
        origin[maxdir-dir-1] = myGH->sp2xyz[gdata.dim - 1][dir];
      }
  
      /* set the direction vector */
      memset (direction, 0, sizeof (direction));
      direction[dir_i] = direction[gdata.dim + dir_j] = 1;
  
      mapping = Hyperslab_GlobalMappingByIndex (GH, vindex, 2,
                                                direction, origin, extent,
                                                NULL, -1, NULL, hsize);
      if (mapping < 0)
      {
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Failed to define hyperslab mapping for variable '%s'",
                    fullname);
        continue;
      }
      {
        // Sanity check
        // (if this fails, we requested an insane number of grid points)
        int imax = INT_MAX;
        int d;
        for (d=0; d<2; ++d) {
          assert (hsize[d] >= 0 && hsize[d] <= imax);
          if (hsize[d] > 0) imax /= hsize[d];
        }
        assert (CCTK_VarTypeSize (gdata.vartype) <= imax);
      }
      total_hsize = hsize[0] * hsize[1] * CCTK_VarTypeSize (htype);
      assert (total_hsize >= 0);
      if (total_hsize <= 0)
      {
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Selected hyperslab has zero size for variable '%s' "
                    "direction %d", fullname, dir);
        Hyperslab_FreeMapping (mapping);
        continue;
      }
  
      /* allocate hyperslab buffer */
      hdata = myproc == 0 ? malloc (total_hsize) : NULL;
  
      /* get the hyperslab */
      int i = Hyperslab_Get (GH, mapping, 0, vindex, 0, htype, hdata);
      if (i)
      {
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Failed to extract hyperslab for variable '%s'", fullname);
      }
  
      /* release the mapping structure */
      Hyperslab_FreeMapping (mapping);
  
      /* and dump the data to file */
      if (myproc == 0)
      {
        if (i == 0)
        {
          /* refine the 2D slice if requested */
          void *outdata;
          outdata = refinement_factor > 1 ? RefineData (hsize, hdata) : hdata;
          if (! outdata)
          {
            CCTK_WARN (1, "Failed to refine 2D array");
            outdata = hdata;
          }
          WriteData (GH, vindex, alias, gdata.dim, dir, min, max, hsize, outdata);
          if (outdata != hdata)
          {
            free (outdata);
          }
        }
  
        /* clean up */
        free (hdata);
      }
    } /* end of looping through xyz directions */

  }
  else if (CCTK_EQUALS(gridpoints, "interpolate"))
  {
    
    int interpolator;
    int options_table;
    int coord_handle;
    
    CCTK_REAL * restrict coordsx;
    CCTK_REAL * restrict coordsy;
    CCTK_REAL * restrict coordsz;
    CCTK_POINTER_TO_CONST coords[3];
    CCTK_INT * restrict inputs;
    CCTK_INT * restrict output_types;
    CCTK_POINTER * restrict outputs;
    
    const int nvars = 1;
    int npoints;
    
    int n;
    int i, j;
    int ierr;
    
    CCTK_REAL *outdata;
    
    
    
    interpolator = CCTK_InterpHandle (interpolator_name);
    assert (interpolator >= 0);
    
    options_table = Util_TableCreateFromString (interpolator_options);
    assert (options_table >= 0);
    
    coord_handle = CCTK_CoordSystemHandle (interpolator_coordinates);
    assert (coord_handle >= 0);
    
    
    
    {
      ierr = Util_TableSetInt (options_table, 1, "want_global_mode");
      assert (! ierr);
    }
    
    
    
    /* 2D Arrays */
    npoints = myproc==0 ? array2d_npoints_i * array2d_npoints_j : 0;
    
    coordsx = malloc (npoints * sizeof * coordsx);
    assert (npoints==0 || coordsx);
    coordsy = malloc (npoints * sizeof * coordsy);
    assert (npoints==0 || coordsy);
    coordsz = malloc (npoints * sizeof * coordsz);
    assert (npoints==0 || coordsz);
    coords[0] = coordsx;
    coords[1] = coordsy;
    coords[2] = coordsz;
    
    if (myproc==0) {
      n = 0;
      for (j=0; j<array2d_npoints_j; ++j) {
        for (i=0; i<array2d_npoints_i; ++i) {
          assert (n <= npoints);
          coordsx[n] = array2d_x0 + i * array2d_dx_i + j * array2d_dx_j;
          coordsy[n] = array2d_y0 + i * array2d_dy_i + j * array2d_dy_j;
          coordsz[n] = array2d_z0 + i * array2d_dz_i + j * array2d_dz_j;
          ++n;
        }
      }
      assert (n == npoints);
    }
    
    inputs = malloc (nvars * sizeof * inputs);
    assert (inputs);
    
    for (n=0; n<nvars; ++n) {
      inputs[n] = vindex;
      if (inputs[n] < 0) {
        inputs[n] = -1;
      }
    }
    
    output_types = malloc (nvars * sizeof * output_types);
    assert (output_types);
    
    for (n=0; n<nvars; ++n) {
      output_types[n] = CCTK_VARIABLE_REAL;
    }
    
    outputs = malloc (nvars * sizeof * outputs);
    assert (outputs);
    
    outdata = malloc (npoints * sizeof * outdata);
    assert (npoints==0 || outdata);
    
    for (n=0; n<nvars; ++n) {
      outputs[n] = outdata;
    }
    
    ierr = CCTK_InterpGridArrays
      (GH, 3, interpolator, options_table, coord_handle,
       npoints, CCTK_VARIABLE_REAL, coords,
       nvars, inputs,
       nvars, output_types, outputs);
    assert (! ierr);
    
    
    
    /* Multiply values by r */
    if (myproc==0) {
      for (n=0; n<npoints; ++n) {
        const CCTK_REAL x = coordsx[n];
        const CCTK_REAL y = coordsy[n];
        const CCTK_REAL z = coordsz[n];
        const CCTK_REAL r = sqrt (pow (x, 2) + pow (y, 2) + pow (z, 2));
        outdata[n] *= r;
      }
    }
    
    

    free (coordsx);
    free (coordsy);
    free (coordsz);
    free (inputs);
    free (output_types);
    free (outputs);
    
    
    
    ierr = Util_TableDestroy (options_table);
    assert (! ierr);
    
    
    
    if (myproc==0) {
      const int dir = 0;
      hsize[0] = array2d_npoints_i;
      hsize[1] = array2d_npoints_j;
      
      assert (CCTK_VarTypeI(vindex) == CCTK_VARIABLE_REAL);
      
      WriteData (GH, vindex, alias, gdata.dim, dir,
                 min, max, hsize, outdata);
    }
    
    free (outdata);
    
  }
  else
  {
    CCTK_WARN (CCTK_WARN_ABORT, "unknown setting for parameter IOJpeg::gridpoints");
  }

  free (fullname);

  return (0);
}


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
/*@@
   @routine    RefineData
   @date       Thu 4 November 2004
   @author     Thomas Radke
   @desc
               Refines a given 2D input array by a certain factor
   @enddesc

   @returntype void *
   @returndesc
               a pointer to the allocated refined output array (must be freed
               by the caller), or NULL in case of an error
   @endreturndesc
 @@*/
static void *RefineData (CCTK_INT input_dims[2], const void *input_data)
{
  const CCTK_REAL coord_origin[2] = {0.0, 0.0};
  CCTK_REAL coord_delta[2];
  const CCTK_INT array_type_codes[2] = {CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL};
  CCTK_REAL *interp_coords[2], *refined_data;
  int i, j, interp_handle, table_handle, N_interp_points;
  CCTK_INT output_dims[2];
  DECLARE_CCTK_PARAMETERS


  interp_handle = CCTK_InterpHandle ("uniform cartesian");
  if (interp_handle < 0)
  {
    CCTK_WARN (1, "Couldn't get handle for interpolation operator 'uniform "
               "cartesian'. Did you forget to activate a thorn providing "
               "CCTK_InterpLocalUniform() ?");
    return (NULL);
  }

  table_handle = Util_TableCreateFromString ("order = 1");
  coord_delta[0] = coord_delta[1] = refinement_factor;

  /* leave one grid point at the boundary as interpolation stencil */
  output_dims[0] = refinement_factor * (input_dims[0] - 1) - 2;
  output_dims[1] = refinement_factor * (input_dims[1] - 1) - 2;
  N_interp_points = output_dims[0] * output_dims[1];

  interp_coords[0] = malloc (2 * N_interp_points * sizeof (CCTK_REAL));
  interp_coords[1] = interp_coords[0] + N_interp_points;
  for (j = 0; j < output_dims[1]; j++)
  {
    for (i = 0; i < output_dims[0]; i++)
    {
      *interp_coords[0]++ = i + 1;
      *interp_coords[1]++ = j + 1;
    }
  }
  interp_coords[0] -= N_interp_points;
  interp_coords[1] -= N_interp_points;

  refined_data = malloc (N_interp_points * sizeof (CCTK_REAL));
  if (CCTK_InterpLocalUniform (2, interp_handle, table_handle,
                               coord_origin, coord_delta, N_interp_points,
                               CCTK_VARIABLE_REAL,
                               (const void *const *) interp_coords,
                               1, input_dims, array_type_codes,
                               &input_data, 1, array_type_codes,
                               (void *const *) &refined_data) != 0)
  {
    CCTK_WARN (1, "Failed to interpolate 2D array");
    free (refined_data);
    refined_data = NULL;
  }
  else
  {
    input_dims[0] = output_dims[0];
    input_dims[1] = output_dims[1];
  }

  free (interp_coords[0]);
  Util_TableDestroy (table_handle);

  return (refined_data);
}


/*@@
   @routine    WriteData
   @date       Thu 18 April 2002
   @author     Thomas Radke
   @desc
               Writes the given hyperslab into a Jpeg output file.
   @enddesc
 @@*/
static void WriteData (const cGH *GH, int vindex, const char *alias, int dim,
                       int dir, CCTK_REAL min, CCTK_REAL max,
                       const CCTK_INT hsize[2], const void *hdata)
{
  FILE *file;
  unsigned char *dataout;
  const ioJpegGH *myGH;
  ioAdvertisedFileDesc advertised_file;
  char *filename, *tmpfilename, *fullname;
  char slicename[30];
  const char *extensions[] = {"xy", "xz", "yz"};
  DECLARE_CCTK_PARAMETERS


  myGH = (const ioJpegGH *) CCTK_GHExtension (GH, "IOJpeg");

  /* allocate the RGB image buffer */
  dataout = (unsigned char *) malloc (3 * hsize[0] * hsize[1]);

  AutoColorDataSlice (hsize[0], hsize[1], hdata, dataout, min, max,
                      colormap_bias, colormap_factor);

  /* open the file */
  if (dim == 2)
  {
    strcpy (slicename, "2D");
  }
  else
  {
    /* give the slice origin as range [1 .. n] */
    sprintf (slicename, "%s_[%d]", extensions[dir], myGH->sp2xyz[dim-1][dir]);
  }

  filename = malloc (strlen (myGH->out_dir) + strlen (alias) +
                     sizeof (slicename) + 20);

  if (CCTK_Equals (mode, "remove"))
  {
    sprintf (filename, "%s%s_%s.jpeg", myGH->out_dir, alias, slicename);
    tmpfilename = malloc (strlen (filename) + 5);
    sprintf (tmpfilename, "%s.tmp", filename);
  }
  else
  {
    sprintf (filename, "%s%s_%s.it_%d.jpeg", myGH->out_dir, alias, slicename,
             GH->cctk_iteration);
    tmpfilename = NULL;
  }

  /* Write a JPEG file to be advertised to a temporary file first
     and rename it later.
     This fixes a racing problem with HTTPD when running with pthreads support:
     a file could be downloaded (through the HTTPD thread) while is was still
     being written to (in the main simulation thread).
     Now the JPEG is written to a temporary file first and then (atomically)
     renamed. */
  file = fopen (tmpfilename ? tmpfilename : filename, "w");
  if (file)
  {
    /* write the data */
    WriteJPEGToFileRGB (hsize[0], hsize[1], dataout, colormap_quality, file);

    /* close the file */
    fclose (file);

    /* in "remove" mode: rename and advertise the file for downloading */
    if (tmpfilename)
    {
      if (rename (tmpfilename, filename))
      {
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Cannot rename temporary output file '%s' into '%s'",
                    tmpfilename, filename);
      }
      else
      {
        fullname = CCTK_FullName (vindex);
        advertised_file.slice = slicename;
        advertised_file.thorn = CCTK_THORNSTRING;
        advertised_file.varname = fullname;
        advertised_file.description = "Jpegs of slices";
        advertised_file.mimetype = "image/jpeg";

        IOUtil_AdvertiseFile (GH, filename, &advertised_file);

        free (fullname);
      }
    }
  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Cannot open IOJpeg %s output file '%s'",
                tmpfilename ? "temporary" : "",
                tmpfilename ? tmpfilename : filename);
  }

  /* clean up */
  free (dataout);
  free (tmpfilename);
  free (filename);
}
