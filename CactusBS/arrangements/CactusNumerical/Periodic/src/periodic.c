#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "cctk.h"
#include "cctk_Parameters.h"
#include "Slab.h"



CCTK_INT
BndPeriodicVI (CCTK_POINTER_TO_CONST _GH,
               CCTK_INT size,
               CCTK_INT const * restrict const stencil,
               CCTK_INT const do_periodic[3],
               CCTK_INT const vi)
{
  cGH const * restrict const cctkGH = _GH;
  DECLARE_CCTK_PARAMETERS;
  
  cGroup group;
  cGroupDynamicData data;
  void * restrict varptr;
  struct xferinfo * restrict xferinfo;
  int global_bbox[6];
  int global_lbnd[3], global_ubnd[3];
  int fake_bbox[6];
  int gi;
  int dir, face;
  int d;
  int ierr;
  
  /* Check arguments */
  assert (cctkGH);
  assert (stencil);
  assert (vi>=0 && vi<CCTK_NumVars());
  
  /* Get and check group info */
  gi = CCTK_GroupIndexFromVarI (vi);
  assert (gi>=0 && gi<CCTK_NumGroups());
  
  ierr = CCTK_GroupData (gi, &group);
  assert (!ierr);
  assert (group.grouptype == CCTK_GF);
  assert (group.disttype == CCTK_DISTRIB_DEFAULT);
  int const vartypesize = CCTK_VarTypeSize(group.vartype);
  assert (vartypesize>0);

  if(group.dim != size) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "The group \"%s\" has dimension %d, but the given stencil has "
               "size %d.", CCTK_GroupNameFromVarI(vi), group.dim, (int)size);
  }
  
  ierr = CCTK_GroupDynamicData (cctkGH, gi, &data);
  assert (!ierr);
  
  varptr = CCTK_VarDataPtrI (cctkGH, 0, vi);
  assert (varptr);
  
  {
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
    
    for (d=0; d<3; ++d) {
      fake_bbox[2*d  ] = data.lbnd[d] == global_lbnd[d];
      fake_bbox[2*d+1] = data.ubnd[d] == global_ubnd[d];
    }
  }
  
  if (poison_boundaries) {
    /* poison destination grid points */
    
    for (dir=0; dir<group.dim; ++dir) {
      if (dir<3 && do_periodic[dir]) {
        assert (stencil[dir] >= 0);
        for (face=0; face<2; ++face) {
          if (cctkGH->cctk_bbox[2*dir+face]) {
            
            int imin[3], imax[3];
            for (d=0; d<group.dim; ++d) {
              imin[d] = 0;
              imax[d] = data.lsh[d];
            }
            for (d=group.dim; d<3; ++d) {
              imin[d] = 0;
              imax[d] = 1;
            }
            if (face==0) {
              imax[dir] = imin[dir] + stencil[dir];
            } else {
              imin[dir] = imax[dir] - stencil[dir];
            }
            
#pragma omp parallel for
            for (int k=imin[2]; k<imax[2]; ++k) {
              for (int j=imin[1]; j<imax[1]; ++j) {
                int const ind3d = imin[0] + data.ash[0] * (j + data.ash[1] * k);
                memset ((char*)varptr + ind3d*vartypesize,
                        poison_value, (imax[0]-imin[0])*vartypesize);
              }
            }
            
          }
        }
      }
    }
    
  }
  
  /* Allocate slab transfer description */
  xferinfo = malloc(group.dim * sizeof *xferinfo);
  assert (xferinfo);
  
  for (dir=0; dir<group.dim; ++dir) {
    if (dir<3 && do_periodic[dir]) {
      
      assert (stencil[dir] >= 0);
      
      if (data.gsh[dir] < 2*stencil[dir]+1) {
        return dir+1;
      }
      
      /* Loop over faces */
      for (face=0; face<2; ++face) {
        
        if (global_bbox[2*dir+face]) {
          
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
          
          {
            int const domain_size = data.gsh[dir] - 2*stencil[dir];
            int step_size = stencil[dir];
            int num_steps = 1;
            int step;
            int source;
            
            assert (domain_size > 0);
            if (domain_size < step_size) {
              /* TODO: this could be made more efficient by taking
                 larger steps */
              step_size = 1;
              num_steps = stencil[dir];
            }
            assert (num_steps*step_size == stencil[dir]);
            
            for (step=0; step<num_steps; ++step) {
              
              if (face==0) {
                /* Fill in lower face */
                source = data.gsh[dir] - 2*stencil[dir] + step*step_size;
                while (source < stencil[dir]) {
                  source += domain_size;
                }
                assert (source+step_size <= data.gsh[dir]-stencil[dir]);
                xferinfo[dir].src.off = source;
                xferinfo[dir].src.len = step_size;
                xferinfo[dir].dst.off = step*step_size;
                xferinfo[dir].dst.len = step_size;
              } else {
                /* Fill in upper face */
                source = stencil[dir] + step*step_size;
                while (source + step_size > data.gsh[dir] - stencil[dir]) {
                  source -= domain_size;
                }
                assert (source >= stencil[dir]);
                xferinfo[dir].src.off = source;
                xferinfo[dir].src.len = step_size;
                xferinfo[dir].dst.off
                  = data.gsh[dir] - stencil[dir] + step*step_size;
                xferinfo[dir].dst.len = step_size;
              }
              
              ierr = Slab_Transfer
                (cctkGH, group.dim, xferinfo, -1,
                 group.vartype, varptr, group.vartype, varptr);
              assert (!ierr);
              
            } /* for step */
          }
          
        } /* if bbox */
        
      } /* for f */
      
    } /* if dir */
  } /* for dir */
  
  free (xferinfo);
  
  if (check_boundaries) {
    /* check destination grid points for poison */
    
    char poison[vartypesize];
    memset (poison, poison_value, vartypesize);
    
    for (dir=0; dir<group.dim; ++dir) {
      if (dir<3 && do_periodic[dir]) {
        assert (stencil[dir] >= 0);
        for (face=0; face<2; ++face) {
          if (cctkGH->cctk_bbox[2*dir+face]) {
            
            int imin[3], imax[3];
            for (d=0; d<group.dim; ++d) {
              imin[d] = 0;
              imax[d] = data.lsh[d];
            }
            for (d=group.dim; d<3; ++d) {
              imin[d] = 0;
              imax[d] = 1;
            }
            if (face==0) {
              imax[dir] = imin[dir] + stencil[dir];
            } else {
              imin[dir] = imax[dir] - stencil[dir];
            }
            
            int nerrors = 0;
#pragma omp parallel for reduction(+: nerrors)
            for (int k=imin[2]; k<imax[2]; ++k) {
              for (int j=imin[1]; j<imax[1]; ++j) {
                for (int i=imin[0]; i<imax[0]; ++i) {
                  int const ind3d = i + data.ash[0] * (j + data.ash[1] * k);
                  if (! memcmp((char*)varptr + ind3d*vartypesize, poison,
                               vartypesize))
                  {
                    ++nerrors;
                  }
                }
              }
            }
            if (nerrors > 0) {
              char *const fullname = CCTK_FullName(vi);
              CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                         "Found poison on symmetry boundary of variable \"%s\" direction %d face %d after applying periodicity boundary condition",
                         fullname, dir, face);
              free (fullname);
            }
            
          }
        }
      }
    }
    
  }
  
  return 0;
}



CCTK_INT
BndPeriodicVN (CCTK_POINTER_TO_CONST cctkGH,
               CCTK_INT size,
               CCTK_INT const * restrict const stencil,
               CCTK_INT const do_periodic[3],
               char const * restrict const vn)
{
  int const vi = CCTK_VarIndex (vn);
  assert (vi>=0 && vi<CCTK_NumVars());
  return BndPeriodicVI (cctkGH, size, stencil, do_periodic, vi);
}



CCTK_INT
BndPeriodicGI (CCTK_POINTER_TO_CONST cctkGH,
               CCTK_INT size,
               CCTK_INT const * restrict const stencil,
               CCTK_INT const do_periodic[3],
               CCTK_INT const gi)
{
  int v1, nv;
  int vi;
  int ierr;
  assert (gi>=0 && gi<CCTK_NumGroups());
  nv = CCTK_NumVarsInGroupI(gi);
  assert (nv>=0);
  if (nv>0) {
    v1 = CCTK_FirstVarIndexI(gi);
    assert (v1>=0 && v1<CCTK_NumVars());
    for (vi=v1; vi<v1+nv; ++vi) {
      ierr = BndPeriodicVI (cctkGH, size, stencil, do_periodic, vi);
      if (ierr) return ierr;
    }
  }
  return 0;
}



CCTK_INT
BndPeriodicGN (CCTK_POINTER_TO_CONST cctkGH,
               CCTK_INT size,
               CCTK_INT const * restrict const stencil,
               CCTK_INT const do_periodic[3],
               char const * restrict const gn)
{
  int const gi = CCTK_GroupIndex (gn);
  assert (gi>=0 && gi<CCTK_NumGroups());
  return BndPeriodicGI (cctkGH, size, stencil, do_periodic, gi);
}

void
Periodic_RegisterBC (cGH * restrict const cctkGH)
{
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT handle;
  CCTK_INT faces[6];
  CCTK_INT width[6];
  CCTK_INT ierr;
  CCTK_INT is_internal[6];
  CCTK_INT is_staggered[6];
  CCTK_INT shiftout[6];

  faces[0] = faces[1] = periodic || periodic_x;
  faces[2] = faces[3] = periodic || periodic_y;
  faces[4] = faces[5] = periodic || periodic_z;

  ierr = GetBoundarySpecification
    (6, width, is_internal, is_staggered, shiftout);
  if (ierr < 0)
  {
    CCTK_WARN (0, "Could not get the boundary specification");
  }

  handle = SymmetryRegister ("periodic");
  if (handle < 0) {
    CCTK_WARN (0, "Could not register periodicity boundary condition");
  }
  
  ierr = SymmetryRegisterGrid (cctkGH, handle, faces, width);
  if (ierr < 0) {
    CCTK_WARN (0, "Could not register the periodic boundaries -- probably some other thorn has already registered the same boundary faces for a different symmetry");
  }
}

void
Periodic_ApplyBC (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT do_periodic[3];
  do_periodic[0] = periodic || periodic_x;
  do_periodic[1] = periodic || periodic_y;
  do_periodic[2] = periodic || periodic_z;

  char * fullname;
  int vi;
  int dim;
  int * restrict stencili;
  CCTK_INT * restrict stencil;
  int ierr;

  int nvars;
  CCTK_INT * restrict indices;
  CCTK_INT * restrict faces;
  CCTK_INT * restrict widths;
  CCTK_INT * restrict tables;

  assert (cctkGH);

  nvars = Boundary_SelectedGVs (cctkGH, 0, 0, 0, 0, 0, 0);
  assert (nvars>=0);

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

  for (int i=0; i<nvars; ++i) {
    vi = indices[i];
    assert (vi>=0 && vi<CCTK_NumVars());

    assert (widths[i] >= 0);

    dim = CCTK_GroupDimFromVarI (vi);
    assert (dim>=0);

    stencili = malloc (dim * sizeof *stencili);
    assert (stencili);
    ierr = CCTK_GroupnghostzonesVI (cctkGH, dim, stencili, vi);
    assert (!ierr);
    stencil = malloc (dim * sizeof *stencil);
    assert (stencil);
    for (int d=0; d<dim; ++d) {
      stencil[d] = stencili[d];
    }

    if (verbose) {
      fullname = CCTK_FullName(vi);
      assert (fullname);
      CCTK_VInfo (CCTK_THORNSTRING,
                  "Applying periodicity boundary conditions to \"%s\"",
                  fullname);
      free (fullname);
    }

    ierr = BndPeriodicVI (cctkGH, dim, stencil, do_periodic,  vi);

    if(ierr > 0 && ierr < 4) {
      int dir = ierr-1;

      int gi;
      cGroup group;
      cGroupDynamicData data;

      gi = CCTK_GroupIndexFromVarI (vi);
      ierr = CCTK_GroupData (gi, &group);
      assert (!ierr);

      ierr = CCTK_GroupDynamicData (cctkGH, gi, &data);
      assert (!ierr);

      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "The group \"%s\" has in the %c-direction only %d grid points.  "
                    "This is not large enough for a periodic boundary that is %d grid points wide.  "
                    "The group needs to have at least %d grid points in that direction.",
                    CCTK_GroupNameFromVarI(vi), "xyz"[dir], data.gsh[dir],
                    (int)stencil[dir],
                    (int)(2*stencil[dir]+1));
    }

    assert (!ierr);

    free (stencili);
    free (stencil);
  }

  free (indices);
  free (faces);
  free (widths);
  free (tables);
}
