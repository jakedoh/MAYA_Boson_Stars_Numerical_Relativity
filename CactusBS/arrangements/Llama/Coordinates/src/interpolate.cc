
/* Copyright 2013 Peter Diener, Nils Dorband, Roland Haas, Ian Hinder,
Christian Ott, Denis Pollney, Thomas Radke, Christian Reisswig, Erik
Schnetter, Barry Wardell and Burkhard Zink

This file is part of Llama.

Llama is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 2 of the License, or (at your
option) any later version.

Llama is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with Llama.  If not, see <http://www.gnu.org/licenses/>. */

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <vector>

#include <cctk.h>
#include <util_Table.h>

#include "patchsystem.hh"
#include "tensortypes.h"



// Calculate the number of elements in an array
#define ARRSIZE(T) (sizeof(T) / sizeof(*(T)))



namespace Coordinates {
  
  using namespace std;
  
  
  
  // These can be increased if necessary, but not all necessary tensor
  // types may be defined.
  #define MAX_RANK                2 /* maximum tensor rank */
  #define MAX_TIME_LEVELS         3 /* maximum number of time levels */
  #define MAX_SPATIAL_DERIV_ORDER 2 /* maximum spatial derivative order */
  #define MAX_TIME_DERIV_ORDER    2 /* maximum time derivative order */
  
  
  extern "C"
  CCTK_INT
  Coordinates_SymmetryInterpolate (CCTK_POINTER_TO_CONST const cctkGH_,
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
    // Check arguments
    assert (N_dims == ndirs);
    assert (N_interp_points >= 0);
    assert (interp_coords_type >= 0);
    assert (interp_coords);
    for (int d=0; d<ndirs; ++d) {
      assert (N_interp_points == 0 or interp_coords[d]);
    }
    assert (N_output_arrays >= 0);
    
    // Coordinates must be CCTK_REAL
    assert (interp_coords_type == CCTK_VARIABLE_REAL);
    assert (output_array_types);
    for (int m=0; m<N_output_arrays; ++m) {
      assert (output_array_types[m] == CCTK_VARIABLE_REAL);
    }
    
    
    
    // Claim faces
    for (int f=0; f<ndirs*nfaces; ++f) {
      assert (faces & (1 << f));
    }
    CCTK_INT new_faces = faces;
    for (int f=0; f<ndirs*nfaces; ++f) {
      new_faces &= ~ (1 << f);
    }
    
    
    
    // Convert to patch-local coordinates
    CCTK_POINTER new_interp_coords[ndirs];
    CCTK_INT * new_interp_maps;
    vector<CCTK_REAL> new_interp_coords_arrays[ndirs];
    vector<CCTK_INT> new_interp_maps_array;
    for (int d=0; d<ndirs; ++d) {
      new_interp_coords_arrays[d].resize (N_interp_points);
      new_interp_coords[d] = & new_interp_coords_arrays[d].front();
    }
    new_interp_maps_array.resize (N_interp_points);
    new_interp_maps = & new_interp_maps_array.front();
    
    // Calculate basis transformation coefficients
    // (Note that these are in Fortran array order)
    vector<CCTK_REAL> dadx[ndirs][ndirs];
    vector<CCTK_REAL> ddadxdx[ndirs][ndirs][ndirs];
    CCTK_POINTER dadx_ptrs[ndirs*ndirs];
    CCTK_POINTER ddadxdx_ptrs[ndirs*ndirs*(ndirs+1)/2];
    for (int d=0; d<ndirs; ++d) {
      for (int e=0; e<ndirs; ++e) {
        dadx[e][d].resize (N_interp_points);
        for (int f=0; f<ndirs; ++f) {
          ddadxdx[f][e][d].resize (N_interp_points);
        }
      }
    }
  #if 0
    CCTK_REAL const poison = -424242.0;
    for (int d=0; d<ndirs; ++d) {
      for (int e=0; e<ndirs; ++e) {
        dadx[e][d].assign (N_interp_points, poison);
        for (int f=0; f<ndirs; ++f) {
          ddadxdx[f][e][d].assign (N_interp_points, poison);
        }
      }
    }
  #endif
    {
      size_t cnt = 0;
      for (int d=0; d<ndirs; ++d) {
        for (int e=0; e<ndirs; ++e) {
          dadx_ptrs[d*ndirs+e] = & dadx[e][d].front();
          for (int f=0; f<ndirs; ++f) {
            if (e <= f) {
              assert (cnt < ARRSIZE (ddadxdx_ptrs));
              ddadxdx_ptrs[cnt] = & ddadxdx[f][e][d].front();
              ++ cnt;
            }
          }
        }
      }
      assert (cnt == ARRSIZE (ddadxdx_ptrs));
    }
    
    {
      CCTK_INT const ierr =
        MultiPatch_GlobalToLocal (cctkGH_, N_dims,
                                  N_interp_points,
                                  interp_coords,
                                  new_interp_maps,
                                  new_interp_coords,
                                  dadx_ptrs,
                                  ddadxdx_ptrs);
      assert (not ierr);
    }
    for (int d=0; d<ndirs; ++d) {
      for (int e=0; e<ndirs; ++e) {
        for (int f=0; f<ndirs; ++f) {
          if (e > f) {
            ddadxdx[f][e][d] = ddadxdx[e][f][d];
          }
        }
      }
    }
    
    
    // Add map numbers to parameter table
    {
      int const icnt =
        Util_TableSetIntArray (param_table_handle, N_interp_points,
                               new_interp_maps,
                               "source_map");
      assert (icnt >= 0);
    }
    
      
    
    // Allocate new output arrays
    vector <CCTK_POINTER> new_output_arrays (N_output_arrays);
    vector <vector <CCTK_REAL> > new_output_arrays_arrays (N_output_arrays);
    for (int m=0; m<N_output_arrays; ++m) {
      assert (output_array_types[m] == CCTK_VARIABLE_REAL);
      new_output_arrays_arrays.at(m).resize (N_interp_points);
      new_output_arrays.at(m) = & new_output_arrays_arrays.at(m).front();
    }
    
    
    
    // Recursive interpolation call
    CCTK_INT const iret =
      SymmetryInterpolateFaces
      (cctkGH_, N_dims,
       local_interp_handle, param_table_handle, coord_system_handle,
       N_interp_points, interp_coords_type, new_interp_coords,
       N_input_arrays, input_array_indices,
       N_output_arrays, output_array_types, & new_output_arrays.front(),
       new_faces);
    
    
    
    // Remove map numbers from parameter table
    {
      int const ierr = Util_TableDeleteKey (param_table_handle, "source_map");
      assert (not ierr);
    }
    
    
    
    // Convert tensor types to the global tensor basis
    
    // Find output variable indices
    vector<CCTK_INT> input_array_time_levels (N_input_arrays);
    {
      int const icnt =
        Util_TableGetIntArray
        (param_table_handle, N_input_arrays,
         & input_array_time_levels.front(), "input_array_time_levels");
      if (icnt == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
        for (int m = 0; m < N_input_arrays; ++ m) {
          input_array_time_levels.at(m) = 0; // time level is 0
        }
      } else {
        assert (icnt == N_input_arrays);
      }
    }
    
    vector<CCTK_INT> operand_indices (N_output_arrays);
    {
      int const icnt =
        Util_TableGetIntArray
        (param_table_handle, N_output_arrays,
         & operand_indices.front(), "operand_indices");
      if (icnt == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
        assert (N_output_arrays == N_input_arrays);
        for (int m = 0; m < N_output_arrays; ++ m) {
          operand_indices.at(m)= m; // output index equals input index
        }
      } else {
        assert (icnt == N_output_arrays);
      }
    }
    
    vector<CCTK_INT> operation_codes (N_output_arrays);
    {
      int const icnt =
        Util_TableGetIntArray
        (param_table_handle, N_output_arrays,
         & operation_codes.front(), "operation_codes");
      if (icnt == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
        assert (N_output_arrays == N_input_arrays); // why?
        for (int m = 0; m < N_output_arrays; ++ m) {
          operation_codes.at(m) = 0; // do not take spatial derivatives
        }
      } else {
        assert (icnt == N_output_arrays);
      }
      for (int m = 0; m < N_output_arrays; ++ m) {
        assert (operation_codes.at(m) >= 0);
      }
    }
    
    vector<CCTK_INT> time_deriv_order (N_output_arrays);
    {
      int const icnt =
        Util_TableGetIntArray
        (param_table_handle, N_output_arrays,
         & time_deriv_order.front(), "time_deriv_order");
      if (icnt == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
        for (int m = 0; m < N_output_arrays; ++ m) {
          time_deriv_order.at(m) = 0; // do not take time derivatives
        }
      } else {
        assert (icnt == N_output_arrays);
      }
      for (int m = 0; m < N_output_arrays; ++ m) {
        assert (time_deriv_order.at(m) >= 0 and
                time_deriv_order.at(m) <= MAX_TIME_DERIV_ORDER);
      }
    }
    
    vector<CCTK_INT> output_array_indices (N_output_arrays);
    for (int m = 0; m < N_output_arrays; ++ m) {
      assert (operand_indices.at(m) >= 0 and
              operand_indices.at(m) < N_input_arrays);
      output_array_indices.at(m) = input_array_indices[operand_indices.at(m)];
      assert (output_array_indices.at(m) == -1 or
              (output_array_indices.at(m) >= 0 and
               output_array_indices.at(m) < CCTK_NumVars()));
    }
    
    
    
    // Map Cactus variables to tensor objects
    
    // [MAX_TIME_DERIV_ORDER+1][MAX_SPATIAL_DERIV_ORDER+1][MAX_TIME_LEVELS+1]
    // [NumVars]
    vector <vector <vector <vector <struct tensor const *> > > > thetensor;
    vector <vector <vector <vector <vector <int>         > > > > thevars;
    thetensor.resize (MAX_TIME_DERIV_ORDER+1);
    thevars  .resize (MAX_TIME_DERIV_ORDER+1);
    for (int p=0; p<=MAX_TIME_DERIV_ORDER; ++p) {
      thetensor.at(p).resize (MAX_SPATIAL_DERIV_ORDER+1);
      thevars  .at(p).resize (MAX_SPATIAL_DERIV_ORDER+1);
      for (int q=0; q<=MAX_SPATIAL_DERIV_ORDER; ++q) {
        thetensor.at(p).at(q).resize (MAX_TIME_LEVELS+1);
        thevars  .at(p).at(q).resize (MAX_TIME_LEVELS+1);
        for (int tl=0; tl<=MAX_TIME_LEVELS; ++tl) {
          thetensor.at(p).at(q).at(tl).resize (CCTK_NumVars());
          thevars  .at(p).at(q).at(tl).resize (CCTK_NumVars());
          for (int n=0; n<CCTK_NumVars(); ++n) {
            thetensor.at(p).at(q).at(tl).at(n) = NULL;
          }
        }
      }
    }
    
    // Map output arrays to the base Cactus variable (i.e., the first
    // component in the tensor objects)
    vector<int> thevar (N_output_arrays, -1);
    for (int m=0; m<N_output_arrays; ++m) {
      if (output_array_indices.at(m) != -1) {
        
        // Get some variable information
        int const vi = output_array_indices.at(m);
        assert (vi>=0 and vi<CCTK_NumVars());
        int const gi = CCTK_GroupIndexFromVarI (vi);
        assert (gi>=0 and gi<CCTK_NumGroups());
        
        // Find the tensor type
        struct tensor const * const tensortype = &TT_scalar;
        assert (tensortype);
        int const var = 0;
        
        int const basevar = vi;
        
        // Determine the component
        int const component = tensortype->comps[var];
        
        // Translate into indices
        vector<int> indices (MAX_RANK+1);
        TT_component2indices (tensortype->dim, tensortype->rank, component,
                              & indices.front());
        
        // Take time derivative order (i.e., time derivatives) into account
        int const num_time_derivs = time_deriv_order.at(m);
        assert (num_time_derivs>=0 and num_time_derivs<=MAX_TIME_DERIV_ORDER);
        
        // Take operation code (i.e., spatial derivatives) into account
        int const num_derivs = TT_derivcode2num_derivs (operation_codes.at(m));
        assert (num_derivs>=0 and num_derivs<=MAX_SPATIAL_DERIV_ORDER);
        
        // Get the time level
        int const time_level =
          input_array_time_levels.at(operand_indices.at(m));
        assert (time_level>=0 and time_level<=MAX_TIME_LEVELS);
        
        // Find the derivative of the tensor type
        struct tensor const * const derivtensortype =
          TT_derivative (tensortype, num_derivs);
        
        // Find the additional indices
        {
          assert (derivtensortype->rank>=0 and derivtensortype->rank<=MAX_RANK);
          assert (derivtensortype->rank >= tensortype->rank);
          int const code = operation_codes.at(m);
          TT_derivcode2derivs
            (code, num_derivs, & indices.at(tensortype->rank));
          for (int r=tensortype->rank; r<derivtensortype->rank; ++r) {
            assert (indices.at(r)>=0 and indices.at(r)<derivtensortype->dim);
          }
        }
        
        // Re-calculate component
        int const derivcomponent =
          TT_indices2component (derivtensortype->dim, derivtensortype->rank,
                                & indices.front());
        
        int const derivvar = derivtensortype->vars[derivcomponent];
        assert (derivvar>=0 and derivvar<derivtensortype->nvars);
        thevar.at(m) = derivvar;
        
        // Create or cross-check the tensor object
        if (not thetensor.at(num_time_derivs).at(num_derivs).at(time_level).at(basevar)) {
          thetensor.at(num_time_derivs).at(num_derivs).at(time_level).at(basevar) = derivtensortype;
          assert (thevars.at(num_time_derivs).at(num_derivs).at(time_level).at(basevar).empty());
          thevars.at(num_time_derivs).at(num_derivs).at(time_level).at(basevar).resize (derivtensortype->nvars);
          for (int i=0; i<derivtensortype->nvars; ++i) {
            thevars.at(num_time_derivs).at(num_derivs).at(time_level).at(basevar).at(i) = -1;
          }
        }
        assert (thetensor.at(num_time_derivs).at(num_derivs).at(time_level).at(basevar) == derivtensortype);
        // This does not hold if the caller requests the same
        // interpolation to be done into different output arrays.
        // This may happen e.g. when CarpetInterp needs to
        // differentiate in time.  This is arguably a performance bug
        // in CarpetInterp.  (See whether this goes away now.)  (It
        // should!)
        assert (thevars.at(num_time_derivs).at(num_derivs).at(time_level).at(basevar).at(derivvar) == -1);
        thevars.at(num_time_derivs).at(num_derivs).at(time_level).at(basevar).at(derivvar) = m;
        
      } // if != -1
    } // for m
    
    
    
    
    /*CCTK_INT const ierr2 =
      MultiPatch_LocalToGlobal (cctkGH_, N_dims,
                                N_interp_points,
                                new_interp_maps,
                                new_interp_coords,
                                NULL,
                                NULL,
                                NULL,
                                dadx_ptrs,
                                NULL,
                                ddadxdx_ptrs,
                                NULL);*/
    
  #if 0
    for (int n=0; n<N_interp_points; ++n) {
      for (int d=0; d<ndirs; ++d) {
        for (int e=0; e<ndirs; ++e) {
          assert (dadx[e][d].at(n) != poison);
          for (int f=0; f<ndirs; ++f) {
            assert (ddadxdz[f][e][d].at(n) != poison);
          }
        }
      }
    }
  #endif
    
    
    
    // Loop over all output arrays
    for (int m=0; m<N_output_arrays; ++m) {
      if (output_array_indices.at(m) != -1) {
        
        // Take time derivative order (i.e., time derivatives) into account
        int const num_time_derivs = time_deriv_order.at(m);
        assert (num_time_derivs>=0 and num_time_derivs<=MAX_TIME_DERIV_ORDER);
        
        // Take operation code (i.e., spatial derivatives) into account
        int const num_derivs = TT_derivcode2num_derivs (operation_codes.at(m));
        assert (num_derivs>=0 and num_derivs<=MAX_SPATIAL_DERIV_ORDER);
        
        // Get the time level
        int const time_level =
          input_array_time_levels.at(operand_indices.at(m));
        assert (time_level>=0 and time_level<=MAX_TIME_LEVELS);
        
        // Get the tensor type
        int const vi = output_array_indices.at(m);
        int const basevar = vi;
        assert (basevar>=0 and basevar<=CCTK_NumVars());
        struct tensor const * const tensortype =
          thetensor.at(num_time_derivs).at(num_derivs).at(time_level).at(basevar);
        assert (tensortype);
        vector<int> const & vars =
          thevars.at(num_time_derivs).at(num_derivs).at(time_level).at(basevar);
        int const var = thevar.at(m);
        assert (var>=0 and var<tensortype->nvars);
        
        // Transform into indices
        int const component = tensortype->comps[var];
        assert (tensortype->rank>=0 and tensortype->rank<=MAX_RANK);
        vector<int> indices(MAX_RANK+1);
        TT_component2indices (tensortype->dim, tensortype->rank, component,
                              & indices.front());
        
        
        
        switch (num_derivs) {
          
        case 0: {
          // T
          assert (tensortype == &TT_scalar);
          
          // do nothing
          memcpy (output_arrays[m], new_output_arrays[m],
                  N_interp_points * sizeof (CCTK_REAL));
          break;
        }
          
        case 1: {
          // T,i = T,a dx^a/dx^i
          assert (tensortype == &TT_vector);
          
          int const i = indices[0];
          
          CCTK_REAL * restrict const T_i =
            static_cast <CCTK_REAL *> (output_arrays[m]);
          
          for (int n=0; n<N_interp_points; ++n) {
            T_i[n] = 0;
          }
          
          for (int a=0; a<ndirs; ++a) {
            
            int myindices[1];
            myindices[0] = a;
            int const mycomponent =
              TT_indices2component
              (tensortype->dim, tensortype->rank, myindices);
            int const myvar = tensortype->vars[mycomponent];
            assert (myvar >= 0);
            int const mym = vars[myvar];
            assert (mym >= 0);
            
            vector<CCTK_REAL> const & dadx_ai = dadx[a][i];
            CCTK_REAL const * restrict const T_a =
              static_cast <CCTK_REAL const *> (new_output_arrays[mym]);
            
            for (int n=0; n<N_interp_points; ++n) {
              T_i[n] += T_a[n] * dadx_ai[n];
            }
            
          } // for a
          break;
        }
          
        case 2: {
          // T,(kl)
          //    = (T,c dx^c/dx^k),l
          //    = T,cl dx^c/dx^k + T,c ddx^c/(dx^k dx^l)
          //    = T,cd dx^c/dx^k dx^d/dx^l + T,c ddx^c/(dx^k dx^l)
          assert (tensortype == &TT_symmtensor);
          
          // Get the tensor type for the first derivative
          int const num_derivs1 = num_derivs - 1;
          struct tensor const * const tensortype1 =
            thetensor.at(num_time_derivs).at(num_derivs1).at(time_level).at(basevar);
          assert (tensortype1);
          vector<int> const & vars1 =
            thevars.at(num_time_derivs).at(num_derivs1).at(time_level).at(basevar);
          
          int const k = indices[0];
          int const l = indices[1];
          
          CCTK_REAL * restrict const T_kl =
            static_cast <CCTK_REAL *> (output_arrays[m]);
          
          for (int n=0; n<N_interp_points; ++n) {
            T_kl[n] = 0;
          }
          
          for (int c=0; c<ndirs; ++c) {
            for (int d=0; d<ndirs; ++d) {
              
              int myindices[2];
              myindices[0] = c;
              myindices[1] = d;
              int const mycomponent =
                TT_indices2component
                (tensortype->dim, tensortype->rank, myindices);
              int const myvar = tensortype->vars[mycomponent];
              assert (myvar >= 0);
              int const mym = vars[myvar];
              assert (mym >= 0);
              
              CCTK_REAL const * restrict const T_cd =
                static_cast <CCTK_REAL const *> (new_output_arrays[mym]);
              
              vector<CCTK_REAL> const & dadx_ck = dadx[c][k];
              vector<CCTK_REAL> const & dadx_dl = dadx[d][l];
              
              for (int n=0; n<N_interp_points; ++n) {
                T_kl[n] += T_cd[n] * dadx_ck[n] * dadx_dl[n];
              }
              
            } // for d
            
            int myindices1[1];
            myindices1[0] = c;
            
            int const mycomponent1 =
              TT_indices2component
              (tensortype1->dim, tensortype1->rank, myindices1);
            int const myvar1 = tensortype1->vars[mycomponent1];
            assert (myvar1 >= 0);
            int const mym1 = vars1[myvar1];
            assert (mym1 >= 0);
            
            CCTK_REAL const * restrict const T_c =
              static_cast <CCTK_REAL const *> (new_output_arrays[mym1]);
            
            vector<CCTK_REAL> const & ddadxdx_ckl = ddadxdx[c][k][l];
            
            for (int n=0; n<N_interp_points; ++n) {
              T_kl[n] += T_c[n] * ddadxdx_ckl[n];
            }
            
          } // for c
          break;
        }
          
        default:
          assert (0);
        }
        
      } // if != -1
    } // for m
    
    
    
    return iret;
  }
  
} // namespace Interpolate
