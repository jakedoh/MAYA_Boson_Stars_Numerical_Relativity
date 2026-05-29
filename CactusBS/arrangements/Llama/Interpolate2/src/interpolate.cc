
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
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <list>
#include <vector>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <util_Table.h>

#include <carpet.hh>
#include <loopcontrol.h>

#include <carpetinterp2.hh>

#ifndef DECLARE_CCTK_ARGUMENTS_CHECKED
#define DECLARE_CCTK_ARGUMENTS_CHECKED(func) DECLARE_CCTK_ARGUMENTS
#endif



namespace Interpolate2 {
  
  using namespace std;
  using namespace Carpet;
  using namespace CarpetInterp2;
  
  
  
  struct jacobian_t {
    CCTK_REAL J[3][3];          // J[a][b]
    CCTK_REAL det() const
    {
      return fabs(+ J[0][0]*(+ J[1][1]*J[2][2]
                             - J[1][2]*J[2][1])
                  + J[0][1]*(+ J[1][2]*J[2][0]
                             - J[1][0]*J[2][2])
                  + J[0][2]*(+ J[1][0]*J[2][1]
                             - J[1][1]*J[2][0]));
    }
  };
  
  struct scatter_setup_t {
    int m;                      // map
    int c;                      // component
    vector<int> indices;        // 3D indices of grid points
    int npoints;
  };
  
  template <typename FASTERP>
  struct interp2_setup_gen_t {
    list<scatter_setup_t> scatter_setups;
    FASTERP * fasterp_setup;
    // Jacobians for transforming local components
    vector<jacobian_t> jacobian_th; // J[a][b] = d[there]^a / d[here]^b
    vector<jacobian_t> jacobian_ht; // J[a][b] = d[here]^a / d[there]^b
    int npoints;
    
    interp2_setup_gen_t () : fasterp_setup (NULL), npoints (0) {}
    ~interp2_setup_gen_t () { if (fasterp_setup) delete fasterp_setup; }
    
    // Forbid copying
  private:
    interp2_setup_gen_t (interp2_setup_gen_t const &);
    interp2_setup_gen_t& operator= (interp2_setup_gen_t const &);
  };
  
  typedef interp2_setup_gen_t<fasterp_setup_t> interp2_setup_t;
  
  typedef interp2_setup_gen_t<fasterp_eno2_setup_t> interp2_setup_eno2_t;
  
  
  
  // LAPACK wrapper to invert generic 3x3 matrices
  extern "C" CCTK_FCALL void CCTK_FNAME(IP2_calc_inv3)
    (CCTK_REAL const * restrict g3, CCTK_REAL * restrict gu3);
  
  
  
  // Global variable containing the interpolation setups.
  // This is the standard Lagrange interpolator
  vector<interp2_setup_t *> interp2_setups;
  
  // This is a 1st-order Lagrange matter interpolator.
  vector<interp2_setup_t *> interp2_setups_matter_1st;
  
  // This is a 2nd-order ENO matter interpolator.
  vector<interp2_setup_eno2_t *> interp2_setups_matter_2nd;
  
  
  
  // Calculate Jacobians of the coordinate transformations
  template <typename FASTERP>
  static void
  setup_coordinate_transformation (CCTK_ARGUMENTS,
                                   interp2_setup_gen_t<FASTERP> & interp2_setup)
  {
    // Do nothing if the coordinates have already been set up
    if (not interp2_setup.jacobian_ht.empty() or interp2_setup.npoints == 0) {
      return;
    }
    
    // Allocate memory for the Jacobians
    interp2_setup.jacobian_ht.resize (interp2_setup.npoints);
    interp2_setup.jacobian_th.resize (interp2_setup.npoints);
    
    // Loop over all components and maps
    list<scatter_setup_t>::const_iterator scatter_setup_it =
      interp2_setup.scatter_setups.begin();
    int pos = 0;
    BEGIN_MAP_LOOP (cctkGH, CCTK_GF) {
      BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
        assert (scatter_setup_it != interp2_setup.scatter_setups.end());
        scatter_setup_t const & scatter_setup = *scatter_setup_it;
        assert (scatter_setup.m == Carpet::map);
        assert (scatter_setup.c == Carpet::component);
        
        {
          DECLARE_CCTK_ARGUMENTS_CHECKED(Interpolate2ApplyBC);
          
          // Obtain global coordinates for all grid points
          vector<CCTK_REAL> globalcoords[3];
          CCTK_REAL * globalcoords_ptr[3];
          for (int d=0; d<3; ++d) {
            globalcoords[d].resize (scatter_setup.npoints);
            globalcoords_ptr[d] = &globalcoords[d].front();
          }
          for (int i=0; i<scatter_setup.npoints; ++i) {
            int const ind3d = scatter_setup.indices.AT(i);
            globalcoords[0].AT(i) = x[ind3d];
            globalcoords[1].AT(i) = y[ind3d];
            globalcoords[2].AT(i) = z[ind3d];
          }
          
          // Allocate memory for patch numbers, local coordinates, and
          // the "other" Jacobian of these grid points (i.e., where
          // they are interpolated from)
          vector<CCTK_INT> patch (scatter_setup.npoints);
          CCTK_INT * const patch_ptr = &patch.front();
          CCTK_REAL * localcoords_ptr[3];
          vector<CCTK_REAL> localcoords[3];
          for (int d=0; d<3; ++d) {
            localcoords[d].resize (scatter_setup.npoints);
            localcoords_ptr[d] = &localcoords[d].front();
          }
          vector<CCTK_REAL> dadx[3][3]; // dadx[d][e] = J.J[d][e]
          CCTK_REAL * dadx_ptr[9];
          for (int d=0; d<3; ++d) {
            for (int e=0; e<3; ++e) {
              dadx[e][d].resize (scatter_setup.npoints);
              dadx_ptr[d*3+e] = &dadx[e][d].front();
            }
          }
          
          // Calculate patch numbers, local coordinates, and the
          // "other" Jacobian of these grid points
          int const ierr = MultiPatch_GlobalToLocal
            (cctkGH, 3, scatter_setup.npoints,
             globalcoords_ptr,
             patch_ptr, localcoords_ptr, dadx_ptr, NULL);
          assert (not ierr);
          
          // Set up pointers to the "local" Jacobians, which are
          // stored in grid functions
          CCTK_REAL const * restrict Jh_ptr[3][3];
          // Jh[a][b] = d[local]^a/d[global]^b
          Jh_ptr[0][0] = J11;
          Jh_ptr[0][1] = J12;
          Jh_ptr[0][2] = J13;
          Jh_ptr[1][0] = J21;
          Jh_ptr[1][1] = J22;
          Jh_ptr[1][2] = J23;
          Jh_ptr[2][0] = J31;
          Jh_ptr[2][1] = J32;
          Jh_ptr[2][2] = J33;
          
          // Calculate the combined Jacobian for all grid points
          for (int i=0; i<scatter_setup.npoints; ++i) {
            int const ind3d = scatter_setup.indices.AT(i);
            
            // reference to the combined Jacobians
            jacobian_t & Jht = interp2_setup.jacobian_ht.AT(pos+i);
            jacobian_t & Jth = interp2_setup.jacobian_th.AT(pos+i);
            
            // the "here" (local) Jacobian for this grid point
            jacobian_t Jh;
            for (int d=0; d<3; ++d) {
              for (int e=0; e<3; ++e) {
                Jh.J[d][e] = Jh_ptr[d][e][ind3d];
              }
            }
            
            // the "there" (other) Jacobian for this grid point --
            // actually, the inverse of what we need later
            jacobian_t Jti;
            for (int d=0; d<3; ++d) {
              for (int e=0; e<3; ++e) {
                Jti.J[d][e] = dadx[d][e][i];
              }
            }
            
            // the "there" Jacobian
            jacobian_t Jt;
            
            // Calculate the "there" Jacobian by calling LAPACK to
            // invert the inverse Jacobian
            CCTK_FNAME(IP2_calc_inv3) (&Jti.J[0][0], &Jt.J[0][0]);
            
            // These are our definitions and their transformation
            // rules:
            
            // Jti^a_b = d[there]^a/d[global]^b
            // Jt^a_b = d[global]^a/d[there]^b
            
            // Jh^a_b = d[here]^a/d[global]^b
            
            // Jht^a_b = d[here]^a/d[there]^b
            //         = d[here]^a/d[global]^c  d[global]^c/d[there]^b
            //         = Jh^a_c Jt^c_b
            
            // Calculate the above
            for (int d=0; d<3; ++d) {
              for (int e=0; e<3; ++e) {
                Jht.J[d][e] = 0;
                for (int f=0; f<3; ++f) {
                  Jht.J[d][e] += Jh.J[d][f] * Jt.J[f][e];
                }
              }
            }
            
            // Invert this to obtain the inverse combined Jacobian
            CCTK_FNAME(IP2_calc_inv3) (&Jht.J[0][0], &Jth.J[0][0]);
            
          } // for i
          
        } // DECLARE_CCTK_ARGUMENTS_CHECKED(Interpolate2ApplyBC)
        
        ++ scatter_setup_it;
        pos += scatter_setup.npoints;
      } END_LOCAL_COMPONENT_LOOP;
    } END_MAP_LOOP;
    assert (scatter_setup_it == interp2_setup.scatter_setups.end());
    assert (pos == interp2_setup.npoints);
  }
  
  
  
  // Transform the coordinates in the interpolated grid points
  template <typename FASTERP>
  static void
  transform_coordinates (CCTK_ARGUMENTS,
                         int const nvars,
                         vector<CCTK_INT> const & indices,
                         vector<CCTK_REAL *> const & values,
                         interp2_setup_gen_t<FASTERP> & interp2_setup)
  {
    // Loop over all interpolated variables
    for (int n = 0; n < nvars; ++ n) {
      
      // Get some information for these variables
      int const vi = indices.AT(n);
      assert (vi>=0 and vi<CCTK_NumVars());
      int const gi = CCTK_GroupIndexFromVarI (vi);
      assert (gi>=0 and gi<CCTK_NumGroups());
      int const v0 = CCTK_FirstVarIndexI (gi);
      assert (v0>=0 and v0<=vi);
      
      cGroup groupdata;
      {
        int const ierr = CCTK_GroupData (gi, &groupdata);
        assert (not ierr);
      }
      
      char jacobian[100];
      {
        int const ierr = Util_TableGetString
          (groupdata.tagstable, sizeof jacobian, jacobian, "jacobian");
        if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
          strcpy (jacobian, "");
        } else if (ierr<0) {
          char * const fullname = CCTK_FullName (vi);
          CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "Error while getting entry \"jacobian\" from tags table of variable \"%s\"",
                      fullname);
          free (fullname);
        }
      }
      
      // If the grid variable does not have a Jacobian in its tags, do
      // nothing
      if (strcmp (jacobian, "") != 0) {
        
        // Set up the coordinate transformation
        setup_coordinate_transformation (CCTK_PASS_CTOC, interp2_setup);
        
        // Find the tensor type
        char tensortypealias[100];
        // Only scalars, vectors, and co-vectors are supported yet, as
        // this is what hydro needs.  Support for other types could
        // easily be added.
        enum tensortype_t { tt_unknown, tt_scalar, tt_d, tt_u };
        tensortype_t tensortype;
        {
          int const ierr = Util_TableGetString
            (groupdata.tagstable,
             sizeof tensortypealias, tensortypealias, "tensortypealias");
          if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
            tensortype = tt_unknown;
          } else if (ierr<0) {
            CCTK_VERROR("Could not obtain 'tensortypealias' tag: %d", ierr);
          } else {
            if (CCTK_EQUALS (tensortypealias, "scalar")) {
              tensortype = tt_scalar;
            } else if (CCTK_EQUALS (tensortypealias, "d")) {
              tensortype = tt_d;
            } else if (CCTK_EQUALS (tensortypealias, "u")) {
              tensortype = tt_u;
            } else {
              tensortype = tt_unknown;
            }
          }
        }
        
        CCTK_REAL tensorweight = NAN;
        {
           int const ierr = Util_TableGetReal
            (groupdata.tagstable,
             &tensorweight, "tensorweight");
          if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
            tensorweight = 0;
          } else if (ierr<0) {
            CCTK_VERROR("Could not obtain 'tensorweight' tag: %d", ierr);
          } else {
            // nothing to do anymore
          }
        }
        
        switch (tensortype) {
          
        case tt_scalar: {
          // do nothing, if tensorweight == 0
          if (tensorweight == 0.0) break;
          for (int i=0; i<interp2_setup.npoints; ++i) {
             jacobian_t const & J = interp2_setup.jacobian_th.AT(i);
             values.at(n)[i] *= pow(J.det(), tensorweight);
          }
          break;
        }
          
        case tt_d: {
          assert (groupdata.numvars == 3);
          if (vi == v0) {
            
            // find all components
            int comp[3];
            for (int d=0; d<3; ++d) {
              for (int nn=0; nn<nvars; ++nn) {
                if (indices.AT(nn) == v0+d) {
                  comp[d] = nn;
                  goto found_d;
                }
              }
              CCTK_VERROR("Could not find x,y,z components for group '%s' among the interpolated variables",
                          CCTK_GroupNameFromVarI(v0));
            found_d:;
            }
            
            // get variable pointers
            assert (groupdata.vartype == CCTK_VARIABLE_REAL);
            CCTK_REAL * restrict ptrs[3];
            for (int d=0; d<3; ++d) {
              ptrs[d] = values.at(comp[d]);
              assert (interp2_setup.npoints==0 or ptrs[d]);
            }
            
            // transform all interpolated points
            for (int i=0; i<interp2_setup.npoints; ++i) {
              CCTK_REAL there[3];
              for (int d=0; d<3; ++d) {
                there[d] = ptrs[d][i];
              }
              jacobian_t const & J = interp2_setup.jacobian_th.AT(i);
              CCTK_REAL here[3];
              for (int d=0; d<3; ++d) {
                here[d] = 0;
                for (int e=0; e<3; ++e) {
                  // [here]_d = d[there]^e/d[here]^d [there]_e
                  here[d] += J.J[e][d] * there[e];
                }
              }
              jacobian_t const & Jdet = interp2_setup.jacobian_th.AT(i);
              const CCTK_REAL det =
                tensorweight == 0 ? 1.0 : pow(Jdet.det(), tensorweight);
              for (int d=0; d<3; ++d) {
                ptrs[d][i] = here[d] * det;
              }
            }
            
          } // if vi==v0
          break;
        } // case tt_d
          
        case tt_u: {
          assert (groupdata.numvars == 3);
          if (vi == v0) {
            
            // find all components
            int comp[3];
            for (int d=0; d<3; ++d) {
              for (int nn=0; nn<nvars; ++nn) {
                if (indices.AT(nn) == v0+d) {
                  comp[d] = nn;
                  goto found_u;
                }
              }
              CCTK_VERROR("Could not find x,y,z components for group '%s' among the interpolated variables",
                          CCTK_GroupNameFromVarI(v0));
            found_u:;
            }
            
            // get variable pointers
            assert (groupdata.vartype == CCTK_VARIABLE_REAL);
            CCTK_REAL * restrict ptrs[3];
            for (int d=0; d<3; ++d) {
              ptrs[d] = values.at(comp[d]);
              assert (interp2_setup.npoints==0 or ptrs[d]);
            }
            
            // transform all interpolated points
            for (int i=0; i<interp2_setup.npoints; ++i) {
              CCTK_REAL there[3];
              for (int d=0; d<3; ++d) {
                there[d] = ptrs[d][i];
              }
              jacobian_t const & J = interp2_setup.jacobian_ht.AT(i);
              CCTK_REAL here[3];
              for (int d=0; d<3; ++d) {
                here[d] = 0;
                for (int e=0; e<3; ++e) {
                  // [here]^d = d[here]^d/d[there]^e [there]^e
                  here[d] += J.J[d][e] * there[e];
                }
              }
              jacobian_t const & Jdet = interp2_setup.jacobian_th.AT(i);
              const CCTK_REAL det =
                tensorweight == 0 ? 1.0 : pow(Jdet.det(), tensorweight);
              for (int d=0; d<3; ++d) {
                ptrs[d][i] = here[d] * det;
              }
            }
            
          } // if vi==v0
          break;
        } // case tt_u
          
        default: {
          char * const fullname = CCTK_FullName (vi);
          CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "Unknown or unsupported tensor type alias \"%s\" for variable \"%s\"",
                      tensortypealias, fullname);
          free (fullname);
        }
          
        } // switch tensortype
        
      } // if jacobian != ""
      
    } // for n
  }
  
  
  
  extern "C"
  void
  Interpolate2ApplyBC (CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_PARAMETERS;

    assert (cctkGH);
    
    assert (cctkGH->cctk_dim == 3);
    
    
    
    int nvars = Boundary_SelectedGVs (cctkGH, 0, 0, 0, 0, 0, 0);
    assert (nvars>=0);
    if (verbose) CCTK_VInfo (CCTK_THORNSTRING,
                             "Interpolating to %d variables", nvars);
    if (nvars == 0) return;
    
    vector<CCTK_INT> indices (nvars);
    vector<CCTK_INT> faces   (nvars);
    vector<CCTK_INT> widths  (nvars);
    vector<CCTK_INT> tables  (nvars);
    
    {
      int const iret =
        Boundary_SelectedGVs
        (cctkGH, nvars,
         & indices.front(), & faces.front(), & widths.front(), & tables.front(),
         0);
      assert (iret == nvars);
    }
    
    // Initially, we assume no matter variables, i.e. we use the standard interpolator!
    int nvars_matter = 0;
    vector<CCTK_INT> indices_matter (0);
    
    // Check all variables and get interpolator
    for (int i=0; i<nvars; ++i) {
      
      int const vi = indices[i];
      assert (vi>=0 and vi<CCTK_NumVars());
      
      assert (widths[i] >= 0);
      
      int const gi = CCTK_GroupIndexFromVarI (vi);
      assert (gi>=0 and gi<CCTK_NumGroups());
      
      cGroup group;
      int const ierr = CCTK_GroupData (gi, &group);
      assert (not ierr);
      
      assert (group.grouptype == CCTK_GF);
      assert (group.vartype == CCTK_VARIABLE_REAL);
      assert (group.disttype == CCTK_DISTRIB_DEFAULT);
      
      // If requested, get interpolator from variables and sort variable indices accordingly
      if (interpolator_order_matter >= 0)
      {
        char interpolator[100];
        int const ierr = Util_TableGetString
          (group.tagstable, sizeof interpolator, interpolator, "interpolator");
        if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
          strcpy (interpolator, "");
        } else if (ierr<0) {
          char * const fullname = CCTK_FullName (vi);
          CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "Error while getting entry \"interpolator\" from tags table of variable \"%s\"",
                      fullname);
          free (fullname);
        }
        if (CCTK_EQUALS(interpolator, "matter"))
        {
           // remove var from standard Lagrange interpolation
           // and add to matter interpolator!
           --nvars;
           ++nvars_matter;
           indices_matter.push_back(indices[i]);
           indices.erase(indices.begin()+i);
           --i;
        }
      }
      
    }
    
    
    
    // Work only on the coarsest grid
    if (Carpet::reflevel > 0) return;
    
    
    // Check, if standard Lagrange interpolator has
    // something to do, or if it is the matter interpolator
    // that needs to do all the work!
    if (nvars != 0)
    {
    
    // Adjust size of interpolation setups
    interp2_setups.resize (Carpet::maxreflevels, NULL);
    
    
    
    assert (is_level_mode());
    if (not interp2_setups.AT(Carpet::reflevel) or
        interp2_setups.AT(Carpet::reflevel)->fasterp_setup->outofdate()) {
      if (interp2_setups.AT(Carpet::reflevel))
        delete interp2_setups.AT(Carpet::reflevel);
      interp2_setups.AT(Carpet::reflevel) = new interp2_setup_t;
      CCTK_VInfo (CCTK_THORNSTRING,
                  "Setting up interpolation for level %d", Carpet::reflevel);
    }
    interp2_setup_t & interp2_setup = * interp2_setups.AT(Carpet::reflevel);
    
    if (not interp2_setup.fasterp_setup) {
      if (verbose) CCTK_INFO ("Preparing setting up interpolation");
      
      assert (interp2_setup.scatter_setups.empty());
      assert (interp2_setup.npoints == 0);
      
      // Count points
      if (verbose) CCTK_INFO ("Counting grid points");
      
      BEGIN_LOCAL_MAP_LOOP (cctkGH, CCTK_GF) {
        BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
          DECLARE_CCTK_ARGUMENTS_CHECKED(Interpolate2ApplyBC);
          
          int npoints = 0;
#pragma omp parallel reduction (+: npoints)
          LC_LOOP3 (count_points,
                    i, j, k,
                    0, 0, 0,
                    cctk_lsh[0], cctk_lsh[1], cctk_lsh[2],
                    cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
          {
            int const ind3d = CCTK_GFINDEX3D (cctkGH, i, j, k);
            if (Sn[ind3d] >= 0) {
              ++ npoints;
            }
          } LC_ENDLOOP3 (count_points);
          
          scatter_setup_t scatter_setup;
          scatter_setup.m = Carpet::map;
          scatter_setup.c = Carpet::component;
          scatter_setup.npoints = npoints;
          
          interp2_setup.scatter_setups.push_back (scatter_setup);
          interp2_setup.npoints += scatter_setup.npoints;
          
        } END_LOCAL_COMPONENT_LOOP;
      } END_LOCAL_MAP_LOOP;
      
      // Allocate storage for coordinates and values
      fasterp_glocs_t locations (interp2_setup.npoints);
      
      // Collect coordinates
      if (verbose) CCTK_INFO ("Collecting coordinates");
      
      list<scatter_setup_t>::iterator scatter_setup_it =
        interp2_setup.scatter_setups.begin();
      int pos = 0;
      BEGIN_LOCAL_MAP_LOOP (cctkGH, CCTK_GF) {
        BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
          assert (scatter_setup_it != interp2_setup.scatter_setups.end());
          scatter_setup_t & scatter_setup = *scatter_setup_it;
          assert (scatter_setup.m == Carpet::map);
          assert (scatter_setup.c == Carpet::component);
          
          DECLARE_CCTK_ARGUMENTS_CHECKED(Interpolate2ApplyBC);
          
          int const di = 1;
          assert(di ==
                 CCTK_GFINDEX3D(cctkGH, 1,0,0) - CCTK_GFINDEX3D(cctkGH, 0,0,0));
          int const dj =
            CCTK_GFINDEX3D(cctkGH, 0,1,0) - CCTK_GFINDEX3D(cctkGH, 0,0,0);
          
          scatter_setup.indices.reserve (scatter_setup.npoints);
          
          // This loop is not parallel
          //#pragma omp parallel reduction (+: pos)
          LC_LOOP3 (collect_coordinates,
                    i, j, k,
                    0, 0, 0,
                    cctk_lsh[0], cctk_lsh[1], cctk_lsh[2],
                    cctk_ash[0], cctk_ash[1], cctk_ash[2])
          {
            int const ind3d = CCTK_GFINDEX3D (cctkGH, i, j, k);
            if (Sn[ind3d] >= 0) {
              locations.coords[0].AT(pos) = x[ind3d];
              locations.coords[1].AT(pos) = y[ind3d];
              locations.coords[2].AT(pos) = z[ind3d];
              if (shift_edges) {
                bool const ilo = cctk_lbnd[0] + i == 0;
                bool const jlo = cctk_lbnd[1] + j == 0;
                bool const klo = cctk_lbnd[2] + k == 0;
                bool const ihi = cctk_lbnd[0] + i == cctk_gsh[0]-1;
                bool const jhi = cctk_lbnd[1] + j == cctk_gsh[1]-1;
                bool const khi = cctk_lbnd[2] + k == cctk_gsh[2]-1;
                if (ilo && (jlo || jhi || klo || khi)) {
                  locations.coords[0].AT(pos) = x[ind3d+di];
                  locations.coords[1].AT(pos) = y[ind3d+di];
                  locations.coords[2].AT(pos) = z[ind3d+di];
                } else if (ihi && (jlo || jhi || klo || khi)) {
                  locations.coords[0].AT(pos) = x[ind3d-di];
                  locations.coords[1].AT(pos) = y[ind3d-di];
                  locations.coords[2].AT(pos) = z[ind3d-di];
                } else if (jlo && (klo || khi)) {
                  locations.coords[0].AT(pos) = x[ind3d+dj];
                  locations.coords[1].AT(pos) = y[ind3d+dj];
                  locations.coords[2].AT(pos) = z[ind3d+dj];
                } else if (jhi && (klo || khi)) {
                  locations.coords[0].AT(pos) = x[ind3d-dj];
                  locations.coords[1].AT(pos) = y[ind3d-dj];
                  locations.coords[2].AT(pos) = z[ind3d-dj];
                }
              } // if shift_edges
              scatter_setup.indices.push_back (ind3d);
              ++ pos;
            }
          } LC_ENDLOOP3 (collect_coordinates);
          
          assert (int(scatter_setup.indices.size()) == scatter_setup.npoints);

          ++ scatter_setup_it;
        } END_LOCAL_COMPONENT_LOOP;
      } END_LOCAL_MAP_LOOP;
      assert (scatter_setup_it == interp2_setup.scatter_setups.end());
      assert (pos == interp2_setup.npoints);
      
      // Create the interpolation setup
      if (verbose) {
        CCTK_VInfo (CCTK_THORNSTRING,
                    "Creating interpolation setup for %d grid points",
                    interp2_setup.npoints);
      }
      
      assert (not interp2_setup.fasterp_setup);
      interp2_setup.fasterp_setup =
        new fasterp_setup_t (cctkGH, locations, interpolator_order);
      
    } // if not interp2_setup.fasterp_setup
    
    
    
    // Poison boundaries, so that they are not accidentally used as
    // interpolation source
    if (verbose) CCTK_INFO ("Poisoning inter-patch boundaries");
    
    {
      list<scatter_setup_t>::const_iterator scatter_setup_it =
        interp2_setup.scatter_setups.begin();
      BEGIN_LOCAL_MAP_LOOP (cctkGH, CCTK_GF) {
        BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
          assert (scatter_setup_it != interp2_setup.scatter_setups.end());
          scatter_setup_t const & scatter_setup = *scatter_setup_it;
          assert (scatter_setup.m == Carpet::map);
          assert (scatter_setup.c == Carpet::component);
          
          DECLARE_CCTK_ARGUMENTS_CHECKED(Interpolate2ApplyBC);
          
          for (int n = 0; n < nvars; ++ n) {
            
            CCTK_REAL * restrict const var =
              static_cast <CCTK_REAL *>
              (CCTK_VarDataPtrI (cctkGH, 0, indices.AT(n)));
            assert (var);
            
#pragma omp parallel for
            for (int i=0; i<int(scatter_setup.indices.size()); ++i) {
              var[scatter_setup.indices.AT(i)] = NAN; // poison;
            }
            
#if 0
#pragma omp parallel
            LC_LOOP3 (check_poison_before,
                      i, j, k,
                      0, 0, 0,
                      cctk_lsh[0], cctk_lsh[1], cctk_lsh[2],
                      cctk_ash[0], cctk_ash[1], cctk_ash[2])
            {
              int const ind3d = CCTK_GFINDEX3D (cctkGH, i, j, k);
              if (Sn[ind3d] >= 0) {
                assert(CCTK_isnan(var[ind3d]) || fabs(var[ind3d]) >= 1.0e+10);
              } else {
                assert(!CCTK_isnan(var[ind3d]) && fabs(var[ind3d]) <= 1.0e+100);
              }
            } LC_ENDLOOP3 (check_poison_before);
#endif
            
          } // for n
          
          ++ scatter_setup_it;
        } END_LOCAL_COMPONENT_LOOP;
      } END_LOCAL_MAP_LOOP;
      assert (scatter_setup_it == interp2_setup.scatter_setups.end());
    }
    
    
    
    // Interpolate
    if (verbose) CCTK_INFO ("Interpolating");
    
    vector<int> varinds (nvars);
    for (int n=0; n<nvars; ++n) {
      varinds.AT(n) = indices.AT(n);
    }
    
    vector<vector<CCTK_REAL> > valuess (nvars);
    vector<CCTK_REAL *> values (nvars);
    for (int n=0; n<nvars; ++n) {
      valuess.AT(n).resize (interp2_setup.npoints, NAN); // poison;
      values.AT(n) = & valuess.AT(n).front();
    }
    
    interp2_setup.fasterp_setup->interpolate (cctkGH, varinds, values);
    
    
    
    // Transform coordinates in result
    if (verbose) CCTK_INFO ("Transforming coordinates");
    
    transform_coordinates (CCTK_PASS_CTOC, nvars, indices, values, interp2_setup);
    
    
    
    // Write back result
    if (verbose) CCTK_INFO ("Writing back interpolation result");
    
    {
      list<scatter_setup_t>::const_iterator scatter_setup_it =
        interp2_setup.scatter_setups.begin();
      int pos = 0;
      BEGIN_LOCAL_MAP_LOOP (cctkGH, CCTK_GF) {
        BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
          assert (scatter_setup_it != interp2_setup.scatter_setups.end());
          scatter_setup_t const & scatter_setup = *scatter_setup_it;
          assert (scatter_setup.m == Carpet::map);
          assert (scatter_setup.c == Carpet::component);
          
          DECLARE_CCTK_ARGUMENTS_CHECKED(Interpolate2ApplyBC);
          
          for (int n = 0; n < nvars; ++ n) {
            
            if (verbose) {
              char * restrict const fullname = CCTK_FullName (indices.AT(n));
              const int patch = MultiPatch_GetMap(cctkGH);
              CCTK_VInfo (CCTK_THORNSTRING,
                          "Filling variable #%d %s patch %d",
                          n, fullname, patch);
              free (fullname);
            }
            CCTK_REAL * restrict const var =
              static_cast <CCTK_REAL *>
              (CCTK_VarDataPtrI (cctkGH, 0, indices.AT(n)));
            assert (var);
            
            if (interpolate_zero) {
              // for debugging only
              
#pragma omp parallel for
              for (int i=0; i<scatter_setup.npoints; ++i) {
#ifdef DEBUG
                assert (Sn[scatter_setup.indices.AT(i)] >= 0);
#endif
                var[scatter_setup.indices.AT(i)] = 0.0;
              }
              
            } else {
              // regular case
              
#pragma omp parallel for
              for (int i=0; i<scatter_setup.npoints; ++i) {
#ifdef DEBUG
                assert (Sn[scatter_setup.indices.AT(i)] >= 0);
#endif
                var[scatter_setup.indices.AT(i)] = valuess.AT(n).AT(pos + i);
              }
              
            }
            
          } // for n
          
          ++ scatter_setup_it;
          pos += scatter_setup.npoints;
        } END_LOCAL_COMPONENT_LOOP;
      } END_LOCAL_MAP_LOOP;
      assert (scatter_setup_it == interp2_setup.scatter_setups.end());
      assert (pos == interp2_setup.npoints);
    }
    
    } // nvars != 0
    
    
    // Treat matter variables differently...?
    if (interpolator_order_matter < 0) {
       if (verbose) CCTK_INFO ("Done.");
       return;
    }
    
    // ...no? Then let's do special matter interpolation!
    
    /////////////////////////////////////////////////////////////////////////////////////
    //
    //  Standard Lagrange interpolation ends.
    //  Below, we exclusively interpolate matter variables with another interp_setup!
    //
    ////////////////////////////////////////////////////////////////////////////////////
    
    
    
    if (verbose) CCTK_VInfo (CCTK_THORNSTRING,
                             "Interpolating to %d matter variables", nvars_matter);
    if (nvars_matter == 0) {
       if (verbose) CCTK_INFO ("Done.");
       return;
    }
    
    //
    // 1st-order Lagrange interpolation!
    //
    
    if (interpolator_order_matter == 1)
    {
    
    typedef interp2_setup_t interp2_setup_matter_t;
    
    // Adjust size of interpolation setups
    interp2_setups_matter_1st.resize (Carpet::maxreflevels, NULL);
    
    
    
    assert (is_level_mode());
    if (not interp2_setups_matter_1st.AT(Carpet::reflevel) or
        interp2_setups_matter_1st.AT(Carpet::reflevel)->fasterp_setup->outofdate()) {
      if (interp2_setups_matter_1st.AT(Carpet::reflevel))
        delete interp2_setups_matter_1st.AT(Carpet::reflevel);
      interp2_setups_matter_1st.AT(Carpet::reflevel) = new interp2_setup_matter_t;
      CCTK_VInfo (CCTK_THORNSTRING,
                  "Setting up matter interpolation for level %d", Carpet::reflevel);
    }
    interp2_setup_matter_t & interp2_setup_matter = * interp2_setups_matter_1st.AT(Carpet::reflevel);
    
    if (not interp2_setup_matter.fasterp_setup) {
      if (verbose) CCTK_INFO ("Preparing setting up matter interpolation");
      
      assert (interp2_setup_matter.scatter_setups.empty());
      assert (interp2_setup_matter.npoints == 0);
      
      // Count points
      if (verbose) CCTK_INFO ("Counting grid points for matter interpolation");
      
      BEGIN_LOCAL_MAP_LOOP (cctkGH, CCTK_GF) {
        BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
          DECLARE_CCTK_ARGUMENTS_CHECKED(Interpolate2ApplyBC);
          
          int npoints = 0;
#pragma omp parallel reduction (+: npoints)
          LC_LOOP3 (count_points,
                    i, j, k,
                    0, 0, 0,
                    cctk_lsh[0], cctk_lsh[1], cctk_lsh[2],
                    cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
          {
            int const ind3d = CCTK_GFINDEX3D (cctkGH, i, j, k);
            if (Sn[ind3d] >= 0) {
              ++ npoints;
            }
          } LC_ENDLOOP3 (count_points);
          
          scatter_setup_t scatter_setup;
          scatter_setup.m = Carpet::map;
          scatter_setup.c = Carpet::component;
          scatter_setup.npoints = npoints;
          
          interp2_setup_matter.scatter_setups.push_back (scatter_setup);
          interp2_setup_matter.npoints += scatter_setup.npoints;
          
        } END_LOCAL_COMPONENT_LOOP;
      } END_LOCAL_MAP_LOOP;
      
      // Allocate storage for coordinates and values
      fasterp_glocs_t locations (interp2_setup_matter.npoints);
      
      // Collect coordinates
      if (verbose) CCTK_INFO ("Collecting coordinates for matter interpolation");
      
      list<scatter_setup_t>::iterator scatter_setup_it =
        interp2_setup_matter.scatter_setups.begin();
      int pos = 0;
      BEGIN_LOCAL_MAP_LOOP (cctkGH, CCTK_GF) {
        BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
          assert (scatter_setup_it != interp2_setup_matter.scatter_setups.end());
          scatter_setup_t & scatter_setup = *scatter_setup_it;
          assert (scatter_setup.m == Carpet::map);
          assert (scatter_setup.c == Carpet::component);
          
          DECLARE_CCTK_ARGUMENTS_CHECKED(Interpolate2ApplyBC);
          
          scatter_setup.indices.reserve (scatter_setup.npoints);
          
          // This loop is not parallel
          //#pragma omp parallel reduction (+: pos)
          LC_LOOP3 (collect_coordinates,
                    i, j, k,
                    0, 0, 0,
                    cctk_lsh[0], cctk_lsh[1], cctk_lsh[2],
                    cctk_ash[0], cctk_ash[1], cctk_ash[2])
          {
            int const ind3d = CCTK_GFINDEX3D (cctkGH, i, j, k);
            if (Sn[ind3d] >= 0) {
              locations.coords[0].AT(pos) = x[ind3d];
              locations.coords[1].AT(pos) = y[ind3d];
              locations.coords[2].AT(pos) = z[ind3d];
              scatter_setup.indices.push_back (ind3d);
              ++ pos;
            }
          } LC_ENDLOOP3 (collect_coordinates);
          
          assert (int(scatter_setup.indices.size()) == scatter_setup.npoints);
          
          ++ scatter_setup_it;
        } END_LOCAL_COMPONENT_LOOP;
      } END_LOCAL_MAP_LOOP;
      assert (scatter_setup_it == interp2_setup_matter.scatter_setups.end());
      assert (pos == interp2_setup_matter.npoints);
      
      // Create the interpolation setup
      if (verbose) {
        CCTK_VInfo (CCTK_THORNSTRING,
                    "Creating matter interpolation setup for %d grid points",
                    interp2_setup_matter.npoints);
      }
      
      assert (not interp2_setup_matter.fasterp_setup);
      interp2_setup_matter.fasterp_setup =
        new fasterp_setup_t (cctkGH, locations, interpolator_order_matter);
      
    } // if not interp2_setup_matter.fasterp_setup
    
    
    
    // Poison boundaries, so that they are not accidentally used as
    // interpolation source
    if (verbose) CCTK_INFO ("Poisoning inter-patch boundaries (matter interpolator)");
    
    {
      list<scatter_setup_t>::const_iterator scatter_setup_it =
        interp2_setup_matter.scatter_setups.begin();
      BEGIN_LOCAL_MAP_LOOP (cctkGH, CCTK_GF) {
        BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
          assert (scatter_setup_it != interp2_setup_matter.scatter_setups.end());
          scatter_setup_t const & scatter_setup = *scatter_setup_it;
          assert (scatter_setup.m == Carpet::map);
          assert (scatter_setup.c == Carpet::component);
          
          DECLARE_CCTK_ARGUMENTS_CHECKED(Interpolate2ApplyBC);
          
          for (int n = 0; n < nvars_matter; ++ n) {
            
            CCTK_REAL * restrict const var =
              static_cast <CCTK_REAL *>
              (CCTK_VarDataPtrI (cctkGH, 0, indices_matter.AT(n)));
            assert (var);
            
#pragma omp parallel for
            for (int i=0; i<int(scatter_setup.indices.size()); ++i) {
              var[scatter_setup.indices.AT(i)] = poison;
            }
            
          } // for n
          
          ++ scatter_setup_it;
        } END_LOCAL_COMPONENT_LOOP;
      } END_LOCAL_MAP_LOOP;
      assert (scatter_setup_it == interp2_setup_matter.scatter_setups.end());
    }
    
    
    
    // Interpolate
    if (verbose) CCTK_INFO ("Interpolating with matter interpolator");
    
    vector<int> varinds (nvars_matter);
    for (int n=0; n<nvars_matter; ++n) {
      varinds.AT(n) = indices_matter.AT(n);
    }
    
    vector<vector<CCTK_REAL> > valuess (nvars_matter);
    vector<CCTK_REAL *> values (nvars_matter);
    for (int n=0; n<nvars_matter; ++n) {
      valuess.AT(n).resize (interp2_setup_matter.npoints);
#ifndef NDEBUG
#pragma omp parallel for
      for (int i=0; i<int(valuess.AT(n).size()); ++i) {
        valuess.AT(n).AT(i) = poison;
      }
#endif
      values.AT(n) = & valuess.AT(n).front();
    }
    
    interp2_setup_matter.fasterp_setup->interpolate (cctkGH, varinds, values);
    
    
    
    // Transform coordinates in result
    if (verbose) CCTK_INFO ("Transforming coordinates of matter variables");
    
    transform_coordinates (CCTK_PASS_CTOC, nvars_matter, indices_matter, values, interp2_setup_matter);
    
    
    
    // Write back result
    if (verbose) CCTK_INFO ("Writing back matter interpolation result");
    
    {
      list<scatter_setup_t>::const_iterator scatter_setup_it =
        interp2_setup_matter.scatter_setups.begin();
      int pos = 0;
      BEGIN_LOCAL_MAP_LOOP (cctkGH, CCTK_GF) {
        BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
          assert (scatter_setup_it != interp2_setup_matter.scatter_setups.end());
          scatter_setup_t const & scatter_setup = *scatter_setup_it;
          assert (scatter_setup.m == Carpet::map);
          assert (scatter_setup.c == Carpet::component);
          
          DECLARE_CCTK_ARGUMENTS_CHECKED(Interpolate2ApplyBC);
          
          for (int n = 0; n < nvars_matter; ++ n) {
            
            if (verbose) {
              char * restrict const fullname = CCTK_FullName (indices_matter.AT(n));
              const int patch = MultiPatch_GetMap(cctkGH);
              CCTK_VInfo (CCTK_THORNSTRING,
                          "Filling matter variable #%d %s patch %d",
                          n, fullname, patch);
              free (fullname);
            }
            CCTK_REAL * restrict const var =
              static_cast <CCTK_REAL *>
              (CCTK_VarDataPtrI (cctkGH, 0, indices_matter.AT(n)));
            assert (var);
            
            if (interpolate_zero) {
              // for debugging only
              
#pragma omp parallel for
              for (int i=0; i<scatter_setup.npoints; ++i) {
#ifdef DEBUG
                assert (Sn[scatter_setup.indices.AT(i)] >= 0);
#endif
                var[scatter_setup.indices.AT(i)] = 0.0;
              }
              
            } else {
              // regular case
              
#pragma omp parallel for
              for (int i=0; i<scatter_setup.npoints; ++i) {
#ifdef DEBUG
                assert (Sn[scatter_setup.indices.AT(i)] >= 0);
#endif
                var[scatter_setup.indices.AT(i)] = valuess.AT(n).AT(pos + i);
              }
              
            }
            
          } // for n
          
          ++ scatter_setup_it;
          pos += scatter_setup.npoints;
        } END_LOCAL_COMPONENT_LOOP;
      } END_LOCAL_MAP_LOOP;
      assert (scatter_setup_it == interp2_setup_matter.scatter_setups.end());
      assert (pos == interp2_setup_matter.npoints);
    }
    
    } // end 1st-order Lagrange
    
    
    //
    // 2nd-order ENO interpolation!
    //
    
    else if (interpolator_order_matter == 2)
    {
    
    typedef interp2_setup_eno2_t interp2_setup_matter_t;
    
    // Adjust size of interpolation setups
    interp2_setups_matter_2nd.resize (Carpet::maxreflevels, NULL);
    
    
    
    assert (is_level_mode());
    if (not interp2_setups_matter_2nd.AT(Carpet::reflevel) or
        interp2_setups_matter_2nd.AT(Carpet::reflevel)->fasterp_setup->outofdate()) {
      if (interp2_setups_matter_2nd.AT(Carpet::reflevel))
        delete interp2_setups_matter_2nd.AT(Carpet::reflevel);
      interp2_setups_matter_2nd.AT(Carpet::reflevel) = new interp2_setup_matter_t;
      CCTK_VInfo (CCTK_THORNSTRING,
                  "Setting up matter interpolation for level %d", Carpet::reflevel);
    }
    interp2_setup_matter_t & interp2_setup_matter = * interp2_setups_matter_2nd.AT(Carpet::reflevel);
    
    if (not interp2_setup_matter.fasterp_setup) {
      if (verbose) CCTK_INFO ("Preparing setting up matter interpolation");
      
      assert (interp2_setup_matter.scatter_setups.empty());
      assert (interp2_setup_matter.npoints == 0);
      
      // Count points
      if (verbose) CCTK_INFO ("Counting grid points for matter interpolation");
      
      BEGIN_LOCAL_MAP_LOOP (cctkGH, CCTK_GF) {
        BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
          DECLARE_CCTK_ARGUMENTS_CHECKED(Interpolate2ApplyBC);
          
          int npoints = 0;
#pragma omp parallel reduction (+: npoints)
          LC_LOOP3 (count_points,
                    i, j, k,
                    0, 0, 0,
                    cctk_lsh[0], cctk_lsh[1], cctk_lsh[2],
                    cctk_ash[0], cctk_ash[1], cctk_ash[2])
          {
            int const ind3d = CCTK_GFINDEX3D (cctkGH, i, j, k);
            if (Sn[ind3d] >= 0) {
              ++ npoints;
            }
          } LC_ENDLOOP3 (count_points);
          
          scatter_setup_t scatter_setup;
          scatter_setup.m = Carpet::map;
          scatter_setup.c = Carpet::component;
          scatter_setup.npoints = npoints;
          
          interp2_setup_matter.scatter_setups.push_back (scatter_setup);
          interp2_setup_matter.npoints += scatter_setup.npoints;
          
        } END_LOCAL_COMPONENT_LOOP;
      } END_LOCAL_MAP_LOOP;
      
      // Allocate storage for coordinates and values
      fasterp_glocs_t locations (interp2_setup_matter.npoints);
      
      // Collect coordinates
      if (verbose) CCTK_INFO ("Collecting coordinates for matter interpolation");
      
      list<scatter_setup_t>::iterator scatter_setup_it =
        interp2_setup_matter.scatter_setups.begin();
      int pos = 0;
      BEGIN_LOCAL_MAP_LOOP (cctkGH, CCTK_GF) {
        BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
          assert (scatter_setup_it != interp2_setup_matter.scatter_setups.end());
          scatter_setup_t & scatter_setup = *scatter_setup_it;
          assert (scatter_setup.m == Carpet::map);
          assert (scatter_setup.c == Carpet::component);
          
          DECLARE_CCTK_ARGUMENTS_CHECKED(Interpolate2ApplyBC);
          
          scatter_setup.indices.reserve (scatter_setup.npoints);
          
          // This loop is not parallel
          //#pragma omp parallel reduction (+: pos)
          LC_LOOP3 (collect_coordinates,
                    i, j, k,
                    0, 0, 0,
                    cctk_lsh[0], cctk_lsh[1], cctk_lsh[2],
                    cctk_ash[0], cctk_ash[1], cctk_ash[2])
          {
            int const ind3d = CCTK_GFINDEX3D (cctkGH, i, j, k);
            if (Sn[ind3d] >= 0) {
              locations.coords[0].AT(pos) = x[ind3d];
              locations.coords[1].AT(pos) = y[ind3d];
              locations.coords[2].AT(pos) = z[ind3d];
              scatter_setup.indices.push_back (ind3d);
              ++ pos;
            }
          } LC_ENDLOOP3 (collect_coordinates);
          
          assert (int(scatter_setup.indices.size()) == scatter_setup.npoints);
          
          ++ scatter_setup_it;
        } END_LOCAL_COMPONENT_LOOP;
      } END_LOCAL_MAP_LOOP;
      assert (scatter_setup_it == interp2_setup_matter.scatter_setups.end());
      assert (pos == interp2_setup_matter.npoints);
      
      // Create the interpolation setup
      if (verbose) {
        CCTK_VInfo (CCTK_THORNSTRING,
                    "Creating matter interpolation setup for %d grid points",
                    interp2_setup_matter.npoints);
      }
      
      assert (not interp2_setup_matter.fasterp_setup);
      interp2_setup_matter.fasterp_setup =
        new fasterp_eno2_setup_t (cctkGH, locations, interpolator_order_matter);
      
    } // if not interp2_setup_matter.fasterp_setup
    
    
    
    // Poison boundaries, so that they are not accidentally used as
    // interpolation source
    if (verbose) CCTK_INFO ("Poisoning inter-patch boundaries (matter interpolator)");
    
    {
      list<scatter_setup_t>::const_iterator scatter_setup_it =
        interp2_setup_matter.scatter_setups.begin();
      BEGIN_LOCAL_MAP_LOOP (cctkGH, CCTK_GF) {
        BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
          assert (scatter_setup_it != interp2_setup_matter.scatter_setups.end());
          scatter_setup_t const & scatter_setup = *scatter_setup_it;
          assert (scatter_setup.m == Carpet::map);
          assert (scatter_setup.c == Carpet::component);
          
          DECLARE_CCTK_ARGUMENTS_CHECKED(Interpolate2ApplyBC);
          
          for (int n = 0; n < nvars_matter; ++ n) {
            
            CCTK_REAL * restrict const var =
              static_cast <CCTK_REAL *>
              (CCTK_VarDataPtrI (cctkGH, 0, indices_matter.AT(n)));
            assert (var);
            
#pragma omp parallel for
            for (int i=0; i<int(scatter_setup.indices.size()); ++i) {
              var[scatter_setup.indices.AT(i)] = poison;
            }
            
          } // for n
          
          ++ scatter_setup_it;
        } END_LOCAL_COMPONENT_LOOP;
      } END_LOCAL_MAP_LOOP;
      assert (scatter_setup_it == interp2_setup_matter.scatter_setups.end());
    }
    
    
    
    // Interpolate
    if (verbose) CCTK_INFO ("Interpolating with matter interpolator");
    
    vector<int> varinds (nvars_matter);
    for (int n=0; n<nvars_matter; ++n) {
      varinds.AT(n) = indices_matter.AT(n);
    }
    
    vector<vector<CCTK_REAL> > valuess (nvars_matter);
    vector<CCTK_REAL *> values (nvars_matter);
    for (int n=0; n<nvars_matter; ++n) {
      valuess.AT(n).resize (interp2_setup_matter.npoints);
#ifndef NDEBUG
#pragma omp parallel for
      for (int i=0; i<int(valuess.AT(n).size()); ++i) {
        valuess.AT(n).AT(i) = poison;
      }
#endif
      values.AT(n) = & valuess.AT(n).front();
    }
    
    interp2_setup_matter.fasterp_setup->interpolate (cctkGH, varinds, values);
    
    
    
    // Transform coordinates in result
    if (verbose) CCTK_INFO ("Transforming coordinates of matter variables");
    
    transform_coordinates (CCTK_PASS_CTOC, nvars_matter, indices_matter, values, interp2_setup_matter);
    
    
    
    // Write back result
    if (verbose) CCTK_INFO ("Writing back matter interpolation result");
    
    {
      list<scatter_setup_t>::const_iterator scatter_setup_it =
        interp2_setup_matter.scatter_setups.begin();
      int pos = 0;
      BEGIN_LOCAL_MAP_LOOP (cctkGH, CCTK_GF) {
        BEGIN_LOCAL_COMPONENT_LOOP (cctkGH, CCTK_GF) {
          assert (scatter_setup_it != interp2_setup_matter.scatter_setups.end());
          scatter_setup_t const & scatter_setup = *scatter_setup_it;
          assert (scatter_setup.m == Carpet::map);
          assert (scatter_setup.c == Carpet::component);
          
          DECLARE_CCTK_ARGUMENTS_CHECKED(Interpolate2ApplyBC);
          
          for (int n = 0; n < nvars_matter; ++ n) {
            
            if (verbose) {
              char * restrict const fullname = CCTK_FullName (indices_matter.AT(n));
              const int patch = MultiPatch_GetMap(cctkGH);
              CCTK_VInfo (CCTK_THORNSTRING,
                          "Filling matter variable #%d %s patch %d",
                          n, fullname, patch);
              free (fullname);
            }
            CCTK_REAL * restrict const var =
              static_cast <CCTK_REAL *>
              (CCTK_VarDataPtrI (cctkGH, 0, indices_matter.AT(n)));
            assert (var);
            
            if (interpolate_zero) {
              // for debugging only
              
#pragma omp parallel for
              for (int i=0; i<scatter_setup.npoints; ++i) {
#ifdef DEBUG
                assert (Sn[scatter_setup.indices.AT(i)] >= 0);
#endif
                var[scatter_setup.indices.AT(i)] = 0.0;
              }
              
            } else {
              // regular case
              
#pragma omp parallel for
              for (int i=0; i<scatter_setup.npoints; ++i) {
#ifdef DEBUG
                assert (Sn[scatter_setup.indices.AT(i)] >= 0);
#endif
                var[scatter_setup.indices.AT(i)] = valuess.AT(n).AT(pos + i);
              }
              
            }
            
          } // for n
          
          ++ scatter_setup_it;
          pos += scatter_setup.npoints;
        } END_LOCAL_COMPONENT_LOOP;
      } END_LOCAL_MAP_LOOP;
      assert (scatter_setup_it == interp2_setup_matter.scatter_setups.end());
      assert (pos == interp2_setup_matter.npoints);
    }
    
    } // end 2nd-order ENO
    
    else 
    {
       CCTK_WARN(0, "Matter interpolation order not supported!");
    }
    
    if (verbose) CCTK_INFO ("Done with matter interpolation.");
  }
  
  
} // namespace Interpolate2
