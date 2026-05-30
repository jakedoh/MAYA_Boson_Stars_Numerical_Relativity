
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
#include <iostream>

#include <cctk.h>
#include <cctk_Parameters.h>

#include "coordinates.hh"
#include "patchsystem.hh"


namespace Coordinates {

  using namespace std;
  
  int npatches;
  
  vector <patch_description_t> patch_descriptions;
  
  extern "C"
  int
  Coordinates_ChoosePatchSystem ()
  {
    DECLARE_CCTK_PARAMETERS

      if (CCTK_Equals(coordinate_system, "cartesian")) {

        // Undo shift in grid point locations imposed in GetDomainSpecification
        if (stagger_patch_boundaries != stagger_outer_boundaries)
          CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "Staggering of patch and outer boundaries differs. Expect CCTK_DELTA_SPACE to be wrong.");
        int const stagger_correction = !!stagger_patch_boundaries;

        npatches = 1;
        patch_descriptions.resize (npatches);
        for (int patch = 0; patch < npatches; ++ patch) {
          assert (ndirs == 3);
          patch_descriptions.at(patch).xmin[0] = patch_xmin;
          patch_descriptions.at(patch).xmin[1] = patch_ymin;
          patch_descriptions.at(patch).xmin[2] = patch_zmin;
          patch_descriptions.at(patch).xmax[0] = patch_xmax;
          patch_descriptions.at(patch).xmax[1] = patch_ymax;
          patch_descriptions.at(patch).xmax[2] = patch_zmax;
          patch_descriptions.at(patch).ncells[0] = ncells_x-stagger_correction;
          patch_descriptions.at(patch).ncells[1] = ncells_y-stagger_correction;
          patch_descriptions.at(patch).ncells[2] = ncells_z-stagger_correction;
          for (int dir = 0; dir < ndirs; ++ dir) {
            for (int face = 0; face < nfaces; ++ face) {
              patch_descriptions.at(patch).is_outer_boundary[dir][face] = true;
            }
          }
        }

      } else if (CCTK_Equals(coordinate_system, "twopatchcartesian")) {

        // Undo shift in grid point locations imposed in GetDomainSpecification
        if (stagger_patch_boundaries != stagger_outer_boundaries)
          CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "Staggering of patch and outer boundaries differs. Expect CCTK_DELTA_SPACE to be wrong.");
        int const stagger_correction = !!stagger_patch_boundaries;

        npatches = 2;
        patch_descriptions.resize (npatches);
        assert (ndirs == 3);
        patch_descriptions.at(0).xmin[0] = patch_one_xmin;
        patch_descriptions.at(0).xmin[1] = patch_one_ymin;
        patch_descriptions.at(0).xmin[2] = patch_one_zmin;
        patch_descriptions.at(0).xmax[0] = patch_one_xmax;
        patch_descriptions.at(0).xmax[1] = patch_one_ymax;
        patch_descriptions.at(0).xmax[2] = patch_one_zmax;
        patch_descriptions.at(0).ncells[0] = patch_one_ncells_x-stagger_correction;
        patch_descriptions.at(0).ncells[1] = patch_one_ncells_y-stagger_correction;
        patch_descriptions.at(0).ncells[2] = patch_one_ncells_z-stagger_correction;
        patch_descriptions.at(1).xmin[0] = patch_two_xmin;
        patch_descriptions.at(1).xmin[1] = patch_two_ymin;
        patch_descriptions.at(1).xmin[2] = patch_two_zmin;
        patch_descriptions.at(1).xmax[0] = patch_two_xmax;
        patch_descriptions.at(1).xmax[1] = patch_two_ymax;
        patch_descriptions.at(1).xmax[2] = patch_two_zmax;
        patch_descriptions.at(1).ncells[0] = patch_two_ncells_x-stagger_correction;
        patch_descriptions.at(1).ncells[1] = patch_two_ncells_y-stagger_correction;
        patch_descriptions.at(1).ncells[2] = patch_two_ncells_z-stagger_correction;
        for (int patch = 0; patch < npatches; ++ patch) {
          for (int dir = 0; dir < ndirs; ++ dir) {
            for (int face = 0; face < nfaces; ++ face) {
              patch_descriptions.at(patch).is_outer_boundary[dir][face] = true;
            }
          }
        }
        patch_descriptions.at(0).is_outer_boundary[0][1] = false;
        patch_descriptions.at(1).is_outer_boundary[0][0] = false;

      } else if (CCTK_Equals(coordinate_system, "TwoPatchDistorted")) {

        // Undo shift in grid point locations imposed in GetDomainSpecification
        if (stagger_patch_boundaries != stagger_outer_boundaries)
          CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "Staggering of patch and outer boundaries differs. Expect CCTK_DELTA_SPACE to be wrong.");
        int const stagger_correction = !!stagger_patch_boundaries;

        npatches = 2;
        patch_descriptions.resize (npatches);
        assert (ndirs == 3);

        for (int patch = 0; patch < npatches; ++ patch) {
          for (int dir = 0; dir-1 < ndirs; ++ dir) {
            for (int face = 0; face < nfaces; ++ face) {
              patch_descriptions.at(patch).is_outer_boundary[dir][face] = true;
            }
          }
        }

	patch_descriptions.at(TPD_LEFT).is_outer_boundary[0][1] = false;
	patch_descriptions.at(TPD_RIGHT).is_outer_boundary[0][1] = false;
	
        // Left patch (Cartesian)
	
	patch_descriptions.at(TPD_LEFT).xmin[0] = patch_one_xmin;
	patch_descriptions.at(TPD_LEFT).xmax[0] = patch_one_xmax;
	patch_descriptions.at(TPD_LEFT).xmin[1] = patch_one_ymin;
	patch_descriptions.at(TPD_LEFT).xmax[1] = patch_one_ymax;
	patch_descriptions.at(TPD_LEFT).xmin[2] = patch_one_zmin;
	patch_descriptions.at(TPD_LEFT).xmax[2] = patch_one_zmax;

	patch_descriptions.at(TPD_LEFT).ncells[0] = patch_one_ncells_x-stagger_correction;
	patch_descriptions.at(TPD_LEFT).ncells[1] = patch_one_ncells_y-stagger_correction;
	patch_descriptions.at(TPD_LEFT).ncells[2] = patch_one_ncells_z-stagger_correction;

        // Right patch (distorted in x direction)
	
	patch_descriptions.at(TPD_RIGHT).xmin[0] = 1;
	patch_descriptions.at(TPD_RIGHT).xmax[0] = 2;
	patch_descriptions.at(TPD_RIGHT).xmin[1] = patch_two_ymin;
	patch_descriptions.at(TPD_RIGHT).xmax[1] = patch_two_ymax;
	patch_descriptions.at(TPD_RIGHT).xmin[2] = patch_two_zmin;
	patch_descriptions.at(TPD_RIGHT).xmax[2] = patch_two_zmax;

	patch_descriptions.at(TPD_RIGHT).ncells[0] = patch_two_ncells_x-stagger_correction;
	patch_descriptions.at(TPD_RIGHT).ncells[1] = patch_two_ncells_y-stagger_correction;
	patch_descriptions.at(TPD_RIGHT).ncells[2] = patch_two_ncells_z-stagger_correction;

      } else if (CCTK_Equals(coordinate_system, "Thornburg04")) {

        if (CCTK_EQUALS(symmetry, "full")) {
           npatches = 7;
        } else if (CCTK_EQUALS(symmetry, "+z bitant")) {
           npatches = 6;
        }
        
        patch_descriptions.resize (npatches);

        CCTK_REAL cube_boundary_location = sphere_inner_radius;
	
	// check if sphere_inner_radius is an integer-multiple of h_cartesian!
	if (fabs(floor((2*cube_boundary_location) / h_cartesian + 0.5) 
	         - (2*cube_boundary_location) / h_cartesian) > 1e-8)
	   CCTK_WARN(0, "sphere_inner_radius is not a half-integer-multiple of h_cartesian!");
        // check if outer radius is outside of corner of inner cube (we repeat
        // this test again in SetGlobalCoords when we have the boundary
        // information around)
        if (sqrt(3.)*cube_boundary_location > sphere_outer_radius) {
	  CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "inner cube intersects sphere_outer_radius. You need at least sphere_outer_radius = %g + overlap_size!",
                     sqrt(3.)*cube_boundary_location);
        }
	
        CCTK_INT cube_ncells =
          lrint((2*cube_boundary_location) / h_cartesian
                - (stagger_patch_boundaries ? 1.0 : 0.0));
        CCTK_INT angular_ncells = n_angular;
      
	CCTK_REAL Rmin = stretch_rmin_1;
	CCTK_REAL Rmax = stretch_rmax_1;
      
        CCTK_REAL h0 = h_radial;
        CCTK_REAL h1 = h_radial_1;
      
        if (h1 < 0)
          h1 = h0;

        // r_physical is outer edge of last cell
        CCTK_REAL R0=sphere_inner_radius, R=R0, Rl=R0, r_physical=R0;
        CCTK_REAL Rstart = R0;
        CCTK_INT radial_ncells;
        // 1e-8 is to counteract roundoff
        for (radial_ncells=0; r_physical<sphere_outer_radius-1e-8; ++radial_ncells)
          {
            if (radial_stretch)
	      {
		Rl += h0;
#include "thornburg04_local_to_global.hh"
	      }
	    else
	      R += h0;

            r_physical = R;
          }
        // Undo shift in grid point locations imposed in GetDomainSpecification
        if (stagger_patch_boundaries != stagger_outer_boundaries)
          CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "Staggering of patch and outer boundaries differs. Expect CCTK_DELTA_SPACE to be wrong.");
        int const stagger_correction = !!stagger_patch_boundaries;
      
        // Central Cartesian patch
        
        for (int dir = 0; dir < ndirs; ++ dir) {
          patch_descriptions.at(CENTRAL_CUBE).xmin[dir] = -sphere_inner_radius;
          patch_descriptions.at(CENTRAL_CUBE).xmax[dir] = +sphere_inner_radius;
          patch_descriptions.at(CENTRAL_CUBE).ncells[dir] = cube_ncells;
          for (int face = 0; face < nfaces; ++ face) {
            patch_descriptions.at(CENTRAL_CUBE).is_outer_boundary[dir][face] = false;
            patch_descriptions.at(CENTRAL_CUBE).bndry_type[dir][face] = 0;
          }
          // The central Cartesian patch's symmetry direction coincide with the actual coordinate directions!
          patch_descriptions.at(CENTRAL_CUBE).sym_dir[dir][0] = dir;
          patch_descriptions.at(CENTRAL_CUBE).sym_dir[dir][1] = -1;  // the upper faces are not compatible with symmetries
        }
        
        if (CCTK_EQUALS(symmetry, "+z bitant")) {
           
           if ((stagger_patch_boundaries || stagger_outer_boundaries) && additional_symmetry_size) {
              CCTK_WARN(0, "If patch/outer boundaries are staggered (e.g. for cell-centered AMR), then additional_symmetry_size must be 0! Otherwise it must be 1!");
           } else if (!stagger_patch_boundaries && !stagger_outer_boundaries && !additional_symmetry_size) {
              CCTK_WARN(0, "If patch/outer boundaries are not staggered (e.g. for vertex-centered AMR), then additional_symmetry_size must be 1! Otherwise it must be 0!");
           }
           
           patch_descriptions.at(CENTRAL_CUBE).xmin[2] = 0.0;
           patch_descriptions.at(CENTRAL_CUBE).ncells[2] = lrint((cube_boundary_location) / h_cartesian
                                                         - (stagger_patch_boundaries ? 1.0 : 0.0));
           patch_descriptions.at(CENTRAL_CUBE).is_outer_boundary[2][0] = true;
           patch_descriptions.at(CENTRAL_CUBE).bndry_type[2][0] = 1;  // 1 = symmetry boundary, i.e. a boundary where symmetry conditions can be applied
        }
      
        // Angular patches
      
        for (int patch = 1; patch < npatches; ++ patch) {
          patch_descriptions.at(patch).xmin[0] = -0.25*PI;
          patch_descriptions.at(patch).xmin[1] = -0.25*PI;
          patch_descriptions.at(patch).xmin[2] = sphere_inner_radius;
          patch_descriptions.at(patch).xmax[0] = +0.25*PI;
          patch_descriptions.at(patch).xmax[1] = +0.25*PI;
	  patch_descriptions.at(patch).xmax[2] = sphere_inner_radius
	    + radial_ncells*h_radial;
	  //patch_descriptions.at(patch).xmax[2] = r_physical;
          patch_descriptions.at(patch).ncells[0] = angular_ncells;
          patch_descriptions.at(patch).ncells[1] = angular_ncells;
          patch_descriptions.at(patch).ncells[2] = radial_ncells-stagger_correction;
          for (int dir = 0; dir < ndirs-1; ++ dir) {
            for (int face = 0; face < nfaces; ++ face) {
              patch_descriptions.at(patch).is_outer_boundary[dir][face]=false;
              patch_descriptions.at(patch).bndry_type[dir][face]=0;
              patch_descriptions.at(1).sym_dir[0][0] = -1;  // set invalid symmetry direction
            }
          }
          patch_descriptions.at(patch).is_outer_boundary[ndirs-1][0] = false;
          patch_descriptions.at(patch).is_outer_boundary[ndirs-1][1] = true;
        }
        
        if (CCTK_EQUALS(symmetry, "+z bitant")) {
           
           patch_descriptions.at(1).xmin[0] = 0.0;
           patch_descriptions.at(2).xmax[0] = 0.0;
           patch_descriptions.at(3).xmin[0] = 0.0;
           patch_descriptions.at(4).xmax[0] = 0.0;
           
           patch_descriptions.at(1).ncells[0] = (angular_ncells + stagger_correction)/2 - stagger_correction;
           patch_descriptions.at(2).ncells[0] = (angular_ncells + stagger_correction)/2 - stagger_correction;
           patch_descriptions.at(3).ncells[0] = (angular_ncells + stagger_correction)/2 - stagger_correction;
           patch_descriptions.at(4).ncells[0] = (angular_ncells + stagger_correction)/2 - stagger_correction;
           
           patch_descriptions.at(1).is_outer_boundary[0][0] = true;
           patch_descriptions.at(2).is_outer_boundary[0][1] = true;
           patch_descriptions.at(3).is_outer_boundary[0][0] = true;
           patch_descriptions.at(4).is_outer_boundary[0][1] = true;
           
           // bndry_type = 1 is symmetry boundary, i.e. a boundary where symmetry conditions can be applied
           patch_descriptions.at(1).bndry_type[0][0] = 1;
           patch_descriptions.at(2).bndry_type[0][1] = 1;
           patch_descriptions.at(3).bndry_type[0][0] = 1;
           patch_descriptions.at(4).bndry_type[0][1] = 1;
           
           // set symmetry direction
           patch_descriptions.at(1).sym_dir[0][0] = 2;
           patch_descriptions.at(2).sym_dir[0][1] = 2;
           patch_descriptions.at(3).sym_dir[0][0] = 2;
           patch_descriptions.at(4).sym_dir[0][1] = 2;
           
        }

      }  else if (CCTK_Equals(coordinate_system, "Thornburg04nc")) {

        npatches = 6;
        patch_descriptions.resize (npatches);

        CCTK_INT angular_ncells = n_angular;
      
	CCTK_REAL Rmin = stretch_rmin_1;
	CCTK_REAL Rmax = stretch_rmax_1;
      
        CCTK_REAL h0 = h_radial;
        CCTK_REAL h1 = h_radial_1;
      
        if (h1 < 0)
          h1 = h0;

        CCTK_REAL R0=sphere_inner_radius, R=R0, Rl=R0, r_physical=R0;
        CCTK_REAL Rstart = R0;
        CCTK_INT radial_ncells, stagger_correction = 0;
        // 1e-8 is to counteract roundoff
        for (radial_ncells=0; r_physical<sphere_outer_radius-1e-8; ++radial_ncells)
        {
          Rl += h0;
#include "thornburg04_local_to_global.hh"
          r_physical = R;
        }
        // undo shift in grid point locations imposed in getdomainspecification
        if (stagger_patch_boundaries && stagger_outer_boundaries)
          stagger_correction = 1;
        else if (stagger_patch_boundaries != stagger_outer_boundaries)
          CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "Staggering of patch and outer boundaries differs. Expect cctk_delta_space to be wrong.");
              
        for (int patch = 0; patch < npatches; ++ patch) {
          patch_descriptions.at(patch).xmin[0] = -0.25*PI;
          patch_descriptions.at(patch).xmin[1] = -0.25*PI;
          patch_descriptions.at(patch).xmin[2] = sphere_inner_radius;
          patch_descriptions.at(patch).xmax[0] = +0.25*PI;
          patch_descriptions.at(patch).xmax[1] = +0.25*PI;
	  patch_descriptions.at(patch).xmax[2] = sphere_inner_radius
	    + radial_ncells*h_radial;
	  //patch_descriptions.at(patch).xmax[2] = r_physical;
          patch_descriptions.at(patch).ncells[0] = angular_ncells;
          patch_descriptions.at(patch).ncells[1] = angular_ncells;
          patch_descriptions.at(patch).ncells[2] = radial_ncells-stagger_correction;
          for (int dir = 0; dir < ndirs-1; ++ dir) {
            for (int face = 0; face < nfaces; ++ face) {
              patch_descriptions.at(patch).is_outer_boundary[dir][face]=false;
            }
          }
          patch_descriptions.at(patch).is_outer_boundary[ndirs-1][0] = true;
          patch_descriptions.at(patch).is_outer_boundary[ndirs-1][1] = true;
        }


      } else if (CCTK_Equals(coordinate_system, "Thornburg13")) {
        npatches = 13;
        patch_descriptions.resize (npatches);

        CCTK_REAL cube_boundary_location = sphere_inner_radius;
	
	// check if sphere_inner_radius is an integer-multiple of h_cartesian!
	if (fabs(floor((2*cube_boundary_location) / h_cartesian + 0.5) 
	         - (2*cube_boundary_location) / h_cartesian) > 1e-8)
	   CCTK_WARN(0, "sphere_inner_radius is not an integer-multiple of h_cartesian!");
	
	// check if sphere_medium_radius is an integer-multiple of h_radial_inner!
	if (fabs(floor((sphere_medium_radius-sphere_inner_radius) / h_radial_inner + 0.5)
	         - (sphere_medium_radius-sphere_inner_radius) / h_radial_inner) > 1e-8)
	   CCTK_WARN(0, "(sphere_medium_radius-sphere_inner_radius) is not an integer-multiple of h_radial_inner!");
	
        CCTK_INT cube_ncells =
          lrint((2*cube_boundary_location) / h_cartesian
                - (stagger_patch_boundaries ? 1.0 : 0.0));
        
        CCTK_INT angular_ncells_inner = n_angular_inner;
        CCTK_INT angular_ncells_outer = n_angular_outer;
        
        // set up radial grid for inner 6-patches
        CCTK_REAL h_inner = h_radial_inner;

        CCTK_REAL R0=sphere_inner_radius, R=R0, Rl=R0, r_physical=R0;
        CCTK_INT radial_ncells_inner, stagger_correction = 0;
        // 1e-8 is to counteract roundoff
        for (radial_ncells_inner=0; r_physical<sphere_medium_radius-1e-8; ++radial_ncells_inner)
          {
            R += h_inner;

            r_physical = R;
          }
        
        
        // Set up radial grid for outer 6-patches.
        CCTK_REAL Rmin = stretch_rmin_1;
	CCTK_REAL Rmax = stretch_rmax_1;
      
        CCTK_REAL h0 = h_radial_outer;
        CCTK_REAL h1 = h_radial_1;
        
        // make sure inner spherical grid ends where sphere_medium_radius starts!
        //cout << sphere_medium_radius << ", " << (sphere_inner_radius + radial_ncells_inner*h_radial_inner) << endl;
        assert(fabs(sphere_medium_radius - (sphere_inner_radius + radial_ncells_inner*h_radial_inner)) < 1e-8);
        
        R0=sphere_medium_radius, R=R0, Rl=R0, r_physical=R0;
        CCTK_REAL Rstart = R0;
        CCTK_INT radial_ncells_outer;
        // 1e-8 is to counteract roundoff
        for (radial_ncells_outer=0; r_physical<sphere_outer_radius-1e-8; ++radial_ncells_outer)
          {
            if (radial_stretch)
	      {
		Rl += h0;
#include "thornburg04_local_to_global.hh"
	      }
	    else
	      R += h0;

            r_physical = R;
            assert(r_physical >= R0);
          }
        if (stagger_patch_boundaries && stagger_outer_boundaries)
          stagger_correction = 1;
        else if (stagger_patch_boundaries != stagger_outer_boundaries)
          CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "Staggering of patch and outer boundaries differs. Expect CCTK_DELTA_SPACE to be wrong.");
      
        // Central Cartesian patch
      
        for (int dir = 0; dir < ndirs; ++ dir) {
          patch_descriptions.at(CENTRAL13_CUBE).xmin[dir] = -sphere_inner_radius;
          patch_descriptions.at(CENTRAL13_CUBE).xmax[dir] = +sphere_inner_radius;
          patch_descriptions.at(CENTRAL13_CUBE).ncells[dir] = cube_ncells;
          for (int face = 0; face < nfaces; ++ face) {
            patch_descriptions.at(CENTRAL13_CUBE).is_outer_boundary[dir][face] =
              false;
          }
        }
      
        // Angular patches
      
        for (int patch = 1; patch < npatches; ++ patch) {
          patch_descriptions.at(patch).xmin[0] = -0.25*PI;
          patch_descriptions.at(patch).xmin[1] = -0.25*PI;
          if (patch<7) 
            patch_descriptions.at(patch).xmin[2] = sphere_inner_radius;
          else 
            patch_descriptions.at(patch).xmin[2] = sphere_medium_radius;
          patch_descriptions.at(patch).xmax[0] = +0.25*PI;
          patch_descriptions.at(patch).xmax[1] = +0.25*PI;
          if (patch<7) {
	    patch_descriptions.at(patch).xmax[2] = sphere_medium_radius;
            patch_descriptions.at(patch).ncells[0] = angular_ncells_inner;
            patch_descriptions.at(patch).ncells[1] = angular_ncells_inner;
          } else {
	    patch_descriptions.at(patch).xmax[2] = sphere_medium_radius
	      + radial_ncells_outer*h_radial_outer;
            patch_descriptions.at(patch).ncells[0] = angular_ncells_outer;
            patch_descriptions.at(patch).ncells[1] = angular_ncells_outer;
          }
	  //patch_descriptions.at(patch).xmax[2] = r_physical;
          if (patch<7)
            patch_descriptions.at(patch).ncells[2] = radial_ncells_inner-stagger_correction;
          else 
            patch_descriptions.at(patch).ncells[2] = radial_ncells_outer-stagger_correction;
          for (int dir = 0; dir < ndirs-1; ++ dir) {
            for (int face = 0; face < nfaces; ++ face) {
              patch_descriptions.at(patch).is_outer_boundary[dir][face]=false;
            }
          }
          if (patch>6) {
            patch_descriptions.at(patch).is_outer_boundary[ndirs-1][0] = false;
            patch_descriptions.at(patch).is_outer_boundary[ndirs-1][1] = true;
          }
        }
        
      }  else if (CCTK_Equals(coordinate_system, "Thornburg13nc")) {

        // needs to be implemented!


      } else if (CCTK_Equals(coordinate_system, "CylinderInBox")) {
        
        npatches = 2;
        patch_descriptions.resize (npatches);
        
        CCTK_INT const box_ncells =
          floor(2 * box_radius / h_cartesian
                - (stagger_outer_boundaries ? 1.0 : 0.0)
                + 0.5);
        CCTK_INT const radial_ncells =
          floor ((transition_radius - cylinder_radius) / h_radial
                 - (stagger_outer_boundaries ? 0.5 : 0.0)
                 - (stagger_patch_boundaries ? 0.5 : 0.0)
                 + 0.5);
        CCTK_INT const angular_ncells =
          n_angular - (stagger_patch_boundaries ? 1 : 0);
        
        // Cartesian box
        for (int dir = 0; dir < ndirs; ++ dir) {
          patch_descriptions.at(CIB_BOX).xmin[dir] = -box_radius;
          patch_descriptions.at(CIB_BOX).xmax[dir] = +box_radius;
          patch_descriptions.at(CIB_BOX).ncells[dir] = box_ncells;
          for (int face = 0; face < nfaces; ++ face) {
            patch_descriptions.at(CIB_BOX).is_outer_boundary[dir][face] = true;
          }
        }
        
        // Cylinder
        patch_descriptions.at(CIB_CYLINDER).xmin[0] = cylinder_radius;
        patch_descriptions.at(CIB_CYLINDER).xmin[1] = -PI;
        patch_descriptions.at(CIB_CYLINDER).xmin[2] = -box_radius;
        patch_descriptions.at(CIB_CYLINDER).xmax[0] = transition_radius;
        patch_descriptions.at(CIB_CYLINDER).xmax[1] = +PI;
        patch_descriptions.at(CIB_CYLINDER).xmax[2] = +box_radius;
        patch_descriptions.at(CIB_CYLINDER).ncells[0] = radial_ncells;
        patch_descriptions.at(CIB_CYLINDER).ncells[1] = angular_ncells;
        patch_descriptions.at(CIB_CYLINDER).ncells[2] = box_ncells;
        for (int dir = 0; dir < ndirs; ++ dir) {
          for (int face = 0; face < nfaces; ++ face) {
            patch_descriptions.at(CIB_CYLINDER).is_outer_boundary[dir][face] =
              false;
          }
        }
        patch_descriptions.at(CIB_CYLINDER).is_outer_boundary[0][0] = true;
        patch_descriptions.at(CIB_CYLINDER).is_outer_boundary[2][0] = true;
        patch_descriptions.at(CIB_CYLINDER).is_outer_boundary[2][1] = true;
        
      } else if (CCTK_Equals(coordinate_system, "Sphere+Column")) {
        npatches = 3;
        patch_descriptions.resize (npatches);
        assert (ndirs == 3);

        CCTK_REAL Rmin = stretch_rmin_1;
        CCTK_REAL Rmax = stretch_rmax_1;

        CCTK_REAL h0 = h_radial;
        CCTK_REAL h1 = h_radial_1;

        if (h1 < 0)
          h1 = h0;

        CCTK_REAL R0=sphere_inner_radius, R=R0, Rl=R0, r_physical=R0;
        CCTK_REAL Rstart=R0;
        CCTK_INT radial_ncells;
        for (radial_ncells=0; r_physical<sphere_outer_radius; ++radial_ncells)
          {
            if (radial_stretch)
              {
                Rl += h0;
#include "thornburg04_local_to_global.hh"
              }
            else
              R += h0;

            r_physical = R;
          }
        --radial_ncells;

        for (int patch = 0; patch < npatches; ++ patch) {
          for (int dir = 0; dir-1 < ndirs; ++ dir) {
            for (int face = 0; face < nfaces; ++ face) {
              patch_descriptions.at(patch).is_outer_boundary[dir][face] = false;
            }
          }
        }

        for (int patch = 0; patch < npatches; ++ patch) {
          for (int face = 0; face < nfaces; ++ face) {
            patch_descriptions.at(patch).is_outer_boundary[2][face] = true;
          }
        }

        // Spherical patch

        patch_descriptions.at(SPC_SPHERE).xmin[0] = theta_min/180.0*PI;
        patch_descriptions.at(SPC_SPHERE).xmin[1] = -PI;
        patch_descriptions.at(SPC_SPHERE).xmin[2] = sphere_inner_radius;
        patch_descriptions.at(SPC_SPHERE).xmax[0] = PI - theta_min/180.0*PI;
        patch_descriptions.at(SPC_SPHERE).xmax[1] = PI;
        patch_descriptions.at(SPC_SPHERE).xmax[2] = sphere_inner_radius
            + radial_ncells*h_radial;
        patch_descriptions.at(SPC_SPHERE).ncells[0] = n_angular_theta;
        patch_descriptions.at(SPC_SPHERE).ncells[1] = n_angular_phi;
        patch_descriptions.at(SPC_SPHERE).ncells[2] = radial_ncells;

        // Column patches covering the poles.

        for (int patch = 1; patch < npatches; ++ patch) {
          patch_descriptions.at(patch).xmin[0] = -1.0;
          patch_descriptions.at(patch).xmin[1] = -1.0;
          patch_descriptions.at(patch).xmin[2] = sphere_inner_radius;
          patch_descriptions.at(patch).xmax[0] = 1.0;
          patch_descriptions.at(patch).xmax[1] = 1.0;
          patch_descriptions.at(patch).xmax[2] = sphere_inner_radius
              + radial_ncells*h_radial;
        patch_descriptions.at(patch).ncells[0] = n_xy;
        patch_descriptions.at(patch).ncells[1] = n_xy;
        patch_descriptions.at(patch).ncells[2] = radial_ncells;
        }
        

      } else if (CCTK_Equals(coordinate_system, "Cylinder+Column")) {
        npatches = 2;
        patch_descriptions.resize (npatches);
        assert (ndirs == 3);

        for (int patch = 0; patch < npatches; ++ patch) {
          for (int dir = 0; dir-1 < ndirs; ++ dir) {
            for (int face = 0; face < nfaces; ++ face) {
              patch_descriptions.at(patch).is_outer_boundary[dir][face] = false;
            }
          }
        }

	patch_descriptions.at(0).is_outer_boundary[0][1] = true;
	for (int face = 0; face < nfaces; ++ face)
	  patch_descriptions.at(0).is_outer_boundary[2][face] = true;
	
	for (int face = 0; face < nfaces; ++ face)
	  patch_descriptions.at(1).is_outer_boundary[2][face] = true;

        // Cylindrical patch
	
	patch_descriptions.at(CC_CYLINDER).xmin[0] = cylinder_inner_radius;
	patch_descriptions.at(CC_CYLINDER).xmax[0] = cylinder_outer_radius;
	patch_descriptions.at(CC_CYLINDER).xmin[1] = -PI;
	patch_descriptions.at(CC_CYLINDER).xmax[1] = PI;
	patch_descriptions.at(CC_CYLINDER).xmin[2] = cylinder_zmin;
	patch_descriptions.at(CC_CYLINDER).xmax[2] = cylinder_zmax;

	patch_descriptions.at(CC_CYLINDER).ncells[0] = (cylinder_outer_radius - cylinder_inner_radius) 
	  / h_radial;
	patch_descriptions.at(CC_CYLINDER).ncells[1] = n_angular_phi;
	patch_descriptions.at(CC_CYLINDER).ncells[2] = (cylinder_zmax - cylinder_zmin) / h_z;

        // Column patch covering the center.

	patch_descriptions.at(CC_COLUMN).xmin[0] = -cylinder_inner_radius;
	patch_descriptions.at(CC_COLUMN).xmax[0] = cylinder_inner_radius;
	patch_descriptions.at(CC_COLUMN).xmin[1] = -cylinder_inner_radius;
	patch_descriptions.at(CC_COLUMN).xmax[1] = cylinder_inner_radius;
	patch_descriptions.at(CC_COLUMN).xmin[2] = cylinder_zmin;
	patch_descriptions.at(CC_COLUMN).xmax[2] = cylinder_zmax;

	patch_descriptions.at(CC_COLUMN).ncells[0] = 2 * cylinder_inner_radius / h_cartesian;
	patch_descriptions.at(CC_COLUMN).ncells[1] = patch_descriptions.at(CC_COLUMN).ncells[0];
	patch_descriptions.at(CC_COLUMN).ncells[2] = patch_descriptions.at(CC_CYLINDER).ncells[2];

      } else {
      
        CCTK_WARN (CCTK_WARN_ABORT, "unknown patch system");

      }

    return 0;
  }
  
} // namespace Coordinates
