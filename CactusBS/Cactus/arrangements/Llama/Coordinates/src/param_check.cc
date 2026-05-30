
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

#include <cctk.h>
#include <cctk_Parameters.h>
#include <cctk_Arguments.h>

#ifndef DECLARE_CCTK_ARGUMENTS_CHECKED
#define DECLARE_CCTK_ARGUMENTS_CHECKED(func) DECLARE_CCTK_ARGUMENTS
#endif

namespace Coordinates
{
  extern "C"
  void
  Coordinates_ParamCheck(CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS_CHECKED(Coordinates_ParamCheck);
    DECLARE_CCTK_PARAMETERS;
    
    if (CCTK_EQUALS (coordinate_system, "Thornburg04") || CCTK_EQUALS (coordinate_system, "Thornburg04nc") || CCTK_EQUALS (coordinate_system, "Thornburg13")) {
      
      if (sphere_outer_radius == 0.)
        CCTK_PARAMWARN("Please specify a non-zero outer boundary radius.");
      
      if (sphere_inner_radius == 0.)
        CCTK_PARAMWARN("Please specify a non-zero sphere inner radius.");
      
      if (sphere_inner_radius > sphere_outer_radius)
        CCTK_VParamWarn(CCTK_THORNSTRING,
                        "sphere_inner_radius (%g) should be less than or equal to sphere_outer_radius (%g).",
                        double(sphere_inner_radius), double(sphere_outer_radius));

      if ((radial_stretch) && 
	  ((h_radial_1 < 0) || 
	   (stretch_rmin_1 >= stretch_rmax_1) || 
	   (stretch_rmin_1 < 0) ||
	   (stretch_rmax_1 < 0)))

	    CCTK_VParamWarn(CCTK_THORNSTRING,
                            "Incorrect specification of stretching region 1: h_radial_1 = %g, stretch_rmin_1 = %g, stretch_rmax_1 = %g",
                            double(h_radial_1), double(stretch_rmin_1), double(stretch_rmax_1));

      CCTK_INT const angular_ncells = n_angular;
      if (CCTK_EQUALS(symmetry, "+z bitant")) {
        if (!(stagger_patch_boundaries || stagger_outer_boundaries) && angular_ncells % 2 != 0)
          CCTK_PARAMWARN("For bitant symmetry with NON-STAGGERED boundaries, angular_ncells must be an EVEN number!");
        if ((stagger_patch_boundaries || stagger_outer_boundaries) && (angular_ncells+1) % 2 != 0)
          CCTK_PARAMWARN("For bitant symmetry with STAGGERED boundaries, angular_ncells must be an ODD number!");
      }


      if (CCTK_EQUALS(symmetry, "+z bitant") && !CCTK_IsImplementationActive("CoordinatesSymmetry")) {
         CCTK_WARN(1, "For +z bitant mode, the thorn 'CoordinatesSymmetry' is recommended to be active!");
      }

    }
    if (CCTK_EQUALS (coordinate_system, "Sphere+Column")) {
      if ((radial_stretch) && 
	  ((h_radial_1 < 0) || 
	   (stretch_rmin_1 >= stretch_rmax_1) || 
	   (stretch_rmin_1 < 0) ||
	   (stretch_rmax_1 < 0)))
	    CCTK_VParamWarn(CCTK_THORNSTRING,
                            "Incorrect specification of stretching region 1: h_radial_1 = %g, stretch_rmin_1 = %g, stretch_rmax_1 = %g",
                            double(h_radial_1), double(stretch_rmin_1), double(stretch_rmax_1));
	
    }
    if (CCTK_EQUALS (coordinate_system, "Thornburg13")) {
       
       if (sphere_medium_radius > sphere_outer_radius || sphere_medium_radius < sphere_inner_radius)
          CCTK_VParamWarn(CCTK_THORNSTRING,
                          "sphere medium radius (%g) must be between sphere inner (%g) and sphere outer radius (%g)!",
                          double(sphere_medium_radius),
                          double(sphere_inner_radius),
                          double(sphere_outer_radius));
       
       if ((radial_stretch) && 
           stretch_rmin_1 < sphere_medium_radius)
           CCTK_VParamWarn(CCTK_THORNSTRING,
                           "Stretching region (%g) must be outside of sphere medium radius (%g).",
                           double(stretch_rmin_1), double(sphere_medium_radius));
    }
  }
}
