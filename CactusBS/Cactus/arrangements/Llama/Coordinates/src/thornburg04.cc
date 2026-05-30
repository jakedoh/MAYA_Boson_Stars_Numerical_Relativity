
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

/*
 * Sets up a 7-patch multigrid, with coordinates based on those described
 * in Thornburg (2004a,b).
 *
 * Assumes a central cubical grid, surrounded by 6 patches. The central
 * grid overlaps the spherical patches, and is set up with usual
 * cartesian coordinates. The spherical patches use coordinates:
 *
 * +/-X : (a,b,c) == (nu,ph,rho)
 *   nu = rotation about local y = arctan(z/x)
 *   ph = rotation about local z = arctan(y/x)
 *   rho = radius = +/-x*sqrt(1 + tan^2(ph) + tan^2(nu))
 *
 * +/-Y : (a,b,c) == (mu,ph,rho)
 *   mu = rotation about local x = arctan(z/y)
 *   ph = rotation about local z = arctan(x/y)
 *   rho = radius = +/-y*sqrt(1 + tan^2(ph) + tan^2(nu))
 *
 * +/-Z : (a,b,c) == (mu,nu,rho)
 *   mu = rotation about local x = arctan(y/z)
 *   nu = rotation about local y = arctan(x/z)
 *   rho = radius = +/-z*sqrt(1 + tan^2(mu) + tan^2(nu))
 *
 * The local angular coordinates have a range [-pi/4; +pi/4], the
 * local radial coordinate ranges from rmin to rmax, same as the
 * global r coordinate.  The local Cartesian coordinates are also
 * identical to the global Cartesian coordinates.
 *
 */
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <time.h>

#include <cctk.h>
#include <cctk_Parameters.h>
#include <cctk_Arguments.h>

#include "coordinates.hh"
#include "patchsystem.hh"

#define EPSILON 1e-10

#define REFLEVEL (GetRefinementLevel(cctkGH))

#ifndef DECLARE_CCTK_ARGUMENTS_CHECKED
#define DECLARE_CCTK_ARGUMENTS_CHECKED(func) DECLARE_CCTK_ARGUMENTS
#endif


namespace Coordinates
{
  
  using namespace std;
  
  // Determine the extent of the intermediate coordinate system.  This
  // is the local coordinate system of the coarsest grid.
  static
  void
  get_intermediate_coordinates (CCTK_INT const map,
                                CCTK_INT const size,
                                CCTK_REAL * const exterior_min,
                                CCTK_REAL * const exterior_max,
                                CCTK_REAL * const h)
  {
    /*CCTK_REAL physical_min[3];
    CCTK_REAL physical_max[3];
    CCTK_INT ncells[3];
    for (int d=0; d<3; ++d) {
      physical_min[d] = patch_descriptions.at(map).xmin[d];
      physical_max[d] = patch_descriptions.at(map).xmax[d];
      ncells[d] = patch_descriptions.at(map).ncells[d];
      h[d] = (physical_max[d] - physical_min[d]) / ncells[d];
    }
    
    // Calculate the exterior boundaries
    CCTK_REAL interior_min[3];
    CCTK_REAL interior_max[3];
    CCTK_INT const ierr = MultiPatch_ConvertFromPhysicalBoundary
      (map, 3,
       physical_min, physical_max,
       interior_min, interior_max,
       exterior_min, exterior_max,
       h);
    assert (not ierr);*/

    CCTK_REAL physical_min[3];
    CCTK_REAL physical_max[3];
    CCTK_REAL interior_min[3];
    CCTK_REAL interior_max[3];
    CCTK_INT const ierr = Coordinates_GetDomainSpecification (map,
                                               size,
                                               physical_min,
                                               physical_max,
                                               interior_min,
                                               interior_max,
                                               exterior_min,
                                               exterior_max,
                                               h);

    assert (not ierr);

  }

  static
  void
  get_local_coordinates (cGH const * const cctkGH,
                         CCTK_REAL const * const x0,
                         CCTK_REAL const * const dx,
                         int       const * const ipos,
                         CCTK_REAL       * const a)
  {
    CCTK_REAL a0[3], da[3];
    for (int d=0; d<3; ++d) {
      // grid spacing on this refinement level
      da[d] = dx[d] / cctkGH->cctk_levfac[d];
      // offset of this refinement level
      a0[d] =
        x0[d] + da[d] * cctkGH->cctk_levoff[d] / cctkGH->cctk_levoffdenom[d];
      // coordinate of this grid point
      a[d] = a0[d] + da[d] * (cctkGH->cctk_lbnd[d] + ipos[d]);
    }
  }



  extern "C"
  void
  Coordinates_SetGlobalCoords_Thornburg04(CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS_CHECKED(Coordinates_SetGlobalCoords_Thornburg04);
    DECLARE_CCTK_PARAMETERS;
    
    CCTK_INT const map = MultiPatch_GetMap (cctkGH);
    assert (map >= 0);
    
    if (verbose) {
      CCTK_INFO("Setting up global coordinates for thornburg04.");
    }

    // Determine the extent of the intermediate (see comment above)
    // coordinate system
    CCTK_REAL exterior_min[3];
    CCTK_REAL exterior_max[3];
    CCTK_REAL h[3];
    get_intermediate_coordinates
      (map, 3,
       exterior_min, exterior_max,
       h);

    const CCTK_REAL Rmin = stretch_rmin_1;
    const CCTK_REAL Rmax = stretch_rmax_1;

    const CCTK_REAL h0 = h_radial;
    const CCTK_REAL h1 = ( h_radial_1<0 ? h0 : h_radial_1 );
    const CCTK_REAL Rstart = sphere_inner_radius;
    
    const CCTK_INT do_radial_stretch = radial_stretch \
      && (abs(h1 - h0) > EPSILON);

    switch (map) {
      
    case CENTRAL_CUBE: {
      if (sqrt (pow(exterior_max[0], 2) + pow(exterior_max[1], 2) + pow(exterior_max[2], 2)) >= sphere_outer_radius)
        {
          CCTK_VWarn 
            (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
              "edge of inner cube intersects sphere_outer_radius. Minimum sphere_outer_radius is %g!",
              sqrt (pow(exterior_max[0], 2) + pow(exterior_max[1], 2) + pow(exterior_max[2], 2)));
        }

      for (int k=0; k<cctk_lsh[2]; ++k)
        for (int j=0; j<cctk_lsh[1]; ++j)
          for (int i=0; i<cctk_lsh[0]; ++i)
            {
              const int ipos[3] = { i, j, k };
              const int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
              CCTK_REAL a[3];
              get_local_coordinates (cctkGH, exterior_min, h, ipos, a);
            
              x[ijk] = a[0];
              y[ijk] = a[1];
              z[ijk] = a[2];
              r[ijk] = sqrt (pow(x[ijk],2) + pow(y[ijk],2) + pow(z[ijk],2));
            }
      break;
    }
      
      /*
       * sphere +x : (a,b,c) == (nu,ph,rho)
       *   nu = rotation about local y = arctan(z/x)
       *   ph = rotation about local z = arctan(y/x)
       *   rho = radius = x*sqrt(1 + tan^2(ph) + tan^2(nu))
       */
    case SPHERE_PLUS_X: {
      for (int k=0; k<cctk_lsh[2]; ++k)
        for (int j=0; j<cctk_lsh[1]; ++j)
          for (int i=0; i<cctk_lsh[0]; ++i)
            {
              const int ipos[3] = { i, j, k };
              const int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
              CCTK_REAL a[3];
              get_local_coordinates (cctkGH, exterior_min, h, ipos, a);
            
              CCTK_REAL const nu  = a[0];
              CCTK_REAL const ph  = a[1];
              CCTK_REAL const rho = a[2];
            
              CCTK_REAL const tan_nu = tan(nu);
              CCTK_REAL const tan_ph = tan(ph);

              CCTK_REAL R;
	      if (do_radial_stretch)
		{
		  CCTK_REAL Rl=rho;
#include "include/thornburg04_local_to_global.hh"
		}
	      else
		R = rho;

              x[ijk] = R/sqrt(1.0 + tan_nu*tan_nu + tan_ph*tan_ph);
              y[ijk] = x[ijk] * tan_ph;
              z[ijk] = x[ijk] * tan_nu;
              r[ijk] = R;
            }
      break;
    }
      
      /*
       * sphere -x : (a,b,c) == (nu,ph,rho)
       *   nu = rotation about y = arctan(z/x)
       *   ph = rotation about z = arctan(y/x)
       *   rho = radius = -x*sqrt(1 + tan^2(ph) + tan^2(nu))
       */
    case SPHERE_MINUS_X: {
      for (int k=0; k<cctk_lsh[2]; ++k)
        for (int j=0; j<cctk_lsh[1]; ++j)
          for (int i=0; i<cctk_lsh[0]; ++i)
            {
              const int ipos[3] = { i, j, k };
              const int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
              CCTK_REAL a[3];
              get_local_coordinates (cctkGH, exterior_min, h, ipos, a);
            
              CCTK_REAL const nu  = a[0];
              CCTK_REAL const ph  = a[1];
              CCTK_REAL const rho = a[2];
            
              CCTK_REAL const tan_nu = tan(nu);
              CCTK_REAL const tan_ph = tan(ph);

              CCTK_REAL R;
	      if (do_radial_stretch)
		{
		  CCTK_REAL Rl=rho;
#include "include/thornburg04_local_to_global.hh"
		}
	      else
		R = rho;

              x[ijk] = -R/sqrt(1.0 + tan_nu*tan_nu + tan_ph*tan_ph);
              y[ijk] = x[ijk] * tan_ph;
              z[ijk] = x[ijk] * tan_nu;
              r[ijk] = R;
            }
      break;
    }
      
      /*
       * sphere +y : (a,b,c) == (mu,ph,rho)
       *   mu = rotation about local x = arctan(z/y)
       *   ph = rotation about local z = arctan(x/y)
       *   rho = radius = y*sqrt(1 + tan^2(ph) + tan^2(nu))
       */
    case SPHERE_PLUS_Y: {
      for (int k=0; k<cctk_lsh[2]; ++k)
        for (int j=0; j<cctk_lsh[1]; ++j)
          for (int i=0; i<cctk_lsh[0]; ++i)
            {
              const int ipos[3] = { i, j, k };
              const int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
              CCTK_REAL a[3];
              get_local_coordinates (cctkGH, exterior_min, h, ipos, a);
            
              CCTK_REAL const mu  = a[0];
              CCTK_REAL const ph  = a[1];
              CCTK_REAL const rho = a[2];
            
              CCTK_REAL const tan_mu = tan(mu);
              CCTK_REAL const tan_ph = tan(ph);

              CCTK_REAL R;
	      if (do_radial_stretch)
		{
		  CCTK_REAL Rl=rho;
#include "include/thornburg04_local_to_global.hh"
		}
	      else
		R = rho;

              y[ijk] = R/sqrt(1.0 + tan_mu*tan_mu + tan_ph*tan_ph);
              x[ijk] = y[ijk] * tan_ph;
              z[ijk] = y[ijk] * tan_mu;
              r[ijk] = R;
            }
      break;
    }
      
      /*
       * sphere +y : (a,b,c) == (mu,ph,rho)
       *   mu = rotation about local x = arctan(z/y)
       *   ph = rotation about local z = arctan(x/y)
       *   rho = radius = y*sqrt(1 + tan^2(ph) + tan^2(nu))
       */
    case SPHERE_MINUS_Y: {
      for (int k=0; k<cctk_lsh[2]; ++k)
        for (int j=0; j<cctk_lsh[1]; ++j)
          for (int i=0; i<cctk_lsh[0]; ++i)
            {
              const int ipos[3] = { i, j, k };
              const int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
              CCTK_REAL a[3];
              get_local_coordinates (cctkGH, exterior_min, h, ipos, a);
            
              CCTK_REAL const mu  = a[0];
              CCTK_REAL const ph  = a[1];
              CCTK_REAL const rho = a[2];
            
              CCTK_REAL const tan_mu = tan(mu);
              CCTK_REAL const tan_ph = tan(ph);

              CCTK_REAL R;
	      if (do_radial_stretch)
		{
		  CCTK_REAL Rl=rho;
#include "include/thornburg04_local_to_global.hh"
		}
	      else
		R = rho;
          
              y[ijk] = -R/sqrt(1.0 + tan_mu*tan_mu + tan_ph*tan_ph);
              x[ijk] = y[ijk] * tan_ph;
              z[ijk] = y[ijk] * tan_mu;
              r[ijk] = R;
            }
      break;
    }
      
      /*
       * sphere +z : (a,b,c) == (mu,nu,rho)
       *   mu = rotation about local x = arctan(y/z)
       *   nu = rotation about local y = arctan(x/z)
       *   rho = radius = z*sqrt(1 + tan^2(mu) + tan^2(nu))
       */
    case SPHERE_PLUS_Z: {
      for (int k=0; k<cctk_lsh[2]; ++k)
        for (int j=0; j<cctk_lsh[1]; ++j)
          for (int i=0; i<cctk_lsh[0]; ++i)
            {
              const int ipos[3] = { i, j, k };
              const int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
              CCTK_REAL a[3];
              get_local_coordinates (cctkGH, exterior_min, h, ipos, a);
            
              CCTK_REAL const mu  = a[0];
              CCTK_REAL const nu  = a[1];
              CCTK_REAL const rho = a[2];
            
              CCTK_REAL const tan_mu = tan(mu);
              CCTK_REAL const tan_nu = tan(nu);

              CCTK_REAL R;
	      if (do_radial_stretch)
		{
		  CCTK_REAL Rl=rho;
#include "include/thornburg04_local_to_global.hh"
		}
	      else
		R = rho;

              z[ijk] = R/sqrt(1.0 + tan_mu*tan_mu + tan_nu*tan_nu);
              x[ijk] = z[ijk] * tan_nu;
              y[ijk] = z[ijk] * tan_mu;
              r[ijk] = R;
            }
      break;
    }
      
      /*
       * sphere +z : (a,b,c) == (mu,nu,rho)
       *   mu = rotation about local x = arctan(y/z)
       *   nu = rotation about local y = arctan(x/z)
       *   rho = radius = z*sqrt(1 + tan^2(mu) + tan^2(nu))
       */
    case SPHERE_MINUS_Z: {
      for (int k=0; k<cctk_lsh[2]; ++k)
        for (int j=0; j<cctk_lsh[1]; ++j)
          for (int i=0; i<cctk_lsh[0]; ++i)
            {
              const int ipos[3] = { i, j, k };
              const int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
              CCTK_REAL a[3];
              get_local_coordinates (cctkGH, exterior_min, h, ipos, a);
            
              CCTK_REAL const mu  = a[0];
              CCTK_REAL const nu  = a[1];
              CCTK_REAL const rho = a[2];
            
              CCTK_REAL const tan_mu = tan(mu);
              CCTK_REAL const tan_nu = tan(nu);

              CCTK_REAL R;
	      if (do_radial_stretch)
		{
		  CCTK_REAL Rl=rho;
#include "include/thornburg04_local_to_global.hh"
		}
	      else
		R = rho;

              z[ijk] = -R/sqrt(1.0 + tan_mu*tan_mu + tan_nu*tan_nu);
              x[ijk] = z[ijk] * tan_nu;
              y[ijk] = z[ijk] * tan_mu;
              r[ijk] = R;
            }
      break;
    }
      
    default:
      CCTK_WARN(0, "Patch number must be 0..6");
      
    }
  }

  extern "C"
  void
  Coordinates_SetJacobian_Thornburg04(CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS_CHECKED(Coordinates_SetJacobian_Thornburg04);
    DECLARE_CCTK_PARAMETERS;

    using namespace std;

    int i, j, k;

    if (verbose) {
      CCTK_INFO("Setting local coordinates: Thornburg04");
    }

    *general_coordinates = 1;
    *interpolate_boundary_points = CCTK_EQUALS(symmetry,"full") ? 1 : 0;
    *jacobian_state = 1;
    *jacobian_derivative_state = 1;

    const CCTK_REAL Rmin = stretch_rmin_1;
    const CCTK_REAL Rmax = stretch_rmax_1;

    const CCTK_REAL h0 = h_radial;
    const CCTK_REAL h1 = ( h_radial_1<0 ? h0 : h_radial_1 );
    const CCTK_REAL Rstart = sphere_inner_radius;

    const CCTK_INT do_radial_stretch = radial_stretch \
      && (abs(h1 - h0) > EPSILON);

    const CCTK_INT patch = Coordinates_GetMap(cctkGH);

    switch (patch)
      {
      case CENTRAL_CUBE:
        {
          for (k=0; k<cctk_lsh[2]; ++k)
            for (j=0; j<cctk_lsh[1]; ++j)
              for (i=0; i<cctk_lsh[0]; ++i)
                {
                  const int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
  
                  J11[ijk] = 1.0;
                  J12[ijk] = 0.0;
                  J13[ijk] = 0.0;
                  J21[ijk] = 0.0;
                  J22[ijk] = 1.0;
                  J23[ijk] = 0.0;
                  J31[ijk] = 0.0;
                  J32[ijk] = 0.0;
                  J33[ijk] = 1.0;


                  dJ111[ijk] = 0.0;
                  dJ112[ijk] = 0.0;
                  dJ113[ijk] = 0.0;
                  dJ122[ijk] = 0.0;
                  dJ123[ijk] = 0.0;
                  dJ133[ijk] = 0.0;
          
                  dJ211[ijk] = 0.0;
                  dJ212[ijk] = 0.0;
                  dJ213[ijk] = 0.0;
                  dJ222[ijk] = 0.0;
                  dJ223[ijk] = 0.0;
                  dJ233[ijk] = 0.0;

                  dJ311[ijk] = 0.0;
                  dJ312[ijk] = 0.0;
                  dJ313[ijk] = 0.0;
                  dJ322[ijk] = 0.0;
                  dJ323[ijk] = 0.0;
                  dJ333[ijk] = 0.0;
                }
        }
        break;

      case SPHERE_PLUS_X:
        {
          for (k=0; k<cctk_lsh[2]; ++k)
            for (j=0; j<cctk_lsh[1]; ++j)
              for (i=0; i<cctk_lsh[0]; ++i)
                {
                  const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
                  CCTK_REAL xp = x[ijk];
                  CCTK_REAL yp = y[ijk];
                  CCTK_REAL zp = z[ijk];
                  CCTK_REAL J[3][3];
                  CCTK_REAL dJ[3][3][3];
		  
		  if (do_radial_stretch)
		    {
#include "include/thornburg04_x_plus_J.hh"
#include "include/thornburg04_x_plus_dJ.hh"
		    }
		  else
		    {
#include "include/thornburg04_nostretch_x_plus_J.hh"
#include "include/thornburg04_nostretch_x_plus_dJ.hh"
		    }

                  J11[ijk] = J[0][0];
                  J12[ijk] = J[0][1];
                  J13[ijk] = J[0][2];
                  J21[ijk] = J[1][0];
                  J22[ijk] = J[1][1];
                  J23[ijk] = J[1][2];
                  J31[ijk] = J[2][0];
                  J32[ijk] = J[2][1];
                  J33[ijk] = J[2][2];
		  
                  dJ111[ijk] = dJ[0][0][0];
                  dJ112[ijk] = dJ[0][0][1];
                  dJ113[ijk] = dJ[0][0][2];
                  dJ122[ijk] = dJ[0][1][1];
                  dJ123[ijk] = dJ[0][1][2];
                  dJ133[ijk] = dJ[0][2][2];

                  dJ211[ijk] = dJ[1][0][0];
                  dJ212[ijk] = dJ[1][0][1];
                  dJ213[ijk] = dJ[1][0][2];
                  dJ222[ijk] = dJ[1][1][1];
                  dJ223[ijk] = dJ[1][1][2];
                  dJ233[ijk] = dJ[1][2][2];

                  dJ311[ijk] = dJ[2][0][0];
                  dJ312[ijk] = dJ[2][0][1];
                  dJ313[ijk] = dJ[2][0][2];
                  dJ322[ijk] = dJ[2][1][1];
                  dJ323[ijk] = dJ[2][1][2];
                  dJ333[ijk] = dJ[2][2][2];
                }
        }
        break;

      case SPHERE_MINUS_X:
        {
          for (k=0; k<cctk_lsh[2]; ++k)
            for (j=0; j<cctk_lsh[1]; ++j)
              for (i=0; i<cctk_lsh[0]; ++i)
                {
                  const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
                  CCTK_REAL xp = x[ijk];
                  CCTK_REAL yp = y[ijk];
                  CCTK_REAL zp = z[ijk];
                  CCTK_REAL J[3][3];
                  CCTK_REAL dJ[3][3][3];

		  if (do_radial_stretch)
		    {
#include "include/thornburg04_x_plus_J.hh"
#include "include/thornburg04_x_plus_dJ.hh"
		    }
		  else
		    {
#include "include/thornburg04_nostretch_x_plus_J.hh"
#include "include/thornburg04_nostretch_x_plus_dJ.hh"
		    }

                  J11[ijk] = J[0][0];
                  J12[ijk] = J[0][1];
                  J13[ijk] = J[0][2];
                  J21[ijk] = J[1][0];
                  J22[ijk] = J[1][1];
                  J23[ijk] = J[1][2];
                  J31[ijk] = J[2][0];
                  J32[ijk] = J[2][1];
                  J33[ijk] = J[2][2];
		  
                  dJ111[ijk] = dJ[0][0][0];
                  dJ112[ijk] = dJ[0][0][1];
                  dJ113[ijk] = dJ[0][0][2];
                  dJ122[ijk] = dJ[0][1][1];
                  dJ123[ijk] = dJ[0][1][2];
                  dJ133[ijk] = dJ[0][2][2];

                  dJ211[ijk] = dJ[1][0][0];
                  dJ212[ijk] = dJ[1][0][1];
                  dJ213[ijk] = dJ[1][0][2];
                  dJ222[ijk] = dJ[1][1][1];
                  dJ223[ijk] = dJ[1][1][2];
                  dJ233[ijk] = dJ[1][2][2];

                  dJ311[ijk] = dJ[2][0][0];
                  dJ312[ijk] = dJ[2][0][1];
                  dJ313[ijk] = dJ[2][0][2];
                  dJ322[ijk] = dJ[2][1][1];
                  dJ323[ijk] = dJ[2][1][2];
                  dJ333[ijk] = dJ[2][2][2];
                }
        }
        break;

      case SPHERE_PLUS_Y:
        {
          for (k=0; k<cctk_lsh[2]; ++k)
            for (j=0; j<cctk_lsh[1]; ++j)
              for (i=0; i<cctk_lsh[0]; ++i)
                {
                  const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
                  CCTK_REAL xp = x[ijk];
                  CCTK_REAL yp = y[ijk];
                  CCTK_REAL zp = z[ijk];
                  CCTK_REAL J[3][3];
                  CCTK_REAL dJ[3][3][3];

		  if (do_radial_stretch)
		    {
#include "include/thornburg04_y_plus_J.hh"
#include "include/thornburg04_y_plus_dJ.hh"
		    }
		  else
		    {
#include "include/thornburg04_nostretch_y_plus_J.hh"
#include "include/thornburg04_nostretch_y_plus_dJ.hh"
		    }

                  J11[ijk] = J[0][0];
                  J12[ijk] = J[0][1];
                  J13[ijk] = J[0][2];
                  J21[ijk] = J[1][0];
                  J22[ijk] = J[1][1];
                  J23[ijk] = J[1][2];
                  J31[ijk] = J[2][0];
                  J32[ijk] = J[2][1];
                  J33[ijk] = J[2][2];
		  
                  dJ111[ijk] = dJ[0][0][0];
                  dJ112[ijk] = dJ[0][0][1];
                  dJ113[ijk] = dJ[0][0][2];
                  dJ122[ijk] = dJ[0][1][1];
                  dJ123[ijk] = dJ[0][1][2];
                  dJ133[ijk] = dJ[0][2][2];

                  dJ211[ijk] = dJ[1][0][0];
                  dJ212[ijk] = dJ[1][0][1];
                  dJ213[ijk] = dJ[1][0][2];
                  dJ222[ijk] = dJ[1][1][1];
                  dJ223[ijk] = dJ[1][1][2];
                  dJ233[ijk] = dJ[1][2][2];

                  dJ311[ijk] = dJ[2][0][0];
                  dJ312[ijk] = dJ[2][0][1];
                  dJ313[ijk] = dJ[2][0][2];
                  dJ322[ijk] = dJ[2][1][1];
                  dJ323[ijk] = dJ[2][1][2];
                  dJ333[ijk] = dJ[2][2][2];
		  
                }
        }
        break;

      case SPHERE_MINUS_Y:
        {
          for (k=0; k<cctk_lsh[2]; ++k)
            for (j=0; j<cctk_lsh[1]; ++j)
              for (i=0; i<cctk_lsh[0]; ++i)
                {
                  const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
                  CCTK_REAL xp = x[ijk];
                  CCTK_REAL yp = y[ijk];
                  CCTK_REAL zp = z[ijk];
                  CCTK_REAL J[3][3];
                  CCTK_REAL dJ[3][3][3];

		  if (do_radial_stretch)
		    {
#include "include/thornburg04_y_plus_J.hh"
#include "include/thornburg04_y_plus_dJ.hh"
		    }
		  else
		    {
#include "include/thornburg04_nostretch_y_plus_J.hh"
#include "include/thornburg04_nostretch_y_plus_dJ.hh"
		    }

                  J11[ijk] = J[0][0];
                  J12[ijk] = J[0][1];
                  J13[ijk] = J[0][2];
                  J21[ijk] = J[1][0];
                  J22[ijk] = J[1][1];
                  J23[ijk] = J[1][2];
                  J31[ijk] = J[2][0];
                  J32[ijk] = J[2][1];
                  J33[ijk] = J[2][2];
		  
                  dJ111[ijk] = dJ[0][0][0];
                  dJ112[ijk] = dJ[0][0][1];
                  dJ113[ijk] = dJ[0][0][2];
                  dJ122[ijk] = dJ[0][1][1];
                  dJ123[ijk] = dJ[0][1][2];
                  dJ133[ijk] = dJ[0][2][2];

                  dJ211[ijk] = dJ[1][0][0];
                  dJ212[ijk] = dJ[1][0][1];
                  dJ213[ijk] = dJ[1][0][2];
                  dJ222[ijk] = dJ[1][1][1];
                  dJ223[ijk] = dJ[1][1][2];
                  dJ233[ijk] = dJ[1][2][2];

                  dJ311[ijk] = dJ[2][0][0];
                  dJ312[ijk] = dJ[2][0][1];
                  dJ313[ijk] = dJ[2][0][2];
                  dJ322[ijk] = dJ[2][1][1];
                  dJ323[ijk] = dJ[2][1][2];
                  dJ333[ijk] = dJ[2][2][2];
                }
        }
        break;

      case SPHERE_PLUS_Z:
        {
          for (k=0; k<cctk_lsh[2]; ++k)
            for (j=0; j<cctk_lsh[1]; ++j)
              for (i=0; i<cctk_lsh[0]; ++i)
                {
                  const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
                  CCTK_REAL xp = x[ijk];
                  CCTK_REAL yp = y[ijk];
                  CCTK_REAL zp = z[ijk];
                  CCTK_REAL J[3][3];
                  CCTK_REAL dJ[3][3][3];

		  if (do_radial_stretch)
		    {
#include "include/thornburg04_z_plus_J.hh"
#include "include/thornburg04_z_plus_dJ.hh"
		    }
		  else
		    {
#include "include/thornburg04_nostretch_z_plus_J.hh"
#include "include/thornburg04_nostretch_z_plus_dJ.hh"
		    }

                  J11[ijk] = J[0][0];
                  J12[ijk] = J[0][1];
                  J13[ijk] = J[0][2];
                  J21[ijk] = J[1][0];
                  J22[ijk] = J[1][1];
                  J23[ijk] = J[1][2];
                  J31[ijk] = J[2][0];
                  J32[ijk] = J[2][1];
                  J33[ijk] = J[2][2];
		  
                  dJ111[ijk] = dJ[0][0][0];
                  dJ112[ijk] = dJ[0][0][1];
                  dJ113[ijk] = dJ[0][0][2];
                  dJ122[ijk] = dJ[0][1][1];
                  dJ123[ijk] = dJ[0][1][2];
                  dJ133[ijk] = dJ[0][2][2];

                  dJ211[ijk] = dJ[1][0][0];
                  dJ212[ijk] = dJ[1][0][1];
                  dJ213[ijk] = dJ[1][0][2];
                  dJ222[ijk] = dJ[1][1][1];
                  dJ223[ijk] = dJ[1][1][2];
                  dJ233[ijk] = dJ[1][2][2];

                  dJ311[ijk] = dJ[2][0][0];
                  dJ312[ijk] = dJ[2][0][1];
                  dJ313[ijk] = dJ[2][0][2];
                  dJ322[ijk] = dJ[2][1][1];
                  dJ323[ijk] = dJ[2][1][2];
                  dJ333[ijk] = dJ[2][2][2];
                }
        }
        break;

      case SPHERE_MINUS_Z:
        {
          for (k=0; k<cctk_lsh[2]; ++k)
            for (j=0; j<cctk_lsh[1]; ++j)
              for (i=0; i<cctk_lsh[0]; ++i)
                {
                  const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
                  CCTK_REAL xp = x[ijk];
                  CCTK_REAL yp = y[ijk];
                  CCTK_REAL zp = z[ijk];
                  CCTK_REAL J[3][3];
                  CCTK_REAL dJ[3][3][3];

		  if (do_radial_stretch)
		    {
#include "include/thornburg04_z_plus_J.hh"
#include "include/thornburg04_z_plus_dJ.hh"
		    }
		  else
		    {
#include "include/thornburg04_nostretch_z_plus_J.hh"
#include "include/thornburg04_nostretch_z_plus_dJ.hh"
		    }

                  J11[ijk] = J[0][0];
                  J12[ijk] = J[0][1];
                  J13[ijk] = J[0][2];
                  J21[ijk] = J[1][0];
                  J22[ijk] = J[1][1];
                  J23[ijk] = J[1][2];
                  J31[ijk] = J[2][0];
                  J32[ijk] = J[2][1];
                  J33[ijk] = J[2][2];
		  
                  dJ111[ijk] = dJ[0][0][0];
                  dJ112[ijk] = dJ[0][0][1];
                  dJ113[ijk] = dJ[0][0][2];
                  dJ122[ijk] = dJ[0][1][1];
                  dJ123[ijk] = dJ[0][1][2];
                  dJ133[ijk] = dJ[0][2][2];

                  dJ211[ijk] = dJ[1][0][0];
                  dJ212[ijk] = dJ[1][0][1];
                  dJ213[ijk] = dJ[1][0][2];
                  dJ222[ijk] = dJ[1][1][1];
                  dJ223[ijk] = dJ[1][1][2];
                  dJ233[ijk] = dJ[1][2][2];

                  dJ311[ijk] = dJ[2][0][0];
                  dJ312[ijk] = dJ[2][0][1];
                  dJ313[ijk] = dJ[2][0][2];
                  dJ322[ijk] = dJ[2][1][1];
                  dJ323[ijk] = dJ[2][1][2];
                  dJ333[ijk] = dJ[2][2][2];
                }
        }
        break;

      default:
        CCTK_WARN(0, "Patch number must be 0..6");

      }
    
    return;
  }
  
  
  
  extern "C"
  CCTK_INT
  global_to_local_Thornburg04
  (CCTK_REAL const x, CCTK_REAL const y, CCTK_REAL const z,
   CCTK_REAL localcoord[ndirs])
  {
    DECLARE_CCTK_PARAMETERS;

    const CCTK_REAL rp2 = x*x + y*y + z*z;
    
    if (rp2 < pow (sphere_inner_radius, 2))
      {
        localcoord[0] = x;
        localcoord[1] = y;
        localcoord[2] = z;
	
        return CENTRAL_CUBE;
      }
    else
      {
        const CCTK_REAL Rmin = stretch_rmin_1;
        const CCTK_REAL Rmax = stretch_rmax_1;

        const CCTK_REAL h0 = h_radial;
	const CCTK_REAL h1 = ( h_radial_1<0 ? h0 : h_radial_1 );
        const CCTK_REAL Rstart = sphere_inner_radius;
        
	const CCTK_INT do_radial_stretch = radial_stretch 
	  && (abs(h1 - h0) > EPSILON);


        const CCTK_REAL ax = fabs(x);
        const CCTK_REAL ay = fabs(y);
        const CCTK_REAL az = fabs(z);

        const CCTK_REAL Rg = sqrt(ax*ax + ay*ay + az*az);

        CCTK_REAL R;
	if (do_radial_stretch)
	  {
#include "include/thornburg04_global_to_local.hh"
	  }
	else
	  R = Rg;

        localcoord[2] = R;

        /*   +/-x spheres   */
        if ((ax >= ay) && (ax >= az))
          {
#ifdef USE_ATAN2
            localcoord[0] = atan2(x >= 0 ? z : -z, ax);
            localcoord[1] = atan2(x >= 0 ? y : -y, ax);
#else
            assert(x);
            localcoord[0] = atan(z/x);
            localcoord[1] = atan(y/x);
#endif
	    
            return x > 0 ? SPHERE_PLUS_X : SPHERE_MINUS_X;
          }
	
        /*   +/-y spheres   */
        else if ((ay >= ax) && (ay >= az))
          {
#ifdef USE_ATAN2
            localcoord[0] = atan2(y >= 0 ? z : -z, ay);
            localcoord[1] = atan2(y >= 0 ? x : -x, ay);
#else
            assert(y);
            localcoord[0] = atan(z/y);
            localcoord[1] = atan(x/y);
#endif
	    
            return y > 0 ? SPHERE_PLUS_Y : SPHERE_MINUS_Y;
          }

        /*   +/-z spheres   */
        else
          {
#ifdef USE_ATAN2
            localcoord[0] = atan2(z >= 0 ? y : -y, az);
            localcoord[1] = atan2(z >= 0 ? x : -x, az);
#else
            assert(z);
            localcoord[0] = atan(y/z);
            localcoord[1] = atan(x/z);
#endif
	    
            return z > 0 ? SPHERE_PLUS_Z : SPHERE_MINUS_Z;
          }
      }

    CCTK_WARN(0, "Thornburg04: Provided coordinates outside of patches");
    return -1;
  }

  extern "C"
  void
  local_to_global_Thornburg04
  (CCTK_INT patch, CCTK_REAL a, CCTK_REAL b, CCTK_REAL c,
   CCTK_REAL globalcoord[ndirs])
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert(!radial_stretch);
    
    switch (patch) {
    case CENTRAL_CUBE: {
      globalcoord[0] = a;
      globalcoord[1] = b;
      globalcoord[2] = c;
      break;
    }
    case SPHERE_PLUS_X:
    case SPHERE_MINUS_X: {
      CCTK_REAL const dir = patch==SPHERE_PLUS_X ? +1.0 : -1.0;
      CCTK_REAL const tan_a = tan(a);
      CCTK_REAL const tan_b = tan(b);
      CCTK_REAL const R = c;
      // assert(fabs(tan_a) <= 1.0);
      // assert(fabs(tan_b) <= 1.0);
      assert(R > 0.0);
      CCTK_REAL const x = dir * R / sqrt(1.0 + pow(tan_a,2) + pow(tan_b,2));
      CCTK_REAL const z = x * tan(a);
      CCTK_REAL const y = x * tan(b);
      globalcoord[0] = x;
      globalcoord[1] = y;
      globalcoord[2] = z;
      break;
    }
    case SPHERE_PLUS_Y:
    case SPHERE_MINUS_Y: {
      CCTK_REAL const dir = patch==SPHERE_PLUS_Y ? +1.0 : -1.0;
      CCTK_REAL const tan_a = tan(a);
      CCTK_REAL const tan_b = tan(b);
      CCTK_REAL const R = c;
      // assert(fabs(tan_a) <= 1.0);
      // assert(fabs(tan_b) <= 1.0);
      assert(R > 0.0);
      CCTK_REAL const y = dir * R / sqrt(1.0 + pow(tan_a,2) + pow(tan_b,2));
      CCTK_REAL const z = y * tan(a);
      CCTK_REAL const x = y * tan(b);
      globalcoord[0] = x;
      globalcoord[1] = y;
      globalcoord[2] = z;
      break;
    }
    case SPHERE_PLUS_Z:
    case SPHERE_MINUS_Z: {
      CCTK_REAL const dir = patch==SPHERE_PLUS_Z ? +1.0 : -1.0;
      CCTK_REAL const tan_a = tan(a);
      CCTK_REAL const tan_b = tan(b);
      CCTK_REAL const R = c;
      // assert(fabs(tan_a) <= 1.0);
      // assert(fabs(tan_b) <= 1.0);
      assert(R > 0.0);
      CCTK_REAL const z = dir * R / sqrt(1.0 + pow(tan_a,2) + pow(tan_b,2));
      CCTK_REAL const y = z * tan(a);
      CCTK_REAL const x = z * tan(b);
      globalcoord[0] = x;
      globalcoord[1] = y;
      globalcoord[2] = z;
      break;
    }
    default:
      CCTK_WARN(0, "Thornburg04: Provided patch number is incorrect");
    }
  }
    
  extern "C"
  void
  dadx_Thornburg04
  (CCTK_INT const patch,
   CCTK_REAL const xp, CCTK_REAL const yp, CCTK_REAL const zp,
   CCTK_REAL J[ndirs][ndirs])
  {
    DECLARE_CCTK_PARAMETERS;

    const CCTK_REAL Rmin = stretch_rmin_1;
    const CCTK_REAL Rmax = stretch_rmax_1;

    const CCTK_REAL h0 = h_radial;
    const CCTK_REAL h1 = ( h_radial_1<0 ? h0 : h_radial_1 );
    const CCTK_REAL Rstart = sphere_inner_radius;
    
    const CCTK_INT do_radial_stretch = radial_stretch 
      && (abs(h1 - h0) > EPSILON);
    
    switch (patch)
      {
      case CENTRAL_CUBE:
        {
          for (int i=0; i < 3; ++i)
            for (int j=0; j < 3; ++j)
              J[j][i] = (i == j);
        }
        break;
      
      case SPHERE_PLUS_X:
        {
	  if (do_radial_stretch)
	    {
#include "include/thornburg04_x_plus_J.hh"
	    }
	  else
	    {
#include "include/thornburg04_nostretch_x_plus_J.hh"
	    }

        }
        break;
      
      case SPHERE_MINUS_X:
        {
	  if (do_radial_stretch)
	    {
#include "include/thornburg04_x_plus_J.hh"
	    }
	  else
	    {
#include "include/thornburg04_nostretch_x_plus_J.hh"
	    }
        }
        break;
      
      case SPHERE_PLUS_Y:
        {
	  if (do_radial_stretch)
	    {
#include "include/thornburg04_y_plus_J.hh"
	    }
	  else
	    {
#include "include/thornburg04_nostretch_y_plus_J.hh"
	    }
        }
        break;
      
      case SPHERE_MINUS_Y:
        {
	  if (do_radial_stretch)
	    {
#include "include/thornburg04_y_plus_J.hh"
	    }
	  else
	    {
#include "include/thornburg04_nostretch_y_plus_J.hh"
	    }
        }
        break;
        
      case SPHERE_PLUS_Z:
        {
	  if (do_radial_stretch)
	    {
#include "include/thornburg04_z_plus_J.hh"
	    }
	  else
	    {
#include "include/thornburg04_nostretch_z_plus_J.hh"
	    }
        }
        break;
      
      case SPHERE_MINUS_Z:
        {
	  if (do_radial_stretch)
	    {
#include "include/thornburg04_z_plus_J.hh"
	    }
	  else
	    {
#include "include/thornburg04_nostretch_z_plus_J.hh"
	    }
        }
        break;
      
      default:
        {
          CCTK_WARN(0,"Patch number must be 0..6");
        }
      }
  } 
  
  extern "C"
  void
  ddadxdx_Thornburg04
  (CCTK_INT const patch,
   CCTK_REAL const xp, CCTK_REAL const yp, CCTK_REAL const zp,
   CCTK_REAL dJ[ndirs][ndirs][ndirs])
  {
    DECLARE_CCTK_PARAMETERS;

    const CCTK_REAL Rmin = stretch_rmin_1;
    const CCTK_REAL Rmax = stretch_rmax_1;

    const CCTK_REAL h0 = h_radial;
    const CCTK_REAL h1 = ( h_radial_1<0 ? h0 : h_radial_1 );
    const CCTK_REAL Rstart = sphere_inner_radius;

    const CCTK_INT do_radial_stretch = radial_stretch 
      && (abs(h1 - h0) > EPSILON);


    switch (patch)
      {
      case CENTRAL_CUBE:
        {
          for (int i=0; i < 3; ++i)
            for (int j=0; j < 3; ++j)
              for (int k=0; k < 3; ++k)
                dJ[k][j][i] = 0;
        }
        break;
      
      case SPHERE_PLUS_X:
        {
	  if (do_radial_stretch)
	    {
#include "include/thornburg04_x_plus_dJ.hh"
	    }
	  else
	    {
#include "include/thornburg04_nostretch_x_plus_dJ.hh"
	    }

        }
        break;
      
      case SPHERE_MINUS_X:
        {
	  if (do_radial_stretch)
	    {
#include "include/thornburg04_x_plus_dJ.hh"
	    }
	  else
	    {
#include "include/thornburg04_nostretch_x_plus_dJ.hh"
	    }
        }
        break;
      
      case SPHERE_PLUS_Y:
        {
	  if (do_radial_stretch)
	    {
#include "include/thornburg04_y_plus_dJ.hh"
	    }
	  else
	    {
#include "include/thornburg04_nostretch_y_plus_dJ.hh"
	    }
        }
        break;
      
      case SPHERE_MINUS_Y:
        {
	  if (do_radial_stretch)
	    {
#include "include/thornburg04_y_plus_dJ.hh"
	    }
	  else
	    {
#include "include/thornburg04_nostretch_y_plus_dJ.hh"
	    }
        }
        break;
        
      case SPHERE_PLUS_Z:
        {
	  if (do_radial_stretch)
	    {
#include "include/thornburg04_z_plus_dJ.hh"
	    }
	  else
	    {
#include "include/thornburg04_nostretch_z_plus_dJ.hh"
	    }
        }
        break;
      
      case SPHERE_MINUS_Z:
        {
	  if (do_radial_stretch)
	    {
#include "include/thornburg04_z_plus_dJ.hh"
	    }
	  else
	    {
#include "include/thornburg04_nostretch_z_plus_dJ.hh"
	    }
        }
        break;
      
      default:
        {
          CCTK_WARN(0,"Patch number must be 0..6");
        }
      }
  }
  
  
  // This computes the overlap of a Cartesian cell with a sphere of radius r
  // and is needed for the computation of the volume form
  double calculate_alpha(const double x0, const double y0, const double z0, 
                         const double xc, const double yc, const double zc, 
                         const double dx, const double dy, const double dz, 
                         const double r, const int nMonteCarloParticles,
                         const int seed)
  {
     const int n = nMonteCarloParticles; // number of random numbers used for calculation of the volume fraction
     int n_points_in = 0;

     double const r2 = pow(r,2);
     double x, y, z, rp2;

     unsigned int rand_state = (unsigned int)seed;

     for (int i=0; i<n; ++i) {
        x = xc+dx*(0.5-1.0*rand_r(&rand_state)/(double)RAND_MAX);
        y = yc+dy*(0.5-1.0*rand_r(&rand_state)/(double)RAND_MAX);
        z = zc+dz*(0.5-1.0*rand_r(&rand_state)/(double)RAND_MAX);
        rp2 = pow(x-x0,2)+pow(y-y0,2)+pow(z-z0,2);
        if (rp2 < r2)
           ++n_points_in;   
     }
     return (double)n_points_in/(double)n;
  } 
  
  // Set weight function
  extern "C"
  void Coordinates_SetVolumeForm_Thornburg04(CCTK_ARGUMENTS)
  {
     DECLARE_CCTK_ARGUMENTS_CHECKED(Coordinates_SetVolumeForm_Thornburg04)
     DECLARE_CCTK_PARAMETERS
     
     CCTK_INT const map = MultiPatch_GetMap (cctkGH);
     assert (map >= 0);
     
     if (!store_volume_form)
     {
        *volume_form_state = 0;
        return;
     }
     
     CCTK_INFO("Setting up volume form for thornburg04.");
     
     *volume_form_state = 1;
     
     // Determine the extent of the intermediate (see comment above)
     // coordinate system
     CCTK_REAL exterior_min[3];
     CCTK_REAL exterior_max[3];
     CCTK_REAL h[3];
     get_intermediate_coordinates
        (map, 3,
         exterior_min, exterior_max,
         h);
     
     switch (map)
     {
        case CENTRAL_CUBE:
        {
        
           // Compute volume form on coarse grid only once initially
           if (REFLEVEL == 0 && cctk_iteration != 0)
              return;
        
           // handle central patch separately
           // TODO: Assuming the coarse grid does not change, we can 
           //       compute the overlap of cells (which is expensive) only initially.
           //       We should not recompute after recovery, since the computation is (pseudo) stochastic!!!
           #pragma omp parallel for
           for (int k=0; k<cctk_lsh[2]; ++k)
              for (int j=0; j<cctk_lsh[1]; ++j)
                 for (int i=0; i<cctk_lsh[0]; ++i)
                 {
                    const int ipos[3] = { i, j, k };
                    const int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
                    CCTK_REAL a[3];
                    get_local_coordinates (cctkGH, exterior_min, h, ipos, a);
            
                    // coordinates of current cell center
                    const CCTK_REAL xc[3] = { x[ijk], y[ijk], z[ijk] };
                    
                    // get radii of vertices of a cell and
                    // count how many of them are inside sphere_inner_radius
                    const CCTK_REAL rv[8] = { sqrt(pow(xc[0]-h[0]/2, 2.) + pow(xc[1]-h[1]/2, 2.) + pow(xc[2]-h[2]/2, 2.)),
                                              sqrt(pow(xc[0]-h[0]/2, 2.) + pow(xc[1]-h[1]/2, 2.) + pow(xc[2]+h[2]/2, 2.)),
                                              sqrt(pow(xc[0]-h[0]/2, 2.) + pow(xc[1]+h[1]/2, 2.) + pow(xc[2]+h[2]/2, 2.)),
                                              sqrt(pow(xc[0]+h[0]/2, 2.) + pow(xc[1]+h[1]/2, 2.) + pow(xc[2]+h[2]/2, 2.)),
                                              sqrt(pow(xc[0]+h[0]/2, 2.) + pow(xc[1]-h[1]/2, 2.) + pow(xc[2]-h[2]/2, 2.)),
                                              sqrt(pow(xc[0]+h[0]/2, 2.) + pow(xc[1]+h[1]/2, 2.) + pow(xc[2]-h[2]/2, 2.)),
                                              sqrt(pow(xc[0]+h[0]/2, 2.) + pow(xc[1]-h[1]/2, 2.) + pow(xc[2]+h[2]/2, 2.)),
                                              sqrt(pow(xc[0]-h[0]/2, 2.) + pow(xc[1]+h[1]/2, 2.) + pow(xc[2]-h[2]/2, 2.)) };
                    
                    int n_vertices_in = 0;
                    int n_vertices_out = 0;
                    
                    for (int n=0; n < 8; ++n)
                    {
                       if (rv[n] > sphere_inner_radius)
                          ++n_vertices_out;
                       else
                          ++n_vertices_in;
                    }
                    
                    if (n_vertices_out == 0 && n_vertices_in > 0)
                       // cell is completely inside
                       // set volume form to product of _coarse_ grid spacing (Carpet will rescale this appropriately during reduction for finer levels!)
                       volume_form[ijk] = cctk_delta_space[0]*cctk_delta_space[1]*cctk_delta_space[2];
                    else if (n_vertices_in == 0 && n_vertices_out > 0)
                       // cell is completely outside
                       volume_form[ijk] = 0;
                    else {
                       // cell is overlapping with spherical grid
                       const CCTK_REAL alpha = calculate_alpha(0,0,0, xc[0],xc[1],xc[2], h[0],h[1],h[2], sphere_inner_radius, nMonteCarloParticles, MonteCarloSeed);
                       volume_form[ijk] = alpha * cctk_delta_space[0]*cctk_delta_space[1]*cctk_delta_space[2];
                    }
                 }
           
           break;
        }
        default:
        {
           // angular patches are easy to handle
           for (int k=0; k<cctk_lsh[2]; ++k)
              for (int j=0; j<cctk_lsh[1]; ++j)
                 for (int i=0; i<cctk_lsh[0]; ++i)
                 {
                    const int ipos[3] = { i, j, k };
                    const int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
                    CCTK_REAL a[3];
                    get_local_coordinates (cctkGH, exterior_min, h, ipos, a);
            
                    CCTK_REAL const rho = a[2];
                    
                    if (a[0] < -PI/4 || a[0] > PI/4 
                     || a[1] < -PI/4 || a[1] > PI/4
                     || rho < sphere_inner_radius || rho > sphere_outer_radius)
                    {
                       // set volume form to zero since we are outside of the nominal domain
                       volume_form[ijk] = 0;
                    }
                    else
                    {
                       // set volume form to deterimant of Jacobian
                       const CCTK_REAL det =  fabs((  J11[ijk] * J22[ijk] * J33[ijk]
                            			    + J12[ijk] * J23[ijk] * J31[ijk]
                            			    + J13[ijk] * J21[ijk] * J32[ijk]
                            			    - J11[ijk] * J23[ijk] * J32[ijk]
                            			    - J12[ijk] * J21[ijk] * J33[ijk]
                            			    - J13[ijk] * J22[ijk] * J31[ijk]));
                       
                       volume_form[ijk] = h[0]*h[1]*h[2]/det;
                    }
                 }
           break;
        }
     }
  }
  
}
