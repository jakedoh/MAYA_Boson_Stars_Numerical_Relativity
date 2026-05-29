
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
 * Sets up a 3-patch multipatch grid using spherical coordinates but with the
 * polar directions (from theta=0 to theta=theta_min) covered by column
 * patches. Based on the corresponding patch system in the old multipatch
 * infrastructure by Burkhard Zink.
 *
 * On the spherical grid:
 *
 * (a,b,c) == (theta,phi,r)
 * X = r*sin(theta)*cos(phi)
 * Y = r*sin(theta)*sin(phi)
 * Z = r*cos(theta)  
 *
 * Polar +/- patches:
 * (a,b,c) == (a,b,r)
 * Define:
 * alpha = sin(theta_min)
 * root = sqrt(1-alpha^2*(a^2+b^2))
 * sign = +1/-1 depending on pole

 * X = alpha*a*r
 * Y = alpha*b*r
 * Z = sign*r*root

 * The local coordinate ranges on the spherical patch are:
 * theta: [theta_min,Pi-theta_min]
 * phi:   [-Pi,Pi]
 * r:     [rmin,rmax]

 * The local coordinate ranges on the column patches are:
 * a: [-1,1]
 * b: [-1,1]
 * r: [rmin,rmax]
 *
 */
#include <cassert>
#include <cmath>
#include <iostream>

#include <cctk.h>
#include <cctk_Parameters.h>
#include <cctk_Arguments.h>

#include "coordinates.hh"
#include "patchsystem.hh"

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
  Coordinates_SetGlobalCoords_SphereColumn(CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS_CHECKED(Coordinates_SetGlobalCoords_SphereColumn);
    DECLARE_CCTK_PARAMETERS;
    
    CCTK_INT const map = MultiPatch_GetMap (cctkGH);
    assert (map >= 0);
    
    CCTK_INFO("Setting up global coordinates for sphere+column.");

    // Determine the extent of the intermediate (see comment above)
    // coordinate system
    CCTK_REAL exterior_min[3];
    CCTK_REAL exterior_max[3];
    CCTK_REAL h[3];
    get_intermediate_coordinates
      (map, 3,
       exterior_min, exterior_max,
       h);

    CCTK_REAL Rmin = stretch_rmin_1;
    CCTK_REAL Rmax = stretch_rmax_1;
    CCTK_REAL Rstart = sphere_inner_radius;

    CCTK_REAL h0 = h_radial;
    CCTK_REAL h1 = h_radial_1;

    CCTK_REAL alpha = sin(theta_min/180.0*PI);

    switch (map) {
      
    case SPC_SPHERE: {
#pragma omp parallel for collapse(3)
      for (int k=0; k<cctk_lsh[2]; ++k)
        for (int j=0; j<cctk_lsh[1]; ++j)
          for (int i=0; i<cctk_lsh[0]; ++i)
            {
              const int ipos[3] = { i, j, k };
              const int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
              CCTK_REAL a[3];
              get_local_coordinates (cctkGH, exterior_min, h, ipos, a);

              CCTK_REAL const sinphi = sin(a[1]);
              CCTK_REAL const cosphi = cos(a[1]);
              CCTK_REAL const sintheta = sin(a[0]);
              CCTK_REAL const costheta = cos(a[0]);

              CCTK_REAL const rho = a[2];

              CCTK_REAL R;
              if (radial_stretch)
                {
                  CCTK_REAL Rl=rho;
#include "include/thornburg04_local_to_global.hh"
                }
              else
                R = rho;

              x[ijk] = R*cosphi*sintheta;
              y[ijk] = R*sinphi*sintheta;
              z[ijk] = R*costheta;
              r[ijk] = R;
//              printf("(i,j,k) = %d,%d,%d\n",i,j,k);
//              printf("(x,y,z,r) = %f,%f,%f,%f\n",x[ijk],y[ijk],z[ijk],r[ijk]);
            }
      break;
    }
      
    case SPC_COLUMN_PLUS_Z: {
#pragma omp parallel for collapse(3)
      for (int k=0; k<cctk_lsh[2]; ++k)
        for (int j=0; j<cctk_lsh[1]; ++j)
          for (int i=0; i<cctk_lsh[0]; ++i)
            {
              const int ipos[3] = { i, j, k };
              const int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
              CCTK_REAL a[3];
              get_local_coordinates (cctkGH, exterior_min, h, ipos, a);
            
              CCTK_REAL const root = sqrt(1.0-pow(alpha,2)*(pow(a[0],2)+pow(a[1],2)));
              CCTK_REAL const rho = a[2];

              CCTK_REAL R;
	      if (radial_stretch)
		{
		  CCTK_REAL Rl=rho;
#include "include/thornburg04_local_to_global.hh"
		}
	      else
		R = rho;

              x[ijk] = alpha*a[0]*R;
              y[ijk] = alpha*a[1]*R;
              z[ijk] = R*root;
              r[ijk] = R;
//              printf("(i,j,k) = %d,%d,%d\n",i,j,k);
//              printf("(x,y,z,r) = %f,%f,%f,%f\n",x[ijk],y[ijk],z[ijk],r[ijk]);
            }
      break;
    }
      
    case SPC_COLUMN_MINUS_Z: {
#pragma omp parallel for collapse(3)
      for (int k=0; k<cctk_lsh[2]; ++k)
        for (int j=0; j<cctk_lsh[1]; ++j)
          for (int i=0; i<cctk_lsh[0]; ++i)
            {
              const int ipos[3] = { i, j, k };
              const int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
              CCTK_REAL a[3];
              get_local_coordinates (cctkGH, exterior_min, h, ipos, a);
            
              CCTK_REAL const root = sqrt(1.0-pow(alpha,2)*(pow(a[0],2)+pow(a[1],2)));
              CCTK_REAL const rho = a[2];

              CCTK_REAL R;
	      if (radial_stretch)
		{
		  CCTK_REAL Rl=rho;
#include "include/thornburg04_local_to_global.hh"
		}
	      else
		R = rho;

              x[ijk] = alpha*a[0]*R;
              y[ijk] = alpha*a[1]*R;
              z[ijk] = -R*root;
              r[ijk] = R;
//              printf("(i,j,k) = %d,%d,%d\n",i,j,k);
//              printf("(x,y,z,r) = %f,%f,%f,%f\n",x[ijk],y[ijk],z[ijk],r[ijk]);
            }
      break;
    }
      
    default:
      CCTK_WARN(0, "Patch number must be 0..2");
      
    }
  }

  extern "C"
  void
  Coordinates_SetJacobian_SphereColumn(CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS_CHECKED(Coordinates_SetJacobian_SphereColumn);
    DECLARE_CCTK_PARAMETERS;

    using namespace std;

    int i, j, k;

    CCTK_INFO("Setting jacobian: sphere+column");

    *general_coordinates = 1;
    *interpolate_boundary_points = 0; /* not sure about this - safer to say no */
    *jacobian_state = 1;
    *jacobian_derivative_state = 1;

    CCTK_REAL Rmin = stretch_rmin_1;
    CCTK_REAL Rmax = stretch_rmax_1;

    CCTK_REAL h0 = h_radial;
    CCTK_REAL h1 = h_radial_1;

    const CCTK_INT patch = Coordinates_GetMap(cctkGH);

    CCTK_REAL alpha = sin(theta_min/180.0*PI);

    switch (patch)
      {
      case SPC_SPHERE:
        {
          for (k=0; k<cctk_lsh[2]; ++k)
            for (j=0; j<cctk_lsh[1]; ++j)
              for (i=0; i<cctk_lsh[0]; ++i)
                {
                  const int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

                  const CCTK_REAL xp = x[ijk];
                  const CCTK_REAL yp = y[ijk];
                  const CCTK_REAL zp = z[ijk];
                  CCTK_REAL J[3][3];
                  CCTK_REAL dJ[3][3][3];

                  if (radial_stretch)
                    {
#include "include/spherecolumn_sphere_J.hh"
#include "include/spherecolumn_sphere_dJ.hh"
                    }
                  else
                    {
#include "include/spherecolumn_nostretch_sphere_J.hh"
#include "include/spherecolumn_nostretch_sphere_dJ.hh"
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

      case SPC_COLUMN_PLUS_Z:
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
		  
		  if (radial_stretch)
		    {
#include "include/spherecolumn_column_J.hh"
#include "include/spherecolumn_column_dJ.hh"
		    }
		  else
		    {
#include "include/spherecolumn_nostretch_column_J.hh"
#include "include/spherecolumn_nostretch_column_dJ.hh"
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

      case SPC_COLUMN_MINUS_Z:
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

		  if (radial_stretch)
		    {
#include "include/spherecolumn_column_J.hh"
#include "include/spherecolumn_column_dJ.hh"
		    }
		  else
		    {
#include "include/spherecolumn_nostretch_column_J.hh"
#include "include/spherecolumn_nostretch_column_dJ.hh"
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
        CCTK_WARN(0, "Patch number must be 0..2");

      }
    
    return;
  }
  
  
  
  extern "C"
  CCTK_INT
  global_to_local_SphereColumn
  (CCTK_REAL const x, CCTK_REAL const y, CCTK_REAL const z,
   CCTK_REAL localcoord[ndirs])
  {
    DECLARE_CCTK_PARAMETERS;

    CCTK_REAL Rmin = stretch_rmin_1;
    CCTK_REAL Rmax = stretch_rmax_1;
    CCTK_REAL Rstart = sphere_inner_radius;

    CCTK_REAL h0 = h_radial;
    CCTK_REAL h1 = h_radial_1;

    const CCTK_REAL rp2 = x*x + y*y + z*z;
    
    const CCTK_REAL theta = atan2(sqrt(x*x+y*y),z);
    const CCTK_REAL alpha = sin(theta_min/180.0*PI);

    const CCTK_REAL Rg = sqrt(rp2);

    CCTK_REAL R;
    if (radial_stretch)
      {
#include "include/thornburg04_global_to_local.hh"
      }
    else
      R = Rg;

    /* SPC_COLUMN_PLUS_Z */
    if (theta < theta_min/180.0*PI)
      {
        localcoord[0] = x/(alpha*R);
        localcoord[1] = y/(alpha*R);
        localcoord[2] = R;

//        printf("SPC_COLUMN_PLUS_Z: (%f,%f,%f)\n",localcoord[0],localcoord[1],localcoord[2]);
        return SPC_COLUMN_PLUS_Z;
      }
      /* SPC_SPHERE */
      else if (theta_min/180.0*PI<=theta && theta<=PI-theta_min/180.0*PI)
      { 
        localcoord[0] = theta;
        localcoord[1] = atan2(y,x);
        localcoord[2] = R;
	    
//        printf("SPC_SPHERE: (%f,%f,%f)\n",localcoord[0],localcoord[1],localcoord[2]);
        return SPC_SPHERE;
      }
      /* SPC_COLUMN_MINUS_Z */
      else
      {
        localcoord[0] = x/(alpha*R);
        localcoord[1] = y/(alpha*R);
        localcoord[2] = R;
	    
//        printf("SPC_COLUMN_MINUS_Z: (%f,%f,%f)\n",localcoord[0],localcoord[1],localcoord[2]);
        return SPC_COLUMN_MINUS_Z;
      }

    CCTK_WARN(0, "SphereColumn: Provided coordinates outside of patches");
    return -1;
  }
    
  extern "C"
  void
  dadx_SphereColumn
  (CCTK_INT const patch,
   CCTK_REAL const xp, CCTK_REAL const yp, CCTK_REAL const zp,
   CCTK_REAL J[ndirs][ndirs])
  {
    DECLARE_CCTK_PARAMETERS;

    CCTK_REAL Rmin = stretch_rmin_1;
    CCTK_REAL Rmax = stretch_rmax_1;

    CCTK_REAL h0 = h_radial;
    CCTK_REAL h1 = h_radial_1;
    
    const CCTK_REAL alpha = sin(theta_min/180.0*PI);

    switch (patch)
      {
      case SPC_SPHERE:
        {
          if (radial_stretch)
            {
#include "include/spherecolumn_sphere_J.hh"
            }
          else
            {
#include "include/spherecolumn_nostretch_sphere_J.hh"
            }

        }
        break;
      
      case SPC_COLUMN_PLUS_Z:
        {
	  if (radial_stretch)
	    {
#include "include/spherecolumn_column_J.hh"
	    }
	  else
	    {
#include "include/spherecolumn_nostretch_column_J.hh"
	    }

        }
        break;
      
      case SPC_COLUMN_MINUS_Z:
        {
	  if (radial_stretch)
	    {
#include "include/spherecolumn_column_J.hh"
	    }
	  else
	    {
#include "include/spherecolumn_nostretch_column_J.hh"
	    }
        }
        break;
      
      default:
        {
          CCTK_WARN(0,"Patch number must be 0..2");
        }
      }
  } 
  
  extern "C"
  void
  ddadxdx_SphereColumn
  (CCTK_INT const patch,
   CCTK_REAL const xp, CCTK_REAL const yp, CCTK_REAL const zp,
   CCTK_REAL dJ[ndirs][ndirs][ndirs])
  {
    DECLARE_CCTK_PARAMETERS;

    CCTK_REAL Rmin = stretch_rmin_1;
    CCTK_REAL Rmax = stretch_rmax_1;

    CCTK_REAL h0 = h_radial;
    CCTK_REAL h1 = h_radial_1;

    const CCTK_REAL alpha = sin(theta_min/180.0*PI);

    if (h1 < 0)
      h1 = h0;

    switch (patch)
      {
      
      case SPC_SPHERE:
        {
	  if (radial_stretch)
	    {
#include "include/spherecolumn_sphere_dJ.hh"
	    }
	  else
	    {
#include "include/spherecolumn_nostretch_sphere_dJ.hh"
	    }

        }
        break;
      
      case SPC_COLUMN_PLUS_Z:
        {
	  if (radial_stretch)
	    {
#include "include/spherecolumn_column_dJ.hh"
	    }
	  else
	    {
#include "include/spherecolumn_nostretch_column_dJ.hh"
	    }
        }
        break;
      
      case SPC_COLUMN_MINUS_Z:
        {
	  if (radial_stretch)
	    {
#include "include/spherecolumn_column_dJ.hh"
	    }
	  else
	    {
#include "include/spherecolumn_nostretch_column_dJ.hh"
	    }
        }
        break;
      
      default:
        {
          CCTK_WARN(0,"Patch number must be 0..2");
        }
      }
  }
}
