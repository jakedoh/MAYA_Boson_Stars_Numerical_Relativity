
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

#include <cctk.h>
#include <cctk_Parameters.h>

#include <carpet.hh>

#include "patchsystem.hh"
#include "coordinates.hh"


namespace Coordinates {

  using namespace std;
  
  
  extern "C"
  CCTK_INT
  Coordinates_GlobalToLocal (CCTK_POINTER_TO_CONST const cctkGH_,
                             CCTK_INT const ndims,
                             CCTK_INT const npoints,
                             CCTK_POINTER_TO_CONST const globalcoords,
                             CCTK_INT * patch,
                             CCTK_POINTER const localcoords,
                             CCTK_POINTER const dadx_,
                             CCTK_POINTER const ddadxdx_)
  {
    DECLARE_CCTK_PARAMETERS
    (void)cctkGH_; // pacify compiler
    
    assert(ndims == ndirs);

    CCTK_REAL* const * const       abc = (CCTK_REAL**) localcoords;
    const CCTK_REAL* const * const xyz = (CCTK_REAL**) globalcoords;
    CCTK_REAL* const* const        dadx = (CCTK_REAL**) dadx_;
    CCTK_REAL* const* const        ddadxdx = (CCTK_REAL**) ddadxdx_;
    
    
    
    if (CCTK_Equals(coordinate_system, "cartesian"))
    {
      // a Cartesian one-patch system global coords = local coords
#pragma omp parallel for
      for (int ijk=0; ijk < npoints; ++ijk)
      {
        patch[ijk] = 0;
        if (abc)
        {
          for (int d=0; d<ndirs; ++d)
          {
            abc[d][ijk] = xyz[d][ijk];
          }
        }
        if (dadx)
        {
          for (int d=0; d<ndirs*ndirs; ++d)
          {
            dadx[d][ijk] = (d/ndirs == d%ndirs);
          }
        }
        if (ddadxdx)
        {
          for (int d=0; d<ndirs*ndirs*(ndirs+1)/2; ++d)
          {
            ddadxdx[d][ijk] = 0;
          }
        }
      }
    }
    else if (CCTK_Equals(coordinate_system, "twopatchcartesian"))
    {
#pragma omp parallel for
      for (int ijk=0; ijk < npoints; ++ijk)
      {
        CCTK_REAL localcoord[ndirs];
        patch[ijk] = global_to_local_TwoPatchCartesian
          (xyz[0][ijk], xyz[1][ijk], xyz[2][ijk], localcoord);
        if (abc) {
           for (int d=0; d<ndirs; ++d)
           {
             abc[d][ijk] = localcoord[d];
           }
        }
        if (dadx)
        {
          CCTK_REAL J[ndirs][ndirs];
          dadx_TwoPatchCartesian
            (patch[ijk], xyz[0][ijk], xyz[1][ijk], xyz[2][ijk], J);
          for (int d=0; d<ndirs; ++d) {
            for (int e=0; e<ndirs; ++e) {
              dadx[d*ndirs+e][ijk] = J[e][d];
            }
          }
        }
        if (ddadxdx)
        {
          CCTK_REAL dJ[ndirs][ndirs][ndirs];
          ddadxdx_TwoPatchCartesian
            (patch[ijk], xyz[0][ijk], xyz[1][ijk], xyz[2][ijk], dJ);
          for (int d=0; d<ndirs; ++d) {
            for (int e=0; e<ndirs; ++e) {
              for (int f=e; f<ndirs; ++f) {
                const int idx =
                  d*ndirs*(ndirs+1)/2 + ((ndirs-1)*e - e*(e-1)/2) + f;
                ddadxdx[idx][ijk] = dJ[f][e][d];
              }
            }
          }
        }
      }
    }
    else if (CCTK_Equals(coordinate_system, "TwoPatchDistorted"))
    {
#pragma omp parallel for
      for (int ijk=0; ijk < npoints; ++ijk)
      {
        CCTK_REAL localcoord[ndirs];
        patch[ijk] = global_to_local_TwoPatchDistorted
          (xyz[0][ijk], xyz[1][ijk], xyz[2][ijk], localcoord);
        if (abc)
        {
          for (int d=0; d<ndirs; ++d)
          {
            abc[d][ijk] = localcoord[d];
          }
        }
        if (dadx)
        {
          CCTK_REAL J[ndirs][ndirs];
          dadx_TwoPatchDistorted
            (patch[ijk], xyz[0][ijk], xyz[1][ijk], xyz[2][ijk], J);
          for (int d=0; d<ndirs; ++d) {
            for (int e=0; e<ndirs; ++e) {
              dadx[d*ndirs+e][ijk] = J[e][d];
            }
          }
        }
        if (ddadxdx)
        {
          CCTK_REAL dJ[ndirs][ndirs][ndirs];
          ddadxdx_TwoPatchDistorted
            (patch[ijk], xyz[0][ijk], xyz[1][ijk], xyz[2][ijk], dJ);
          for (int d=0; d<ndirs; ++d) {
            for (int e=0; e<ndirs; ++e) {
              for (int f=e; f<ndirs; ++f) {
                const int idx =
                  d*ndirs*(ndirs+1)/2 + ((ndirs-1)*e - e*(e-1)/2) + f;
                ddadxdx[idx][ijk] = dJ[f][e][d];
              }
            }
          }
        }
      }
    }
    else if (CCTK_Equals(coordinate_system, "Thornburg04"))
    {
#pragma omp parallel for
      for (int ijk=0; ijk < npoints; ++ijk)
      {
        CCTK_REAL localcoord[ndirs];
        patch[ijk] = global_to_local_Thornburg04
          (xyz[0][ijk], xyz[1][ijk], xyz[2][ijk], localcoord);
        if (abc)
        {
          for (int d=0; d<ndirs; ++d)
          {
            abc[d][ijk] = localcoord[d];
          }
        }
        if (dadx)
        {
          CCTK_REAL J[ndirs][ndirs];
          dadx_Thornburg04
            (patch[ijk], xyz[0][ijk], xyz[1][ijk], xyz[2][ijk], J);
          for (int d=0; d<ndirs; ++d) {
            for (int e=0; e<ndirs; ++e) {
              dadx[d*ndirs+e][ijk] = J[e][d];
            }
          }
        }
        if (ddadxdx)
        {
          CCTK_REAL dJ[ndirs][ndirs][ndirs];
          ddadxdx_Thornburg04
            (patch[ijk], xyz[0][ijk], xyz[1][ijk], xyz[2][ijk], dJ);
          for (int d=0; d<ndirs; ++d) {
            for (int e=0; e<ndirs; ++e) {
              for (int f=e; f<ndirs; ++f) {
                const int idx =
                  d*ndirs*(ndirs+1)/2 + ((ndirs-1)*e - e*(e-1)/2) + f;
                ddadxdx[idx][ijk] = dJ[f][e][d];
              }
            }
          }
        }
      }
    }
    else if (CCTK_Equals(coordinate_system, "Thornburg13"))
    {
#pragma omp parallel for
      for (int ijk=0; ijk < npoints; ++ijk)
      {
        CCTK_REAL localcoord[ndirs];
        patch[ijk] = global_to_local_Thornburg13
          (xyz[0][ijk], xyz[1][ijk], xyz[2][ijk], localcoord);
        if (abc)
        {
          for (int d=0; d<ndirs; ++d)
          {
            abc[d][ijk] = localcoord[d];
          }
        }
        if (dadx)
        {
          CCTK_REAL J[ndirs][ndirs];
          dadx_Thornburg13
            (patch[ijk], xyz[0][ijk], xyz[1][ijk], xyz[2][ijk], J);
          for (int d=0; d<ndirs; ++d) {
            for (int e=0; e<ndirs; ++e) {
              dadx[d*ndirs+e][ijk] = J[e][d];
            }
          }
        }
        if (ddadxdx)
        {
          CCTK_REAL dJ[ndirs][ndirs][ndirs];
          ddadxdx_Thornburg13
            (patch[ijk], xyz[0][ijk], xyz[1][ijk], xyz[2][ijk], dJ);
          for (int d=0; d<ndirs; ++d) {
            for (int e=0; e<ndirs; ++e) {
              for (int f=e; f<ndirs; ++f) {
                const int idx =
                  d*ndirs*(ndirs+1)/2 + ((ndirs-1)*e - e*(e-1)/2) + f;
                ddadxdx[idx][ijk] = dJ[f][e][d];
              }
            }
          }
        }
      }
    }

    else if (CCTK_Equals(coordinate_system, "Thornburg04nc"))
    {
#pragma omp parallel for
      for (int ijk=0; ijk < npoints; ++ijk)
      {
        CCTK_REAL localcoord[ndirs];
        patch[ijk] = global_to_local_Thornburg04nc
          (xyz[0][ijk], xyz[1][ijk], xyz[2][ijk], localcoord);
        if (abc)
        {
          for (int d=0; d<ndirs; ++d)
          {
            abc[d][ijk] = localcoord[d];
          }
        }
        if (dadx)
        {
          CCTK_REAL J[ndirs][ndirs];
          dadx_Thornburg04nc
            (patch[ijk], xyz[0][ijk], xyz[1][ijk], xyz[2][ijk], J);
          for (int d=0; d<ndirs; ++d) {
            for (int e=0; e<ndirs; ++e) {
              dadx[d*ndirs+e][ijk] = J[e][d];
            }
          }
        }
        if (ddadxdx)
        {
          CCTK_REAL dJ[ndirs][ndirs][ndirs];
          ddadxdx_Thornburg04nc
            (patch[ijk], xyz[0][ijk], xyz[1][ijk], xyz[2][ijk], dJ);
          for (int d=0; d<ndirs; ++d) {
            for (int e=0; e<ndirs; ++e) {
              for (int f=e; f<ndirs; ++f) {
                const int idx =
                  d*ndirs*(ndirs+1)/2 + ((ndirs-1)*e - e*(e-1)/2) + f;
                ddadxdx[idx][ijk] = dJ[f][e][d];
              }
            }
          }
        }
      }
    }
    else if (CCTK_Equals(coordinate_system, "CylinderInBox"))
    {
#pragma omp parallel for
      for (int ijk=0; ijk < npoints; ++ijk)
      {
        CCTK_REAL localcoord[ndirs];
        patch[ijk] = global_to_local_CylinderInBox
          (xyz[0][ijk], xyz[1][ijk], xyz[2][ijk], localcoord);
        if (abc)
        {
          for (int d=0; d<ndirs; ++d)
          {
            abc[d][ijk] = localcoord[d];
          }
        }
        if (dadx)
        {
          CCTK_REAL J[ndirs][ndirs];
          dadx_CylinderInBox
            (patch[ijk], xyz[0][ijk], xyz[1][ijk], xyz[2][ijk], J);
          for (int d=0; d<ndirs; ++d) {
            for (int e=0; e<ndirs; ++e) {
              dadx[d*ndirs+e][ijk] = J[e][d];
            }
          }
        }
        if (ddadxdx)
        {
          CCTK_REAL dJ[ndirs][ndirs][ndirs];
          ddadxdx_CylinderInBox
            (patch[ijk], xyz[0][ijk], xyz[1][ijk], xyz[2][ijk], dJ);
          for (int d=0; d<ndirs; ++d) {
            for (int e=0; e<ndirs; ++e) {
              for (int f=e; f<ndirs; ++f) {
                const int idx =
                  d*ndirs*(ndirs+1)/2 + ((ndirs-1)*e - e*(e-1)/2) + f;
                ddadxdx[idx][ijk] = dJ[f][e][d];
              }
            }
          }
        }
      }
    }
    else if (CCTK_Equals(coordinate_system, "Sphere+Column"))
    {
#pragma omp parallel for
      for (int ijk=0; ijk < npoints; ++ijk)
      {
        CCTK_REAL localcoord[ndirs];
        patch[ijk] = global_to_local_SphereColumn
          (xyz[0][ijk], xyz[1][ijk], xyz[2][ijk], localcoord);
        if (abc)
        {
          for (int d=0; d<ndirs; ++d)
          {
            abc[d][ijk] = localcoord[d];
          }
        }
        if (dadx)
        {
          CCTK_REAL J[ndirs][ndirs];
          dadx_SphereColumn
            (patch[ijk], xyz[0][ijk], xyz[1][ijk], xyz[2][ijk], J);
          for (int d=0; d<ndirs; ++d) {
            for (int e=0; e<ndirs; ++e) {
              dadx[d*ndirs+e][ijk] = J[e][d];
            }
          }
        }
        if (ddadxdx)
        {
          CCTK_REAL dJ[ndirs][ndirs][ndirs];
          ddadxdx_SphereColumn
            (patch[ijk], xyz[0][ijk], xyz[1][ijk], xyz[2][ijk], dJ);
          for (int d=0; d<ndirs; ++d) {
            for (int e=0; e<ndirs; ++e) {
              for (int f=e; f<ndirs; ++f) {
                const int idx =
                  d*ndirs*(ndirs+1)/2 + ((ndirs-1)*e - e*(e-1)/2) + f;
                ddadxdx[idx][ijk] = dJ[f][e][d];
              }
            }
          }
        }
      }
    }
    else if (CCTK_Equals(coordinate_system, "Cylinder+Column"))
    {
#pragma omp parallel for
      for (int ijk=0; ijk < npoints; ++ijk)
      {
        CCTK_REAL localcoord[ndirs];
        patch[ijk] = global_to_local_CylinderColumn
          (xyz[0][ijk], xyz[1][ijk], xyz[2][ijk], localcoord);
        if (abc)
        {
          for (int d=0; d<ndirs; ++d)
          {
            abc[d][ijk] = localcoord[d];
          }
        }
        if (dadx)
        {
          CCTK_REAL J[ndirs][ndirs];
          dadx_CylinderColumn
            (patch[ijk], xyz[0][ijk], xyz[1][ijk], xyz[2][ijk], J);
          for (int d=0; d<ndirs; ++d) {
            for (int e=0; e<ndirs; ++e) {
              dadx[d*ndirs+e][ijk] = J[e][d];
            }
          }
        }
        if (ddadxdx)
        {
          CCTK_REAL dJ[ndirs][ndirs][ndirs];
          ddadxdx_CylinderColumn
            (patch[ijk], xyz[0][ijk], xyz[1][ijk], xyz[2][ijk], dJ);
          for (int d=0; d<ndirs; ++d) {
            for (int e=0; e<ndirs; ++e) {
              for (int f=e; f<ndirs; ++f) {
                const int idx =
                  d*ndirs*(ndirs+1)/2 + ((ndirs-1)*e - e*(e-1)/2) + f;
                ddadxdx[idx][ijk] = dJ[f][e][d];
              }
            }
          }
        }
      }
    }
    else
    {
      
      CCTK_WARN (CCTK_WARN_ABORT, "unknown patch system");
      
    }

    return 0;
  }
  
} // namespace Coordinates
