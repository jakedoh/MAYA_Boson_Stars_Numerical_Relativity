
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
  Coordinates_LocalToGlobal (CCTK_POINTER_TO_CONST const cctkGH_,
                             CCTK_INT const ndims,
                             CCTK_INT const npoints,
                             CCTK_INT const* const patch,
                             CCTK_POINTER_TO_CONST const localcoords,
                             CCTK_POINTER const globalcoords,
                             CCTK_POINTER const dxda,
                             CCTK_POINTER const det_dxda,
                             CCTK_POINTER const dadx,
                             CCTK_POINTER const ddxdada,
                             CCTK_POINTER const ddadxdx,
                             CCTK_POINTER const dddxdadada)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert(ndims == 3);
    
    assert(not dxda);           // not supported
    assert(not det_dxda);       // not supported
    assert(not ddxdada);        // not supported
    assert(not dddxdadada);     // not supported
    
    (void)cctkGH_; // pacify compiler
    CCTK_REAL* const* const xyz = static_cast<CCTK_REAL* const*>(globalcoords);
    CCTK_REAL const* const* const abc =
      static_cast<CCTK_REAL const* const*>(localcoords);
    CCTK_REAL* const* const J =
      static_cast<CCTK_REAL* const*>(dadx);
    CCTK_REAL* const* const dJ =
      static_cast<CCTK_REAL* const*>(ddadxdx);
    
    if (CCTK_Equals(coordinate_system, "Cartesian")) {
      // a Cartesian one-patch system global coords = local coords

      if (globalcoords) {
#pragma omp parallel for collapse(2)
        for (int d=0; d<ndirs; ++d) {
          for (int ijk=0; ijk < npoints; ++ijk) {
            assert(patch[ijk] == 0);
            xyz[d][ijk] = abc[d][ijk];
          }
        }
      }
      
      if (dadx) {
#pragma omp parallel for collapse(3)
        for (int a=0; a<ndirs; ++a) {
          for (int b=0; b<ndirs; ++b) {
            for (int ijk=0; ijk < npoints; ++ijk) {
              J[a*ndirs+b][ijk] = a==b;
            }
          }
        }
      }

      if (ddadxdx) {
#pragma omp parallel for collapse(2)
        for (int d=0; d<ndirs*ndirs*(ndirs+1)/2; ++d) {
          for (int ijk=0; ijk < npoints; ++ijk) {
            dJ[d][ijk] = 0;
          }
        }
      }
      
    } else if (CCTK_Equals(coordinate_system, "TwoPatchCartesian")) {
      assert(not dadx);         // not supported
      assert(not ddadxdx);      // not supported
      
      if (globalcoords) {
#pragma omp parallel for
        for (int ijk=0; ijk < npoints; ++ijk) {
          assert(patch[ijk] == 0 or patch[ijk] == 1);
          for (int d=0; d<ndirs; ++d) {
            xyz[d][ijk] = abc[d][ijk];
          }
        }
      }
      
    } else if (CCTK_Equals(coordinate_system, "TwoPatchDistorted")) {
      assert(not ddadxdx);      // not supported
      
      if (globalcoords) {
#pragma omp parallel for
        for (int ijk=0; ijk < npoints; ++ijk) {
          CCTK_REAL globalcoord[ndirs];
          local_to_global_TwoPatchDistorted
            (patch[ijk], abc[0][ijk], abc[1][ijk], abc[2][ijk], globalcoord);
          for (int d=0; d<ndirs; ++d) {
            xyz[d][ijk] = globalcoord[d];
          }
        }
      }
      
      if (dadx) {
        assert(globalcoords);
#pragma omp parallel for
        for (int ijk=0; ijk < npoints; ++ijk) {
          CCTK_REAL JJ[ndirs][ndirs];
          dadx_TwoPatchDistorted
            (patch[ijk], xyz[0][ijk], xyz[1][ijk], xyz[2][ijk], JJ);
          for (int a=0; a<ndirs; ++a) {
            for (int b=0; b<ndirs; ++b) {
              J[a*ndirs+b][ijk] = JJ[a][b];
            }
          }
        }
      }
      
    } else if (CCTK_Equals(coordinate_system, "Thornburg04")) {
      assert(not ddadxdx);      // not supported
      
      if (globalcoords) {
#pragma omp parallel for
        for (int ijk=0; ijk < npoints; ++ijk) {
          CCTK_REAL globalcoord[3];
          local_to_global_Thornburg04
            (patch[ijk], abc[0][ijk], abc[1][ijk], abc[2][ijk], globalcoord);
          for (int d=0; d<ndirs; ++d) {
            xyz[d][ijk] = globalcoord[d];
          }
        }
      }

      if (dadx) {
        assert(globalcoords);
#pragma omp parallel for
        for (int ijk=0; ijk < npoints; ++ijk) {
          CCTK_REAL JJ[ndirs][ndirs];
          dadx_Thornburg04
            (patch[ijk], xyz[0][ijk], xyz[1][ijk], xyz[2][ijk], JJ);
          for (int a=0; a<ndirs; ++a) {
            for (int b=0; b<ndirs; ++b) {
              J[a*ndirs+b][ijk] = JJ[a][b];
            }
          }
        }
      }
      
    } else if (CCTK_Equals(coordinate_system, "Thornburg04nc")) {
      
      if (globalcoords) {
#pragma omp parallel for
        for (int ijk=0; ijk < npoints; ++ijk) {
          CCTK_REAL globalcoord[3];
          local_to_global_Thornburg04nc
            (patch[ijk], abc[0][ijk], abc[1][ijk], abc[2][ijk], globalcoord);
          for (int d=0; d<ndirs; ++d) {
            xyz[d][ijk] = globalcoord[d];
          }
        }
      }

      if (dadx) {
        assert(globalcoords);
#pragma omp parallel for
        for (int ijk=0; ijk < npoints; ++ijk) {
          CCTK_REAL JJ[ndirs][ndirs];
          dadx_Thornburg04nc
            (patch[ijk], xyz[0][ijk], xyz[1][ijk], xyz[2][ijk], JJ);
          for (int a=0; a<ndirs; ++a) {
            for (int b=0; b<ndirs; ++b) {
              J[a*ndirs+b][ijk] = JJ[a][b];
            }
          }
        }
      }

      if (ddadxdx) {
#pragma omp parallel for
        for (int ijk=0; ijk < npoints; ++ijk) {
          CCTK_REAL dJJ[ndirs][ndirs][ndirs];
          ddadxdx_Thornburg04nc
            (patch[ijk], xyz[0][ijk], xyz[1][ijk], xyz[2][ijk], dJJ);
          for (int d=0; d<ndirs; ++d) {
            for (int e=0; e<ndirs; ++e) {
              for (int f=e; f<ndirs; ++f) {
                const int idx =
                  d*ndirs*(ndirs+1)/2 + ((ndirs-1)*e - e*(e-1)/2) + f;
                dJ[idx][ijk] = dJJ[f][e][d];
              }
            }
          }
        }
      }
      
    } else {
      // transformation not yet implemented for this patch system, or
      // unknown patch system
      
      if (globalcoords) {
        CCTK_ERROR("LocalToGlobal: globalcoords calculation not implemented yet!");
      }
      return 1;
      
    }
    
    return 0;
  }
  
} // namespace Coordinates
