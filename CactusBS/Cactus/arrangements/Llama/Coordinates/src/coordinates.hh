
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

#ifndef COORDINATES_HH
#define COORDINATES_HH

#include "patchsystem.hh"



namespace Coordinates {

  // Good old ratio
  CCTK_REAL const PI = 3.14159265358979323846264338327950288419716939937510582097494;
  
  // TwoPatchDistorted

  enum { TPD_LEFT = 0,
	 TPD_RIGHT = 1
  };


  // Thornburg04

  enum { CENTRAL_CUBE   = 0,
         SPHERE_PLUS_X  = 1,
         SPHERE_MINUS_X = 2,
         SPHERE_PLUS_Y  = 3,
         SPHERE_MINUS_Y = 4,
         SPHERE_PLUS_Z  = 5,
         SPHERE_MINUS_Z = 6
  };

  // Thornburg13

  enum { CENTRAL13_CUBE   = 0,
         SPHERE13_PLUS_X  = 1,
         SPHERE13_MINUS_X = 2,
         SPHERE13_PLUS_Y  = 3,
         SPHERE13_MINUS_Y = 4,
         SPHERE13_PLUS_Z  = 5,
         SPHERE13_MINUS_Z = 6,
         SPHERE13_PLUS_X_OUTER  = 7,
         SPHERE13_MINUS_X_OUTER = 8,
         SPHERE13_PLUS_Y_OUTER  = 9,
         SPHERE13_MINUS_Y_OUTER = 10,
         SPHERE13_PLUS_Z_OUTER  = 11,
         SPHERE13_MINUS_Z_OUTER = 12
  };

  // Thornburg04nc

  enum { NC_SPHERE_PLUS_X  = 0,
         NC_SPHERE_MINUS_X = 1,
         NC_SPHERE_PLUS_Y  = 2,
         NC_SPHERE_MINUS_Y = 3,
         NC_SPHERE_PLUS_Z  = 4,
         NC_SPHERE_MINUS_Z = 5
  };

  // A cylinder in a cube

  enum { CIB_BOX      = 0,
         CIB_CYLINDER = 1
  };

  // A sphere with columns covering the poles

  enum { SPC_SPHERE         = 0,
         SPC_COLUMN_PLUS_Z  = 1,
         SPC_COLUMN_MINUS_Z = 2
  };

  // A cylinder with a column covering the center

  enum { CC_CYLINDER        = 0,
         CC_COLUMN          = 1
  };


  extern "C"
  CCTK_INT
  global_to_local_TwoPatchCartesian
  (CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL localcoord[ndirs]);
  
  extern "C"
  void
  dadx_TwoPatchCartesian
  (CCTK_INT patch, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z,
   CCTK_REAL J[ndirs][ndirs]);
  
  extern "C"
  void
  ddadxdx_TwoPatchCartesian
  (CCTK_INT patch, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z,
   CCTK_REAL dJ[ndirs][ndirs][ndirs]);


  
  extern "C"
  CCTK_INT
  global_to_local_TwoPatchDistorted
  (CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL localcoord[ndirs]);
  
  extern "C"
  void
  local_to_global_TwoPatchDistorted
  (CCTK_INT patch, CCTK_REAL a, CCTK_REAL b, CCTK_REAL c,
   CCTK_REAL globalcoord[ndirs]);
  
  extern "C"
  void
  dadx_TwoPatchDistorted
  (CCTK_INT patch, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z,
   CCTK_REAL J[ndirs][ndirs]);
  
  extern "C"
  void
  ddadxdx_TwoPatchDistorted
  (CCTK_INT patch, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z,
   CCTK_REAL dJ[ndirs][ndirs][ndirs]);


  
  
  extern "C"
  CCTK_INT
  global_to_local_Thornburg04
  (CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL localcoord[ndirs]);
  
  extern "C"
  void
  local_to_global_Thornburg04
  (CCTK_INT patch, CCTK_REAL a, CCTK_REAL b, CCTK_REAL c,
   CCTK_REAL globalcoord[ndirs]);
  
  extern "C"
  void
  dadx_Thornburg04
  (CCTK_INT patch, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z,
   CCTK_REAL J[ndirs][ndirs]);
  
  extern "C"
  void
  ddadxdx_Thornburg04
  (CCTK_INT patch, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z,
   CCTK_REAL dJ[ndirs][ndirs][ndirs]);



  extern "C"
  CCTK_INT
  global_to_local_Thornburg13
  (CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL localcoord[ndirs]);
  
  extern "C"
  void
  dadx_Thornburg13
  (CCTK_INT patch, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z,
   CCTK_REAL J[ndirs][ndirs]);
  
  extern "C"
  void
  ddadxdx_Thornburg13
  (CCTK_INT patch, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z,
   CCTK_REAL dJ[ndirs][ndirs][ndirs]);




  extern "C"
  CCTK_INT
  global_to_local_Thornburg04nc
  (CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL localcoord[ndirs]);
  
  extern "C"
  void
  local_to_global_Thornburg04nc
  (CCTK_INT patch, CCTK_REAL a, CCTK_REAL b, CCTK_REAL c,
   CCTK_REAL globalcoord[ndirs]);

  extern "C"
  extern "C"
  void
  dadx_Thornburg04nc
  (CCTK_INT patch, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z,
   CCTK_REAL J[ndirs][ndirs]);
  
  extern "C"
  void
  ddadxdx_Thornburg04nc
  (CCTK_INT patch, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z,
   CCTK_REAL dJ[ndirs][ndirs][ndirs]);



  extern "C"
  CCTK_INT
  global_to_local_CylinderInBox
  (CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL localcoord[ndirs]);
  
  extern "C"
  void
  dadx_CylinderInBox
  (CCTK_INT patch, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z,
   CCTK_REAL J[ndirs][ndirs]);
  
  extern "C"
  void
  ddadxdx_CylinderInBox
  (CCTK_INT patch, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z,
   CCTK_REAL dJ[ndirs][ndirs][ndirs]);



  extern "C"
  CCTK_INT
  global_to_local_SphereColumn
  (CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL localcoord[ndirs]);
  
  extern "C"
  void
  dadx_SphereColumn
  (CCTK_INT patch, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z,
   CCTK_REAL J[ndirs][ndirs]);
  
  extern "C"
  void
  ddadxdx_SphereColumn
  (CCTK_INT patch, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z,
   CCTK_REAL dJ[ndirs][ndirs][ndirs]);



  extern "C"
  CCTK_INT
  global_to_local_CylinderColumn
  (CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL localcoord[ndirs]);
  
  extern "C"
  void
  dadx_CylinderColumn
  (CCTK_INT patch, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z,
   CCTK_REAL J[ndirs][ndirs]);
  
  extern "C"
  void
  ddadxdx_CylinderColumn
  (CCTK_INT patch, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z,
   CCTK_REAL dJ[ndirs][ndirs][ndirs]);



} // namespace

#endif
