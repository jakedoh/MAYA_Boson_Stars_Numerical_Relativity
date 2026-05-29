
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
  Coordinates_SetJacobian_Cartesian(CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS_CHECKED(Coordinates_SetJacobian_Cartesian);
    DECLARE_CCTK_PARAMETERS;

    int i, j, k, ijk;

    CCTK_INFO("Calculating 'cartesian' jacobian components.");

    *general_coordinates = 0;
    *interpolate_boundary_points = 0;
    *jacobian_state = 1;
    *jacobian_derivative_state = 1;

    for (k=0; k<cctk_lsh[2]; ++k)
      for (j=0; j<cctk_lsh[1]; ++j)
        for (i=0; i<cctk_lsh[0]; ++i)
          {
            ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
  
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

    return;
  }

} // namespace Coordinates
